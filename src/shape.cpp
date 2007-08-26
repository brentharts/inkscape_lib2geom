#include "shape.h"
#include "utils.h"
#include "sweep.h"

#include <iostream>
#include <algorithm>

namespace Geom {

//Utility funcs
//yes, xor is !=, but I'm pretty sure this is safer in the event of strange bools
bool logical_xor (bool a, bool b) { return (a || b) && !(a && b); }

bool disjoint(Path const & a, Path const & b) {
    return !contains(a, b.initialPoint()) && !contains(b, a.initialPoint());
}

template<typename T>
void append(T &a, T const &b) {
    a.insert(a.end(), b.begin(), b.end());
}

void first_false(std::vector<std::vector<bool> > visited, unsigned &i, unsigned &j) {
    for(i = 0, j = 0; i < visited.size(); i++) {
        std::vector<bool>::iterator unvisited = std::find(visited[i].begin(), visited[i].end(), false);
        if(unvisited != visited[i].end()) {
            j = unvisited - visited[i].begin();
            break;
        }
    }
}

unsigned find_crossing(Crossings const &cr, Crossing x, unsigned i) {
    return std::lower_bound(cr.begin(), cr.end(), x, CrossingOrder(i)) - cr.begin();
}

/* This function handles boolean ops on shapes.  The first parameter is a bool
 * which determines its behavior in each combination of these cases.  For proper
 * fill information and noncrossing behavior, the fill data of the regions must
 * be correct.  The boolean parameter determines whether the operation is a
 * union or a subtraction.  Reversed paths represent inverse regions, where
 * everything is included in the fill except for the insides.
 *
 * Here is a chart of the behavior under various circumstances:
 * 
 * rev = false (union)
 *            A
 *        F         H
 * F  A+B -> F  A-B -> H
 *B
 * H  B-A -> H  AxB -> H
 *
 * rev = true (intersect)
 *            A
 *        F         H
 * F  AxB -> F  B-A -> F
 *B
 * H  A-B -> F  A+B -> H
 *
 * F/H = Fill outer / Hole outer
 * A/B specify operands
 * + = union, - = subtraction, x = intersection
 * -> read as "produces"
 *
 * The operation of this function isn't very complicated.  It just traverses the crossings, and uses the crossing direction to decide whether the next segment should be from A or from B.
 */
Shape shape_boolean(bool rev, Shape const & a, Shape const & b, CrossingSet const & crs) {
    const Regions ac = a.content, bc = b.content;

    //Keep track of which crossings we've hit.
    std::vector<std::vector<bool> > visited;
    for(unsigned i = 0; i < crs.size(); i++)
        visited.push_back(std::vector<bool>(crs[i].size(), false));

    //Traverse the crossings, creating chunks
    Regions chunks;
    while(true) {
        unsigned i, j;
        first_false(visited, i, j);
        if(i == visited.size()) break;

        Path res;
        do {
            Crossing cur = crs[i][j];
            visited[i][j] = true;
            //get indices of the dual:
            unsigned io = cur.getOther(i), jo = find_crossing(crs[io], cur, io);
            if(jo < visited[io].size()) visited[io][jo] = true;
            
            if(logical_xor(cur.dir, rev)) {
                if(i >= ac.size()) { i = io; j = jo; }
                j++;
                if(j >= crs[i].size()) j = 0;
                Crossing next = crs[i][j];
                ac[next.a].boundary.appendPortionTo(res, cur.ta, next.ta);
            } else {
                if(i < ac.size()) { i = io; j = jo; }
                j++;
                if(j >= crs[i].size()) j = 0;
                Crossing next = crs[i][j];
                bc[next.b - ac.size()].boundary.appendPortionTo(res, cur.tb, next.tb);
            }
        } while (!visited[i][j]);
        chunks.push_back(Region(res));
    }
    
    //If true, then the hole must be inside the other to be included
    bool const on_sub = logical_xor(a.fill, b.fill);
    bool const a_mode = logical_xor(logical_xor(!rev, a.fill), on_sub),
               b_mode = logical_xor(logical_xor(!rev, b.fill), on_sub);
    //Handle unintersecting portions
    for(unsigned i = 0; i < crs.size(); i++) {
        if(crs[i].size() == 0) {
            Region r;
            bool mode;
            if(i < ac.size()) {
                r = ac[i];
                mode = a_mode;
            } else {
                r = bc[i - ac.size()];
                mode = b_mode;
            }
            
            if(logical_xor(r.fill, i < ac.size() ? a.fill : b.fill)) {
                //is an inner
                Point exemplar = r.boundary.initialPoint();
                Regions const & others = i < ac.size() ? bc : ac;
                for(unsigned j = 0; j < others.size(); j++) {
                    if(others[j].contains(exemplar)) {
                        //contained in other
                        if(mode) chunks.push_back(r);
                        goto skip;
                    }
                }
            }
            //disjoint
            if(!mode) chunks.push_back(r);
            skip: (void)0;
        }
    }
    
    return Shape(chunks);
}

//Returns a vector of crossings, such that those associated with B are in the range [a.size(), a.size() + b.size())
CrossingSet crossings_between(Shape const &a, Shape const &b) { 
    CrossingSet results(a.content.size() + b.content.size(), Crossings());
    
    //TODO: sweep
    for(unsigned i = 0; i < a.content.size(); i++) {
        for(unsigned jx = 0; jx < b.content.size(); jx++) {
            unsigned j = jx + a.content.size();
            Crossings cr = crossings(a.content[i].getBoundary(), b.content[jx].getBoundary());
            for(unsigned k = 0; k < cr.size(); k++) { cr[k].a = i; cr[k].b = j; }
            //Sort & add I crossings
            sort_crossings(cr, i);
            Crossings n(results[i].size() + cr.size());
            std::merge(results[i].begin(), results[i].end(), cr.begin(), cr.end(), n.begin(), CrossingOrder(i));
            results[i] = n;
            //Sort & add J crossings
            sort_crossings(cr, j);
            n.resize(results[j].size() + cr.size());
            std::merge(results[j].begin(), results[j].end(), cr.begin(), cr.end(), n.begin(), CrossingOrder(j));
            results[j] = n;
        }
    }
    return results;
}

Shape shape_boolean(bool rev, Shape const & a, Shape const & b) {
    CrossingSet crs = crossings_between(a, b);
    
    /* for(unsigned i = 0; i < crs.size(); i++) {
        std::cout << i << "\n";
        for(unsigned j = 0; j < crs[i].size(); j++) {
            std::cout << " " << crs[i][j].a << " " << crs[i][j].b << " :" << crs[i][j].ta << " to " << crs[i][j].tb << "\n";
        }
    } */
    
    return shape_boolean(rev, a, b, crs);
}

Shape shape_boolean(Shape const &a, Shape const &b, unsigned flags) {
    flags &= 15;
    if(flags > SHAPE_UNION) return shape_boolean(a, b, ~flags).inverse();
    //trivial cases
    switch(flags) {
        case SHAPE_NULL:         return Shape();
        case SHAPE_INTERSECT:    return shape_boolean(true, a, b);
        case SHAPE_SUBTRACT_A_B: return shape_boolean(true, a, b.inverse());
        case SHAPE_IDENTITY_A:   return a;
        case SHAPE_SUBTRACT_B_A: return shape_boolean(true, b, a.inverse());
        case SHAPE_IDENTITY_B:   return b;
        case SHAPE_EXCLUSION: {
            Shape res = shape_boolean(true, a, b.inverse());
            append(res.content, shape_boolean(true, b, a.inverse()).content);
            return res;
        }
        case SHAPE_UNION:        return shape_boolean(false, a, b);
    }
}

int paths_winding(std::vector<Path> const &ps, Point p) {
    int ret;
    for(unsigned i = 0; i < ps.size(); i++)
        ret += winding(ps[i], p);
    return ret;
}

//sanitize
//We have two phases, one for each winding direction.
Shape sanitize_paths(std::vector<Path> const &ps, bool evenodd) {
    CrossingSet crs = crossings_among(ps);
    
    //two-phase process - g
    Regions chunks;
    for(bool phase = 0; phase < 2; phase++) {
        
        //Keep track of which crossings we've hit.
        std::vector<std::vector<bool> > visited;
        for(unsigned i = 0; i < crs.size(); i++)
            visited.push_back(std::vector<bool>(crs[i].size(), false));

        while(true) {
            unsigned i, j;
            first_false(visited, i, j);
            if(i == visited.size()) break;
            
            bool use = paths_winding(ps, ps[i].initialPoint()) % 2 == 1;
            Path res;
            do {
                Crossing cur = crs[i][j];
                visited[i][j] = true;
                
                //get indices of the dual:
                i = cur.getOther(i), j = find_crossing(crs[i], cur, i);
                if(j < visited[i].size()) visited[i][j] = true;
                
                if(logical_xor(phase, cur.dir)) {
                    // forwards
                    j++;
                    if(j >= crs[i].size()) j = 0;
                } else {
                    // backwards
                    if(j == 0) j = crs[i].size() - 1; else j--;
                }
                if(use) {
                    Crossing next = crs[i][j];
                    ps[i].appendPortionTo(res, cur.ta, next.ta);
                }
            } while(!visited[i][j]);
            
            if(use) {
                chunks.push_back(Region(res, true));
            }
        }
    }
    return Shape(chunks);
}

Shape Shape::operator*(Matrix const &m) const {
    Shape ret;
    for(unsigned i = 0; i < content.size(); i++)
        ret.content.push_back(content[i] * m);
    ret.fill = fill;
    return ret;
}

// inverse is a boolean not
Shape Shape::inverse() const {
    Shape ret;
    for(unsigned i = 0; i < content.size(); i++)
        ret.content.push_back(content[i].inverse());
    ret.fill = !fill;
    return ret;
}

struct ContainmentOrder {
    std::vector<Region> const *rs;
    explicit ContainmentOrder(std::vector<Region> const *r) : rs(r) {}
    bool operator()(unsigned a, unsigned b) const { return (*rs)[b].contains((*rs)[a]); }
};

bool Shape::contains(Point const &p) const {
    std::vector<unsigned> ixs;
    //TODO: sweep
    for(unsigned i = 0; i < content.size(); i++) {
        if(content[i].boundary.initialPoint() == p) continue; //bad hack for inside_invariants
        if(content[i].contains(p)) ixs.push_back(i);
    }
    if(ixs.size() == 0) return !fill;
    return content[*min_element(ixs.begin(), ixs.end(), ContainmentOrder(&content))].isFill();
}

bool Shape::inside_invariants() const {  //semi-slow & easy to violate
    for(unsigned i = 0; i < content.size(); i++)
        if( logical_xor(content[i].isFill(), contains(content[i].boundary.initialPoint())) ) return false;
    return true;
}
bool Shape::region_invariants() const { //semi-slow
    for(unsigned i = 0; i < content.size(); i++)
        if(!content[i].invariants()) return false;
    return true;
}
bool Shape::cross_invariants() const { //slow
    CrossingSet crs = crossings_among(paths_from_regions(content));
    for(unsigned i = 0; i < crs.size(); i++)
        if(!crs[i].empty()) return false;
    return true;
}

bool Shape::invariants() const {
    return inside_invariants() && region_invariants() && cross_invariants();
}

}