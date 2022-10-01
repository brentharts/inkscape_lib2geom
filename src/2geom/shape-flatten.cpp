/** @file Flatten algorithm – uncross and regularize a path-vector bounding a shape.
 */
/*
 * Authors:
 *   Rafał Siejakowski <rs@rs-math.net>
 *
 * Copyright 2022 the Authors
 *
 * This library is free software; you can redistribute it and/or
 * modify it either under the terms of the GNU Lesser General Public
 * License version 2.1 as published by the Free Software Foundation
 * (the "LGPL") or, at your option, under the terms of the Mozilla
 * Public License Version 1.1 (the "MPL"). If you do not alter this
 * notice, a recipient may use your version of this file under either
 * the MPL or the LGPL.
 *
 * You should have received a copy of the LGPL along with this library
 * in the file COPYING-LGPL-2.1; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * You should have received a copy of the MPL along with this library
 * in the file COPYING-MPL-1.1
 *
 * The contents of this file are subject to the Mozilla Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * This software is distributed on an "AS IS" basis, WITHOUT WARRANTY
 * OF ANY KIND, either express or implied. See the LGPL or the MPL for
 * the specific language governing rights and limitations.
 */

#include <set>
#include <vector>

#include <2geom/math-utils.h>
#include <2geom/pathvector.h>
#include <2geom/sweeper.h>
#include <2geom/shape.h>

#include "planar-graph.h"

#define VERBOSE 0

#if VERBOSE // Extremely verbose mode
    #include <iostream>
    #include <2geom/svg-path-writer.h>
    #define LOG_DEBUG(s) (std::cout << s << std::endl)
#else
    #define LOG_DEBUG(s)
#endif

namespace Geom {

namespace { // anonymous namespace containing primitives needed by the flatten algorithm.

/** @brief An edge label for use with PlanarGraph class in the flatten algorithm.
 *
 * The label keeps track of winding numbers to the left and to the right of the edge.
 * We store only the winding number to our right, whereas the winding number on the
 * left differs from it by the "weight" of the edge.
 */
class FlattenLabel
{
private:
    /// Winding number of the input path-vector about points immediately to the right of this edge.
    int _winding_right = 0;
    /** Coefficient of this edge in the integral singular 1-chain given by the
     * input path-vector. The initial value of 1 can only be changed by the
     * PlanarGraph callbacks. */
    int _weight = 1;

public:
    PathExtrema extrema_x; ///< Extrema in the x-axis direction.
    /// Indices of paths in the input path-vector containing this edge as a subarc.
    std::set<unsigned> path_indices;
    bool has_windings = false; ///< Whether the winding numbers have been calculated.
    bool processed    = false; ///< Whether the edge was included as a part of the boundary or rejected.
    bool active_start = false; ///< Whether the edge is the starting edge of the current walk.

    /// Create a new label for a portion of an input path with the given index.
    FlattenLabel(unsigned origin_index) { path_indices.insert(origin_index); }
    FlattenLabel(FlattenLabel const &) = default;
    FlattenLabel &operator=(FlattenLabel const &) = default;

    // Getters and setters for winding numbers on either side.
    int windingOnRight() const { return _winding_right; }
    int windingOnLeft() const { return _winding_right + _weight; }
    int windingOnRL(bool on_right) const { return on_right ? windingOnRight() : windingOnLeft(); }
    void setWindingOnRight(int wr)
    {
        _winding_right = wr;
        has_windings = true;
    }
    void setWindingOnLeft(int wl)
    {
        _winding_right = wl - _weight;
        has_windings = true;
    }
    void setWindingRL(bool on_right, int winding)
    {
        on_right ? setWindingOnRight(winding) : setWindingOnLeft(winding);
    }
    int getWeight() const { return _weight; }

    /** @brief Determine whether the edge with this label is a part of the conventional
     *         boundary of the shape (see FillRule::CONVENTIONAL).
    *
    * @param fr The fill rule which defines the shape.
    * @return +1 if the edge, with its current orientation, is a part of the boundary;
    *         -1 if the edge is part of the boundary but should be reversed as per the convention;
    *          0 if the edge is not part of the boundary.
    */
    int8_t isBoundary(FillRule fr)
    {
        if (!_weight) {
            return 0; // We contribute zero to the 1-chain.
        }
        switch (fr)
        {
            case FillRule::CONVENTIONAL:
                return sgn(_weight);

            case FillRule::EVEN_ODD:
                if (_weight % 2 == 0) {
                    return 0; // Even weight never separates inside from outside.
                }
                // If the weight is odd, we must return +1 when _winding_right is even
                // and -1 when _winding_right is odd; we do it branchlessly.
                return 1 - ((_winding_right & 0x1) * 2);

            case FillRule::NONZERO: {
                auto const wl = windingOnLeft();
                if (wl && !_winding_right) {
                    return 1;
                } else if (_winding_right && !wl) {
                    return -1;
                }
                return 0;
            }
        }
        return 0;
    }

    /// Add another path's indices to our set of path indices.
    void unionIndices(FlattenLabel const &other)
    {
        path_indices.insert(other.path_indices.begin(), other.path_indices.end());
    }

    /// PlanarGraph callbacks.
    void onReverse()
    {
        setWindingOnRight(windingOnLeft()); // Swap right and left.
        _weight = -_weight; // Keep the 1-chain unchanged.
    }
    void onMergeWith(FlattenLabel const &other)
    {
        unionIndices(other);
        _weight += other._weight; // Combine like terms in the 1-chain.
    }
    void onDetach() { _weight = 0; }
};

using FlattenGraph = PlanarGraph<FlattenLabel>;

/** @brief A generator class yielding pieces of a PathVector between self-intersections.
 *
 * Construct an instance by passing a reference to a path-vector and the desired precision.
 * A range-for loop iteration over the generator object will then yield pairs (std::pair)
 * of the form (index, portion), where `portion` is a piece of a path from the path-vector
 * and index is the original path's index in the path-vector. However, components without
 * any intersections are not chopped or included in the generated values; you must handle
 * them separately.
 */
class PathVectorCutter
{
private:
    PathVector const &_subject;
    std::vector< std::pair<Point, PathVectorTime> > _cuts;
    size_t _length_estimate = 0;

public:

    /** @brief Construct a path-vector cutter for the given path-vector.
     *
     * @param to_cut The path-vector, consisting of closed paths only, which is to
     *               be cut up at its self-intersections.
     * @param precision Numerical precision for the self-intersection search.
     *
     * The instance stores a reference to the path-vector (first argument),
     * so the owner must guarantee that the reference remains valid for as
     * long as the generator is used.
     */
    PathVectorCutter(PathVector const &to_cut, Coord precision)
        : _subject{to_cut}
    {
        // Find self-intersections in the path-vector and for each of them,
        // store both crossing times, unless they are identical.
        auto const xings = _subject.intersectSelf(precision);
        _length_estimate = 2 * xings.size(); // This is quite generous.
        _cuts.reserve(_length_estimate);
        for (auto const &xing : xings) {
            _cuts.emplace_back(std::make_pair(xing.point(), xing.first));
            if (xing.second != xing.first) {
                _cuts.emplace_back(std::make_pair(xing.point(), xing.second));
            }
        }

        // Sort by the PathVector times of the cuts.
        std::sort(_cuts.begin(), _cuts.end(), [](std::pair<Point, PathVectorTime> const &a,
                                                 std::pair<Point, PathVectorTime> const &b) -> bool {
                                                    return a.second < b.second;
                                              });
        LOG_DEBUG("PVC: A total of " << _cuts.size() << " cut locations, (precision = " << precision << ").");
    }

    /// Iterator for going over the generator object.
    class Iterator
    {
    private:
        using InternalIt = std::vector< std::pair<Point, PathVectorTime> >::iterator;
        InternalIt _it;
        InternalIt _first_on_path;
        InternalIt const _end;
        PathVector const &_pv;

        Iterator(InternalIt iterator, InternalIt fop, InternalIt end, PathVector const& pv)
            : _it{iterator}
            , _first_on_path{fop}
            , _end{end}
            , _pv{pv}
        {}

    public:
        Iterator operator++()
        {
            auto const old_path = _it->second.path_index;
            if (++_it != _end && _it->second.path_index != old_path) {
                _first_on_path = _it;
            }
            return *this;
        }

        bool operator==(Iterator const &other) const { return _it == other._it && &_pv == &other._pv; }
        bool operator!=(Iterator const &other) const { return _it != other._it || &_pv != &other._pv; }

        std::pair<unsigned, Path> operator*() const
        {
            auto const index = _it->second.path_index;
            auto next = std::next(_it);

            // Handle the wraparound to the first intersection on the path.
            bool is_last_on_path = (next == _end) || (next->second.path_index != index);
            auto const &from = _it->second;
            auto const until_it = is_last_on_path ? _first_on_path : next;
            auto const &until = until_it->second;

            Path path; ///< The resulting path portion.

            auto const START = PathVectorTime(index, PathTime(0, 0.0));
            auto const END = PathVectorTime(index, PathTime(_pv[index].size() - 1, 1.0));
            if (is_last_on_path && _it == until_it && (from == START || from == END)) {
                // Special case: only one cut location, occurring at the start/end of
                // a closed path. In this case, we open up the path instead of cutting it.
                path = _pv[index];
                auto const &closing_seg = path.closingSegment();
                if (closing_seg.isDegenerate()) { // Simply open up the path
                    path.close(false);
                } else { // Opening up deletes the closing segment; we readd it manually.
                    auto real_closing_segment = closing_seg;
                    path.close(false);
                    path.append(real_closing_segment);
                }
            } else {
                path = _pv[index].portion(from.asPathTime(), until.asPathTime(), is_last_on_path);
                // Snap the path's ends to the exact intersection locations.
                path.setInitial(_it->first);
                if (!path.empty()) {
                    path.setFinal(until_it->first);
                }
            }
            return std::make_pair(index, std::move(path));
        }
        friend class PathVectorCutter;
    };

    /** Return an approximate length of the generator. */
    size_t getLengthEstimate() { return _length_estimate; }

    /** Get iterators to the beginning and past the end of the generator's range. */
    Iterator begin() { return Iterator(_cuts.begin(), _cuts.begin(), _cuts.end(), _subject); }
    Iterator end()
    {
        auto const end = _cuts.end();
        return Iterator(end, end, end, _subject);
    }
 };

/** @brief Return a PathVector containing only closed paths by either filtering
 * out or closing open paths. The paths in the returned path-vector are cleaned
 * of degenerate segments.
 *
 * @param pv The source path-vector.
 * @param close_opens Whether open paths should be closed.
 * @return A new path-vector in which open paths have been either closed or removed.
 */
PathVector clean_paths(PathVector const &pv, bool close_opens)
{
    PathVector result;
    for (auto const &path : pv) {
        if (path.closed()) {
            result.push_back(path.withoutDegenerateCurves());
        } else if (close_opens) {
            auto to_close = path.withoutDegenerateCurves();
            to_close.close();
            result.push_back(std::move(to_close));
        }
    }
    return result;
}

/** @brief Return the algebraic sum of winding numbers of a subset
 *         of a path-vector's paths about the specified point.
 *
 * @param pt The point about which the winding number will be calculated.
 * @param pv The path-vector from which to source the paths.
 * @param blacklist A set of indices of paths that will be excluded from the summation.
 * @return The algebraic sum of winding numbers of paths of pv about pt,
 *         excluding the paths whose indices are in the blacklist.
 */
int partial_winding_number(Point const &pt, PathVector const &pv, std::set<unsigned> const &blacklist)
{
    int result = 0;
    for (unsigned i = 0; i < pv.size(); i++) {
        // TODO: After C++20, replace with !blacklist.contains(i)
        if (std::find(blacklist.begin(), blacklist.end(), i) == blacklist.end()) {
            result += pv.at(i).winding(pt);
        }
    }
    return result;
}

/** @page flatten_description Description of the flatten algorithm
 *
 * This is a detailed description of the Shape::flatten algorithm.  The flatten
 * algorithm assembles the boundary components of a 2D shape (orienting them so as
 * to follow the FillRule::CONVENTIONAL fill rule), operating on a pre-regularized
 * \c FlattenGraph.
 *
 * @section flatten_what What does it do?
 * The input of the flatten algorithm is a description of a 2D shape ("fill
 * region") via a PathVector and a fill rule. This description could be quite
 * complicated: there can be intersections between the constituent paths, or even
 * overlaps between them, and the shape can have an arbitrary number of connected
 * components.  It could have holes in it and other pieces inside those holes, to
 * an arbitrary level of nesting.
 *
 * The goal of the flatten algorithm is to find a PathVector which bounds that
 * same fill region, but is greatly simplified. Specifically, we want it to have
 * the following property: if the y-axis points up, then every path in the
 * PathVector has the interior of the shape on its left and the exterior on its
 * right, when travelling in the direction of increasing time coordinate on the
 * path.  Throughout this description, we imagine that the y-axis points up (as in
 * mathematics).
 *
 * The consequences of the aforementioned property are as follows:
 * \li The PathVector satisfies the FillRule::CONVENTIONAL fill rule.
 * \li Every piece of every path in the PathVector separates inside from outside.
 * \li The PathVector has an algebraic (total) winding number of 0 or 1 about
 *     every complementary region in the plane, which means that converting to
 *     FillRule::NONZERO and FillRule::EVEN_ODD is a no-op.
 * \li If the shape's interior is empty, the boundary PathVector is also empty.
 * \li The PathVector does not contain transverse self-intersections.
 *
 * As you can see, the flatten algorithm is quite useful. In particular,
 * performing boolean operations on shapes becomes much easier when their contours
 * are flattened in this way, not even requiring a sweepline approach.
 *
 * @section flatten_how How does it work?
 * The algorithm starts by finding self-intersections in the input path-vector and
 * feeding the path pieces between consecutive self-intersections to the
 * \c FlattenGraph. If there are paths without intersections, they get inserted
 * into the graph as "detached edges" (free loops), not connected to any vertex.
 *
 * Then, we regularize the graph and compute the extrema of the x-coordinate
 * for each of its edges.
 *
 * After this is done, we run the actual sweepline algorithm using the
 * \c FlattenSweep class as the SweepSet (defined in shape-flatten.cpp).
 * The sweepline is vertical and moves from the right towards the left.
 *
 * When the sweepline hits a new valid edge that hasn't been processed yet,
 * the \c FlattenSweep starts walking the graph from that edge, appending to the
 * internally stored output PathVector the correctly oriented boundary
 * component that the encountered edge belongs to (if it actually belongs to any).
 * Visited edges are marked as processed and further sweepline events involving
 * them are ignored.
 *
 * Each walk starts with the determination of the ambient winding number, i.e.,
 * the winding number of the original path-vector around the points "just outside"
 * of an "outer" edge of the component that the sweepline has just hit.  There
 * are two problems here: determining where "just outside" is relative to the
 * encountered edge, and computing that "ambient" winding number.  Since the
 * sweepline is vertical and moves left (towards negative X), it is natural to
 * look at a rightmost point (which is not unique) on the edge.  The exact
 * strategy depends on the nature of this rightmost "witness" point.
 *
 * When the witness point is somewhere in the middle of the path (whether a node
 * or not), then we check if the y-coordinate increases or decreases as the path
 * passes this point; this will tell us whether the rightmost path goes "up" or
 * "down" and, correspondingly, whether it has the "outside" on its right or left.
 *
 * When the rightmost point is a vertex of the graph, then we have no guarantee
 * that the edge is the "outermost" one, because there could be several paths
 * emanating from the vertex roughly "towards the left" but at different angles,
 * like a fan. We need to be careful to pick an outermost edge from such fan.
 * Therefore, we must inspect the vertex (taking advantage of the additional
 * multple tangency resolution already performed by \c FlattenGraph), so that the
 * edge we pick is the first edge encountered when going counter-clockwise around
 * the vertex from an infinite ray emanating from the vertex to the right.  This
 * edge is guaranteed to touch the ambient region, whence the sweepline is
 * arriving. Whether the edge has the ambient region on its left or right depends
 * on the edge's orientation, which the vertex knows about.
 *
 * Note that the "outside" of the component is not necessarily the outside of the
 * shape; there could be a hole in the shape, another piece inside the hole, a
 * hole in a piece in a hole, and so on. Hence, the computation of the "ambient"
 * winding number for the component must take into account the entirety of the
 * input path-vector.
 *
 * When the witness point that the sweepline hits is not an endpoint of the edge,
 * we can simply inspect the edge's label to see which paths in the original
 * path-vector contained arcs that gave rise to this edge. The ambient winding
 * number is then computed as the sum of winding numbers of paths in the original
 * path-vector \b excluding the paths that the edge originates from, about the
 * rightmost "witness" point that the sweepline has just hit.
 *
 * On the other hand, when the witness point is an endpoint, i.e., a vertex of the
 * \c FlattenGraph, we must exclude from the winding calculation all paths that
 * gave rise to all of the edges incident to that vertex, since they all "pass"
 * through the rightmost point (i.e., the vertex). It is of utmost importance to
 * never compute a winding number of a path about a point lying on the path, since
 * the return value of Path::winding() cannot be trusted in this case. (And for a
 * very good reason: the mathematically correct answer is that the winding number
 * of a path about a point on that path is undefined. If you insist on defining
 * it, the most reasonable answer is infinity.)
 *
 * Once this "ambient winding number" is determined, we know the winding numbers
 * of the regions immediately to the left and right of the starting edge for this
 * walk. We then extend this information as we explore the graph to continue
 * assembling the boundary component. This general procedure is inspired by the
 * Livarot library, though the implementation details are quite different.
 *
 * In the process of traversing the component, we must make a decision where to go
 * next at each vertex of the \c FlattenGraph. Since we walk counter-clockwise
 * along the contour of the shape (or clockwise around a hole), always keeping the
 * shape on our left, we make the leftmost possible turn at each vertex, so as to
 * enclose "as little of the shape" as possible. When turning, we ignore any edges
 * that aren't part of the shape's boundary due to not forming a predicate
 * boundary for the given fill rule (predicate boundary is the set of edges such
 * that the fill rule is satisfied on one side but not on the other side).
 *
 * As we traverse the contour, we mark as processed the edges that were picked
 * as pieces of the contour, as well as those that were rejected due to not
 * being predicate boundaries. When we reach a vertex where the leftmost
 * possible turn takes us back onto the starting edge of the walk, then we know
 * that we have completed the circuit and can close the path, adding it to the
 * internally stored output path-vector \c _out.
 *
 * The algorithm finishes when all edges have been processed.
 */
/** @brief A sweep-set which executes the flatten algorithm.
 *
 * The FlattenSweep passes a sweepline through a non-owned FlattenGraph
 * in a regularized state. The graph contains the cut up pieces of a path-vector
 * that were inserted as edges, but may have been modified or glued by
 * FlattenGraph::regularize(). In the course of the sweepline's passage, the
 * FlattenSweep re-orients and assembles regularized boundary components of
 * the shape. See @ref flatten_description for more details.
 */
class FlattenSweep
{
private:
    /// The PlanarGraph on which the flatten algorithm operates.
    FlattenGraph const &_graph;
    /// The input path-vector of the flatten algorithm.
    PathVector const &_original_pathvector;
    FillRule const _fill_rule; ///< The input fill rule of the flatten algorithm.
    PathVector _out; ///< The in-progress flatten output.

public:
    using ItemIterator = FlattenGraph::EdgeConstIterator;

    /** Construct a sweep-set for a given graph with precomputed x-bounds of edges. */
    FlattenSweep(FlattenGraph const &graph, PathVector const &original_pv, FillRule fr)
        : _graph{graph}
        , _original_pathvector{original_pv}
        , _fill_rule{fr}
    {
        assert(graph.isRegularized());
    }

    /** Get the output path-vector. */
    PathVector &&getResult() { return std::move(_out); }

    // === SweepSet API ===
    auto &items() const { return _graph.getEdges(); }

    Interval itemBounds(ItemIterator edge) const
    {
        return Interval(edge->label.extrema_x.min_point[X], edge->label.extrema_x.max_point[X]);
    }

    /// Callback for when the sweepline starts intersecting a new edge of the graph.
    void addActiveItem(ItemIterator edge)
    {
        if (edge->label.processed) {
            return;
        }
        if (edge->detached && !edge->inserted_as_detached) {
            return;
        }

        // Process the boundary component containing this edge.
        _processFromEdge(edge);
    }

    void removeActiveItem(ItemIterator to_remove)
    {
        // No need to do anything when the sweepline stops intersecting an item.
    }
    // ===

private:
    /** @brief Start walking a FlattenGraph from the given edge.
     *
     * The walk will only happen if the edge is a predicate boundary. Of course,
     * the walk can only visit the connected component of the edge in the graph.
     *
     * @param edge The edge where processing will start.
     */
    void _processFromEdge(ItemIterator const &edge)
    {
        auto *initial_edge = &*edge; ///< The edge where component traversal should start.
        /// The rightmost "witness" point where the sweepline hits the path of the edge (conceptually).
        auto const &witness = edge->label.extrema_x.max_point;

        FlattenGraph::Vertex *hit_vertex = nullptr;
        if (!edge->inserted_as_detached) {
            for (FlattenGraph::Vertex *vertex : {edge->start, edge->end}) {
                if (are_near(witness, vertex->point(), _graph.getPrecision())) {
                    hit_vertex = vertex;
                    break;
                }
            }
        }
        LOG_DEBUG("Sweepline hit " << (edge->inserted_as_detached ? "a free loop" : "an ordinary")
                  << " edge with weight = " << initial_edge->label.getWeight()
                  << ", witness point = " << witness << (hit_vertex ? ", at a vertex." : "."));

        /** @brief Whether the initial edge has the ambient region on its left.
         *
         * Imagine an infinite ray emanating from the witness point towards the right.
         * The ambient region is the last region that the intersection point of the sweepline
         * and the ray transited through before hitting the witness point.
         */
        bool ambient_is_on_left;
        if (hit_vertex) {
            // Find the incidence to the "outermost" edge at this vertex that
            // isn't detached or already processed.
            auto const outermost = _combFromRightwardsRay(hit_vertex);
            if (!outermost) {
                LOG_DEBUG("\tAbandoning vertex after failing to find an outgoing edge.");
                return;
            }
            initial_edge = &_graph.getEdge((*outermost)->index);
            ambient_is_on_left = ((*outermost)->sign == FlattenGraph::Incidence::END);
            LOG_DEBUG("\tResolving vertex to walk along the edge with departure azimuth = "
                      << (*outermost)->azimuth << ", weight = " << initial_edge->label.getWeight()
                      << " and d:\n\t\"" << write_svg_path(initial_edge->path) << '\"');

            // Since we're at a vertex, the ambient winding number calculation must exclude
            // all paths in the original path-vector that pass through this vertex. We will
            // therefore combine the indices of these paths with the source path indices
            // stored in the label of the initial edge.
            for (auto const &incidence : hit_vertex->getIncidences()) {
                initial_edge->label.unionIndices(_graph.getEdge(incidence.index).label);
            }
        } else { // We've hit a point somewhere in the middle of the path.
            ambient_is_on_left = _ambientIsOnLeft(*initial_edge);
        }

        // Compute the ambient winding number of the initial edge.
        auto const ambient_winding = partial_winding_number(witness,
                                                            _original_pathvector,
                                                            initial_edge->label.path_indices);
        LOG_DEBUG("Attempting a walk with ambient winding " << ambient_winding << " set to the "
                  << (ambient_is_on_left ? "left" : "right") << " of the edge \""
                  << write_svg_path(initial_edge->path) << "\".");
        initial_edge->label.setWindingRL(!ambient_is_on_left, ambient_winding);

        // Mark the edge as the start of the active walk and commence the walk.
        initial_edge->label.active_start = true;
        _walkFromEdge(*initial_edge);
        initial_edge->label.active_start = false;
    }

    /** @brief Find the first incidence to the vertex encountered when orbiting
     *         counter-clockwise around it, starting at a hypothetical infinite
     *         ray emanating horizontally to the right from the vertex.
     *
     * Viewing the initial unit tangent vectors of departing edges as points on
     * the unit circle, we search for the first such point encountered when
     * following the circle counter-clockwise from (1, 0), with the convention
     * that the y-axis points up. The search skips detached and already processed
     * edges.
     *
     * @param vertex A vertex in a FlattenGraph.
     * @return A valid const-iterator to the first incidence with a greater azimuthal
     *         sort index than the unit vector (1, 0). If no such incidence is found,
     *         an empty optional is returned.
     */
    std::optional<FlattenGraph::IncConstIt> _combFromRightwardsRay(FlattenGraph::Vertex *vertex)
    {
        auto const &incidences = vertex->getIncidences();
        auto const end_iterator = incidences.end();
        // Find the first incidence with a non-negative departure azimuth.
        auto first_nonegative = std::find_if(incidences.begin(), end_iterator,
            [&](FlattenGraph::Incidence const &inc) -> bool
            {
                auto const &e = _graph.getEdge(inc.index);
                if (e.detached || e.label.processed) {
                    return false;
                }
                if (inc.azimuth >= 0.0) {
                    return true;
                }
                return false;
            });
        if (first_nonegative != end_iterator) {
            return first_nonegative;
        }
        // It is possible that all valid edges have negative departure azimuths;
        // in this case, wrap around to -pi.
        auto wrapped = std::find_if(incidences.begin(), end_iterator,
            [&](FlattenGraph::Incidence const &inc) -> bool
            {
                auto const &e = _graph.getEdge(inc.index);
                return !(e.detached || e.label.processed);
            });
        if (wrapped == end_iterator) {
            return std::nullopt;
        }
        return wrapped;
    }

    /** @brief Determines whether the edge colliding with a sweepline
     *         has the ambient region on its left, when looking in the
     *         direction of the edge's orientation, and imagining the
     *         y-axis as pointing up.
     *
     * @param edge An edge which the sweepline has hit from the right.
     * @return True when the ambient region is to the left of the edge,
     *         false when the ambient region is to the right.
     */
    bool _ambientIsOnLeft(FlattenGraph::Edge const &edge)
    {
        auto const &ydir = edge.label.extrema_x.glance_direction_at_max;
        if (ydir < 0.0) {
            LOG_DEBUG("\tPath moving downwards through the witness point.");
            return true;
        } else if (ydir > 0.0) {
            LOG_DEBUG("\tPath moving upwards through the witness point.");
            return false;
        }
        // Corner case: we've hit some kind of an ultra-sharp cusp (parametric self-tangency).
        if (edge.start == edge.end) {
            // The edge is a loop. A regularized PlanarGraph guarantees that every
            // loop is oriented counterclockwise, so the ambient region is on the right.
            LOG_DEBUG("\tThe edge is a pre-oriented loop.");
            return false;
        } else {
            // The ultra-sharp cusp cannot be resolved at order 1.
            // We resort to an expensive signed area calculation.
            auto path_copy = edge.path;
            path_copy.close();
            bool const reversed = _graph.closedPathArea(path_copy) < 0.0;
            LOG_DEBUG("\tDecided that this is " << (reversed ? "a downward" : "an upward")
                      << " glance via an area check.");
            return reversed;
        }
    }

    /** @brief Walk a component of the shape's boundary in a FlattenGraph,
     *         starting from the given initial edge.
     *
     * The algorithm traverses the FlattenGraph, building a representation of
     * the corresponding boundary component of the shape.
     *
     * @param edge The edge where the component traversal will start.
     */
    void _walkFromEdge(FlattenGraph::Edge const &edge)
    {
        edge.label.processed = true; // Set flag first to enable early returns.

        // Check if the edge bounds the shape and whether the shape is on its left.
        auto const sense = edge.label.isBoundary(_fill_rule);
        if (!sense) {
            LOG_DEBUG("\tAborting walk attempt: not predicate boundary.");
            return;
        }

        // OK, so we're going to include this edge in the contour.
        bool const backwards = sense < 0; ///< Whether we're walking backwards on the edge.
        if (edge.inserted_as_detached) {
            // Append an intersection-free path (single loop edge).
            _out.push_back(backwards ? edge.path.reversed() : edge.path);
            return;
        }

        // No more dragging our feet, we have to walk the graph for real!
        // We will assemble a multi-edge boundary component as we go.
        Path component{(backwards ? edge.end : edge.start)->point()};
        if (_assembleContour(component, edge, backwards)) {
            component.close();
            _out.push_back(std::move(component));
        }
    }

    /// The state of an ongoing walk through a FlattenGraph.
    struct WalkState
    {
        size_t index; ///< The index of the edge we have just walked on.
        FlattenGraph::Edge const *edge; ///< The edge just walked.
        FlattenGraph::Vertex *vertex; ///< The vertex that the edge takes us to.
        FlattenGraph::Incidence *arrival; ///< The arrival incidence at the vertex.
        int winding_on_left; ///< The winding of the region to the left of the walk direction.
        bool backwards; ///< Whether the walk direction agrees with the orientation of the edge.
        bool keep_walking; ///< Whether the walk should continue.
        bool success; ///< Whether the walk ended successfully.
    };

    /** @brief Assemble a multi-edge boundary component of the shape starting at the given edge.
     *
     * @param[out] output The path to which the contour pieces will be appended.
     * @param start_edge The initial edge of the contour.
     * @param backwards Whether the initial edge should be traversed backwards relative to its
     *                  time coordinate. Hence, this argument decides the walk direction.
     * @return True when the boundary component has been correctly assembled, false on failure.
     */
    bool _assembleContour(Path &output, FlattenGraph::Edge const &start_edge, bool backwards)
    {
        // Set up the walk state.
        auto state = WalkState{ .keep_walking = true, .success = false };
        output.append(_makeTurn(state, start_edge, backwards));
        assert(start_edge.label.has_windings);
        state.winding_on_left = start_edge.label.windingOnRL(backwards);

        // Walk until the completion of the circuit or an error.
        while (state.keep_walking) {
            _walkStep(output, state);
        }
        return state.success;
    }

    /** @brief Modify a walk state object when we take a new edge.
     *
     * @param[out] state The state of the walk to be modified.
     * @param edge The edge onto which we wish to turn.
     * @param backwards Whether we will be walking the edge backwards in relation
     *                  to its path's orientation.
     * @return The path of the taken edge, reversed if needed, so as to agree with the walk direction.
     */
    Path _makeTurn(WalkState &state, FlattenGraph::Edge const &edge, bool backwards)
    {
        LOG_DEBUG("Taking edge \"" << write_svg_path(edge.path) << "\"\n       in "
                  << (backwards ? "the reverse" : "its natural") << " direction.");
        state.backwards = backwards;
        state.edge = &edge;
        state.index = _graph.getEdgeIndex(edge);
        auto [next_vertex, next_arrival] = _graph.getIncidence(state.index,
            state.backwards ? FlattenGraph::Incidence::START : FlattenGraph::Incidence::END);
        state.keep_walking = next_vertex && next_arrival;
        state.vertex = next_vertex;
        state.arrival = next_arrival;
        return backwards ? edge.path.reversed() : edge.path;
    }

    /** @brief Step one edge forward in a walk through a FlattenGraph.
     *
     * @param[out] output The path to append to.
     * @param[in,out] state The state of the walk.
     */
    void _walkStep(Path &output, WalkState &state)
    {
        // We're at a vertex and need to decide where to turn next.
        // We will iterate through the incidences at the vertex, in the clockwise order,
        // in order to find the leftmost admissible turn from where we've arrived.
        LOG_DEBUG("\tWalk step: at vertex " << state.vertex->point() <<", backtrack azimuth "
                  << state.arrival->azimuth << ".");
        bool made_turn = false;

        /// When a potential turn is skipped, update the winding number to the left of the walk path.
        auto const skip_turn = [&](FlattenGraph::IncConstIt const &turn, FlattenGraph::Edge const &edge) {
            if (turn->sign == FlattenGraph::Incidence::START) {
                state.winding_on_left -= edge.label.getWeight();
            } else {
                state.winding_on_left += edge.label.getWeight();
            }
        };

        // Find the leftmost admissible turn.
        for (auto turn = _graph.nextIncidenceIt(*state.vertex, *state.arrival, true);
               &(*turn) != state.arrival;
                  turn = _graph.nextIncidenceIt(*state.vertex, turn, true))
        {
            auto const &candidate_edge = _graph.getEdge(turn->index);
            if (candidate_edge.label.processed) {
                // This edge was already processed (and classified).
                if (candidate_edge.label.active_start) {
                    LOG_DEBUG("Detected walk start, closing circuit.\n");
                    output.setFinal(state.vertex->point());
                    state.keep_walking = false;
                    state.success = true;
                    return;
                }
                // The edge was seen before (maybe picked in a different component tangent to
                // ours or rejected due to not being a predicate boundary) and it is not the
                // start of the current walk, so we must not connect to it. We skip it and
                // continue our search.
                skip_turn(turn, candidate_edge);
                continue;
            }

            // Compute winding numbers if unknown.
            if (!candidate_edge.label.has_windings) {
                candidate_edge.label.setWindingRL(turn->sign == FlattenGraph::Incidence::END,
                                                  state.winding_on_left);
            }
            candidate_edge.label.processed = true;

            // Check if the candidate edge is a predicate boundary.
            if (auto sense = candidate_edge.label.isBoundary(_fill_rule)) {
                // This edge looks good, so we'll turn onto it!
                made_turn = true;
                // Append the path to output and break out of the turn search.
                auto next_path = _makeTurn(state, candidate_edge, sense < 0);
                next_path.setInitial(output.finalPoint());
                output.append(std::move(next_path));
                break;
            }
            LOG_DEBUG("\tSkipping edge: not a predicate boundary;\n\td = \"" << write_svg_path(candidate_edge.path) << '\"');
            skip_turn(turn, candidate_edge);
        }
        state.keep_walking = made_turn;
    }
};

} // end anonymous namespace

/* Convert a shape with an arbitrary fill rule to the internal representation,
 * i.e., an uncrossed PathVector with the FillRule::CONVENTIONAL fill rule.
 */
PathVector Shape::flatten(PathVector const &pv, FillRule fill_rule, bool close_open_paths, Coord precision)
{
    if (fill_rule == FillRule::CONVENTIONAL) {
        // If the caller claims that the path-vector description of
        // the boundary is already conventional, we trust them and do nothing.
        return pv;
    }

    // Prepare the piece generator and the flattening graph.
    auto cleaned = clean_paths(pv, close_open_paths);
    auto chopper = PathVectorCutter(cleaned, precision);
    auto graph = FlattenGraph(precision);
    graph.reserveEdgeCapacity(chopper.getLengthEstimate());

    // Keep track of which paths contain intersections with other paths.
    auto has_intersections = std::vector<bool>(cleaned.size(), false);

    // Chop the path-vector at self-intersections and feed the pieces into the graph.
    for (auto &&[index, path] : chopper) {
        if (path.empty() || (path.size() == 1 && path.at(0).isDegenerate())) {
            continue;
        }
        graph.insertEdge(path.withoutDegenerateCurves(), FlattenLabel(index));
        has_intersections[index] = true;
    }

    // Feed the pieces without intersections as "detached edges" (free loops).
    for (size_t i = 0; i < cleaned.size(); i++) {
        if (!has_intersections[i]) {
            graph.insertDetached(std::move(cleaned[i]), FlattenLabel(i));
        }
    }

    // Regularize the graph and compute the x-bounds of all edges.
    graph.regularize();
    for (auto &edge : graph.getEdges()) {
        edge.label.extrema_x = edge.path.extrema(X);
    }

    // Assemble the regularized boundary.
    auto flattener = FlattenSweep(graph, pv, fill_rule);
    Sweeper(flattener).process();
    return flattener.getResult();
}

} // namespace Geom

/*
  Local Variables:
  mode:c++
  c-file-style:"stroustrup"
  c-file-offsets:((innamespace . 0)(inline-open . 0)(case-label . +))
  indent-tabs-mode:nil
  fill-column:99
  End:
*/
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
