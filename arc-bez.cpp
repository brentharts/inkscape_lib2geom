#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>

#include <gtk/gtk.h>
#include <cassert>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <vector>
#include "s-basis.h"
#include "interactive-bits.h"
#include "multidim-sbasis.h"
#include "bezier-to-sbasis.h"
#include "sbasis-to-bezier.h"
#include "path-cairo.h"

#include "toy-framework.cpp"

using std::string;
using std::vector;

void draw_cb(cairo_t *cr, multidim_sbasis<2> const &B) {
    Geom::PathBuilder pb;
    subpath_from_sbasis(pb, B, 0.1);
    cairo_path(cr, pb.peek());
}

SBasis curvature(multidim_sbasis<2> & B) {
    multidim_sbasis<2> dB = derivative(B);
    multidim_sbasis<2> ddB = derivative(dB);
    SBasis n = dB[0]*ddB[1] -dB[1]*ddB[0];
    SBasis den = dB[0]*dB[0] + dB[1]*dB[1];
    den = den*den;
    return divide(n*sqrt(den, 4), den, 6);
}

class ArcBez: public Toy {
    virtual void draw(cairo_t *cr, std::ostringstream *notify, int width, int height, bool save) {
        multidim_sbasis<2> B = bezier_to_sbasis<2, 3>(handles.begin());
        draw_cb(cr, B);
        cairo_stroke(cr);
        
        cairo_set_source_rgba (cr, 0.5, 0.5, 0, 0.8);
        multidim_sbasis<2> dB = derivative(B);
        if(0) for(int dim = 0; dim < 2; dim++) {
            multidim_sbasis<2> plot;
            plot[0] = SBasis(width*BezOrd(0.25,0.75));
            plot[1] = BezOrd(height*3/4) - 0.5*dB[dim];
        
            draw_cb(cr, plot);
            cairo_stroke(cr);
        }
        cairo_set_source_rgba (cr, 0.5, 0, 0.5, 0.8);
        if(0){
        multidim_sbasis<2> plot;
        plot[0] = SBasis(width*BezOrd(0.25,0.75));
        plot[1] = derivative((1./height)*(dB[0]*dB[0])
                + (1./height)*(dB[1]*dB[1]));
        std::vector<double> r = roots(plot[1]);
        plot[1] = BezOrd(height*3/4) - plot[1];
        draw_cb(cr, plot);
        cairo_stroke(cr);
        for(int i = 0; i < r.size(); i++) {
                //draw_cross(cr, point_at(B, r[i]));
                //draw_cross(cr, point_at(plot, r[i]));
        }
        }
        cairo_set_source_rgba (cr, 0.25, 0.5, 0, 0.8);
        {
        multidim_sbasis<2> plot;
        plot[0] = SBasis(width*BezOrd(0.25,0.75));
        plot[1] = height*derivative(curvature(B));
        std::vector<double> r = roots(plot[1]);
        plot[1] = BezOrd(height*3/4) - plot[1];
        draw_cb(cr, plot);
        cairo_stroke(cr);
        for(int i = 0; i < r.size(); i++) {
                draw_cross(cr, point_at(B, r[i]));
                draw_cross(cr, point_at(plot, r[i]));
        }
        }
        
        cairo_set_source_rgba (cr, 0., 0.5, 0, 0.8);
        double prev_seg = 0;
        int N = 10;
        for(int subdivi = 0; subdivi < N; subdivi++) {
        double dsubu = 1./N;
        double subu = dsubu*subdivi;
        BezOrd dt(subu, dsubu + subu);
        multidim_sbasis<2> dBp = compose(dB, dt);
        SBasis arc = L2(dBp, 2);
        arc = (1./N)*integral(arc);
        arc = arc - BezOrd(Hat(arc.point_at(0) - prev_seg));
        prev_seg = arc.point_at(1);
        
        multidim_sbasis<2> plot;
        plot[0] = SBasis(width*dt);
        plot[1] = BezOrd(height) - arc;
        
        draw_cb(cr, plot);
        cairo_stroke(cr);
        }
        *notify << "arc length = " << prev_seg << std::endl;
    }
};

int main(int argc, char **argv) {
    for(int i = 0; i < 4; i++)
        handles.push_back(Geom::Point(uniform()*400, uniform()*400));
    
    init(argc, argv, "arc-bez", new ArcBez());

    return 0;
}

/*
  Local Variables:
  mode:c++
  c-file-style:"stroustrup"
  c-file-offsets:((innamespace . 0)(inline-open . 0)(case-label . +))
  indent-tabs-mode:nil
  fill-column:99
  End:
*/
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:encoding=utf-8:textwidth=99 :
