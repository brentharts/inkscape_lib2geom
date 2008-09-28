/*
 * Python bindings for lib2geom
 *
 * Copyright 2007 Aaron Spike <aaron@ekips.org>
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
 *
 */

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "../curve.h"
#include "../path.h"
#include "../pathvector.h"
#include "../sbasis-to-bezier.h"
#include "helpers.h"

#include "../point.h"
#include "../rect.h"
#include "../d2.h"

using namespace boost::python;

struct CurveWrap : Geom::Curve, wrapper<Geom::Curve>
{
    Geom::Point initialPoint() const {return this->get_override("initialPoint")();}
    Geom::Point finalPoint() const {return this->get_override("finalPoint")();}
    bool isDegenerate() const {return this->get_override("isDegenerate")();}
    CurveWrap *duplicate() const {return this->get_override("duplicate")();}
    Geom::Rect boundsFast() const {return this->get_override("boundsFast")();}
    Geom::Rect boundsExact() const {return this->get_override("boundsExact")();}
    Geom::Rect boundsLocal(Geom::Interval i, unsigned deg) const {return this->get_override("boundsLocal")(i,deg);}
    std::vector<double> roots(double v, Geom::Dim2 d) const {return this->get_override("roots")(v,d);}

    int winding(Geom::Point p) const {
        if (override f = this->get_override("winding")) {
            return f(p);
        }
        return Geom::Curve::winding(p);
    }
    int default_winding(Geom::Point p) const { return this->Geom::Curve::winding(p); }

    Geom::Curve *portion(double f, double t) const { return this->get_override("portion")(f,t); }
    Geom::Curve *reverse() const { 
        if (override f = this->get_override("reverse")) {
            return f();
        }
        return Geom::Curve::reverse(); 
    }
    Geom::Curve *default_reverse() const { return this->Geom::Curve::reverse(); }

    Geom::Curve *derivative() const { return this->get_override("derivative")(); }

    void setInitial(Geom::Point v){ this->get_override("setInitial")(v); }
    void setFinal(Geom::Point v){ this->get_override("setFinal")(v); }


    Geom::Curve *transformed(Geom::Matrix const &m) const { return this->get_override("transformed")(m); }

    Geom::Point pointAt(Geom::Coord t) const {
        if (override f = this->get_override("pointAt")) {
            return f(t);
        }
        return Geom::Curve::pointAt(t);
    }
    Geom::Point default_pointAt(Geom::Coord t) { return this->Geom::Curve::pointAt(t); }
    std::vector<Geom::Point> pointAndDerivatives(Geom::Coord t, unsigned n) const {return this->get_override("pointAndDerivatives")();}

    Geom::D2<Geom::SBasis> toSBasis() const {return this->get_override("sbasis")();}
};

void wrap_path()
{
    class_<CurveWrap, boost::noncopyable>("Curve")
        .def("initalPoint", pure_virtual(&Geom::Curve::initialPoint))
        .def("finalPoint", pure_virtual(&Geom::Curve::finalPoint))
        .def("duplicate", pure_virtual(&Geom::Curve::duplicate), return_value_policy<manage_new_object>())
        .def("boundsFast", pure_virtual(&Geom::Curve::boundsFast))
        .def("boundsExact", pure_virtual(&Geom::Curve::boundsExact))
        //.def("pointAt", &Geom::Curve::pointAt, &CurveWrap::default_pointAt)
        //.def("winding", &Geom::Curve::winding, &CurveWrap::default_winding)
        .def("pointAndDerivatives", pure_virtual(&Geom::Curve::pointAndDerivatives))
        .def("toSBasis", pure_virtual(&Geom::Curve::toSBasis))
    ;
    class_<Geom::Path>("Path")
        .def("empty", &Geom::Path::empty)
        .def("closed", &Geom::Path::closed)
        .def("close", &Geom::Path::close)
        .def("boundsFast", &Geom::Path::boundsFast)
        .def("boundsExact", &Geom::Path::boundsExact)
        .def("toPwSb", &Geom::Path::toPwSb)
        .def(self * Geom::Matrix())
        .def("pointAt", &Geom::Path::pointAt)
        .def("valueAt", &Geom::Path::valueAt)
        .def("__call__", &Geom::Path::pointAt)
        .def("roots", &Geom::Path::roots)
        //.def("allNearestPoints", &Geom::Path::allNearestPoints)
        //.def("nearestPoint", &Geom::Path::nearestPoint)
        .def("appendPortionTo", &Geom::Path::appendPortionTo)
        //.def("portion", &Geom::Path::portion)
        .def("reverse", &Geom::Path::reverse)
        //.def("insert", &Geom::Path::insert)
        .def("clear", &Geom::Path::clear)
        //.def("erase", &Geom::Path::erase)
        .def("erase_last", &Geom::Path::erase_last)
        //.def("replace", &Geom::Path::replace)
        .def("start", &Geom::Path::start)
        .def("initialPoint", &Geom::Path::initialPoint)
        .def("finalPoint", &Geom::Path::finalPoint)
        .def("setInitial", &Geom::Path::setInitial)
        .def("setFinal", &Geom::Path::setFinal)
        //.def("append", &Geom::Path::append)
        //.def("appendNew", &Geom::Path::appendNew)
    ;
    def("paths_to_pw",Geom::paths_to_pw);
    class_<std::vector<Geom::Path> >("PathVector")
        .def(vector_indexing_suite<std::vector<Geom::Path> >())
        .def(self * Geom::Matrix())
    ;
    def("reverse_paths_and_order", Geom::reverse_paths_and_order);
    def("bounds_fast", (Geom::Rect (*)(Geom::PathVector const &))&Geom::bounds_fast);
    def("bounds_exact", (Geom::Rect (*)(Geom::PathVector const &))&Geom::bounds_exact);
    def("path_from_piecewise", Geom::path_from_piecewise);
    def("path_from_sbasis", Geom::path_from_sbasis);
    def("cubicbezierpath_from_sbasis", Geom::cubicbezierpath_from_sbasis);
}

/*
  Local Variables:
  mode:c++
  c-file-style:"stroustrup"
  c-file-offsets:((innamespace . 0)(substatement-open . 0))
  indent-tabs-mode:nil
  c-brace-offset:0
  fill-column:99
  End:
  vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4 :
*/
