/** @file A class representing a two-dimensional shape.
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

#include <2geom/exception.h>
#include <2geom/pathvector.h>
#include <2geom/path-sink.h>
#include <2geom/shape.h>

#include "planar-graph.h"

namespace Geom {

/** @brief Calculate an appropriate numerical precision for the given path-vector.
 *
 * @param pv A path-vector for the handling of which we want to determine a precision.
 * @return A small positive number depending on the maximum absolute value of coordinates
 *         occurring in the path-vector.
 */
Coord Shape::_calculateEpsilon(PathVector const &pv)
{
    if (auto const bounds = pv.boundsFast()) {
        auto const &x_interval = (*bounds)[X];
        auto const &y_interval = (*bounds)[Y];
        Coord const absmax_x = std::max(std::abs(x_interval.min()), std::abs(x_interval.max()));
        Coord const absmax_y = std::max(std::abs(y_interval.min()), std::abs(y_interval.max()));
        return PRECISION_COEFFICIENT * std::max(absmax_x, absmax_y);
    }
    return EPSILON;
}

/** Write the contour of a rectangle as a closed path
 *  consisting of four line segments to the output path-vector.
 */
static void write_rect_contour(PathVector &output, Rect const &rect)
{
    PathBuilder builder(output);
    unsigned vertex = 0;
    builder.moveTo(rect.corner(vertex));
    while (++vertex != 4) {
        builder.lineTo(rect.corner(vertex));
    }
    builder.closePath();
}

Shape::Shape(Rect const &rectangle)
{
    if (!rectangle.isFinite()) {
        THROW_UNBOUNDEDREGION("Cannot construct a Shape object from an unbounded rectangle.");
    }
    if (rectangle.hasZeroArea()) {
        return; // Leave the shape in an empty state.
    }

    write_rect_contour(_boundary, rectangle);
    _epsilon = _calculateEpsilon(_boundary);
}

Shape::Shape(Parallelogram const &parallelogram)
{
    auto const from_unit_square = parallelogram.unitSquareTransform();
    if (from_unit_square.det() == 0.0) {
        return; // Leave the shape in an empty state.
    }

    // Start with a unit square and transform ourselves to the desired parallelogram.
    static auto const unit_interval = Interval(0, 1);
    write_rect_contour(_boundary, Rect(unit_interval, unit_interval));
    *this *= from_unit_square;
}

/** Write the contour of a circle as a path consisting of two circular arcs
 *  to the output path-vector.
 */
static void write_circle_contour(PathVector &output, Circle const &circle)
{
    auto const r = circle.radius();
    auto const horiz = Point(r, 0);
    auto const right_point = circle.center() + horiz;
    auto const left_point = circle.center() - horiz;

    PathBuilder builder(output);
    builder.moveTo(right_point);
    builder.arcTo(r, r, 0.0, true, true, left_point);
    builder.arcTo(r, r, 0.0, true, true, right_point);
    builder.closePath();
}

Shape::Shape(Circle const &circle) noexcept
{
    if (circle.isDegenerate()) {
        return; // Leave the shape in an empty state.
    }
    write_circle_contour(_boundary, circle);
    _epsilon = _calculateEpsilon(_boundary);
}

Shape::Shape(Ellipse const &ellipse) noexcept
{
    auto const from_unit_circle = ellipse.unitCircleTransform();
    if (from_unit_circle.det() == 0.0) {
        return; // Leave the shape in an empty state.
    }

    // Start with a unit circle contour and transform it to the desired ellipse.
    write_circle_contour(_boundary, Circle(Point(0, 0), 1));
    *this *= from_unit_circle;
}

Shape &Shape::operator*=(Affine const &tr)
{
    auto const determinant = tr.det();
    if (determinant == 0.0) {
        // We're getting squashed to zero area.
        _boundary.clear();
    } else {
        _boundary *= tr;
        if (determinant < 0.0) { // Orientation has been flipped; we must reverse all paths.
            _boundary.reverse(false);
        }
    }
    _epsilon = _calculateEpsilon(_boundary);
    return *this;
}

Coord Shape::area() const
{
    Coord result = 0.0;
    Point _;
    for (auto const &path : _boundary) {
        double gauss_green;
        centroid(path.toPwSb(), _, gauss_green);
        result -= gauss_green; // The sign convention in PwSb code differs from ours.
    }
    return result;
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
