/**
 * \file
 * \brief Axis-non-aligned rectangle
 *//*
 * Authors:
 *   Sergei Izmailov <sergei.a.izmailov@gmail.com>
 * Copyright 2020 Authors
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
 * in the file COPYING-LGPL-2.1; if not, output to the Free Software
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

#include <2geom/rotated-rect.h>

namespace {

/// Returns: true if line segments AB and CD intersects
/// Precondition: AB and CD must be not collinear
bool non_collinear_segments_intersect(const Geom::Point &A, const Geom::Point &B,
                                      const Geom::Point &C, const Geom::Point &D) {
    double signC = cross(B - A, C - A);
    double signD = cross(B - A, D - A);

    double signA = cross(D - C, A - C);
    double signB = cross(D - C, B - C);

    return signA * signB < 0 && signC * signD < 0;
}
}

namespace Geom {

const Point &RotatedRect::corner(unsigned int i) const {
    return m_corners[i % 4];
}

Coord RotatedRect::diameter() const {
    return distance(corner(0), corner(2));
}

Coord RotatedRect::minExtent() const {
    return std::min(distance(corner(0), corner(1)), distance(corner(1), corner(2)));
}

Coord RotatedRect::maxExtent() const {
    return std::max(distance(corner(0), corner(1)), distance(corner(1), corner(2)));
}

Point RotatedRect::midpoint() const {
    return Geom::middle_point(corner(0), corner(2));
}

Rect RotatedRect::bounds() const {
    Rect result;
    for (int i=0;i<4;i++){
        result.expandTo(corner(i));
    }
    return result;
}

bool RotatedRect::contains(const Point &p) const {
    return  cross(corner(1) - corner(0), p - corner(0)) >= 0 &&
            cross(corner(2) - corner(1), p - corner(1)) >= 0 &&
            cross(corner(3) - corner(2), p - corner(2)) >= 0 &&
            cross(corner(0) - corner(3), p - corner(3)) >= 0;
}

bool RotatedRect::contains(const Geom::Rect &rect) const {
    return  contains(rect.corner(0)) &&
            contains(rect.corner(1)) &&
            contains(rect.corner(2)) &&
            contains(rect.corner(3));
}

bool RotatedRect::intersects(const Geom::Rect &rect) const {
    // Check if a corner is included in other rect
    // otherwise if two rectangles intersect then two or more diagonals intersects
    return rect.contains(corner(0)) ||
           rect.contains(corner(1)) ||
           rect.contains(corner(2)) ||
           rect.contains(corner(3)) ||
           contains(rect.corner(0)) ||
           contains(rect.corner(1)) ||
           contains(rect.corner(2)) ||
           contains(rect.corner(3)) ||
           non_collinear_segments_intersect(rect.corner(0), rect.corner(2), corner(0), corner(2)) ||
           non_collinear_segments_intersect(rect.corner(1), rect.corner(3), corner(1), corner(3));

}


RotatedRect RotatedRect::from_rect_rotate(const Rect &rect, const Rotate &rotate) {
    return RotatedRect::from_rect_rotate(rect, rotate, rect.midpoint());
}

RotatedRect RotatedRect::from_rect_rotate(const Rect &rect, const Rotate &rotate, const Point& rotation_center) {
    auto c0 = rect.corner(0) * rotate;
    auto c1 = rect.corner(1) * rotate;
    auto c2 = rect.corner(2) * rotate;

    auto shift = rotation_center - rotation_center * rotate;

    c0 += shift;
    c1 += shift;
    c2 += shift;

    return Geom::RotatedRect(c0, c1, c2);

}


}