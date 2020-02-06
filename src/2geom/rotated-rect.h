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
#ifndef LIB2GEOM_SEEN_ROTATED_RECT_H
#define LIB2GEOM_SEEN_ROTATED_RECT_H

#include <2geom/rect.h>

namespace Geom{

/**
 * @brief Axis non-aligned, non-empty rectangle.
 * @ingroup Primitives
 *
 */
class RotatedRect {
public:
    /**
     * @brief Returns i'th rectangle corner
     *
     * RotatedRect::corner(i) corresponds to i'th corner of \a rect it was constructed from
     * */
    const Point &corner(unsigned i) const;

    Point midpoint() const;

    Rect bounds() const;

    ///Returns diagonal length
    Coord diameter() const;

    /// Returns shorter rectangle side length
    Coord minExtent() const;

    ///Returns longer rectangle side length
    Coord maxExtent() const;

    bool contains(Point const &p) const;

    bool contains(Rect const &rect) const;

    bool intersects(Rect const &rect) const;

    static RotatedRect from_rect_rotate(Rect const &rect, const Rotate &rotate);
    static RotatedRect from_rect_rotate(Rect const &rect, const Rotate &rotate, const Point & rotation_center);

private:
    RotatedRect(Point const &cw0, Point const &cw1, Point const &cw2) : m_corners{cw0, cw1, cw2, cw0 + (cw2 - cw1)} {}

    const std::array<Point, 4> m_corners;
};

}
#endif //LIB2GEOM_SEEN_ROTATED_RECT_H
