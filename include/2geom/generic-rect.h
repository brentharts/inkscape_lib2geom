/**
 * \file
 * \brief Axis-aligned rectangle
 *//*
 * Authors:
 *   Michael Sloan <mgsloan@gmail.com>
 *   Krzysztof Kosiński <tweenk.pl@gmail.com>
 * Copyright 2007-2011 Authors
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
 * Authors of original rect class:
 *   Lauris Kaplinski <lauris@kaplinski.com>
 *   Nathan Hurst <njh@mail.csse.monash.edu.au>
 *   bulia byak <buliabyak@users.sf.net>
 *   MenTaLguY <mental@rydia.net>
 */

#ifndef LIB2GEOM_SEEN_GENERIC_RECT_H
#define LIB2GEOM_SEEN_GENERIC_RECT_H

#include <limits>
#include <iostream>
#include <optional>
#include <2geom/coord.h>

namespace Geom {

template <typename C>
class GenericOptRect;

/**
 * @brief Axis aligned, non-empty, generic rectangle.
 * @ingroup Primitives
 */
template <typename C>
class GenericRect
    : CoordTraits<C>::RectOps
{
    using CInterval = typename CoordTraits<C>::IntervalType;
    using CPoint = typename CoordTraits<C>::PointType;
    using CRect = typename CoordTraits<C>::RectType;
    using OptCRect = typename CoordTraits<C>::OptRectType;
protected:
    CInterval f[2];
public:
    using D1Value = CInterval;
    using D1Reference = CInterval &;
    using D1ConstReference = CInterval const &;

    /// @name Create rectangles.
    /// @{
    /** @brief Create a rectangle that contains only the point at (0, 0). */
    GenericRect() = default;
    /** @brief Create a rectangle from X and Y intervals. */
    GenericRect(CInterval const &a, CInterval const &b) : f{a, b} {}
    /** @brief Create a rectangle from two points. */
    GenericRect(CPoint const &a, CPoint const &b) : f{{a[X], b[X]}, {a[Y], b[Y]}} {}
    /** @brief Create rectangle from coordinates of two points. */
    GenericRect(C x0, C y0, C x1, C y1) : f{{x0, x1}, {y0, y1}} {}
    /** @brief Create a rectangle from a range of points.
     * The resulting rectangle will contain all points from the range.
     * The return type of iterators must be convertible to Point.
     * The range must not be empty. For possibly empty ranges, see OptRect.
     * @param start Beginning of the range
     * @param end   End of the range
     * @return Rectangle that contains all points from [start, end). */
    template <typename InputIterator>
    static CRect from_range(InputIterator start, InputIterator end) {
        assert(start != end);
        CPoint p1 = *start;
        auto result = CRect(p1, p1);
        for (++start; start != end; ++start) {
            result.expandTo(*start);
        }
        return result;
    }
    /** @brief Create a rectangle from a C-style array of points it should contain. */
    static CRect from_array(CPoint const *c, unsigned n) {
        return GenericRect<C>::from_range(c, c + n);
    }
    /** @brief Create rectangle from origin and dimensions. */
    static CRect from_xywh(C x, C y, C w, C h) {
        return GenericRect<C>::from_xywh(CPoint(x, y), CPoint(w, h));
    }
    /** @brief Create rectangle from origin and dimensions. */
    static CRect from_xywh(CPoint const &xy, CPoint const &wh) {
        return CRect(xy, xy + wh);
    }
    /// Create infinite rectangle.
    static CRect infinite() {
        auto min = std::numeric_limits<C>::min();
        auto max = std::numeric_limits<C>::max();
        return CRect(min, min, max, max);
    }
    /// @}

    /// @name Inspect dimensions.
    /// @{
    CInterval &operator[](unsigned i)             { return f[i]; }
    CInterval const &operator[](unsigned i) const { return f[i]; }
    CInterval &operator[](Dim2 d)                 { return f[d]; }
    CInterval const &operator[](Dim2 d) const     { return f[d]; }

    /** @brief Get the corner of the rectangle with smallest coordinate values.
     * In 2Geom standard coordinate system, this means upper left. */
    CPoint min() const { return CPoint(f[X].min(), f[Y].min()); }
    /** @brief Get the corner of the rectangle with largest coordinate values.
     * In 2Geom standard coordinate system, this means lower right. */
    CPoint max() const { return CPoint(f[X].max(), f[Y].max()); }
    /** @brief Return the n-th corner of the rectangle.
     * Returns corners in the direction of growing angles, starting from
     * the one given by min(). For the standard coordinate system used
     * in 2Geom (+Y downwards), this means clockwise starting from
     * the upper left. */
    CPoint corner(unsigned i) const {
        switch (i % 4) {
            case 0:  return CPoint(f[X].min(), f[Y].min());
            case 1:  return CPoint(f[X].max(), f[Y].min());
            case 2:  return CPoint(f[X].max(), f[Y].max());
            default: return CPoint(f[X].min(), f[Y].max());
        }
    }

    /** @brief Return top coordinate of the rectangle (+Y is downwards). */
    C top() const { return f[Y].min(); }
    /** @brief Return bottom coordinate of the rectangle (+Y is downwards). */
    C bottom() const { return f[Y].max(); }
    /** @brief Return leftmost coordinate of the rectangle (+X is to the right). */
    C left() const { return f[X].min(); }
    /** @brief Return rightmost coordinate of the rectangle (+X is to the right). */
    C right() const { return f[X].max(); }

    /** @brief Get the horizontal extent of the rectangle. */
    C width() const { return f[X].extent(); }
    /** @brief Get the vertical extent of the rectangle. */
    C height() const { return f[Y].extent(); }
    /** @brief Get the ratio of width to height of the rectangle. */
    Coord aspectRatio() const { return static_cast<Coord>(width()) / height(); }

    /** @brief Get rectangle's width and height as a point.
     * @return Point with X coordinate corresponding to the width and the Y coordinate
     *         corresponding to the height of the rectangle. */
    CPoint dimensions() const { return CPoint(f[X].extent(), f[Y].extent()); }
    /** @brief Get the point in the geometric center of the rectangle. */
    CPoint midpoint() const { return CPoint(f[X].middle(), f[Y].middle()); }

    /** @brief Compute the rectangle's area. */
    C area() const { return f[X].extent() * f[Y].extent(); }
    /** @brief Check whether the rectangle has zero area. */
    bool hasZeroArea() const { return f[X].isSingular() || f[Y].isSingular(); }

    /** @brief Get the larger extent (width or height) of the rectangle. */
    C maxExtent() const { return std::max(f[X].extent(), f[Y].extent()); }
    /** @brief Get the smaller extent (width or height) of the rectangle. */
    C minExtent() const { return std::min(f[X].extent(), f[Y].extent()); }

    /** @brief Get rectangle's distance SQUARED away from the given point **/
    C distanceSq(CPoint const &p) const {
        return (p - clamp(p)).lengthSq();
    }

    /** @brief Clamp point to the rectangle. */
    CPoint clamp(CPoint const &p) const {
        return CPoint(f[X].clamp(p[X]), f[Y].clamp(p[Y]));
    }
    /** @brief Get the nearest point on the edge of the rectangle. */
    CPoint nearestEdgePoint(CPoint const &p) const {
        if (!contains(p)) {
            return clamp(p);
        }
        CPoint result = p;
        C cx = f[X].nearestEnd(p[X]);
        C cy = f[Y].nearestEnd(p[Y]);
        if (std::abs(cx - p[X]) <= std::abs(cy - p[Y])) {
            result[X] = cx;
        } else {
            result[Y] = cy;
        }
        return result;
    }
    /// @}

    /// @name Test other rectangles and points for inclusion.
    /// @{
    /** @brief Check whether the rectangles have any common points. */
    bool intersects(GenericRect<C> const &r) const { 
        return f[X].intersects(r[X]) && f[Y].intersects(r[Y]);
    }
    /** @brief Check whether the rectangle includes all points in the given rectangle. */
    bool contains(GenericRect<C> const &r) const { 
        return f[X].contains(r[X]) && f[Y].contains(r[Y]);
    }

    /** @brief Check whether the rectangles have any common points.
     * Empty rectangles will not intersect with any other rectangle. */
    inline bool intersects(OptCRect const &r) const;
    /** @brief Check whether the rectangle includes all points in the given rectangle.
     * Empty rectangles will be contained in any non-empty rectangle. */
    inline bool contains(OptCRect const &r) const;

    /** @brief Check whether the given point is within the rectangle. */
    bool contains(CPoint const &p) const {
        return f[X].contains(p[X]) && f[Y].contains(p[Y]);
    }
    /// @}

    /// @name Modify the rectangle.
    /// @{
    /** @brief Set the minimum X coordinate of the rectangle. */
    void setLeft(C val) {
        f[X].setMin(val);
    }
    /** @brief Set the maximum X coordinate of the rectangle. */
    void setRight(C val) {
        f[X].setMax(val);
    }
    /** @brief Set the minimum Y coordinate of the rectangle. */
    void setTop(C val) {
        f[Y].setMin(val);
    }
    /** @brief Set the maximum Y coordinate of the rectangle. */
    void setBottom(C val) {
        f[Y].setMax(val);
    }
    /** @brief Set the upper left point of the rectangle. */
    void setMin(CPoint const &p) {
        f[X].setMin(p[X]);
        f[Y].setMin(p[Y]);
    }
    /** @brief Set the lower right point of the rectangle. */
    void setMax(CPoint const &p) {
        f[X].setMax(p[X]);
        f[Y].setMax(p[Y]);
    }
    /** @brief Enlarge the rectangle to contain the given point. */
    void expandTo(CPoint const &p)        { 
        f[X].expandTo(p[X]);
        f[Y].expandTo(p[Y]);
    }
    /** @brief Enlarge the rectangle to contain the argument. */
    void unionWith(CRect const &b) { 
        f[X].unionWith(b[X]);
        f[Y].unionWith(b[Y]);
    }
    /** @brief Enlarge the rectangle to contain the argument.
     * Unioning with an empty rectangle results in no changes. */
    void unionWith(OptCRect const &b);

    /** @brief Expand the rectangle in both directions by the specified amount.
     * Note that this is different from scaling. Negative values will shrink the
     * rectangle. If <code>-amount</code> is larger than
     * half of the width, the X interval will contain only the X coordinate
     * of the midpoint; same for height. */
    void expandBy(C amount) {
        expandBy(amount, amount);
    }

    /** @brief Shrink the rectangle in both directions by the specified amount.
     * If the amount of shrinkage exceeds a half of the smaller dimension, the
     * rectangle will be shrunk to a line segment.
    */
    void shrinkBy(C amount) { expandBy(-amount); }

    /** @brief Expand the rectangle in both directions.
     * Note that this is different from scaling. Negative values will shrink the
     * rectangle. If <code>-x</code> is larger than
     * half of the width, the X interval will contain only the X coordinate
     * of the midpoint; same for height. */
    void expandBy(C x, C y) { 
        f[X].expandBy(x);
        f[Y].expandBy(y);
    }

    /** @brief Shrink the rectangle in both directions by possibly different amounts.
     * If the amount of shrinkage exceeds a half of the corresponding dimension, the
     * rectangle will be shrunk to a line segment.
    */
    void shrinkBy(C x, C y) { expandBy(-x, -y); }

    /** @brief Expand the rectangle by the coordinates of the given point.
     * This will expand the width by the X coordinate of the point in both directions
     * and the height by Y coordinate of the point. Negative coordinate values will
     * shrink the rectangle. If <code>-p[X]</code> is larger than half of the width,
     * the X interval will contain only the X coordinate of the midpoint;
     * same for height. */
    void expandBy(CPoint const &p) {
        expandBy(p[X], p[Y]);
    }

    /** @brief Shrink the rectangle by the coordinates of the given point.
     * If the amount of shrinkage exceeds a half of the corresponding dimension, the
     * rectangle will be shrunk to a line segment.
     */
    void shrinkBy(CPoint const &p) {
        shrinkBy(p[X], p[Y]);
    }
    /// @}

    /// @name Return an expanded or shrunk copy of the rectangle.
    /// @{
    /** Return a new rectangle which results from expanding this one by the same
     * amount along both axes.
     */
    CRect expandedBy(Coord amount) const {
        auto copy{*this};
        copy.expandBy(amount);
        return copy;
    }

    /** Return a new rectangle which results from expanding this one by
     * possibly different amounts along both axes.
     */
    CRect expandedBy(Coord x, Coord y) const {
        auto copy{*this};
        copy.expandBy(x, y);
        return copy;
    }

    /** Return a new rectangle which results from expanding this one by
     * possibly different amounts along both axes.
     */
    CRect expandedBy(CPoint const &p) const {
        auto copy{*this};
        copy.expandBy(p);
        return copy;
    }

    /** Return a new rectangle which results from shrinking this one by the same
     * amount along both axes.
     */
    CRect shrunkBy(Coord amount) const {
        auto copy{*this};
        copy.shrinkBy(amount);
        return copy;
    }

    /** Return a new rectangle which results from shrinking this one by
     * possibly different amounts along both axes.
     */
    CRect shrunkBy(Coord x, Coord y) const {
        auto copy{*this};
        copy.shrinkBy(x, y);
        return copy;
    }

    /** Return a new rectangle which results from shrinking this one by
     * possibly different amounts along both axes.
     */
    CRect shrunkBy(CPoint const &p) const {
        auto copy{*this};
        copy.shrinkBy(p);
        return copy;
    }
    /// @}

    /// @name Operators
    /// @{
    /** @brief Offset the rectangle by a vector. */
    GenericRect<C> &operator+=(CPoint const &p) {
        f[X] += p[X];
        f[Y] += p[Y];
        return *this;
    }
    /** @brief Offset the rectangle by the negation of a vector. */
    GenericRect<C> &operator-=(CPoint const &p) {
        f[X] -= p[X];
        f[Y] -= p[Y];
        return *this;
    }
    /** @brief Union two rectangles. */
    GenericRect<C> &operator|=(CRect const &o) {
        unionWith(o);
        return *this;
    }
    GenericRect<C> &operator|=(OptCRect const &o) {
        unionWith(o);
        return *this;
    }
    /** @brief Test for equality of rectangles. */
    bool operator==(CRect const &o) const { return f[X] == o[X] && f[Y] == o[Y]; }
    /// @}
};

/**
 * @brief Axis-aligned generic rectangle that can be empty.
 * @ingroup Primitives
 */
template <typename C>
class GenericOptRect
    : public std::optional<typename CoordTraits<C>::RectType>
    , boost::equality_comparable< typename CoordTraits<C>::OptRectType
    , boost::equality_comparable< typename CoordTraits<C>::OptRectType, typename CoordTraits<C>::RectType
    , boost::orable< typename CoordTraits<C>::OptRectType
    , boost::andable< typename CoordTraits<C>::OptRectType
    , boost::andable< typename CoordTraits<C>::OptRectType, typename CoordTraits<C>::RectType
      >>>>>
{
    using CInterval = typename CoordTraits<C>::IntervalType;
    using OptCInterval = typename CoordTraits<C>::OptIntervalType;
    using CPoint = typename CoordTraits<C>::PointType;
    using CRect = typename CoordTraits<C>::RectType;
    using OptCRect = typename CoordTraits<C>::OptRectType;
    using Base = std::optional<CRect>;
public:
    using D1Value = CInterval;
    using D1Reference = CInterval &;
    using D1ConstReference = CInterval const &;

    /// @name Create potentially empty rectangles.
    /// @{
    GenericOptRect() = default;
    GenericOptRect(GenericRect<C> const &a) : Base(CRect(a)) {}
    GenericOptRect(CPoint const &a, CPoint const &b) : Base(CRect(a, b)) {}
    GenericOptRect(C x0, C y0, C x1, C y1) : Base(CRect(x0, y0, x1, y1)) {}
    /// Creates an empty OptRect when one of the argument intervals is empty.
    GenericOptRect(OptCInterval const &x_int, OptCInterval const &y_int) {
        if (x_int && y_int) {
            *this = CRect(*x_int, *y_int);
        }
        // else, stay empty.
    }

    /** @brief Create a rectangle from a range of points.
     * The resulting rectangle will contain all points from the range.
     * If the range contains no points, the result will be an empty rectangle.
     * The return type of iterators must be convertible to the corresponding
     * point type (Point or IntPoint).
     * @param start Beginning of the range
     * @param end   End of the range
     * @return Rectangle that contains all points from [start, end). */
    template <typename InputIterator>
    static OptCRect from_range(InputIterator start, InputIterator end) {
        OptCRect result;
        for (; start != end; ++start) {
            result.expandTo(*start);
        }
        return result;
    }
    /// @}

    /// @name Check other rectangles and points for inclusion.
    /// @{
    /** @brief Check for emptiness. */
    inline bool empty() const { return !*this; };
    /** @brief Check whether the rectangles have any common points.
     * Empty rectangles will not intersect with any other rectangle. */
    bool intersects(CRect const &r) const { return r.intersects(*this); }
    /** @brief Check whether the rectangle includes all points in the given rectangle.
     * Empty rectangles will be contained in any non-empty rectangle. */
    bool contains(CRect const &r) const { return *this && (*this)->contains(r); }

    /** @brief Check whether the rectangles have any common points.
     * Empty rectangles will not intersect with any other rectangle.
     * Two empty rectangles will not intersect each other. */
    bool intersects(OptCRect const &r) const { return *this && (*this)->intersects(r); }
    /** @brief Check whether the rectangle includes all points in the given rectangle.
     * Empty rectangles will be contained in any non-empty rectangle.
     * An empty rectangle will not contain other empty rectangles. */
    bool contains(OptCRect const &r) const { return *this && (*this)->contains(r); }

    /** @brief Check whether the given point is within the rectangle.
     * An empty rectangle will not contain any points. */
    bool contains(CPoint const &p) const { return *this && (*this)->contains(p); }
    /// @}

    /** @brief Compute the rectangle's area. */
    C area() const { return *this ? (*this)->area() : 0; }
    /** @brief Check whether the rectangle has zero area. */
    bool hasZeroArea() const { return !*this || (*this)->hasZeroArea(); }
    /** @brief Returns an empty optional (testing false) if the rectangle has zero area. */
    OptCRect regularized() const { return hasZeroArea() ? OptCRect() : *this; }

    /// @name Modify the potentially empty rectangle.
    /// @{
    /** @brief Enlarge the rectangle to contain the argument.
     * If this rectangle is empty, after callng this method it will
     * be equal to the argument. */
    void unionWith(CRect const &b) {
        if (*this) {
            (*this)->unionWith(b);
        } else {
            *this = b;
        }
    }
    /** @brief Enlarge the rectangle to contain the argument.
     * Unioning with an empty rectangle results in no changes.
     * If this rectangle is empty, after calling this method it will
     * be equal to the argument. */
    void unionWith(OptCRect const &b) {
        if (b) unionWith(*b);
    }
    /** @brief Leave only the area overlapping with the argument.
     * If the rectangles do not have any points in common, after calling
     * this method the rectangle will be empty. */
    void intersectWith(CRect const &b) {
        if (!*this) return;
        OptCInterval x = (**this)[X] & b[X], y = (**this)[Y] & b[Y];
        if (x && y) {
            *this = CRect(*x, *y);
        } else {
            *this = {};
        }
    }
    /** @brief Leave only the area overlapping with the argument.
     * If the argument is empty or the rectangles do not have any points
     * in common, after calling this method the rectangle will be empty. */
    void intersectWith(OptCRect const &b) {
        if (b) {
            intersectWith(*b);
        } else {
            *this = {};
        }
    }
    /** @brief Create or enlarge the rectangle to contain the given point.
     * If the rectangle is empty, after calling this method it will be non-empty
     * and it will contain only the given point. */
    void expandTo(CPoint const &p) {
        if (*this) {
            (*this)->expandTo(p);
        } else {
            *this = CRect(p, p);
        }
    }
    /// @}
    
    /// @name Operators
    /// @{
    /** @brief Union with @a b */
    GenericOptRect<C> &operator|=(OptCRect const &b) {
        unionWith(b);
        return *this;
    }
    /** @brief Intersect with @a b */
    GenericOptRect<C> &operator&=(CRect const &b) {
        intersectWith(b);
        return *this;
    }
    /** @brief Intersect with @a b */
    GenericOptRect<C> &operator&=(OptCRect const &b) {
        intersectWith(b);
        return *this;
    }
    /** @brief Test for equality.
     * All empty rectangles are equal. */
    bool operator==(OptCRect const &other) const {
        if (!*this != !other) return false;
        return *this ? **this == *other : true;
    }
    bool operator==(CRect const &other) const {
        if (!*this) return false;
        return **this == other;
    }
    /// @}
};

template <typename C>
inline void GenericRect<C>::unionWith(OptCRect const &b) { 
    if (b) {
        unionWith(*b);
    }
}
template <typename C>
inline bool GenericRect<C>::intersects(OptCRect const &r) const {
    return r && intersects(*r);
}
template <typename C>
inline bool GenericRect<C>::contains(OptCRect const &r) const {
    return !r || contains(*r);
}

template <typename C>
inline std::ostream &operator<<(std::ostream &out, GenericRect<C> const &r) {
    return out << "Rect " << r[X] << " x " << r[Y];
}

template <typename C>
inline std::ostream &operator<<(std::ostream &out, GenericOptRect<C> const &r) {
    return r ? (out << *r) : (out << "Rect (empty)");
}

} // namespace Geom

#endif // LIB2GEOM_SEEN_GENERIC_RECT_H

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
