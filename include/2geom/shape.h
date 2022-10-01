/** @file Shape - A planar region with a piecewise-smooth boundary.
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

#ifndef LIB2GEOM_SEEN_SHAPE_H
#define LIB2GEOM_SEEN_SHAPE_H

#include <2geom/affine.h>
#include <2geom/circle.h>
#include <2geom/coord.h>
#include <2geom/parallelogram.h>
#include <2geom/pathvector.h>
#include <2geom/rect.h>

namespace Geom {

/** @brief Determines how a two-dimensional shape is described by a PathVector.
 * @ingroup ShapeOps
 */
enum class FillRule {

    /** @brief Fill rule following the mathematical convention.
     *
     *  Applicable only to a PathVector without transverse self-intersections,
     *  the CONVENTIONAL rule can be summarized as follows:
     *
     *  If a path travels along the x-axis towards positive x, then the points
     *  with y > 0 are inside the shape and points with y < 0 are outside.
     *
     *  Using the mathematics orientation (y-axis pointing up), this means that each
     *  path in the PathVector has the shape "to its left" and the outside "to the
     *  right" when looking in the direction of increasing time coordinate on the path.
     *
     *  Using the computer graphics convention, (y-axis pointing down), the shape is
     *  always "to the right" and the outside is "to the left" when travelling along
     *  each path making up the boundary.
     */
    CONVENTIONAL,

    /** @brief Fill the parts of the plane that have odd winding numbers.
     *
     *  The even-odd fill rule means that the two-dimensional shape consists of
     *  the points of the plane about which the path-vector has an odd winding number.
     */
    EVEN_ODD,

    /** @brief Fill the parts of the plane that have nonzero winding numbers.
     *
     *  Nonzero fill rule means that the two-dimensional shape consists of all
     *  parts of the plane about which the path-vector has nonzero winding numbers.
     */
    NONZERO
};

/** @class Shape
 * @brief A two-dimensional shape in the plane.
 *
 * A Shape object represents a bounded region of the 2D-plane, possibly disconnected
 * and/or containing holes. The contour of a shape is a collection of closed paths.
 * Internally, this boundary is stored as a PathVector for which the FillRule::CONVENTIONAL
 * fill rule is satisfied.
 *
 * It is not specified if the mathematical set represented by a Shape is open or closed.
 * Whether or not the contour is considered part of the shape may depend on the nature of
 * the operation performed on a Shape object.
 *
 * Since constructing a Shape object from a PathVector with self-intersections is an
 * expensive operation, the user is encouraged to hold on to Shape objects to avoid
 * unnecessary work. In contrast, constructing a Shape from 2D primitives such as
 * Circle, Rect, Ellipse, or Parallelogram is cheap.
 *
 * A Shape object can be implicitly converted to a PathVector.
 * Read-only copy-free access to the shape's contour is also possible,
 * via the contour() member function.
 * @ingroup ShapeOps
 */
class Shape
{
private:
    Coord _epsilon = EPSILON; ///< Numerical precision parameter tailored to this particular shape.
    PathVector _boundary;     ///< A regularized path-vector description of the shape's boundary.

    /// Relates precision to the size of a shape.
    inline static Coord const PRECISION_COEFFICIENT = 0x1p-20;

public:
    /** @brief Construct an empty shape. */
    Shape() = default;

    /** @brief Construct a shape from a PathVector via the specified fill rule.
     *
     * @param pathvector A PathVector describing the shape.
     * @param fill_rule  The fill rule to use when deciding which portions of the plane are
     *                   part of the shape and which are not – see Geom::FillRule for details.
     * @param close_open_paths If true, any open paths in the passed path-vector will be closed
     *                         with line segments. If false, open paths will be ignored.
     */
    Shape(PathVector const &pathvector, FillRule fill_rule, bool close_open_paths = true)
        : _epsilon{_calculateEpsilon(pathvector)}
        , _boundary{flatten(pathvector, fill_rule, close_open_paths, _epsilon)}
    {}

    /** @brief Construct a Shape from a rectangle.
     *
     * If the rectangle has zero area, the resulting shape is empty.
     * @throws UnboundedRegion Thrown when the rectangle is unbounded.
     */
    Shape(Rect const &rectangle);

    /** @brief Construct a Shape from a parallelogram.
     *
     * If the parallelogram has zero area, the resulting shape is empty.
     */
    Shape(Parallelogram const &parallelogram);

    /** @brief Construct a Shape from a circle.
     *
     * If the circle is degenerate, the resulting shape will be empty.
     */
    Shape(Circle const &circle) noexcept;

    /** @brief Construct a Shape from an ellipse.
     *
     * If the ellipse is degenerate, the resulting shape will be empty.
     */
    Shape(Ellipse const &circle) noexcept;

    Shape(Shape const &) = default;
    Shape(Shape &&) = default;
    Shape &operator=(Shape const &) = default;
    Shape &operator=(Shape &&) = default;
    ~Shape() = default;

    /** Implicit conversion to a PathVector. */
    operator PathVector() const { return _boundary; }

    /** Read-only access to the contour of the shape. */
    PathVector const &contour() const { return _boundary; }

    /** Get the numerical precision setting for this shape. */
    Coord precision() const { return _epsilon; }

    /** Transform a shape by an affine map. */
    Shape &operator*=(Affine const &tr);
    friend Shape operator*(Shape const &sh, Affine const &tr) { return Shape(sh) *= tr; }

    /** Check whether the shape has an empty interior. */
    bool empty() const { return _boundary.empty(); }

    /** Get an axis-aligned rectangle fully enclosing the shape. */
    OptRect boundsFast() const { return _boundary.boundsFast(); }

    /** Get the smallest possible axis-aligned rectangle fully enclosing the shape. */
    OptRect boundsExact() const {return _boundary.boundsExact(); }

    /** @brief Get the number of the shape's boundary components.
     *
     * A boundary component is a closed path which separates the shape from its exterior.
     * For example, a square has 1 boundary component and an annulus has 2.
     */
    size_t numBoundaryComponents() const { return _boundary.size(); }

    /** Compute the area of the shape. */
    Coord area() const;

    /*     TO DO LIST (2022-10-02):
     *
     * [ ] Binary boolean operations (operators: &, |, -, ^); in-place versions (&=, |=, -=, ^=).
     * [ ] bool Shape::contains(Point const &point) const
     * [ ] bool Shape::intersects(Shape const &other) const [much faster than `operator&`].
     * [ ] bool operator==(Shape const &other) const { return (*this ^ other).empty(); }
     * [ ] Split a shape into connected components.
     * [ ] Fracture: return all bounded connected regions into which a PV cuts the plane.
     * [ ] Cookie-cutting (split the shape into several shapes along a PV).
     * [ ] n-ary boolean operations: union and intersection.
     * [ ] Growing and shrinking by a specified amount (outlining).
     * [ ] Distance from a point to a shape (infimum type; maybe Hausdorff distance too?).
     * [ ] Return the convex hull of a shape as another shape.
     */

    /** @brief Convert any PathVector and a fill rule into a conventional path-vector
     *         description of the boundary of a 2D shape.
     *
     * A path-vector and a fill rule determine a 2D shape in the plane. The flatten
     * algorithm returns a path-vector representing the boundary of this shape in a
     * non-redundant way. All paths in the returned path-vector have the property that
     * they separate the interior of the shape from its exterior. Furthermore, the shape
     * is always on the left and the exterior on the right (assuming y-axis up) when
     * traveling along the paths, in the direction of the increasing time coordinate.
     * (When the y-axis points down, the interior is on the right and exterior on the left).
     *
     * In particular, the returned path-vector doesn't have transverse self-intersections,
     * though it might have non-transverse multiple points.
     * If the shape has an empty interior, an empty path-vector is returned.
     *
     * @param pv A path-vector describing a 2D shape.
     * @param fill_rule Specifies which parts of the plane are part of the shape.
     * @param close_open_paths If true, any open paths in the passed path-vector will
     *                         be closed with line segments before being processed.
     *                         If false, open paths will be ignored.
     * @param precision The numerical epsilon for overlap and crossing detection.
     * @return A path-vector describing the contour of the shape with the FillRule::CONVENTIONAL fill rule.
     *
     * See @ref flatten_description for technical details.
     */
    static PathVector flatten(PathVector const &pv, FillRule fill_rule, bool close_open_paths, Coord precision);

private:
    /** Calculate a numerical precision depending on the size of this shape. */
    static Coord _calculateEpsilon(PathVector const &pv);
};

} // namespace Geom

#endif // LIB2GEOM_SEEN_SHAPE_H

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
