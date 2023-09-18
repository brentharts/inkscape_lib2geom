/* Bezier curve implementation
 *
 * Authors:
 *   MenTaLguY <mental@rydia.net>
 *   Marco Cecchetti <mrcekets at gmail.com>
 *   Krzysztof Kosiński <tweenk.pl@gmail.com>
 *
 * Copyright 2007-2009 Authors
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

#include <2geom/basic-intersection.h>
#include <2geom/bezier-curve.h>
#include <2geom/bezier-utils.h>
#include <2geom/nearest-time.h>
#include <2geom/path-sink.h>
#include <2geom/polynomial.h>

namespace Geom {

/**
 * @class BezierCurve
 * @brief Two-dimensional Bezier curve of arbitrary order.
 *
 * Bezier curves are an expansion of the concept of linear interpolation to n points.
 * Linear segments in 2Geom are in fact Bezier curves of order 1.
 *
 * Let \f$\mathbf{B}_{\mathbf{p}_0\mathbf{p}_1\ldots\mathbf{p}_n}\f$ denote a Bezier curve
 * of order \f$n\f$ defined by the points \f$\mathbf{p}_0, \mathbf{p}_1, \ldots, \mathbf{p}_n\f$.
 * Bezier curve of order 1 is a linear interpolation curve between two points, defined as
 * \f[ \mathbf{B}_{\mathbf{p}_0\mathbf{p}_1}(t) = (1-t)\mathbf{p}_0 + t\mathbf{p}_1 \f]
 * If we now substitute points \f$\mathbf{p_0}\f$ and \f$\mathbf{p_1}\f$ in this definition
 * by linear interpolations, we get the definition of a Bezier curve of order 2, also called
 * a quadratic Bezier curve.
 * \f{align*}{ \mathbf{B}_{\mathbf{p}_0\mathbf{p}_1\mathbf{p}_2}(t)
       &= (1-t) \mathbf{B}_{\mathbf{p}_0\mathbf{p}_1}(t) + t \mathbf{B}_{\mathbf{p}_1\mathbf{p}_2}(t) \\
     \mathbf{B}_{\mathbf{p}_0\mathbf{p}_1\mathbf{p}_2}(t)
       &= (1-t)^2\mathbf{p}_0 + 2(1-t)t\mathbf{p}_1 + t^2\mathbf{p}_2 \f}
 * By substituting points for quadratic Bezier curves in the original definition,
 * we get a Bezier curve of order 3, called a cubic Bezier curve.
 * \f{align*}{ \mathbf{B}_{\mathbf{p}_0\mathbf{p}_1\mathbf{p}_2\mathbf{p}_3}(t)
       &= (1-t) \mathbf{B}_{\mathbf{p}_0\mathbf{p}_1\mathbf{p}_2}(t)
       + t \mathbf{B}_{\mathbf{p}_1\mathbf{p}_2\mathbf{p}_3}(t) \\
     \mathbf{B}_{\mathbf{p}_0\mathbf{p}_1\mathbf{p}_2\mathbf{p}_3}(t)
       &= (1-t)^3\mathbf{p}_0+3(1-t)^2t\mathbf{p}_1+3(1-t)t^2\mathbf{p}_2+t^3\mathbf{p}_3 \f}
 * In general, a Bezier curve or order \f$n\f$ can be recursively defined as
 * \f[ \mathbf{B}_{\mathbf{p}_0\mathbf{p}_1\ldots\mathbf{p}_n}(t)
     = (1-t) \mathbf{B}_{\mathbf{p}_0\mathbf{p}_1\ldots\mathbf{p}_{n-1}}(t)
     + t \mathbf{B}_{\mathbf{p}_1\mathbf{p}_2\ldots\mathbf{p}_n}(t) \f]
 *
 * This substitution can be repeated an arbitrary number of times. To picture this, imagine
 * the evaluation of a point on the curve as follows: first, all control points are joined with
 * straight lines, and a point corresponding to the selected time value is marked on them.
 * Then, the marked points are joined with straight lines and the point corresponding to
 * the time value is marked. This is repeated until only one marked point remains, which is the
 * point at the selected time value.
 *
 * @image html bezier-curve-evaluation.png "Evaluation of the Bezier curve"
 *
 * An important property of the Bezier curves is that their parameters (control points)
 * have an intuitive geometric interpretation. Because of this, they are frequently used
 * in vector graphics editors.
 *
 * Every Bezier curve is contained in its control polygon (the convex polygon composed
 * of its control points). This fact is useful for sweepline algorithms and intersection.
 *
 * @par Implementation notes
 * The order of a Bezier curve is immuable once it has been created. Normally, you should
 * know the order at compile time and use the BezierCurveN template. If you need to determine
 * the order at runtime, use the BezierCurve::create() function. It will create a BezierCurveN
 * for orders 1, 2 and 3 (up to cubic Beziers), so you can later <tt>dynamic_cast</tt>
 * to those types, and for higher orders it will create an instance of BezierCurve.
 *
 * @relates BezierCurveN
 * @ingroup Curves
 */

/**
 * @class BezierCurveN
 * @brief Bezier curve with compile-time specified order.
 *
 * @tparam degree unsigned value indicating the order of the Bezier curve
 *
 * @relates BezierCurve
 * @ingroup Curves
 */


BezierCurve::BezierCurve(std::vector<Point> const &pts)
    : inner(pts)
{
    if (pts.size() < 2) {
        THROW_RANGEERROR("Bezier curve must have at least 2 control points");
    }
}

bool BezierCurve::isDegenerate() const
{
    for (unsigned d = 0; d < 2; ++d) {
        Coord ic = inner[d][0];
        for (unsigned i = 1; i < size(); ++i) {
            if (inner[d][i] != ic)
                return false;
        }
    }
    return true;
}

/** Return false if there are at least 3 distinct control points, true otherwise. */
bool BezierCurve::isLineSegment() const
{
    auto const last_idx = size() - 1;
    if (last_idx == 1) {
        return true;
    }
    auto const start = controlPoint(0);
    auto const end = controlPoint(last_idx);
    for (unsigned i = 1; i < last_idx; ++i) {
        auto const pi = controlPoint(i);
        if (pi != start && pi != end) {
            return false;
        }
    }
    return true;
}

void BezierCurve::expandToTransformed(Rect &bbox, Affine const &transform) const
{
    bbox |= bounds_exact(inner * transform);
}

Coord BezierCurve::length(Coord tolerance) const
{
    switch (order()) {
        case 0:
            return 0.0;
        case 1:
            return distance(initialPoint(), finalPoint());
        case 2: {
            std::vector<Point> pts = controlPoints();
            return bezier_length(pts[0], pts[1], pts[2], tolerance);
        }
        case 3: {
            std::vector<Point> pts = controlPoints();
            return bezier_length(pts[0], pts[1], pts[2], pts[3], tolerance);
        }
        default:
            return bezier_length(controlPoints(), tolerance);
    }
}

std::vector<CurveIntersection> BezierCurve::intersect(Curve const &other, Coord eps) const
{
    std::vector<CurveIntersection> result;

    // in case we encounter an order-1 curve created from a vector
    // or a degenerate elliptical arc
    if (isLineSegment()) {
        LineSegment ls(initialPoint(), finalPoint());
        result = ls.intersect(other);
        return result;
    }

    // here we are sure that this curve is at least a quadratic Bezier
    BezierCurve const *bez = dynamic_cast<BezierCurve const *>(&other);
    if (bez) {
        std::vector<std::pair<double, double>> xs;
        find_intersections(xs, inner, bez->inner, eps);
        for (auto &i : xs) {
            CurveIntersection x(*this, other, i.first, i.second);
            result.push_back(x);
        }
        return result;
    }

    // pass other intersection types to the other curve
    result = other.intersect(*this, eps);
    transpose_in_place(result);
    return result;
}

bool BezierCurve::isNear(Curve const &c, Coord precision) const
{
    if (this == &c)
        return true;

    BezierCurve const *other = dynamic_cast<BezierCurve const *>(&c);
    if (!other)
        return false;

    if (!are_near(inner.at0(), other->inner.at0(), precision))
        return false;
    if (!are_near(inner.at1(), other->inner.at1(), precision))
        return false;

    if (size() == other->size()) {
        for (unsigned i = 1; i < order(); ++i) {
            if (!are_near(inner.point(i), other->inner.point(i), precision)) {
                return false;
            }
        }
        return true;
    } else {
        // Must equalize the degrees before comparing
        BezierCurve elevated_this, elevated_other;
        for (size_t dim = 0; dim < 2; dim++) {
            unsigned const our_degree = inner[dim].degree();
            unsigned const other_degree = other->inner[dim].degree();

            if (our_degree < other_degree) {
                // Elevate our degree
                elevated_this.inner[dim] = inner[dim].elevate_to_degree(other_degree);
                elevated_other.inner[dim] = other->inner[dim];
            } else if (our_degree > other_degree) {
                // Elevate the other's degree
                elevated_this.inner[dim] = inner[dim];
                elevated_other.inner[dim] = other->inner[dim].elevate_to_degree(our_degree);
            } else {
                // Equal degrees: just copy
                elevated_this.inner[dim] = inner[dim];
                elevated_other.inner[dim] = other->inner[dim];
            }
        }
        assert(elevated_other.size() == elevated_this.size());
        return elevated_this.isNear(elevated_other, precision);
    }
}

Curve *BezierCurve::portion(Coord f, Coord t) const
{
    if (f == 0.0 && t == 1.0) {
        return duplicate();
    }
    if (f == 1.0 && t == 0.0) {
        return reverse();
    }
    return new BezierCurve(Geom::portion(inner, f, t));
}

bool BezierCurve::operator==(Curve const &c) const
{
    if (this == &c)
        return true;

    BezierCurve const *other = dynamic_cast<BezierCurve const *>(&c);
    if (!other)
        return false;
    if (size() != other->size())
        return false;

    for (unsigned i = 0; i < size(); ++i) {
        if (controlPoint(i) != other->controlPoint(i))
            return false;
    }
    return true;
}

Coord BezierCurve::nearestTime(Point const &p, Coord from, Coord to) const { return nearest_time(p, inner, from, to); }

void BezierCurve::feed(PathSink &sink, bool moveto_initial) const
{
    if (size() > 4) {
        Curve::feed(sink, moveto_initial);
        return;
    }

    Point ip = controlPoint(0);
    if (moveto_initial) {
        sink.moveTo(ip);
    }
    switch (size()) {
        case 2:
            sink.lineTo(controlPoint(1));
            break;
        case 3:
            sink.quadTo(controlPoint(1), controlPoint(2));
            break;
        case 4:
            sink.curveTo(controlPoint(1), controlPoint(2), controlPoint(3));
            break;
        default:
            // TODO: add a path sink method that accepts a vector of control points
            //       and converts to cubic spline by default
            assert(false);
            break;
    }
}

BezierCurve *BezierCurve::create(std::vector<Point> const &pts)
{
    switch (pts.size()) {
        case 0:
        case 1:
            THROW_LOGICALERROR("BezierCurve::create: too few points in vector");
            return NULL;
        case 2:
            return new LineSegment(pts[0], pts[1]);
        case 3:
            return new QuadraticBezier(pts[0], pts[1], pts[2]);
        case 4:
            return new CubicBezier(pts[0], pts[1], pts[2], pts[3]);
        default:
            return new BezierCurve(pts);
    }
}

/**
 * Computes the times where the radius of curvature of the bezier curve equals the given radius.
 */
std::vector<Coord> BezierCurve::timesWithRadiusOfCurvature(double radius) const
{
    /**
     * The algorithm works as follows:
     * Find the solutions of curvature of curve at t = curvatureValue
     * This is equivalent to (dx*ddy-ddx*dy)/curvatureValue = (dx**2+dy**2)**(3/2)
     * When the left side is positive, taking the square gives
     * ((dx*ddy-ddx*dy)/curvatureValue)**2 - (dx**2+dy**2)**3 = 0
     * This is a polyomial for BezierCurves and can be solved with root finding algos.
     */

    if (this->size() <= 2) {
        return {};
    }

    std::vector<Coord> res;

    auto const dx = Geom::derivative(inner[X]);
    auto const dy = Geom::derivative(inner[Y]);
    auto const ddx = Geom::derivative(dx);
    auto const ddy = Geom::derivative(dy);
    auto const c0 = (dx * ddy - ddx * dy) * radius;
    auto const c1 = dx * dx + dy * dy;
    auto const p = c0 * c0 - c1 * c1 * c1;
    auto const candidates = p.roots();
    // check which candidates have positive nominator
    // as squaring also give negative (spurious) results
    for (Coord const candidate : candidates) {
        if (c0.valueAt(candidate) > 0) {
            res.push_back(candidate);
        }
    }

    return res;
}

// optimized specializations for LineSegment

template <>
Curve *BezierCurveN<1>::derivative() const
{
    double dx = inner[X][1] - inner[X][0], dy = inner[Y][1] - inner[Y][0];
    return new BezierCurveN<1>(Point(dx, dy), Point(dx, dy));
}

template <>
Coord BezierCurveN<1>::nearestTime(Point const &p, Coord from, Coord to) const
{
    using std::swap;

    if (from > to)
        swap(from, to);
    Point ip = pointAt(from);
    Point fp = pointAt(to);
    Point v = fp - ip;
    Coord l2v = L2sq(v);
    if (l2v == 0)
        return 0;
    Coord t = dot(p - ip, v) / l2v;
    if (t <= 0)
        return from;
    else if (t >= 1)
        return to;
    else
        return from + t * (to - from);
}

/* Specialized intersection routine for line segments.
 *
 * NOTE: if the segments overlap in part or in full, the function returns the start and end
 * of the overlapping subsegment as intersections. This behavior is more useful than throwing
 * Geom::InfinitelyManySolutions.
 */
template <>
std::vector<CurveIntersection> BezierCurveN<1>::intersect(Curve const &other, Coord eps) const
{
    std::vector<CurveIntersection> result;

    // only handle intersections with other LineSegments here
    if (!other.isLineSegment()) {
        // pass all other types to the other curve
        result = other.intersect(*this, eps);
        transpose_in_place(result);
        return result;
    }

    Point const u = finalPoint() - initialPoint();
    Point const v = other.initialPoint() - other.finalPoint();
    if (u.isZero() || v.isZero()) {
        return {};
    }
    Coord const uv = u.length() * v.length();
    Coord const u_cross_v = cross(u, v);
    bool const segments_are_parallel = std::abs(u_cross_v) < eps * uv;

    if (segments_are_parallel) {
        // We check if the segments lie on the same line.
        Coord const distance_between_lines = std::abs(cross(u.normalized(), other.initialPoint() - initialPoint()));
        if (distance_between_lines > eps) {
            // Segments are parallel but aren't part of the same line => no intersections.
            return {};
        }
        // The segments are on the same line, so they may overlap in part or in full.
        // Look for the times on this segment's line at which the line passes through
        // the initial and final points of the other segment.
        Coord const ulen_inverse = 1.0 / u.length();
        auto const time_of_passage = [&](Point const &point_on_line) -> Coord {
            return dot(u.normalized(), point_on_line - initialPoint()) * ulen_inverse;
        };
        // Find the range of times on our segment where we travel through the other segment.
        auto time_in_other = Interval(time_of_passage(other.initialPoint()), time_of_passage(other.finalPoint()));
        Coord const eps_utime = eps * ulen_inverse;
        if (time_in_other.min() > 1 + eps_utime || time_in_other.max() < -eps_utime) {
            return {};
        }

        // Create two intersections, one at each end of the overlap interval.
        Coord last_time = infinity();
        for (Coord t : { time_in_other.min(), time_in_other.max() }) {
            t = std::clamp(t, 0.0, 1.0);
            if (t == last_time) {
                continue;
            }
            last_time = t;
            auto const point = pointAt(t);
            Coord const other_t = std::clamp(dot(v.normalized(), other.initialPoint() - point) / v.length(), 0.0, 1.0);
            auto const other_pt = other.pointAt(other_t);
            if (distance(point, other_pt) > eps) {
                continue;
            }
            result.emplace_back(t, other_t, middle_point(point, other_pt));
        }
        return result;
    } else {
        // Segments are not collinear - there is 0 or 1 intersection.
        // This is the typical case so it's important for it to be fast.
        Point const w = other.initialPoint() - initialPoint();
        Coord candidate_time_this = cross(w, v) / u_cross_v;

        /// Filter out candidate times outside of the interval [0, 1] fuzzed so as to accommodate
        /// for epsilon.
        auto const is_outside_seg = [eps](Coord t, Coord len) -> bool {
            Coord const arclen_coord = t * len;
            return arclen_coord < -eps || arclen_coord > len + eps;
        };

        if (is_outside_seg(candidate_time_this, u.length())) {
            return {};
        }

        Coord candidate_time_othr = cross(u, w) / u_cross_v;
        if (is_outside_seg(candidate_time_othr, v.length())) {
            return {};
        }

        // Even if there was some fuzz in the interval, we must now restrict to [0, 1] and verify
        // that the intersection found is still within epsilon precision (the fuzz may have been too
        // large, depending on the angle between the intervals).
        candidate_time_this = std::clamp(candidate_time_this, 0.0, 1.0);
        candidate_time_othr = std::clamp(candidate_time_othr, 0.0, 1.0);

        auto const pt_this = pointAt(candidate_time_this);
        auto const pt_othr = other.pointAt(candidate_time_othr);
        if (distance(pt_this, pt_othr) > eps) {
            return {};
        }
        return { CurveIntersection(candidate_time_this, candidate_time_othr, middle_point(pt_this, pt_othr)) };
    }
}

/** @brief Find intersections of a low-degree Bézier curve with a line segment.
 *
 * Uses algebraic solutions to low-degree polynomial equations which may be faster
 * and more precise than iterative methods.
 *
 * @tparam degree The degree of the Bézier curve; must be 2 or 3.
 * @param curve A Bézier curve of the given degree.
 * @param line A line (but really a segment).
 * @return Intersections between the passed curve and the fundamental segment of the line
 *         (the segment where the time parameter lies in the unit interval).
 */
template <unsigned degree>
static std::vector<CurveIntersection> bezier_line_intersections(BezierCurveN<degree> const &curve, Line const &line)
{
    static_assert(degree == 2 || degree == 3, "bezier_line_intersections<degree>() error: degree must be 2 or 3.");

    auto const length = distance(line.initialPoint(), line.finalPoint());
    if (length == 0) {
        return {};
    }
    std::vector<CurveIntersection> result;

    // Find the isometry mapping the line to the x-axis, taking the initial point to the origin
    // and the final point to (length, 0). Apply this transformation to the Bézier curve and
    // extract the y-coordinate polynomial.
    auto const transform = line.rotationToZero(Y);
    Bezier const y = (curve.fragment() * transform)[Y];
    std::vector<double> roots;

    // Find roots of the polynomial y.
    {
        double const c2 = y[0] + y[2] - 2.0 * y[1];
        double const c1 = y[1] - y[0];
        double const c0 = y[0];

        if constexpr (degree == 2) {
            roots = solve_quadratic(c2, 2.0 * c1, c0);
        } else if constexpr (degree == 3) {
            double const c3 = y[3] - y[0] + 3.0 * (y[1] - y[2]);
            roots = solve_cubic(c3, 3.0 * c2, 3.0 * c1, c0);
        }
    }

    // Filter the roots and assemble intersections.
    for (double root : roots) {
        if (root < 0.0 || root > 1.0) {
            continue;
        }
        Coord x = (curve.pointAt(root) * transform)[X];
        if (x < 0.0 || x > length) {
            continue;
        }
        result.emplace_back(curve, line, root, x / length);
    }
    return result;
}

/* Specialized intersection routine for quadratic Bézier curves. */
template <>
std::vector<CurveIntersection> BezierCurveN<2>::intersect(Curve const &other, Coord eps) const
{
    if (auto other_bezier = dynamic_cast<BezierCurve const *>(&other)) {
        auto const other_degree = other_bezier->order();
        if (other_degree == 1) {
            // Use the exact method to intersect a quadratic Bézier with a line segment.
            auto line = Line(other_bezier->initialPoint(), other_bezier->finalPoint());
            return bezier_line_intersections<2>(*this, line);
        }
        // TODO: implement exact intersection of two quadratic Béziers using the method of resultants.
    }
    return BezierCurve::intersect(other, eps);
}

/* Specialized intersection routine for cubic Bézier curves. */
template <>
std::vector<CurveIntersection> BezierCurveN<3>::intersect(Curve const &other, Coord eps) const
{
    if (auto other_bezier = dynamic_cast<BezierCurve const *>(&other)) {
        if (other_bezier->order() == 1) {
            // Use the exact method to intersect a cubic Bézier with a line segment.
            auto line = Line(other_bezier->initialPoint(), other_bezier->finalPoint());
            return bezier_line_intersections<3>(*this, line);
        }
    }
    return BezierCurve::intersect(other, eps);
}

template <>
int BezierCurveN<1>::winding(Point const &p) const
{
    Point ip = inner.at0(), fp = inner.at1();
    if (p[Y] == std::max(ip[Y], fp[Y]))
        return 0;

    Point v = fp - ip;
    assert(v[Y] != 0);
    Coord t = (p[Y] - ip[Y]) / v[Y];
    assert(t >= 0 && t <= 1);
    Coord xcross = lerp(t, ip[X], fp[X]);
    if (xcross > p[X]) {
        return v[Y] > 0 ? 1 : -1;
    }
    return 0;
}

template <>
void BezierCurveN<1>::feed(PathSink &sink, bool moveto_initial) const
{
    if (moveto_initial) {
        sink.moveTo(controlPoint(0));
    }
    sink.lineTo(controlPoint(1));
}

template <>
void BezierCurveN<2>::feed(PathSink &sink, bool moveto_initial) const
{
    if (moveto_initial) {
        sink.moveTo(controlPoint(0));
    }
    sink.quadTo(controlPoint(1), controlPoint(2));
}

template <>
void BezierCurveN<3>::feed(PathSink &sink, bool moveto_initial) const
{
    if (moveto_initial) {
        sink.moveTo(controlPoint(0));
    }
    sink.curveTo(controlPoint(1), controlPoint(2), controlPoint(3));
}

static void bezier_expand_to_image(Rect &range, Point const &x0, Point const &x1, Point const &x2)
{
    for (auto i : { X, Y }) {
        bezier_expand_to_image(range[i], x0[i], x1[i], x2[i]);
    }
}

static void bezier_expand_to_image(Rect &range, Point const &x0, Point const &x1, Point const &x2, Point const &x3)
{
    for (auto i : { X, Y }) {
        bezier_expand_to_image(range[i], x0[i], x1[i], x2[i], x3[i]);
    }
}

template <>
void BezierCurveN<1>::expandToTransformed(Rect &bbox, Affine const &transform) const
{
    bbox.expandTo(finalPoint() * transform);
}

template <>
void BezierCurveN<2>::expandToTransformed(Rect &bbox, Affine const &transform) const
{
    bezier_expand_to_image(bbox, controlPoint(0) * transform, controlPoint(1) * transform, controlPoint(2) * transform);
}

template <>
void BezierCurveN<3>::expandToTransformed(Rect &bbox, Affine const &transform) const
{
    bezier_expand_to_image(bbox, controlPoint(0) * transform, controlPoint(1) * transform, controlPoint(2) * transform,
                           controlPoint(3) * transform);
}

static Coord bezier_length_internal(std::vector<Point> &v1, Coord tolerance, int level)
{
    /* The Bezier length algorithm used in 2Geom utilizes a simple fact:
     * the Bezier curve is longer than the distance between its endpoints
     * but shorter than the length of the polyline formed by its control
     * points. When the difference between the two values is smaller than the
     * error tolerance, we can be sure that the true value is no further than
     * 0.5 * tolerance from their arithmetic mean. When it's larger, we recursively
     * subdivide the Bezier curve into two parts and add their lengths.
     *
     * We cap the maximum number of subdivisions at 256, which corresponds to 8 levels.
     */
    Coord lower = distance(v1.front(), v1.back());
    Coord upper = 0.0;
    for (size_t i = 0; i < v1.size() - 1; ++i) {
        upper += distance(v1[i], v1[i + 1]);
    }
    if (upper - lower <= 2 * tolerance || level >= 8) {
        return (lower + upper) / 2;
    }


    std::vector<Point> v2 = v1;

    /* Compute the right subdivision directly in v1 and the left one in v2.
     * Explanation of the algorithm used:
     * We have to compute the left and right edges of this triangle in which
     * the top row are the control points of the Bezier curve, and each cell
     * is equal to the arithmetic mean of the cells directly above it
     * to the right and left. This corresponds to subdividing the Bezier curve
     * at time value 0.5: the left edge has the control points of the first
     * portion of the Bezier curve and the right edge - the second one.
     * In the example we subdivide a curve with 5 control points (order 4).
     *
     * Start:
     * 0 1 2 3 4
     *  ? ? ? ?
     *   ? ? ?
     *    ? ?
     *     ?
     * # means we have overwritten the value, ? means we don't know
     * the value yet. Numbers mean the value is at i-th position in the vector.
     *
     * After loop with i==1
     * # 1 2 3 4
     *  0 ? ? ? -> write 0 to v2[1]
     *   ? ? ?
     *    ? ?
     *     ?
     *
     * After loop with i==2
     * # # 2 3 4
     *  # 1 ? ?
     *   0 ? ? -> write 0 to v2[2]
     *    ? ?
     *     ?
     *
     * After loop with i==3
     * # # # 3 4
     *  # # 2 ?
     *   # 1 ?
     *    0 ? -> write 0 to v2[3]
     *     ?
     *
     * After loop with i==4, we have the right edge of the triangle in v1,
     * and we write the last value needed for the left edge in v2[4].
     */

    for (size_t i = 1; i < v1.size(); ++i) {
        for (size_t j = i; j > 0; --j) {
            v1[j - 1] = 0.5 * (v1[j - 1] + v1[j]);
        }
        v2[i] = v1[0];
    }

    return bezier_length_internal(v1, 0.5 * tolerance, level + 1) +
           bezier_length_internal(v2, 0.5 * tolerance, level + 1);
}

/** @brief Compute the length of a bezier curve given by a vector of its control points
 * @relatesalso BezierCurve */
Coord bezier_length(std::vector<Point> const &points, Coord tolerance)
{
    if (points.size() < 2)
        return 0.0;
    std::vector<Point> v1 = points;
    return bezier_length_internal(v1, tolerance, 0);
}

static Coord bezier_length_internal(Point a0, Point a1, Point a2, Coord tolerance, int level)
{
    Coord lower = distance(a0, a2);
    Coord upper = distance(a0, a1) + distance(a1, a2);

    if (upper - lower <= 2 * tolerance || level >= 8) {
        return (lower + upper) / 2;
    }

    Point // Casteljau subdivision
          // b0 = a0,
          // c0 = a2,
        b1 = 0.5 * (a0 + a1),
        c1 = 0.5 * (a1 + a2),
        b2 = 0.5 * (b1 + c1); // == c2
    return bezier_length_internal(a0, b1, b2, 0.5 * tolerance, level + 1) +
           bezier_length_internal(b2, c1, a2, 0.5 * tolerance, level + 1);
}

/** @brief Compute the length of a quadratic bezier curve given by its control points
 * @relatesalso QuadraticBezier */
Coord bezier_length(Point a0, Point a1, Point a2, Coord tolerance)
{
    return bezier_length_internal(a0, a1, a2, tolerance, 0);
}

static Coord bezier_length_internal(Point a0, Point a1, Point a2, Point a3, Coord tolerance, int level)
{
    Coord lower = distance(a0, a3);
    Coord upper = distance(a0, a1) + distance(a1, a2) + distance(a2, a3);

    if (upper - lower <= 2 * tolerance || level >= 8) {
        return (lower + upper) / 2;
    }

    Point // Casteljau subdivision
          // b0 = a0,
          // c0 = a3,
        b1 = 0.5 * (a0 + a1),
        t0 = 0.5 * (a1 + a2), c1 = 0.5 * (a2 + a3), b2 = 0.5 * (b1 + t0), c2 = 0.5 * (t0 + c1),
        b3 = 0.5 * (b2 + c2); // == c3
    return bezier_length_internal(a0, b1, b2, b3, 0.5 * tolerance, level + 1) +
           bezier_length_internal(b3, c2, c1, a3, 0.5 * tolerance, level + 1);
}

/** @brief Compute the length of a cubic bezier curve given by its control points
 * @relatesalso CubicBezier */
Coord bezier_length(Point a0, Point a1, Point a2, Point a3, Coord tolerance)
{
    return bezier_length_internal(a0, a1, a2, a3, tolerance, 0);
}

template <>
Path BezierCurveN<1>::offsetted(double width, double tolerance, bool no_crossing) const
{
    Path ret;
    assert(!isDegenerate());

    // control points difference is used instead of unitTangent,
    // because isDegenerate just checks for equality and unitTangent checks with small tolerance.
    Point const tangent = controlPoint(1) - controlPoint(0);
    Point const normal = rot90(tangent) * width / tangent.length();
    ret.append(LineSegment(this->controlPoint(0) + normal, this->controlPoint(1) + normal));

    return { ret };
}

// todo:
// * check drag tolerance, why always 3
// * why is stiching required?
template <>
Path BezierCurveN<2>::offsetted(double width, double tolerance, bool no_crossing) const
{
    return BezierCurve::offsettedInternal(width, tolerance, no_crossing);
}

template <>
Path BezierCurveN<3>::offsetted(double width, double tolerance, bool no_crossing) const
{
    return BezierCurve::offsettedInternal(width, tolerance, no_crossing);
}

Path BezierCurve::offsettedInternal(double width, double tolerance, bool no_crossing) const
{
    auto splitTimes = timesWithRadiusOfCurvature(width);
    splitTimes.push_back(1.);
    Path ret;
    ret.setStitching(true);
    double stopTime = 1.;
    Coord startTime = 0.;

    auto d = pointAndDerivatives(0, 2);
    double const rCurv0 = std::pow<double>(d[1].lengthSq(), 1.5) / (d[1][X] * d[2][Y] - d[2][X] * d[1][Y]);
    bool const reversed0 = rCurv0 > 0 ? rCurv0 < width : rCurv0 > width;

    d = pointAndDerivatives(1, 2);
    double const rCurv1 = std::pow<double>(d[1].lengthSq(), 1.5) / (d[1][X] * d[2][Y] - d[2][X] * d[1][Y]);
    bool const reversed1 = rCurv1 > 0 ? rCurv1 < width : rCurv1 > width;

    auto const x = inner[X];
    auto const y = inner[Y];
    auto const dx = Geom::derivative(x);
    auto const dy = Geom::derivative(y);
    auto const c1 = dx * dx + dy * dy;
    
    Point const point0 = pointAt(0);
    Point const tangent0 = unitTangentAt(0);
    Point const offset_point0 = point0 + rot90(tangent0) * width;
    std::cout << point0 << " " << tangent0 << " " << offset_point0 << std::endl;
    auto const time0 = nearestTime(offset_point0);
    std::cout << "time " << time0;
    double dist = (offset_point0 - pointAt(time0)).length();
    std::cout << " dist " << dist << " s " << width - tolerance / 10 << std::endl;
    if (dist < width - tolerance / 10) {
        auto const bx = x - x.valueAt(0);
        auto const by = y - y.valueAt(0);
        auto const c0 = bx * bx + by * by;
        auto const c2 = (by * dx - bx * dy) * width;
        auto poly = c0 * c0 * c1 - 4 * c2 * c2;

        auto const candidates = poly.roots();
        startTime = candidates[candidates.size() - 1];
        std::cout << "starttime set " << startTime << std::endl;
    }

    Point const point1 = pointAt(1);
    Point const tangent1 = unitTangentAt(1);
    Point const offset_point1 = point1 + rot90(tangent1) * width;
    std::cout << point1 << " " << tangent1 << " " << offset_point1 << std::endl;
    auto const time1 = nearestTime(offset_point1);
    std::cout << "time " << time1;
    double dist1 = (offset_point1 - pointAt(time1)).length();
    std::cout << " dist " << dist1 << " s " << width - tolerance / 10 << std::endl;
    if (dist1 < width - tolerance / 10) {
        auto const bx = x - x.valueAt(1);
        auto const by = y - y.valueAt(1);
        auto const c0 = bx * bx + by * by;
        auto const c2 = (by * dx - bx * dy) * width;
        auto poly = c0 * c0 * c1 - 4 * c2 * c2;

        auto const candidates = poly.roots();
        for (auto candidate : candidates) {
            std::cout << "stop cand " << candidate << std::endl;    
        }
        stopTime = candidates[0];
        std::cout << "stoptime set " << stopTime << std::endl;
    }

    bool last_segment_cutout = false;
    for (Coord time : splitTimes) {
        // ignore very small steps as it would produce a very small segment only
        if (time < startTime || are_near(time, startTime, 1e-3)) {
            continue;
        }
        if (time > stopTime) {
            time = stopTime;
        }
        if (startTime >= stopTime) {
            break;
        }

        // calculate the curvature in the middle
        // this is required to check if the radius of curvature is smaller or bigger as the width
        // as it decides in which direction the offset tangent should face.
        // The middle point is used as the start and end points might have the same radius of curvature
        // as the offset width which would lead to numerical issues
        auto const d = pointAndDerivatives(0.5 * (startTime + time), 2);
        // todo: denominator zero?
        double const radiusOfCurvature =
            std::pow<double>(d[1].lengthSq(), 1.5) / (d[1][X] * d[2][Y] - d[2][X] * d[1][Y]);
        // direction of the tangent on the offset curve wrt the original curve
        // as it is a precondition for offsettedSimple that this does not change betwenn start and end,
        // we just have to compute it once.
        bool const offset_tangent_reversed =
            radiusOfCurvature > 0 ? radiusOfCurvature < width : radiusOfCurvature > width;

        std::cout << startTime << " " << time << " " << offset_tangent_reversed << std::endl;

        if (no_crossing && offset_tangent_reversed) {
            last_segment_cutout = true;
            startTime = time;
            continue;
        }

        // todo start / endpoint check

        auto const curve_portion = static_cast<BezierCurve *>(portion(startTime, time));
        Path curve_offset_path = curve_portion->offsettedSimple(width, tolerance, offset_tangent_reversed);

        if (!no_crossing || !last_segment_cutout) {
            // the final point of ret is already very close to the first point of curve_offset_path
            // due to numerical inaccuracies, we manually assure that it is the same point.
            if (ret.size() > 0) {
                curve_offset_path.setInitial(ret.finalPoint());
            }
            ret.append(curve_offset_path);
        } else {
            // \todo: remove intersect as it requires too much time and we have a special case.
            auto const intersections = ret.intersect(curve_offset_path);
            if (intersections.size() == 0.) {
                std::cout << "todo no intersections" << std::endl;
                ret.append(curve_offset_path);
            } else {
                auto const intersection = intersections[0];
                std::cout << "intersections " << intersections.size() << " " << intersection.first.asFlatTime() << " "
                          << ret.size() << " " << intersection.second.asFlatTime() << " " << curve_offset_path.size()
                          << std::endl;
                ret = ret.portion(0., intersection.first.asFlatTime());
                std::cout << "i " << ret.initialPoint() << " f" << ret.finalPoint() << std::endl;
                auto bla = curve_offset_path.portion(intersection.second.asFlatTime(), curve_offset_path.size());
                std::cout << "i " << bla.initialPoint() << " f" << bla.finalPoint() << std::endl;
                bla.setInitial(ret.finalPoint());
                ret.append(bla);
            }
            std::cout << "yes" << std::endl;
        }

        startTime = time;
    }

    return ret;
}

Path BezierCurve::offsettedSimple(double width, double tolerance, bool offset_tangent_reversed, size_t levels) const
{
    size_t const sample_n = 10;
    Point fit_points[sample_n + 1];
    Point bezier_points[4];
    Path ret;

    for (size_t i = 0; i <= sample_n; i += 1) {
        Coord time = static_cast<double>(i) / sample_n;
        Point const point = pointAt(time);
        Point const tangent = unitTangentAt(time);
        fit_points[i] = point + rot90(tangent) * width;
    }

    Point offset_tangent0 = unitTangentAt(0);
    Point offset_tangent1 = -unitTangentAt(1.);
    if (offset_tangent_reversed) {
        offset_tangent0 = -offset_tangent0;
        offset_tangent1 = -offset_tangent1;
    }
    // the error between the points and the generated bezier curve does not matter that much when calling
    // bezier_fit_cubic. We are gonna check the error between the curves separately. todo: optimize fitting as this
    // might construct bad estimations sometimes
    bezier_fit_cubic_full(bezier_points, NULL, fit_points, sample_n + 1, offset_tangent0, offset_tangent1,
                          tolerance / 100, 1);
    CubicBezier bez(bezier_points[0], bezier_points[1], bezier_points[2], bezier_points[3]);

    std::vector<std::pair<double, double>> error_times;
    find_collinear_normal(error_times, controlPoints(), bez.controlPoints(), tolerance);
    double max_error = -1;
    Coord max_time = -1;
    for (auto times : error_times) {
        // todo: optimize split point as the find_collinear_normal is numerically bad condition because of the
        // offsetting operation (or other problem)
        Coord const time = times.first;
        Point const point1 = pointAt(time) + rot90(unitTangentAt(time)) * width;
        Point const point2 = bez.pointAt(bez.nearestTime(point1));
        double const error = (point2 - point1).length();
        if (error > max_error) {
            max_error = error;
            max_time = time;
        }
    }

    if (levels > 0 && max_error > tolerance && max_time > 0) {
        levels -= 1;
        auto const curve_portion0 = static_cast<BezierCurve *>(portion(0., max_time));
        auto const curve_offset_path0 =
            curve_portion0->offsettedSimple(width, tolerance, offset_tangent_reversed, levels);
        ret.append(curve_offset_path0);
        auto const curve_portion1 = static_cast<BezierCurve *>(portion(max_time, 1.));
        auto curve_offset_path1 = curve_portion1->offsettedSimple(width, tolerance, offset_tangent_reversed, levels);
        // the final point of curve_offset_path0 is already very close to the first point of curve_offset_path1
        // due to numerical inaccuracies, we manually assure that it is the same point.
        curve_offset_path1.setInitial(curve_offset_path0.finalPoint());
        ret.append(curve_offset_path1);
    } else {
        ret.append(bez);
    }

    return ret;
};

} // end namespace Geom

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
