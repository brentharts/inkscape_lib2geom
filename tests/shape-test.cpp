/** @file
 * @brief Tests for the Shape class.
 */
/*
 * Authors:
 *   Rafał Siejakowski <rs@rs-math.net>
 *
 * Copyright 2022 Authors
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

#include <gtest/gtest.h>
#include <2geom/shape.h>
#include <2geom/svg-path-parser.h>
#include <2geom/svg-path-writer.h>

using namespace Geom;

#define PV(d) (parse_svg_path(d))
#define PTH(d) (PV(d)[0])

class ShapeTest : public testing::Test
{
protected:
    PathVector const _unit_square_positive, _unit_square_negative;
    PathVector const _open_square;
    PathVector const _circle1_positive, _circle1_negative, _circle2_positive, _circle2_negative;
    PathVector const _bowtie, _bowtie_curved, _bowtie_node;
    ShapeTest()
        : _unit_square_positive{PV("M 0 0 H 1 V 1 H 0 Z")}
        , _unit_square_negative{PV("M 0 0 V 1 H 1 V 0 Z")}
        , _open_square{PV("M 0 0 H 1 V 1 H 0")}
        , _circle1_positive{PV("M 1 0 A 1 1 0 1 1 -1 0 A 1 1 0 1 1 1 0 Z")}
        , _circle1_negative{PV("M 1 0 A 1 1 0 1 0 -1 0 A 1 1 0 1 0 1 0 Z")}
        , _circle2_positive{PV("M 2 0 A 2 2 0 1 1 -2 0 A 2 2 0 1 1 2 0 Z")}
        , _circle2_negative{PV("M 2 0 A 2 2 0 1 0 -2 0 A 2 2 0 1 0 2 0 Z")}
        , _bowtie{PV("M 0,0 L 4,2 V 0 L 0,2 Z")} // A polyline with a self-intersection.
        // A curved bow-tie path with a self-intersection @(10,5) between cubic Béziers.
        , _bowtie_curved{PV("M 0,0 V 10 C 10,10 10,0 20,0 V 10 C 10,10 10,0 0,0 Z")}
        // As above, but twice as large and the self-intersection @(20,10) happens at a node.
        , _bowtie_node{PV("M 0,0 V 20 C 0,20 10,20 20,10 25,5 30,0 40,0 V 20 "
                          "C 30,20 25,15 20,10 10,0 0,0 0,0 Z")}
    {}
};

/** Basic tests of the Shape::flatten function. */
TEST_F(ShapeTest, FlattenBasic)
{
    // Check that passing the CONVENTIONAL fill rule is results in a no-op.
    auto square = Shape::flatten(_unit_square_positive, FillRule::CONVENTIONAL, true, EPSILON);
    EXPECT_EQ(square, _unit_square_positive);

    // Ensure that this is the case even when we the path is negatively oriented.
    square = Shape::flatten(_unit_square_negative, FillRule::CONVENTIONAL, true, EPSILON);
    EXPECT_EQ(square, _unit_square_negative);

    // Ensure that CONVENTIONAL is a no-op even for an open path.
    auto dont_close = Shape::flatten(_open_square, FillRule::CONVENTIONAL, true, EPSILON);
    EXPECT_FALSE(dont_close[0].closed());

    // Check that a positively oriented boundary of the unit square is not altered with NONZERO.
    auto nice_square = Shape::flatten(_unit_square_positive, FillRule::NONZERO, true, EPSILON);
    EXPECT_EQ(nice_square, _unit_square_positive);

    // Check that the same is true for EVEN_ODD
    auto nice_square2 = Shape::flatten(_unit_square_positive, FillRule::EVEN_ODD, true, EPSILON);
    EXPECT_EQ(nice_square2, _unit_square_positive);

    // Check that a path consisting of 3 sides of the unit square is correctly closed when requested.
    auto auto_closed = Shape::flatten(_open_square, FillRule::NONZERO, true, EPSILON);
    EXPECT_EQ(auto_closed, _unit_square_positive);

    // Same with EVEN_ODD.
    auto auto_closed2 = Shape::flatten(_open_square, FillRule::EVEN_ODD, true, EPSILON);
    EXPECT_EQ(auto_closed2, _unit_square_positive);

    // Check that an open path is rejected on NONZERO when requested.
    auto reject = Shape::flatten(_open_square, FillRule::NONZERO, false, EPSILON);
    EXPECT_TRUE(reject.empty());

    // Same with EVEN_ODD
    auto reject2 = Shape::flatten(_open_square, FillRule::EVEN_ODD, false, EPSILON);
    EXPECT_TRUE(reject2.empty());
}

/** Test the correct reversal of negatively oriented contours. */
TEST_F(ShapeTest, FlattenReverse)
{
    // Ensure that a negatively oriented square contour is reversed
    auto reverse = Shape::flatten(_unit_square_negative, FillRule::NONZERO, true, EPSILON);
    // TODO: don't compare with a fixed PV since the order of paths is an implementation detail
    EXPECT_EQ(reverse, _unit_square_positive);

    // Same with the EVEN_ODD fill rule.
    auto reverse2 = Shape::flatten(_unit_square_negative, FillRule::EVEN_ODD, true, EPSILON);
    EXPECT_EQ(reverse2, _unit_square_positive);

    // Test with a negatively oriented bigon
    auto const bigon_negative = PV("M 0 0 Q 1 1 2 0 Q 1 -1 0 0 Z");
    auto const bigon_positive = PV("M 0 0 Q 1 -1 2 0 Q 1 1 0 0 Z");
    auto bigon_flat = Shape::flatten(bigon_negative, FillRule::NONZERO, true, EPSILON);
    auto &bigon = bigon_flat[0];
    auto &bigon_pos = bigon_positive[0];
    EXPECT_EQ(bigon.size(), bigon_pos.size());
    ASSERT_EQ(bigon.size(), 2);

    // Same but with EVEN_ODD
    auto bigon2_flat = Shape::flatten(bigon_negative, FillRule::EVEN_ODD, true, EPSILON);
    auto &bigon2 = bigon2_flat[0];
    EXPECT_EQ(bigon2.size(), bigon_pos.size());
    ASSERT_EQ(bigon2.size(), 2);
}

/** Test the area() member function. */
TEST_F(ShapeTest, ShapeArea)
{
    // Circle of radius 1 should have the area approximately pi
    auto circle_shape = Shape(_circle1_positive, FillRule::NONZERO);
    EXPECT_TRUE(are_near(circle_shape.area(), M_PI));

    // The same should be true when the shape is constructed from a negatively oriented path.
    auto circle_neg = Shape(_circle1_negative, FillRule::NONZERO);
    EXPECT_TRUE(are_near(circle_shape.area(), M_PI));

    // Unit square should have area 1
    auto square = Shape(_unit_square_positive, FillRule::EVEN_ODD);
    EXPECT_DOUBLE_EQ(square.area(), 1.0);

    // Check scaling by affine maps.
    auto const stretch_twice = Scale(2.0);
    EXPECT_DOUBLE_EQ((square * stretch_twice).area(), 4.0);

    // Check that orientation-reversing affine transformations don't make the area negative.
    auto const horiz_flip = Scale(-1.0, 1.0);
    EXPECT_TRUE(Affine(horiz_flip).flips());
    EXPECT_DOUBLE_EQ((square * horiz_flip).area(), 1.0);

    // Check that the area of holes is subtracted from the solid area.
    auto square_3x3 = _unit_square_positive * Scale(3) * Translate(-1, -1);
    square_3x3.insert(square_3x3.end(), _unit_square_negative[0]);
    auto square_3x3_with_1x1_hole = Shape(square_3x3, FillRule::NONZERO);
    EXPECT_DOUBLE_EQ(square_3x3_with_1x1_hole.area(), (3 * 3) - (1 * 1));
}

/** Test the Shape constuctor on a pair of concentric circles. */
TEST_F(ShapeTest, SimpleHole)
{
    // Nicely oriented annulus: outer contour oriented positively, hole negatively.
    auto annulus = _circle2_positive;
    annulus.insert(annulus.end(), _circle1_negative[0]);

    // Flattening with either NONZERO or EVEN_ODD should be a no-op
    for (auto fr : {FillRule::NONZERO, FillRule::EVEN_ODD}) {
        auto as_shape = Shape(annulus, fr);
        // We expect that the annulus shape has area (2*2 - 1*1)*pi = 3*pi
        EXPECT_TRUE(are_near(as_shape.area(), 3.0 * M_PI, 1e-5));
    }

    // Now try two positively oriented concentric circles
    auto positive_concentric = _circle2_positive;
    positive_concentric.insert(positive_concentric.end(), _circle1_positive[0]);

    // Flattening with NONZERO should leave us with just the larger circle
    auto positive_concentric_nonzero = Shape(positive_concentric, FillRule::NONZERO);
    EXPECT_TRUE(are_near(positive_concentric_nonzero.area(), 4.0 * M_PI, 1e-5));

    // Flattening with EVEN_ODD should give a shape with a hole, reorienting the hole's boundary.
    auto positive_concentric_even_odd = Shape(positive_concentric, FillRule::EVEN_ODD);
    EXPECT_TRUE(are_near(positive_concentric_even_odd.area(), 3.0 * M_PI, 1e-5));

    // Now try two negatively oriented concentric circles
    auto negative_concentric = _circle2_negative;
    negative_concentric.insert(negative_concentric.end(), _circle1_negative[0]);

    // Flattening with NONZERO should just give the larger circle, but reoriented.
    auto negative_concentric_nonzero = Shape(negative_concentric, FillRule::NONZERO);
    EXPECT_TRUE(are_near(negative_concentric_nonzero.area(), 4.0 * M_PI, 1e-5));

    // Flattening with EVEN_ODD should give a nice annulus
    auto negative_concentric_even_odd = Shape(negative_concentric, FillRule::EVEN_ODD);
    EXPECT_TRUE(are_near(negative_concentric_even_odd.area(), 3.0 * M_PI, 1e-5));

    // Lastly, check the case when the outer circle is negatively oriented and the inner one – positively.
    auto all_wrong = _circle2_negative;
    all_wrong.insert(all_wrong.end(), _circle1_positive[0]);

    // Flattening with NONZERO should give the nice annulus
    auto wrong_way_nonzero = Shape(all_wrong, FillRule::NONZERO);
    EXPECT_TRUE(are_near(wrong_way_nonzero.area(), 3.0 * M_PI, 1e-5));

    // Same for EVEN_ODD
    auto wrong_way_even_odd = Shape(all_wrong, FillRule::EVEN_ODD);
    EXPECT_TRUE(are_near(wrong_way_even_odd.area(), 3.0 * M_PI, 1e-5));
}


/** Test converting a Rect to a Shape. */
TEST_F(ShapeTest, ConstructFromRectangle)
{
    // Convert the unit square to a shape
    auto unit_square = Rect(Point(0, 0), Point(1, 1));
    auto as_shape = Shape(unit_square);
    EXPECT_EQ(as_shape.contour(), _unit_square_positive);

    // Expect that an empty rectangle converts to an empty shape
    auto empty_rect = Rect();
    auto empty_shape = Shape(empty_rect);
    EXPECT_TRUE(empty_shape.empty());

    // Same for a degenerate rectangle
    auto degen = Rect(Interval(42), Interval(1, 2));
    auto degenerate_shape = Shape(degen);
    EXPECT_TRUE(degenerate_shape.empty());

    // Expect that trying to convert an unbounded rectangle throws an exception.
    auto first_quadrant = Rect(Interval(0, infinity()), Interval(0, infinity()));
    EXPECT_ANY_THROW(auto unbounded = Shape(first_quadrant));
}

/** Test constructing a Shape from a Parallelogram. */
TEST_F(ShapeTest, ConstructFromParallelogram)
{
    // Convert a parallelogram to a shape.
    auto my_affine = Scale(2) * Translate(-1, 3) * Rotate(0.25 * M_PI);
    auto para1 = Parallelogram(Rect(Point(0, 0), Point(1, 1)));
    para1 *= my_affine;
    auto as_shape = Shape(para1);
    EXPECT_EQ(as_shape.boundsFast(), para1.bounds());
    EXPECT_DOUBLE_EQ(as_shape.area(), 4);

    // Now try with an orientation-reversing affine transformation
    auto para2 = para1 * Scale(-1, 0.5);
    auto shape_from_flipped = Shape(para2);
    EXPECT_DOUBLE_EQ(shape_from_flipped.area(), 2);

    // Try constructing a shape from a degenerate parallelogram
    auto para3 = Parallelogram(Rect(Point(0, 0), Point(3, 0)));
    auto degen = Shape(para3);
    EXPECT_TRUE(degen.empty());
}

/** Test constructing a Shape from a circle. */
TEST_F(ShapeTest, ConstructFromCircle)
{
    auto unit_circle = Circle(Point(0, 0), 1);
    auto as_shape = Shape(unit_circle);
    EXPECT_TRUE(are_near(as_shape.area(), M_PI));

    auto neg_circle = Circle(Point(0, 0), -1);
    auto neg_shape = Shape(neg_circle);
    EXPECT_TRUE(are_near(neg_shape.area(), M_PI));

    auto degen_circle = Circle(Point(42, 69), 0);
    auto degen_shape = Shape(degen_circle);
    EXPECT_TRUE(degen_shape.empty());
}

/** Test constructing a Shape from an ellipse. */
TEST_F(ShapeTest, ConstructFromEllipse)
{
    auto const ellipse = Ellipse(Point(1, 2), Point(3, 4), M_PI * 0.25);
    auto ellipse_shape1 = Shape(ellipse);
    // Area of an ellipse is π * a * b where a and b are the rays.
    auto const precision = 1e-4; // Precision is a bit junk on this one
    EXPECT_TRUE(are_near(ellipse_shape1.area(), 12 * M_PI, precision));

    // Check flipping the ellipse by an affine with a negative determinant
    auto const fllipse = Ellipse(Circle(Point(0, 0), 1)) * Scale(1, -1);
    auto shape_from_flipped = Shape(fllipse);
    EXPECT_TRUE(are_near(shape_from_flipped.area(), M_PI));

    // Check that a degenerate ellipse produces an empty shape
    auto degen = Ellipse(Point(42, 0), Point(3, 0), 0.0);
    auto degen_shape = Shape(degen);
    EXPECT_TRUE(degen_shape.empty());
}

/** Test construction (hence, `flatten()`) from a PathVector with grazing rectangles. */
TEST_F(ShapeTest, GrazingRectangles)
{
    // A 1x1 square whose side is the middle third of the side of a 3x3 square.
    auto const left_square = PV("M 0,0 V 3 H -3 V 0 Z");
    auto const right_square = PTH("M 1,1 V 2 H 0 V 1 Z");
    PathVector combined{left_square};
    combined.push_back(right_square);
    auto const expected_contour = PV("M 0 1 H 1 V 2 H 0 V 3 H -3 V 0 H 0 Z");

    // Construct a shape consisting of two square contours which touch
    auto test = Shape(combined, FillRule::NONZERO);
    EXPECT_EQ(test.contour(), expected_contour);

    // Two rectangles with complete overlap of edges; one oriented positively,
    // the other negatively.
    PathVector complete_overlap{left_square};
    complete_overlap.push_back(PTH("M 0,0 V 3 H 42 V 0 Z"));
    auto const expected_bbox = Rect(Interval(-3, 42), Interval(0, 3));

    for (auto fr : {FillRule::NONZERO, FillRule::EVEN_ODD}) {
        auto long_rect = Shape(complete_overlap, fr);
        EXPECT_EQ(long_rect.numBoundaryComponents(), 1);
        auto bbox = long_rect.boundsExact();
        ASSERT_TRUE(bbox);
        EXPECT_EQ(*bbox, expected_bbox);
    }

    // Two staggered rectangles with partial-partial ovelap (brick wall).
    auto const brick1 = PTH("M 0,0 H 2 V 1 H 0 Z"); // area 2
    auto const brick2 = PTH("M 1,1 H 3 V 2 H 1 Z"); // area 2
    PathVector wall;
    wall.push_back(brick1);
    wall.push_back(brick2);
    auto tetromino = Shape(wall, FillRule::NONZERO);
    EXPECT_EQ(tetromino.numBoundaryComponents(), 1);
    EXPECT_EQ(tetromino.area(), 4);

    // Two rectangles with a common external vertex (L shape)
    auto const riser = PTH("M 0,1 H 1 V 3 H 0 Z"); // area 2
    PathVector letter_L;
    letter_L.push_back(brick1);
    letter_L.push_back(riser);
    auto L_shape = Shape(letter_L, FillRule::EVEN_ODD);
    EXPECT_EQ(L_shape.area(), 4);
    EXPECT_EQ(L_shape.numBoundaryComponents(), 1);
}

/** Test contour overlap, but on the inside. */
TEST_F(ShapeTest, InternalOverlap)
{
    // Two positive contours with internal overlap
    auto const inside_overlap = PV("M 0,0 H 2 V 3 H 0 Z"
                                   "M 1,1 H 2 V 2 H 1 Z");
    // We expect to get just the larger rectangle when using FillRule::NONZERO
    auto flattened_rect = Shape(inside_overlap, FillRule::NONZERO);
    auto const expected_bounds = Rect(Interval(0, 2), Interval(0, 3));
    EXPECT_EQ(flattened_rect.numBoundaryComponents(), 1);
    auto bounds = flattened_rect.boundsExact();
    ASSERT_TRUE(bounds);
    EXPECT_EQ(*bounds, expected_bounds);
    EXPECT_DOUBLE_EQ(flattened_rect.area(), 6);

    // The same PathVector with EVEN_ODD fill rule should flatten out to a C-shape.
    auto c_shape = Shape(inside_overlap, FillRule::EVEN_ODD);
    EXPECT_EQ(c_shape.numBoundaryComponents(), 1);
    // The c_shape should have the same bounds (the notch doesn't change the bounding rect).
    bounds = c_shape.boundsExact();
    ASSERT_TRUE(bounds);
    EXPECT_EQ(*bounds, expected_bounds);
    // However, the notch reduces the area from 6 to 5:
    EXPECT_DOUBLE_EQ(c_shape.area(), 5);

    // Now for a complete internal overlap: one rectangle occupying the top half of the other.
    auto const half_and_whole = PV("M 0,0 H 1 V 2 H 0 Z"
                                   "M 0,1 H 1 V 2 H 0 Z");
    auto flattened_total = Shape(half_and_whole, FillRule::NONZERO);
    // We expect to get the whole thing (due to NONZERO)
    EXPECT_EQ(flattened_total.numBoundaryComponents(), 1);
    bounds = flattened_total.boundsExact();
    ASSERT_TRUE(bounds);
    EXPECT_EQ(*bounds, Rect(Interval(0, 1), Interval(0, 2)));
    EXPECT_DOUBLE_EQ(flattened_total.area(), 2);

    // With EVEN_ODD, we expect that the top half will be cut out and only the unit square is left in.
    auto flattened_half = Shape(half_and_whole, FillRule::EVEN_ODD);
    EXPECT_EQ(flattened_half.numBoundaryComponents(), 1);
    bounds = flattened_half.boundsExact();
    ASSERT_TRUE(bounds);
    EXPECT_EQ(*bounds, Rect(Interval(0, 1), Interval(0, 1)));
    EXPECT_DOUBLE_EQ(flattened_half.area(), 1);
    EXPECT_EQ(flattened_half.contour().at(0).size(), 4); // the path should have 4 segments

    // Now for a staggered internal overlap.
    auto const staggered_internal = PV("M 0,0 H 2 V 2 H 0 Z" // 2x2 square
                                       "M 1,1 H 3 V 2 H 1 Z"); // 2x1 rectangle
    // With NONZERO, expect a pentomino like this:
    // XXX
    // XX
    auto flattened_extension = Shape(staggered_internal, FillRule::NONZERO);
    EXPECT_EQ(flattened_extension.numBoundaryComponents(), 1);
    bounds = flattened_extension.boundsExact();
    ASSERT_TRUE(bounds);
    EXPECT_EQ(*bounds, Rect(Interval(0, 3), Interval(0, 2)));
    EXPECT_DOUBLE_EQ(flattened_extension.area(), 5);

    // With EVEN_ODD, we expect that the overlap will be cut out:
    // X X
    // XX
    // Due to the the corner adjacency, we allow the contour to
    // consist of either 1 or 2 paths.
    auto flattened_carved = Shape(staggered_internal, FillRule::EVEN_ODD);
    EXPECT_GE(flattened_carved.numBoundaryComponents(), 1);
    EXPECT_LE(flattened_carved.numBoundaryComponents(), 2);
    bounds =  flattened_carved.boundsExact();
    ASSERT_TRUE(bounds);
    EXPECT_EQ(*bounds, Rect(Interval(0, 3), Interval(0, 2)));
    EXPECT_DOUBLE_EQ(flattened_carved.area(), 4);
}

/** Test flattening out several copies of the same contour. */
TEST_F(ShapeTest, TotalOverlap)
{
    // Two identical squares on top of one another!
    auto stacked2 = _unit_square_positive;
    stacked2.push_back(_unit_square_positive[0]);
    auto sq2 = Shape(stacked2, FillRule::NONZERO);
    EXPECT_EQ(sq2.numBoundaryComponents(), 1);
    EXPECT_EQ(sq2.contour().at(0).size(), 4);
    EXPECT_EQ(sq2.area(), 1);

    // With EVEN_ODD, the identical squares should "cancel out" giving an empty shape
    auto cancelled = Shape(stacked2, FillRule::EVEN_ODD);
    EXPECT_TRUE(cancelled.empty());

    // Try with circles for a change
    PathVector circles = Shape(Circle(Point(0, 0), 1));
    circles.push_back(circles[0]);
    auto from_two_circles = Shape(circles, FillRule::NONZERO);
    EXPECT_EQ(from_two_circles.numBoundaryComponents(), 1);
    EXPECT_TRUE(are_near(from_two_circles.area(), M_PI));

    // With EVEN_ODD, the circles should cancel out
    auto nothing = Shape(circles, FillRule::EVEN_ODD);
    EXPECT_TRUE(nothing.empty());

    // And now for 3 identical ellipses:
    PathVector const ellipse = Shape(Ellipse(Point(3, 4), Point(2, 7), 0.8));
    PathVector three_ellipses = ellipse;
    three_ellipses.push_back(ellipse[0]);
    three_ellipses.push_back(ellipse[0]);
    EXPECT_EQ(three_ellipses.size(), 3);

    for (auto fr : {FillRule::NONZERO, FillRule::EVEN_ODD}) {
        auto as_shape = Shape(three_ellipses, fr);
        EXPECT_EQ(as_shape.numBoundaryComponents(), 1);
        EXPECT_EQ(as_shape.contour().curveCount(), ellipse.curveCount());
        EXPECT_TRUE(are_near(as_shape.area(), 14 * M_PI, 1e-4));
    }
}

/** Test figure-eight (infinity symbol, bowtie) paths with self-intersections. */
TEST_F(ShapeTest, Bowties)
{
    for (auto fr : {FillRule::NONZERO, FillRule::EVEN_ODD}) {
        auto bowtie = Shape(_bowtie, fr);
        EXPECT_GE(bowtie.numBoundaryComponents(), 1);
        EXPECT_LE(bowtie.numBoundaryComponents(), 2);
        auto const bbox = bowtie.boundsExact();
        ASSERT_TRUE(bbox);
        EXPECT_EQ(*bbox, Rect(Point(0, 0), Point(4, 2)));
        EXPECT_EQ(bowtie.area(), 4);
    }

    // Meeting paths are Bézier curves
    for (auto fr : {FillRule::NONZERO, FillRule::EVEN_ODD}) {
        auto curved = Shape(_bowtie_curved, fr);
        EXPECT_GE(curved.numBoundaryComponents(), 1);
        EXPECT_LE(curved.numBoundaryComponents(), 2);
        auto const bbox = curved.boundsExact();
        ASSERT_TRUE(bbox);
        EXPECT_EQ(*bbox, Rect(Point(0, 0), Point(20, 10)));
        EXPECT_EQ(curved.area(), 137.5);
    }

    // Bézier bowtie but with a node at the intersection point
    for (auto fr : {FillRule::NONZERO, FillRule::EVEN_ODD}) {
        auto bownode = Shape(_bowtie_node, fr);
        EXPECT_GE(bownode.numBoundaryComponents(), 1);
        EXPECT_LE(bownode.numBoundaryComponents(), 2);
        auto const bbox = bownode.boundsExact();
        ASSERT_TRUE(bbox);
        EXPECT_EQ(*bbox, Rect(Point(0, 0), Point(40, 20)));
    }

    // Bowtie made of elliptical arcs tangent at the intersection!
    auto const bowtie_elliptical = PV("M -2,1 A 2,1 0 0 0 0,0 A 2,1 0 0 1 2,-1 V 1 A 2,1 0 0 1 0,0 A 2,1 0 0 0 -2,-1 V 1 Z");
    for (auto fr : {FillRule::NONZERO, FillRule::EVEN_ODD}) {
        auto bowtie_tangent = Shape(bowtie_elliptical, fr);
        EXPECT_GE(bowtie_tangent.numBoundaryComponents(), 1);
        EXPECT_LE(bowtie_tangent.numBoundaryComponents(), 2);
        auto const bbox = bowtie_tangent.boundsExact();
        ASSERT_TRUE(bbox);
        EXPECT_EQ(*bbox, Rect(Point(-2,-1), Point(2,1)));
        EXPECT_TRUE(are_near(bowtie_tangent.area(), 2 * M_PI, 1e-6));
    }
}

/** Stress-test the vertex resolution mechanism. */
TEST_F(ShapeTest, Pincers)
{
    auto const pincers = PV("M 4,0 L 2,-1 V 1 Z" // Triangular hole
                            "M 4,0 L 0,2.01 V -2.01 Z"); // In a larger triangle.
    auto as_shape = Shape(pincers, FillRule::NONZERO);
    double const expected_area = 0.5 * 4 * 4.02 - 0.5 * 2 * 2;
    EXPECT_GE(as_shape.numBoundaryComponents(), 1);
    EXPECT_LE(as_shape.numBoundaryComponents(), 2);
    EXPECT_EQ(as_shape.area(), expected_area);

    // Two tangent Béziers at the rightmost point (witness at a vertex).
    auto const D_in_a_D = PV("M 1,0 C 1,1 1,2 0,2 V -2 C 1,-2 1,-1 1,0 Z" // Letter D
                             "M 1,0 C 1,2 1,3 0,3 V -3 C 1,-3. 1,-2 1,0 Z"); // Taller letter D
    auto D_xor_D = Shape(D_in_a_D, FillRule::EVEN_ODD);
    auto const area = D_xor_D.area();
    EXPECT_GT(area, 1.8);
    EXPECT_LT(area, 2.0);

    // Ellipse with an elliptical cutout with internally tangent contour
    auto const chicken_egg = PV("M 30,0 A 15,10 0 1 1 0,0 A 15,10 0 1 1 30,0 Z" // large ellipse
                                "M 30,0 A 5,2 0 1 1 20,0 A 5,2 0 1 1 30,0 Z"); // smaller internally tangent ellipse.
    double const cutout_expected_area = M_PI * (15 * 10 - 5 * 2);
    auto cutout = Shape(chicken_egg, FillRule::EVEN_ODD);
    EXPECT_TRUE(are_near(cutout.area(), cutout_expected_area, 0.005));
}

TEST_F(ShapeTest, Teardrops)
{
    auto const large_drop = PV("M 0,0 C -4,2 -4,-2 0,0 Z");
    auto large_drop_shape = Shape(large_drop, FillRule::NONZERO);
    double const expected_large_drop_area = large_drop_shape.area();

    auto const small_drop_shape = large_drop_shape * Scale(0.5);
    double const expected_small_drop_area = 0.25 * expected_large_drop_area;
    EXPECT_DOUBLE_EQ(expected_small_drop_area, small_drop_shape.area());

    // Two internally tangent teardrops, both positively oriented.
    auto const drop_in_drop = PV("M 0,0 C -2,1 -2,-1 0,0 Z"
                                 "M 0,0 C -4,2 -4,-2 0,0 Z");

    // Flattening with NONZERO should discard the smaller teardrop.
    auto ignore_smaller = Shape(drop_in_drop, FillRule::NONZERO);
    EXPECT_EQ(ignore_smaller.numBoundaryComponents(), 1);
    EXPECT_TRUE(are_near(ignore_smaller.area(), expected_large_drop_area, 1e-4));

    // Flattening with EVEN_ODD should stamp a hole with the smaller teardrop.
    auto with_hole = Shape(drop_in_drop, FillRule::EVEN_ODD);
    EXPECT_GE(with_hole.numBoundaryComponents(), 1);
    EXPECT_TRUE(are_near(with_hole.area(), expected_large_drop_area - expected_small_drop_area, 1e-4));

    // Test a teardrop intersecting a 2x2 square centered around the origin.
    auto square_with_tail = large_drop;
    square_with_tail.push_back(PTH("M -1,-1 h 2 v 2 h -2 Z"));
    auto square_tail_shape = Shape(square_with_tail, FillRule::NONZERO);
    EXPECT_EQ(square_tail_shape.numBoundaryComponents(), 1);
    EXPECT_GT(square_tail_shape.area(), 5.5);

    // Cancel out the intersection
    auto square_minus_tear = Shape(square_with_tail, FillRule::EVEN_ODD);
    EXPECT_GE(square_minus_tear.numBoundaryComponents(), 1);
    EXPECT_LE(square_minus_tear.numBoundaryComponents(), 2);
    EXPECT_LE(square_minus_tear.area(), 5.6);

    // Test a teardrop attached to the corner of a 1x1 square
    auto square_with_ear = large_drop;
    square_with_ear.push_back(PTH("M 0,0 H 1 V 1 H 0 Z"));
    auto square_ear_shape = Shape(square_with_ear, FillRule::NONZERO);
    EXPECT_GE(square_ear_shape.numBoundaryComponents(), 1);
    EXPECT_LE(square_ear_shape.numBoundaryComponents(), 2);
    EXPECT_DOUBLE_EQ(square_ear_shape.area(), 1 + expected_large_drop_area);
}
// TODO: More tests for weird and unusual scenarios.

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
