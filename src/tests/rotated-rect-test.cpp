/** @file
 * @brief Unit tests for RotatedRect
 * Uses the Google Testing Framework
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
#include <2geom/coord.h>
#include <2geom/rotated-rect.h>
#include <2geom/transforms.h>

using namespace Geom;

TEST(RotatedRectTest, midpoint) {
    Rect r(-0.5, -0.5, 5.5, 5.5);
    auto center = Point(2.5, 2.5);

    EXPECT_EQ(r.midpoint(), center);
    for (double angle: {0, 1, 25, 45, 90, 135}) {
        auto rotated_rect = RotatedRect::from_rect_rotate(r, Rotate::from_degrees(angle), Point(0, 0));
        auto rotated_center = center * Rotate(angle / 180.0 * M_PI);
        EXPECT_TRUE(Geom::are_near(rotated_rect.midpoint(), rotated_center, 1e-6)) << "Angle = " << angle << " deg";
    }
}

TEST(RotatedRectTest, containsPoint1) {
    Rect r(0, 0, 1, 1);
    auto rotated_rect = RotatedRect::from_rect_rotate(r, Rotate::from_degrees(0), Point(0, 0));
    EXPECT_TRUE(rotated_rect.contains(Point(0, 0)));
    EXPECT_TRUE(rotated_rect.contains(Point(1, 1)));
    EXPECT_TRUE(rotated_rect.contains(Point(0.5, 0.5)));
    EXPECT_FALSE(rotated_rect.contains(Point(1.1, 0.5)));
    EXPECT_FALSE(rotated_rect.contains(Point(0.5, 1.1)));
}

TEST(RotatedRectTest, containsPoint2) {
    Rect r(0, 0, 1, 1);
    auto rotated_rect = RotatedRect::from_rect_rotate(r, Rotate::from_degrees(45), Point(0, 0));
    EXPECT_TRUE(rotated_rect.contains(Point(0, 0)));
    EXPECT_TRUE(rotated_rect.contains(Point(0, 1.2)));
    EXPECT_TRUE(rotated_rect.contains(Point(0.5, 0.9)));
    EXPECT_FALSE(rotated_rect.contains(Point(1, 1)));
    EXPECT_FALSE(rotated_rect.contains(Point(0.1, 0)));
}

TEST(RotatedRectTest, intersects_aligned) {
    Rect r(0, 0, 1, 1);
    auto rotated_rect = RotatedRect::from_rect_rotate(r, Rotate::from_degrees(0));
// point within rect
    EXPECT_TRUE(rotated_rect.intersects(Rect(-1, -1, 2, 2)));
    EXPECT_TRUE(rotated_rect.intersects(Rect(0.1, 0.1, 0.2, 0.2)));
    EXPECT_TRUE(rotated_rect.intersects(Rect(-0.1, -0.1, 0.1, 0.1)));
    EXPECT_FALSE(rotated_rect.intersects(Rect(-0.2, -0.2, -0.1, -0.1)));
    EXPECT_FALSE(rotated_rect.intersects(Rect(1.1, 1.1, 1.2, 1.2)));
// edge intersection
    EXPECT_TRUE(rotated_rect.intersects(Rect(0.5, -0.1, 0.6, 1.2)));
    EXPECT_TRUE(rotated_rect.intersects(Rect(-0.1, 0.5, 1.2, 0.6)));
}

TEST(RotatedRectTest, bounds) {
    auto r = Rect::from_xywh(1.260, 0.547, 8.523, 11.932);
    auto rrect = RotatedRect::from_rect_rotate(r, Rotate::from_degrees(15.59));
    auto bbox = rrect.bounds();
    auto expected_bbox = Rect::from_xywh(-0.186, -0.378, 11.415, 13.783);
    for (int i = 0; i < 4; i++) {
        EXPECT_TRUE(Geom::are_near(bbox.corner(i), expected_bbox.corner(i), 1e-3));
    }
}


class RotatedRectTest : public testing::TestWithParam<std::tuple<Rect/*rect*/, double/*degrees*/, bool/*intersects*/>> {
    void SetUp() override {
        target = Rect::from_xywh(0, 0, 11, 13);
    }

public:
    Rect target;
};

TEST_P(RotatedRectTest, intersects) {
    Rect rect;
    double degrees;
    bool intersects;
    std::tie(rect, degrees, intersects) = GetParam();
    EXPECT_EQ(RotatedRect::from_rect_rotate(rect, Rotate::from_degrees(degrees)).intersects(target), intersects)
        << "ERROR: rect {" << rect << "} rotated by {" << degrees << "} degrees "
        << (!intersects ? "" : "NOT ") << "intersects with {" << target << "} but MUST "
        << (intersects ? "" : "NOT");
}

INSTANTIATE_TEST_CASE_P(intesect_non_aligned, RotatedRectTest, testing::Values(
        std::make_tuple(Rect::from_xywh(10.456, -4.479, 7, 5), 0, true),
        std::make_tuple(Rect::from_xywh(10.456, -4.479, 7, 5), 15, false),
        std::make_tuple(Rect::from_xywh(9.929, 12.313, 7, 5), 93.2, false),
        std::make_tuple(Rect::from_xywh(9.929, 12.313, 7, 5), 91.37, true),
        std::make_tuple(Rect::from_xywh(-1, 4, 13, 3), 0, true),
        std::make_tuple(Rect::from_xywh(4, -2, 3, 16), 0, true),
        std::make_tuple(Rect::from_xywh(-5.113, -3.283, 5.000, 7.000), 11.81, false),
        std::make_tuple(Rect::from_xywh(-5.113, -3.283, 5.000, 7.000), 13.35, true),
        std::make_tuple(Rect::from_xywh(1.260, 0.547, 8.523, 11.932), 15.59, true),
        std::make_tuple(Rect::from_xywh(5.328, 0.404, 11, 2), 28.16, true),
        std::make_tuple(Rect::from_xywh(4.853, 10.691, 11, 2), -30.4, true),
        std::make_tuple(Rect::from_xywh(-4.429, 10.752, 11, 2), 29.7, true),
        std::make_tuple(Rect::from_xywh(-4.538, 0.314, 11, 2), -34.19, true),
        std::make_tuple(Rect::from_xywh(8.398, -3.790, 2, 11), -34, true),
        std::make_tuple(Rect::from_xywh(8.614, 6.163, 2, 11), 30.38, true),
        std::make_tuple(Rect::from_xywh(0.492, 6.904, 2, 11), -37.29, true),
        std::make_tuple(Rect::from_xywh(0.202, -3.148, 2, 11), 31.12, true)
));

