/** @file
 * @brief Unit tests for bezier-curve offsetting.
 */

#include "testing.h"
#include <iostream>

#include <2geom/path.h>
#include <2geom/bezier-curve.h>

using namespace Geom;

TEST(BezierCurveOffsetTest, Linear) {
    auto test_linear = [](LineSegment const &bez, double width, LineSegment const &expected) {
        auto const offset_path = bez.offset(width);
        EXPECT_TRUE(offset_path);
        EXPECT_EQ((*offset_path).size(), 1);
        auto const &actual = (*offset_path).front();
        EXPECT_EQ(actual.degreesOfFreedom(), 4);
        EXPECT_EQ(actual.initialPoint(), expected.initialPoint());
        EXPECT_EQ(actual.finalPoint(), expected.finalPoint());
    };

    double const l2 = 1/sqrt(2.);
    test_linear(
        LineSegment(Point(0, 1), Point(1, 2)),
        1.,
        LineSegment(Point(-l2, 1+l2), Point(1-l2, 2+l2))
    );
    test_linear(
        LineSegment(Point(0, 1), Point(1, 2)),
        -sqrt(2),
        LineSegment(Point(1, 0), Point(2, 1))
    );
    test_linear(
        LineSegment(Point(0, 1), Point(1, 2)),
        1e-10,
        LineSegment(Point(0, 1), Point(1, 2))
    );
    EXPECT_FALSE(LineSegment(Point(0, 1), Point(0, 1)).offset(1.));
}
