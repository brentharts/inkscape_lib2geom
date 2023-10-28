/** @file
 * @brief Unit tests for bezier-curve offsetting.
 */

#include "testing.h"
#include <iostream>

#include <2geom/path.h>
#include <2geom/bezier-curve.h>

using namespace Geom;

TEST(BezierCurveOffsetTest, Linear) {
    auto test_linear = [](LineSegment const &bez, double amount, LineSegment const &expected, double tolerance =1e-4) {
        auto const offset_path = bez.offsetPointwise(amount);
        EXPECT_EQ(offset_path.size(), 1);
        auto const &actual = offset_path.front();
        EXPECT_EQ(actual.degreesOfFreedom(), 4);
        EXPECT_NEAR(actual.initialPoint()[X] - expected.initialPoint()[X], 0., 1e-6);
        EXPECT_NEAR(actual.initialPoint()[Y] - expected.initialPoint()[Y], 0., 1e-6);
        EXPECT_NEAR(actual.finalPoint()[X] - expected.finalPoint()[X], 0., 1e-6);
        EXPECT_NEAR(actual.finalPoint()[Y] - expected.finalPoint()[Y], 0., 1e-6);
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
    
    // tolerance should be irrelevant in the linear case
    double tolerance = 1000000;
    test_linear(
        LineSegment(Point(0, 1), Point(0, 2)),
        1,
        LineSegment(Point(-1, 1), Point(-1, 2)),
        tolerance
    );
}
