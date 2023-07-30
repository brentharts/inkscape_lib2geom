/** @file
 *  @brief Unit tests for conic sections: the xAx class.
 */
 /*
 * Authors:
 *   Rafał Siejakowski <rs@rs-math.net>
 *
 * Copyright 2023 Authors
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

#include "testing.h"

#include <glib.h>
#include <2geom/conicsec.h>

using namespace Geom;

namespace {
double random_coeff() { return g_random_double_range(-100, 100); }
}

TEST(ConicSectionsTest, RandomIntersectionsValid) {
    g_random_set_seed(0xCAFECAFE);
    double constexpr precision = 1e-6;
    for (unsigned _ = 0; _ < 10'000; ++_) {
        xAx conics[] = {
            { random_coeff(), random_coeff(), random_coeff(), random_coeff(), random_coeff(), random_coeff() },
            { random_coeff(), random_coeff(), random_coeff(), random_coeff(), random_coeff(), random_coeff() }
        };
        auto const intersections = intersect(conics[0], conics[1]);
        for (auto const &intersection : intersections) {
            for (unsigned which : {0, 1}) {
                EXPECT_LE(std::abs(conics[which].valueAt(intersection)), precision);
            }
        }
    }
}
