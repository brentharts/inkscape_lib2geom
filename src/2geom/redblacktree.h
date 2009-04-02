/**
 * \file
 * \brief 
 * Implementation of Red-Black Tree as described in 
 * Intorduction to Algorithms. Cormen et al. Mc Grow Hill. 1990. pp 263-280
 * 
 * The intention is to implement interval trees mentioned in the same book, after the red-black.
 * Interval are heavily based on red-black trees (most operations are the same). So, we begin first 
 * with implementing red-black!
 *
 * Authors:
 *      ? <?@?.?>
 * 
 * Copyright 2009-2009 Evangelos Katsikaros
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
 *
 */

#include <vector>
#include <cassert>

#include <2geom/d2.h>

namespace Geom{

class RedBlack{
public:
    RedBlack *left, *right, *parent;
    bool isRed;
    double key; // This will change in the future 
    int data;
    // We'll use 2geom's interval for interval trees. Key wil be the min of the interval
    // and will also be augmented with the max of the subtree (more info will be added later)

    RedBlack(): left(0), right(0), parent(0), isRed(false), key(0.0), data(0) {}
};

class RedBlackTree{
public:
    RedBlack* root;

    RedBlackTree(): root(0) {}

    RedBlack* search(int shape);

    //void insert(RedBlack* z);
    void insert(Rect const &r, int shape);
    void insert(double x_min, int shape);

    void erase(RedBlack* T, int shape);

    void print_tree();
private:
    void left_rotate(RedBlack* x);
    void right_rotate(RedBlack* x);
    void tree_insert(RedBlack* x);
    void inorder_tree_walk(RedBlack* x);
};

}; //close namespace

/*
  Local Variables:
  mode:c++
  c-file-style:"stroustrup"
  c-file-offsets:((innamespace . 0)(inline-open . 0)(case-label . +))
  indent-tabs-mode:nil
  fill-column:99
  End:
*/
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:encoding=utf-8:textwidth=99 :