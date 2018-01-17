#include <iostream>
#include "third_party/kdtree/kdtree.h"

//    A simple kdtree implementation to efficiently compute
//    nearest neighbors in small dimensions.
//    Contact: Ravi Kolluri
//             rkolluri@gmail.com
//
//    Copyright 2006 - Ravi Kolluri.
//
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
//    Commercial use is absolutely prohibited.
//
//    See the accompanying LICENSE file for details.

//  A simple test of the kdtree implementation.
//  TODO: write tests for nearest neighbors,
//  intersect sphere, and points in sphere.

int main(int argc, char **argv) {
  // insert four points
  // and test if the nearest neighbor function works
  double points[12] = {0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3};
  kdtree::kdtree<3, double> mykdtree;
  // don't create a new copy of the points
  mykdtree.initialize(points, 4, false);

  double testpoint[3] = {0.9, 0.9, 0.9};
  double* n1 = mykdtree.nearestNeighbor(testpoint, -1.0);
  if ((n1 - points) == 3) {
    std::cout  << "pass\n ";
  }
}
