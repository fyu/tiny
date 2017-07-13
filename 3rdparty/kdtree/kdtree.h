#ifndef THIRD_PARTY_KDTREE_KDTREE_H__
#define THIRD_PARTY_KDTREE_KDTREE_H__

//
//    A kdtree implementation.
//
//    Contact: Ravi Kolluri
//              rkolluri@gmail.com

//     Copyright 2006 - Ravi Kolluri.


//       This program is free software; you can redistribute it and/or modify
//       it under the terms of the GNU General Public License as published by
//       the Free Software Foundation; either version 2 of the License, or
//       (at your option) any later version.

//       This program is distributed in the hope that it will be useful,
//       but WITHOUT ANY WARRANTY; without even the implied warranty of
//       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//       GNU General Public License for more details.

//       You should have received a copy of the GNU General Public License
//       along with this program; if not, write to the Free Software
//       Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

//     See the accompanying LICENSE file for details.


//  A templated implementation of kdtrees.
//  The template parameter specifies the
//  dimension in which the tree is used.

//  To build the tree call

//  initialize(myreal* points, int numpoints, bool copypoints);

//  The kdtree supports four different functions to lookup
//  nearest neighbors over the input point cloud.

//  1) Returns  the nearest neighbor of the given point.
//     If minDistanceSqr is greater than 0, the nearest
//     neighbor search occurs inside a sphere of squared radius
//     minDistanceSqr.
//
//    myreal* nearestNeighbor(myreal *srcPoint,myreal minDistanceSqr) {
//
//  2) Computes (at most) numNeighbors nearest neighbors of the given point
//     and stores them in heap. Returns the number of neighbors stored
//     in the heap. If minDistanceSqr is greater than 0, the nearest
//     neighbor search occurs inside a sphere of squared radius
//     minDistanceSqr.
//
//   int nearestNeighbors(myreal* srcPoint,myreal **heap,int numNeighbors,
//                        myreal minDistanceSqr);
//
//  3) Calls the function *fptr for each point inside a sphere
//     of specified radius. Use dataptr to (optionally) pass something on to
//     the callback.
//
//      void intersectSphere(myreal *srcPoint,myreal sqrradius,void* dataptr,
//                           bool  (*fptr)(int,void*));
//
//  4) Returns the number of points inside the specified sphere.
//     int pointsInSphere(myreal *srcPoint,myreal sqrradius);


#include "stack.h"

namespace kdtree {

typedef unsigned short     uint16;

//  The template parameter is for the dimension
//  of the input point cloud. Not very useful
//  in high dimensions.
template <int DIM, typename myreal>
class kdtree {
 public:

  kdtree() :
      coords(NULL),
      numPoints(0),
      kdheap(NULL),
      isinitialized(false)
  {}

  //  Function to initialize the kdtree.
  //  copyPoints should be set to true
  //  if you want the kdtree to keep
  //  its own copy of the point cloud.
  //  Set to false to save memory.
  void initialize(myreal *pointCoords,
                  int nPoints, bool copyPoints = true);

  //  find the nearest neighbor of the given point
  myreal* nearestNeighbor(myreal *srcPoint, myreal minDistanceSqr) const;

  //  Find at most numNeighbors neighbors inside
  //  the given sphere.
  //  Returns the number of neighbors
  int nearestNeighbors(myreal* srcPoint, int *heap, int numNeighbors,
                       myreal minDistanceSqr) const;

  //  For each point inside the given sphere, calls fptr.
  void intersectSphere(myreal *srcPoint, myreal sqrradius, void* dataptr,
                       bool  (*fptr)(int, void*)) const;

  //  count the number of points inside the given sphere
  int pointsInSphere(myreal *srcPoint, myreal sqrradius) const;

  ~kdtree() {
    clear();
  }

  //  destroy the heap..
  void clear() {
    if (coords != NULL && pointsStored)
      delete[] coords;
    if (kdheap != NULL)
      delete[] kdheap;
  }

  int pointIndex(myreal *vec) const {
    return (vec-coords)/DIM;
  }

  myreal* getCoords(int index) const {
    return coords+DIM*index;
  }

 private:
  // Size of the stack used for kdtree traversal
  // (at most the depth of the tree).
  static const int kStackSize = 64;

  //  Class definition for the kdtree node.
  class treenode {
   public:
    myreal *pos;
    int plane;  //  the axis chosen for the split..
    treenode() :
        pos(NULL),
        plane(-1)
    {}
  };
  myreal *coords;    //  the coordinates
  int numPoints;     //  number of points
  treenode *kdheap;  //  the kdtree stored as a heap..
  bool isinitialized;

  bool pointsStored;  //  do we own the point array?
  myreal bboxMin[DIM], bboxMax[DIM];  //  bounding box

  myreal DISTANCE(const myreal *vec1,
                  const myreal *vec2) const {
    myreal dist = 0.0;
    for (int i = 0; i < DIM; i++) {
      dist += SQ(vec1[i] - vec2[i]);
    }
    return dist;
  }
  myreal SQ(myreal a) const {
    return a * a;
  }

  myreal Absolute(myreal a) const {
    return ((a) >= 0.0 ? (a) : -(a));
  }

  void swap(treenode** arr, int i, int j) {
    treenode* temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
  }

  //  build a left balanced kdtree and store it as a heap..
  //  should this be a separate function ?
  void balance() {
    //  allocate temporary storage...
    treenode **arr1= new treenode*[numPoints+1];
    treenode **arr2= new treenode*[numPoints+1];
    for (int i = 1;i <= numPoints; i++)
      arr2[i]=kdheap+i;
    balanceSegment(arr1, arr2, 1, 1, numPoints);
    delete[] arr2;

    //  now store the heap...
    int d = 0, j = 1, foo = 1;
    treenode tempNode = kdheap[j];
    for (int i = 1;i <= numPoints; i++) {
      d = arr1[j]-kdheap;
      arr1[j] = NULL;
      if (d != foo) {
        kdheap[j] = kdheap[d];
      } else {
        kdheap[j] = tempNode;
        if (i < numPoints) {
          for (; foo <= numPoints;foo++)
            if (arr1[foo] != NULL)
              break;
          tempNode = kdheap[foo];
          j = foo;
        }
        continue;
      }
      j = d;
    }
    delete[] arr1;
  }


  //  axis denotes the coordinate plane that  is used to do the split..
  void medianSplit(treenode **arr, int start,
                   int end, int median, const int axis) {
    int left = start;
    int right = end;
    //  the quickSelect algorithm for computing the median..
    while (left < right) {
      //  choose the partitioning  element..
      const myreal v = arr[right]->pos[axis];
      int i = left - 1;
      int j = right;
      for (;;) {
        while (arr[++i]->pos[axis] < v);  //  move to the right..
        while (arr[--j]->pos[axis] > v &&  j>left);

        if (i >= j)  //  do we need to swap?
          break;
        swap(arr, i, j);
      }
      swap(arr, i, right);
      if (i >= median)
        right = i - 1;
      if (i <= median)
        left = i + 1;
    }
  }

  void chooseAxis(treenode **arr2, const int start,
                  const int end, int *axis) {
    myreal means[DIM];
    myreal var[DIM];

    //  set the mean and the variance
    for (int dimension = 0; dimension < DIM ; dimension++) {
      means[dimension] = 0.0;
      var[dimension] = 0.0;
    }

    for (int point = start; point <= end; point++) {
      for (int dimension = 0; dimension < DIM; dimension++)
        means[dimension] += arr2[point]->pos[dimension];
    }

    for (int dimension = 0; dimension < DIM; dimension++)
      means[dimension] /= (end - start + 1);

    for (int point = start; point <= end; point++) {
      for (int dimension = 0; dimension < DIM; dimension++) {
        var[dimension] += SQ(arr2[point]->pos[dimension] - means[dimension]);
      }
    }

    //  get the axis..
    *axis = 0;
    myreal maxvar = var[0];
    for (int dimension = 1; dimension < DIM; dimension++) {
      if (var[dimension] > maxvar) {
        maxvar = var[dimension];
        *axis = dimension;
      }
    }
  }

  //  balance a portion of the array
  void balanceSegment(treenode **arr1, treenode **arr2, const int index,
                      const int start, const int end) {
    //  a recursive function , maybe make this iterative later..

    // find the median first...
    int median = 1;

    while (4 * median <= (end - start + 1))
      median += median;

    if ((3 * median) <= (end - start + 1)) {
      median += median;
      median += start-1;
    } else {
      median = end - median + 1;
    }
    // find the axis with the greatest width...
    int axis = 0;
    chooseAxis(arr2, start, end, &axis);
    // we now have the axis and the median.. now split..
    medianSplit(arr2, start, end, median, axis);
    //  add an element to the heap..
    arr1[index] = arr2[median];
    arr1[index]->plane = axis;
    // now recursively balance the children..
    if (median > start) {
      // balance the left child..
      if (start < median-1) {
        const myreal tmp = bboxMax[axis];
        bboxMax[axis] = arr1[index]->pos[axis];
        balanceSegment(arr1, arr2, 2 * index, start, median - 1);
        bboxMax[axis] = tmp;
      } else {
        //  no recursive split..
        arr1[2 * index] = arr2[start];
      }
    }
    if (median < end) {
      // balance the left child..
      if ((median+1) < end) {
        const myreal tmp = bboxMin[axis];
        bboxMin[axis]= arr1[index]->pos[axis];
        balanceSegment(arr1, arr2, 2 * index + 1, median+1, end);
        bboxMin[axis]=tmp;
      } else {
        arr1[2*index+1]= arr2[end];
      }
    }
  }
};

//  Function to initialize the kdtree.
//  copyPoints should be set to true
//  if you want the kdtree to keep
//  its own copy of the point cloud.
//  Set to false to save memory.

template<int DIM, typename myreal>
void kdtree<DIM, myreal>::initialize(myreal *pointCoords,
                                     int nPoints, bool copyPoints) {
  isinitialized = true;
  this->pointsStored = copyPoints;
  this->numPoints = nPoints;
  if (copyPoints) {
    this->coords = new myreal[DIM*nPoints];
    for (int i = 0; i < DIM * nPoints; i++) {
      coords[i] = pointCoords[i];
    }
  } else {
    this->coords = pointCoords;
  }
  this->kdheap= new treenode[nPoints+1];

  for (int i = 1; i < numPoints + 1; i++) {
    kdheap[i].pos = coords+DIM*(i-1);
  }
  // set the bounding box..

  for (int i = 0;i < DIM; i++) {
    bboxMin[i] = pointCoords[i];
    bboxMax[i] = pointCoords[i];
  }

  for (int i = 0; i < nPoints; i++) {
    for (int j = 0;j < DIM; j++) {
      bboxMin[j] = (bboxMin[j] < pointCoords[i*DIM+j]) ?
                   bboxMin[j]:pointCoords[i*DIM+j];
      bboxMax[j] = (bboxMax[j] > pointCoords[i*DIM+j])?
                   bboxMax[j]:pointCoords[i*DIM+j];
    }
  }
  //  build a balanced tree
  balance();
}

template<int DIM, typename myreal>
myreal* kdtree<DIM, myreal>::nearestNeighbor(myreal *srcPoint,
                                             myreal minDistanceSqr) const {
  if (!isinitialized) {
    return NULL;  //  not initialized, quick return
  }

  if (minDistanceSqr <= 0)
    minDistanceSqr = DISTANCE(srcPoint, kdheap[1].pos) + 0.01;

  uint16 UNVISITED = 0, LEFT = 1, RIGHT = 2;
  treenode *minNode = NULL;
  int index = 1;

  Stack<treenode*, 64> nodeStack;
  nodeStack.push(kdheap + index);
  Stack<uint16, 64> stateStack;
  stateStack.push(UNVISITED);
  int visitednodes = 0;

  while (!(nodeStack.empty())) {
    treenode* currNode = nodeStack.top();
    nodeStack.pop();
    uint16 state = stateStack.top();
    stateStack.pop();
    if (state == UNVISITED) {
      // try to see if the left child can be visited...
      visitednodes++;
      int index = currNode-kdheap;
      if ((2*index) <= numPoints) {
        // push this back on to the stack...
        nodeStack.push(currNode);
        const myreal planeDistance = srcPoint[currNode->plane] -
            currNode->pos[currNode->plane];
        // now check which child needs to be pushed on to the stack..
        if (planeDistance <= 0 || (2*index+1) > numPoints) {
          stateStack.push(LEFT);
          // visit the left child first...
          nodeStack.push(kdheap + 2*index);
          stateStack.push(UNVISITED);
        } else {
          stateStack.push(RIGHT);
          nodeStack.push(kdheap+2*index+1);
          stateStack.push(UNVISITED);
        }
        continue;
      }
    }
    //  now for the other cases...
    if (state == LEFT) {
      //  two cases, is this the first child or the second..
      int index = currNode-kdheap;
      const myreal planeDistance = srcPoint[currNode->plane] -
          currNode->pos[currNode->plane];
      if (planeDistance <= 0 || (2 * index + 1) > numPoints) {  //  first child
        if (((2 * index + 1) <= numPoints) &&
            SQ(planeDistance) < minDistanceSqr) {
          nodeStack.push(currNode);
          stateStack.push(RIGHT);
          nodeStack.push(kdheap + 2*index+1);
          stateStack.push(UNVISITED);
          continue;
        }
      }
    }
    if (state == RIGHT) {
      int index = currNode-kdheap;
      const myreal planeDistance = srcPoint[currNode->plane] -
          currNode->pos[currNode->plane];
      if (planeDistance > 0) {  // first child
        if (SQ(planeDistance) < minDistanceSqr) {
          nodeStack.push(currNode);
          stateStack.push(LEFT);
          nodeStack.push(kdheap+2*index);
          stateStack.push(UNVISITED);
          continue;
        }
      }
    }

    //  deal with the node
    double dist = DISTANCE(currNode->pos, srcPoint);
    if (minDistanceSqr < 0 || dist < minDistanceSqr) {
      minNode = currNode;
      minDistanceSqr = dist;
    }
  }

  if (minNode != NULL)
    return minNode->pos;
  else
    return NULL;
}

template<int DIM, typename myreal>
int kdtree<DIM, myreal>::nearestNeighbors(myreal* srcPoint, int *heap,
                                          int numNeighbors,
                                          myreal minDistanceSqr) const {
  if (!isinitialized) {
    return -1;  //  not initialized, quick return
  }

  if (minDistanceSqr  <= 0) {
    minDistanceSqr = SQ(sqrt(DISTANCE(srcPoint, kdheap[1].pos))  +
                        sqrt(DISTANCE(bboxMin, bboxMax)));
  }

  bool buildHeap = false;
  int heapIndex = 0, heapSize = numNeighbors;

  uint16 UNVISITED = 0, LEFT = 1, RIGHT = 2;
  int index = 1;

  Stack<treenode*, 64> nodeStack;
  nodeStack.push(kdheap + index);
  Stack<uint16, 64> stateStack;
  stateStack.push(UNVISITED);
  int visitednodes = 0;

  while (!nodeStack.empty()) {
    treenode* currNode = nodeStack.top();
    nodeStack.pop();
    uint16 state = stateStack.top();
    stateStack.pop();
    if (state == UNVISITED) {
      // try to see if the left child can be visited...
      visitednodes++;
      int index = currNode-kdheap;
      if ((2*index) <= numPoints) {
        const myreal planeDistance = srcPoint[currNode->plane] -
            currNode->pos[currNode->plane];
        // push this back on to the stack...
        // now check which child needs to be pushed on to the stack..
        if (planeDistance <= 0  || (2 * index + 1) > numPoints) {
          nodeStack.push(currNode);
          stateStack.push(LEFT);
          nodeStack.push(kdheap+2*index);
          stateStack.push(UNVISITED);
        } else {
          nodeStack.push(currNode);
          stateStack.push(RIGHT);
          nodeStack.push(kdheap+2*index+1);
          stateStack.push(UNVISITED);
        }
        continue;
      }
    }

    // now for the other cases...
    if (state == LEFT) {
      const myreal planeDistance = srcPoint[currNode->plane] -
          currNode->pos[currNode->plane];
      // two cases, is this the first child or the second..
      if (planeDistance <= 0 || (2*index+1) > numPoints) {  //  first child
        int index = currNode-kdheap;
        if (SQ(planeDistance) <= minDistanceSqr &&
            ((2*index+1) <= numPoints)) {
          nodeStack.push(currNode);
          stateStack.push(RIGHT);
          nodeStack.push(kdheap + 2 * index + 1);
          stateStack.push(UNVISITED);
          continue;
        }
      }
    }

    if (state == RIGHT) {
      int index = currNode-kdheap;
      const myreal planeDistance = srcPoint[currNode->plane] -
          currNode->pos[currNode->plane];
      if (planeDistance > 0) {  // first child
        if (SQ(planeDistance) <= minDistanceSqr) {
          nodeStack.push(currNode);
          stateStack.push(LEFT);
          nodeStack.push(kdheap + 2 * index);
          stateStack.push(UNVISITED);
          continue;
        }
      }
    }

    double dist = DISTANCE(srcPoint, currNode->pos);
    if (dist <= minDistanceSqr && (heapIndex < heapSize)) {
      heap[heapIndex++] = pointIndex(currNode->pos);
    } else {
      if ((dist <= minDistanceSqr)) {
        if (!buildHeap) {
          // the heap structure has not been built yet...
          buildHeap = true;
          for (int rootIndex = heapSize/2-1; rootIndex >= 0; rootIndex--) {
            // build a heap rooted... at heapIndex
            int nodeIndex = rootIndex;
            while (nodeIndex < heapSize/2) {
              myreal rootDistance = DISTANCE(srcPoint,
                                             coords + DIM * heap[nodeIndex]);
              myreal ldistance =
                DISTANCE(srcPoint, coords + DIM * heap[2*(nodeIndex+1)-1]);
              myreal rdistance;
              if ((2 * (nodeIndex + 1)) < heapSize) {
                rdistance = DISTANCE(srcPoint,
                                     coords + DIM * heap[2*(nodeIndex+1)]);
              } else {
                rdistance = 0;
              }
              if ((rootDistance > ldistance) &&
                  (rootDistance > rdistance)) {
                break;
              } else {
                if (ldistance < rdistance) {
                  int tmp = heap[nodeIndex];
                  heap[nodeIndex] = heap[2*(nodeIndex+1)];
                  heap[2*(nodeIndex+1)] = tmp;
                  nodeIndex = 2 * (nodeIndex + 1);
                } else {
                  int tmp = heap[nodeIndex];
                  heap[nodeIndex] = heap[2 * (nodeIndex + 1) - 1];
                  heap[2*(nodeIndex+1)-1] = tmp;
                  nodeIndex = 2 * (nodeIndex+1) - 1;
                }
              }
            }
          }  //  done with building the heap....
        }
        //  now update the heap..
        if (DISTANCE(srcPoint, currNode->pos) <=
            DISTANCE(srcPoint, coords + DIM * heap[0])) {
          heap[0] = pointIndex(currNode->pos);
          //  now update the heap..
          int nodeIndex = 0;
          while (nodeIndex < heapSize/2) {
            myreal rootDistance = DISTANCE(srcPoint,
                                           coords + DIM * heap[nodeIndex]);
            myreal ldistance = DISTANCE(srcPoint,
                                      coords + DIM * heap[2*(nodeIndex+1)-1]);
            myreal rdistance;
            if ((2*(nodeIndex+1)) < heapSize) {
              rdistance = DISTANCE(srcPoint,
                                   coords + DIM * heap[2*(nodeIndex+1)]);
            } else {
              rdistance = 0;
            }
            if (rootDistance > ldistance && rootDistance > rdistance) {
              break;
            } else {
              if (ldistance < rdistance) {
                int tmp = heap[nodeIndex];
                heap[nodeIndex] = heap[2* (nodeIndex + 1)];
                heap[2*(nodeIndex+1)] = tmp;
                nodeIndex = 2 * (nodeIndex + 1);
            } else {
              int tmp = heap[nodeIndex];
              heap[nodeIndex] = heap[2 * (nodeIndex + 1) - 1];
              heap[2*(nodeIndex+1)-1] = tmp;
              nodeIndex = 2 * (nodeIndex + 1) - 1;
            }
            }
          }  //  done with the update..
          minDistanceSqr = DISTANCE(srcPoint, coords + DIM * heap[0]);
        }
      }
    }
  }
  return heapIndex;
}

template<int DIM, typename myreal>
void kdtree<DIM, myreal>::intersectSphere(myreal *srcPoint,
                                          myreal sqrradius, void* dataptr,
                                          bool (*fptr)(int, void*)) const {
  if (!isinitialized) {
    return;  //  not initialized, quick return
  }

  uint16 UNVISITED = 0, LEFT = 1, RIGHT = 2;
  int index = 1;

  Stack<treenode*, 64> nodeStack;
  nodeStack.push(kdheap+index);
  Stack<uint16, 64> stateStack;
  stateStack.push(UNVISITED);
  int visitednodes = 0;

  while (!nodeStack.empty()) {
    treenode* currNode = nodeStack.top();
    nodeStack.pop();
    uint16 state = stateStack.top();
    stateStack.pop();

    if (state == UNVISITED) {
      // try to see if the left child can be visited...
      visitednodes++;
      int index = currNode-kdheap;
      if ((2*index) <= numPoints) {
        // push this back on to the stack...
        nodeStack.push(currNode);
        const myreal planeDistance =
            srcPoint[currNode->plane] - currNode->pos[currNode->plane];
        if (planeDistance <= 0 || (2*index+1) > numPoints) {
          stateStack.push(LEFT);
          nodeStack.push(kdheap + 2 * index);
          stateStack.push(UNVISITED);
        } else {
          stateStack.push(RIGHT);
          nodeStack.push(kdheap+2*index+1);
          stateStack.push(UNVISITED);
        }
        continue;
      }
    }

    // now for the other cases...
    if (state == LEFT) {
      // two cases, is this the first child or the second..
      int index = currNode-kdheap;
      const myreal planeDistance =
          srcPoint[currNode->plane] - currNode->pos[currNode->plane];
      if (planeDistance <= 0 || (2*index+1) > numPoints) {  //  first child
        if (((2*index+1) <= numPoints) &&
            (SQ(planeDistance)) <= sqrradius) {
          nodeStack.push(currNode);
          stateStack.push(RIGHT);
          nodeStack.push(kdheap + 2 * index + 1);
          stateStack.push(UNVISITED);
          continue;
        }
      }
    }

    if (state == RIGHT) {
      int index = currNode-kdheap;
      const myreal planeDistance =
          srcPoint[currNode->plane] - currNode->pos[currNode->plane];
      if (planeDistance > 0) {  // first child
        if (SQ(planeDistance) <= sqrradius) {
          nodeStack.push(currNode);
          stateStack.push(LEFT);
          nodeStack.push(kdheap + 2 * index);
          stateStack.push(UNVISITED);
          continue;
        }
      }
    }

    double dist = DISTANCE(srcPoint, currNode->pos);
    if (dist <= sqrradius) {
      if (!fptr(pointIndex(currNode->pos), dataptr)) {
        while (!nodeStack.empty()) {
          nodeStack.pop();
          stateStack.pop();
        }
        return;
      }
    }
  }
  return;
}

//  count the number of points inside the given sphere
template<int DIM, typename myreal>
int kdtree<DIM, myreal>::pointsInSphere(myreal *srcPoint,
                                        myreal sqrradius) const {
  if (!isinitialized) {
    return 0;  //  not initialized, quick return
  }

  int numInside = 0;
  uint16 UNVISITED = 0, LEFT = 1, RIGHT = 2;
  int index = 1;

  Stack<treenode*, 64> nodeStack;
  nodeStack.push(kdheap + index);
  Stack<uint16, 64> stateStack;
  stateStack.push(UNVISITED);
  int visitednodes = 0;

  while (!nodeStack.empty()) {
    treenode* currNode = nodeStack.top();
    nodeStack.pop();
    uint16 state = stateStack.top();
    stateStack.pop();

    if (state == UNVISITED) {
      // try to see if the left child can be visited...
      visitednodes++;
      int index = currNode - kdheap;
      if ((2*index) <= numPoints) {
        nodeStack.push(currNode);
        const myreal planeDistance = srcPoint[currNode->plane] -
            currNode->pos[currNode->plane];
        if (planeDistance <= 0 || (2*index+1) > numPoints) {
          stateStack.push(LEFT);
          nodeStack.push(kdheap + 2 * index);
          stateStack.push(UNVISITED);
        } else {
          stateStack.push(RIGHT);
          nodeStack.push(kdheap + 2 * index + 1);
          stateStack.push(UNVISITED);
        }
        continue;
      }
    }

    // now for the other cases...
    if (state == LEFT) {
      // two cases, is this the first child or the second..
      int index = currNode-kdheap;
      const myreal planeDistance = srcPoint[currNode->plane] -
          currNode->pos[currNode->plane];
      if (planeDistance <= 0 || (2 * index + 1) > numPoints) {  //  first child
        if (((2 * index + 1) <= numPoints) &&
            (SQ(planeDistance)) <= sqrradius) {
          nodeStack.push(currNode);
          stateStack.push(RIGHT);
          nodeStack.push(kdheap+2*index+1);
          stateStack.push(UNVISITED);
          continue;
        }
      }
    }

    if (state == RIGHT) {
      int index = currNode-kdheap;
      const myreal planeDistance = srcPoint[currNode->plane] -
          currNode->pos[currNode->plane];
      if (planeDistance > 0) {  // first child
        if (SQ(planeDistance) <= sqrradius) {
          nodeStack.push(currNode);
          stateStack.push(LEFT);
          nodeStack.push(kdheap + 2 * index);
          stateStack.push(UNVISITED);
          continue;
        }
      }
    }
    double dist = DISTANCE(srcPoint, currNode->pos);
    if (dist <= sqrradius)
      numInside++;
  }
  return numInside;
}

}  //  namespace kdtree

#endif  // THIRD_PARTY_KDTREE_KDTREE_H__
