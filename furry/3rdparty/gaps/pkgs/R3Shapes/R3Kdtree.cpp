// Source file for R3Kdtree class

#ifndef __R3KDTREE__C__
#define __R3KDTREE__C__




////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "R3Shapes/R3Shapes.h"





////////////////////////////////////////////////////////////////////////
// Constant definitions
////////////////////////////////////////////////////////////////////////

static const int R3kdtree_max_points_per_node = 32;





////////////////////////////////////////////////////////////////////////
// Node class definition
////////////////////////////////////////////////////////////////////////

// Node declaration

template <class PtrType>
class R3KdtreeNode {
public:
  R3KdtreeNode(R3KdtreeNode<PtrType> *parent);

public:
  class R3KdtreeNode<PtrType> *parent;
  class R3KdtreeNode<PtrType> *children[2];
  RNScalar split_coordinate;
  RNDimension split_dimension;
  PtrType points[R3kdtree_max_points_per_node];
  int npoints;
};



// Node constructor

template <class PtrType>
R3KdtreeNode<PtrType>::
R3KdtreeNode(R3KdtreeNode<PtrType> *parent)
  : parent(parent),
    npoints(0)
{
  // Initialize everything
  children[0] = NULL;
  children[1] = NULL;
  split_coordinate = 0;
  split_dimension = 0;
}





////////////////////////////////////////////////////////////////////////
// Public tree-level functions
////////////////////////////////////////////////////////////////////////

template <class PtrType>
R3Kdtree<PtrType>::
R3Kdtree(const R3Box& bbox, int position_offset)
  : bbox(bbox),
    position_offset(position_offset),
    position_callback(NULL),
    position_callback_data(NULL),
    nnodes(0)
{
  // Create root node
  root = new R3KdtreeNode<PtrType>(NULL);
  assert(root);
}



template <class PtrType>
R3Kdtree<PtrType>::
R3Kdtree(const RNArray<PtrType>& points, int position_offset)
  : position_offset(position_offset),
    position_callback(NULL),
    position_callback_data(NULL),
    nnodes(1)
{
  // Create root node
  root = new R3KdtreeNode<PtrType>(NULL);
  assert(root);

  // Determine bounding box 
  bbox = R3null_box;
  for (int i = 0; i < points.NEntries(); i++) {
    bbox.Union(Position(points[i]));
  }

  // Allocate copy of points array (so that it can be sorted)
  PtrType *copy = new PtrType [ points.NEntries() ];
  assert(copy);

  // Copy points
  for (int i = 0; i < points.NEntries(); i++) 
    copy[i] = points[i];

  // Insert points into root
  InsertPoints(root, bbox, copy, points.NEntries());

  // Delete copy of points array
  delete [] copy;
}



template <class PtrType>
R3Kdtree<PtrType>::
R3Kdtree(const RNArray<PtrType>& points, R3Point (*position_callback)(PtrType, void *), void *position_callback_data)
  : position_offset(-1),
    position_callback(position_callback),
    position_callback_data(position_callback_data),
    nnodes(1)
{
  // Create root node
  root = new R3KdtreeNode<PtrType>(NULL);
  assert(root);

  // Determine bounding box 
  bbox = R3null_box;
  for (int i = 0; i < points.NEntries(); i++) {
    bbox.Union(Position(points[i]));
  }

  // Allocate copy of points array (so that it can be sorted)
  PtrType *copy = new PtrType [ points.NEntries() ];
  assert(copy);

  // Copy points
  for (int i = 0; i < points.NEntries(); i++) 
    copy[i] = points[i];

  // Insert points into root
  InsertPoints(root, bbox, copy, points.NEntries());

  // Delete copy of points array
  delete [] copy;

}



template <class PtrType>
R3Kdtree<PtrType>::
~R3Kdtree(void)
{
  // Check root
  if (!root) return;

  // Traverse tree deleting nodes
  RNArray<R3KdtreeNode<PtrType> *> stack;
  stack.InsertTail(root);
  while (!stack.IsEmpty()) {
    R3KdtreeNode<PtrType> *node = stack.Tail();
    stack.RemoveTail();
    if (node->children[0]) stack.Insert(node->children[0]);
    if (node->children[1]) stack.Insert(node->children[1]);
    delete node;
  }
}



template <class PtrType>
const R3Box& R3Kdtree<PtrType>::
BBox(void) const
{
  // Return bounding box of the whole KD tree
  return bbox;
}



template <class PtrType>
int R3Kdtree<PtrType>::
NNodes(void) const
{
  // Return number of nodes
  return nnodes;
}



#if 0

template <class PtrType>
int R3Kdtree<PtrType>::
NEntries(void) const
{
  // Check root
  if (!root) return 0;

  // Traverse tree to count number of points
  int npoints = 0;
  RNArray<R3KdtreeNode<PtrType> *> stack;
  stack.InsertTail(root);
  while (!stack.IsEmpty()) {
    R3KdtreeNode<PtrType> *node = stack.Tail();
    stack.RemoveTail();
    npoints += node->npoints;
    if (node->children[0]) stack.Insert(node->children[0]);
    if (node->children[1]) stack.Insert(node->children[1]);
  }

  // Return total number of points
  return npoints;
}

#endif



////////////////////////////////////////////////////////////////////////
// Finding the closest one point
////////////////////////////////////////////////////////////////////////

template <class PtrType>
void R3Kdtree<PtrType>::
FindClosest(R3KdtreeNode<PtrType> *node, const R3Box& node_box, 
  PtrType query_point, const R3Point& query_position, 
  RNScalar min_distance_squared, RNScalar max_distance_squared,
  int (*IsCompatible)(PtrType, PtrType, void *), void *compatible_data, 
  PtrType& closest_point, RNScalar& closest_distance_squared) const
{
  // Check if node is interior
  if (node->children[0]) {
    assert(node->children[1]);

    // Find and check axial distances from point to node box
    RNLength dx, dy, dz;
    if (RNIsGreater(query_position.X(), node_box.XMax())) dx = query_position.X() - node_box.XMax();
    else if (RNIsLess(query_position.X(), node_box.XMin())) dx = node_box.XMin()- query_position.X();
    else dx = 0.0;
    RNLength dx_squared = dx * dx;
    if (dx_squared >= closest_distance_squared) return;
    if (RNIsGreater(query_position.Y(), node_box.YMax())) dy = query_position.Y() - node_box.YMax();
    else if (RNIsLess(query_position.Y(), node_box.YMin())) dy = node_box.YMin()- query_position.Y();
    else dy = 0.0;
    RNLength dy_squared = dy * dy;
    if (dy_squared >= closest_distance_squared) return;
    if (RNIsGreater(query_position.Z(), node_box.ZMax())) dz = query_position.Z() - node_box.ZMax();
    else if (RNIsLess(query_position.Z(), node_box.ZMin())) dz = node_box.ZMin()- query_position.Z();
    else dz = 0.0;
    RNLength dz_squared = dz * dz;
    if (dz_squared >= closest_distance_squared) return;
    
    // Find and check actual distance from point to node box
    RNLength distance_squared = 0;
    if ((dy == 0.0) && (dz == 0.0)) distance_squared = dx_squared;
    else if ((dx == 0.0) && (dz == 0.0)) distance_squared = dy_squared;
    else if ((dx == 0.0) && (dy == 0.0)) distance_squared = dz_squared;
    else distance_squared = dx_squared + dy_squared + dz_squared;
    if (distance_squared >= closest_distance_squared) return;

    // Compute distance from point to split plane
    RNLength side = query_position[node->split_dimension] - node->split_coordinate;

    // Search children nodes
    if (side <= 0) {
      // Search negative side first
      R3Box child_box(node_box);
      child_box[RN_HI][node->split_dimension] = node->split_coordinate;
      FindClosest(node->children[0], child_box, 
        query_point, query_position, 
        min_distance_squared, max_distance_squared, 
        IsCompatible, compatible_data,
        closest_point, closest_distance_squared);
      if (side*side < closest_distance_squared) {
        R3Box child_box(node_box);
        child_box[RN_LO][node->split_dimension] = node->split_coordinate;
        FindClosest(node->children[1], child_box, 
          query_point, query_position, 
          min_distance_squared, max_distance_squared, 
          IsCompatible, compatible_data,
          closest_point, closest_distance_squared);
      }
    }
    else {
      // Search positive side first
      R3Box child_box(node_box);
      child_box[RN_LO][node->split_dimension] = node->split_coordinate;
      FindClosest(node->children[1], child_box, 
        query_point, query_position, 
        min_distance_squared, max_distance_squared, 
        IsCompatible, compatible_data,
        closest_point, closest_distance_squared);
      if (side*side < closest_distance_squared) {
        R3Box child_box(node_box);
        child_box[RN_HI][node->split_dimension] = node->split_coordinate;
        FindClosest(node->children[0], child_box, 
          query_point, query_position, 
          min_distance_squared, max_distance_squared, 
          IsCompatible, compatible_data,
          closest_point, closest_distance_squared);
      }
    }
  }
  else {
    for (int i = 0; i < node->npoints; i++) {
      PtrType point = node->points[i];
      RNLength distance_squared = R3SquaredDistance(query_position, Position(point));
      if ((distance_squared >= min_distance_squared) && 
         (distance_squared <= closest_distance_squared)) {
        if (!IsCompatible || !query_point || IsCompatible(query_point, point, compatible_data)) {
          closest_distance_squared = distance_squared;
          closest_point = point;
        }
      }
    }
  }
}



template <class PtrType>
PtrType R3Kdtree<PtrType>::
FindClosest(PtrType query_point, 
  RNScalar min_distance, RNScalar max_distance, 
  int (*IsCompatible)(PtrType, PtrType, void *), void *compatible_data, 
  RNScalar *closest_distance) const
{
  // Check root
  if (!root) return NULL;

  // Use squared distances for efficiency
  RNLength min_distance_squared = min_distance * min_distance;
  RNLength max_distance_squared = max_distance * max_distance;

  // Initialize nearest point 
  PtrType nearest_point = NULL;
  RNLength nearest_distance_squared = max_distance_squared;

  // Search nodes recursively
  FindClosest(root, bbox, 
    query_point, Position(query_point),
    min_distance_squared, max_distance_squared, 
    IsCompatible, compatible_data, 
    nearest_point, nearest_distance_squared);

  // Return closest distance
  if (closest_distance) *closest_distance = sqrt(nearest_distance_squared);

  // Return closest point
  return nearest_point;
}



template <class PtrType>
PtrType R3Kdtree<PtrType>::
FindClosest(PtrType query_point, 
  RNScalar min_distance, RNScalar max_distance, 
  RNScalar *closest_distance) const
{
  // Find the closest point
  return FindClosest(query_point, min_distance, max_distance, closest_distance);
}



template <class PtrType>
PtrType R3Kdtree<PtrType>::
FindClosest(const R3Point& query_position, 
  RNScalar min_distance, RNScalar max_distance, 
  RNScalar *closest_distance) const
{
  // Check root
  if (!root) return NULL;

  // Use squared distances for efficiency
  RNLength min_distance_squared = min_distance * min_distance;
  RNLength max_distance_squared = max_distance * max_distance;

  // Initialize nearest point 
  PtrType nearest_point = NULL;
  RNLength nearest_distance_squared = max_distance_squared;

  // Search nodes recursively
  FindClosest(root, bbox, 
    NULL, query_position, 
    min_distance_squared, max_distance_squared, 
    NULL, NULL, 
    nearest_point, nearest_distance_squared);

  // Return closest distance
  if (closest_distance) *closest_distance = sqrt(nearest_distance_squared);

  // Return closest point
  return nearest_point;
}



////////////////////////////////////////////////////////////////////////
// Finding the closest K points
////////////////////////////////////////////////////////////////////////

template <class PtrType>
void R3Kdtree<PtrType>::
FindClosest(R3KdtreeNode<PtrType> *node, const R3Box& node_box, 
  PtrType query_point, const R3Point& query_position, 
  RNScalar min_distance_squared, RNScalar max_distance_squared, int max_points,
  int (*IsCompatible)(PtrType, PtrType, void *), void *compatible_data, 
  RNArray<PtrType>& points, RNLength *distances_squared) const
{
  // Update max distance squared
  if (points.NEntries() == max_points) {
    max_distance_squared = distances_squared[max_points-1];
  }

  // Check if node is interior
  if (node->children[0]) {
    assert(node->children[1]);

    // Find and check axial distances from point to node box
    RNLength dx, dy, dz;
    if (RNIsGreater(query_position.X(), node_box.XMax())) dx = query_position.X() - node_box.XMax();
    else if (RNIsLess(query_position.X(), node_box.XMin())) dx = node_box.XMin()- query_position.X();
    else dx = 0.0;
    RNLength dx_squared = dx * dx;
    if (dx_squared > max_distance_squared) return;
    if (RNIsGreater(query_position.Y(), node_box.YMax())) dy = query_position.Y() - node_box.YMax();
    else if (RNIsLess(query_position.Y(), node_box.YMin())) dy = node_box.YMin()- query_position.Y();
    else dy = 0.0;
    RNLength dy_squared = dy * dy;
    if (dy_squared > max_distance_squared) return;
    if (RNIsGreater(query_position.Z(), node_box.ZMax())) dz = query_position.Z() - node_box.ZMax();
    else if (RNIsLess(query_position.Z(), node_box.ZMin())) dz = node_box.ZMin()- query_position.Z();
    else dz = 0.0;
    RNLength dz_squared = dz * dz;
    if (dz_squared > max_distance_squared) return;
    
    // Find and check actual distance from point to node box
    RNLength distance_squared = 0;
    if ((dy == 0.0) && (dz == 0.0)) distance_squared = dx_squared;
    else if ((dx == 0.0) && (dz == 0.0)) distance_squared = dy_squared;
    else if ((dx == 0.0) && (dy == 0.0)) distance_squared = dz_squared;
    else distance_squared = dx_squared + dy_squared + dz_squared;
    if (distance_squared > max_distance_squared) return;

    // Compute distance from point to split plane
    RNLength side = query_position[node->split_dimension] - node->split_coordinate;

    // Search children nodes
    if ((side <= 0) || (side*side <= max_distance_squared)) {
      // Search negative side 
      R3Box child_box(node_box);
      child_box[RN_HI][node->split_dimension] = node->split_coordinate;
      FindClosest(node->children[0], child_box, query_point, query_position, 
        min_distance_squared, max_distance_squared, max_points, IsCompatible, compatible_data,
        points, distances_squared);
    }
    if ((side >= 0) || (side*side <= max_distance_squared)) {
      R3Box child_box(node_box);
      child_box[RN_LO][node->split_dimension] = node->split_coordinate;
      FindClosest(node->children[1], child_box, query_point, query_position, 
        min_distance_squared, max_distance_squared, max_points, IsCompatible, compatible_data,
        points, distances_squared);
    }
  }
  else {
    // Search points
    for (int i = 0; i < node->npoints; i++) {
      PtrType point = node->points[i];
      RNLength distance_squared = R3SquaredDistance(query_position, Position(point));
      if ((distance_squared >= min_distance_squared) && 
          (distance_squared <= max_distance_squared)) {

        // Check if point is compatible
        if (!IsCompatible || !query_point || IsCompatible(query_point, point, compatible_data)) {

          // Find slot for point (points are sorted by distance)
          int slot = 0;
          while (slot < points.NEntries()) {
            if (distance_squared < distances_squared[slot]) break;
            slot++;
          }
          
          // Insert point and distance into sorted arrays
          if (slot < max_points) {
            int first = points.NEntries();
            if (first >= max_points) first = max_points-1;
            for (int j = first; j > slot; j--) distances_squared[j] = distances_squared[j-1];
            distances_squared[slot] = distance_squared;
            points.InsertKth(point, slot);
            points.Truncate(max_points);
          }
        }
      }
    }
  }
}



template <class PtrType>
int R3Kdtree<PtrType>::
FindClosest(PtrType query_point, RNScalar min_distance, RNScalar max_distance, int max_points, 
  RNArray<PtrType>& points, RNLength *distances) const
{
  // Find closest within some distance
  return FindClosest(query_point, min_distance, max_distance, max_points, NULL, NULL, points, distances);
}




template <class PtrType>
int R3Kdtree<PtrType>::
FindClosest(PtrType query_point, 
  RNScalar min_distance, RNScalar max_distance, int max_points, 
  int (*IsCompatible)(PtrType, PtrType, void *), void *compatible_data, 
  RNArray<PtrType>& points, RNLength *distances) const
{
  // Check root
  if (!root) return 0;

  // Use squared distances for efficiency
  RNLength min_distance_squared = min_distance * min_distance;
  RNLength max_distance_squared = max_distance * max_distance;

  // Allocate temporary array of squared distances to max_points closest points
  RNLength *distances_squared = new RNLength [ max_points ];

  // Search nodes recursively
  FindClosest(root, bbox, 
    query_point, Position(query_point),
    min_distance_squared, max_distance_squared, max_points, 
    IsCompatible, compatible_data,
    points, distances_squared);

  // Update return distances
  if (distances) {
    for (int i = 0; i < points.NEntries(); i++) {
      distances[i] = sqrt(distances_squared[i]);
    }
  }

  // Delete temporary array of squared distances
  delete [] distances_squared;

  // Return number of points
  return points.NEntries();
}



template <class PtrType>
int R3Kdtree<PtrType>::
FindClosest(const R3Point& query_position, RNScalar min_distance, RNScalar max_distance, int max_points, 
  RNArray<PtrType>& points, RNLength *distances) const
{
  // Check root
  if (!root) return 0;

  // Use squared distances for efficiency
  RNLength min_distance_squared = min_distance * min_distance;
  RNLength max_distance_squared = max_distance * max_distance;

  // Allocate temporary array of squared distances to max_points closest points
  RNLength *distances_squared = new RNLength [ max_points ];

  // Search nodes recursively
  FindClosest(root, bbox, 
    NULL, query_position, 
    min_distance_squared, max_distance_squared, max_points, 
    NULL, NULL, 
    points, distances_squared);

  // Update return distances
  if (distances) {
    for (int i = 0; i < points.NEntries(); i++) {
      distances[i] = sqrt(distances_squared[i]);
    }
  }

  // Delete temporary array of squared distances
  delete [] distances_squared;

  // Return number of points
  return points.NEntries();
}



////////////////////////////////////////////////////////////////////////
// Finding all points within a radius
////////////////////////////////////////////////////////////////////////

template <class PtrType>
void R3Kdtree<PtrType>::
FindAll(R3KdtreeNode<PtrType> *node, const R3Box& node_box, 
  PtrType query_point, const R3Point& query_position, 
  RNScalar min_distance_squared, RNScalar max_distance_squared, 
  int (*IsCompatible)(PtrType, PtrType, void *), void *compatible_data, 
  RNArray<PtrType>& points) const
{
  // Check if node is interior
  if (node->children[0]) {
    assert(node->children[1]);

    // Find and check axial distances from point to node box
    RNLength dx, dy, dz;
    if (RNIsGreater(query_position.X(), node_box.XMax())) dx = query_position.X() - node_box.XMax();
    else if (RNIsLess(query_position.X(), node_box.XMin())) dx = node_box.XMin()- query_position.X();
    else dx = 0.0;
    RNLength dx_squared = dx * dx;
    if (dx_squared > max_distance_squared) return;
    if (RNIsGreater(query_position.Y(), node_box.YMax())) dy = query_position.Y() - node_box.YMax();
    else if (RNIsLess(query_position.Y(), node_box.YMin())) dy = node_box.YMin()- query_position.Y();
    else dy = 0.0;
    RNLength dy_squared = dy * dy;
    if (dy_squared > max_distance_squared) return;
    if (RNIsGreater(query_position.Z(), node_box.ZMax())) dz = query_position.Z() - node_box.ZMax();
    else if (RNIsLess(query_position.Z(), node_box.ZMin())) dz = node_box.ZMin()- query_position.Z();
    else dz = 0.0;
    RNLength dz_squared = dz * dz;
    if (dz_squared > max_distance_squared) return;
    
    // Find and check actual distance from point to node box
    RNLength distance_squared = 0;
    if ((dy == 0.0) && (dz == 0.0)) distance_squared = dx_squared;
    else if ((dx == 0.0) && (dz == 0.0)) distance_squared = dy_squared;
    else if ((dx == 0.0) && (dy == 0.0)) distance_squared = dz_squared;
    else distance_squared = dx_squared + dy_squared + dz_squared;
    if (distance_squared > max_distance_squared) return;

    // Compute distance from point to split plane
    RNLength side = query_position[node->split_dimension] - node->split_coordinate;

    // Search children nodes
    if ((side <= 0) || (side*side <= max_distance_squared)) {
      // Search negative side 
      R3Box child_box(node_box);
      child_box[RN_HI][node->split_dimension] = node->split_coordinate;
      FindAll(node->children[0], child_box, query_point, query_position, 
        min_distance_squared, max_distance_squared, 
        IsCompatible, compatible_data, points);
    }
    if ((side >= 0) || (side*side <= max_distance_squared)) {
      R3Box child_box(node_box);
      child_box[RN_LO][node->split_dimension] = node->split_coordinate;
      FindAll(node->children[1], child_box, query_point, query_position, 
        min_distance_squared, max_distance_squared, 
        IsCompatible, compatible_data, points);
    }
  }
  else {
    // Search points
    for (int i = 0; i < node->npoints; i++) {
      PtrType point = node->points[i];
      RNLength distance_squared = R3SquaredDistance(query_position, Position(point));
      if ((distance_squared >= min_distance_squared) && 
          (distance_squared <= max_distance_squared)) {
        if (!IsCompatible || !query_point || IsCompatible(query_point, point, compatible_data)) {
          points.Insert(point);
        }
      }
    }
  }
}



template <class PtrType>
int R3Kdtree<PtrType>::
FindAll(PtrType query_point,
  RNScalar min_distance, RNScalar max_distance, 
  int (*IsCompatible)(PtrType, PtrType, void *), void *compatible_data, 
  RNArray<PtrType>& points) const
{
  // Check root
  if (!root) return 0;

  // Use squared distances for efficiency
  RNLength min_distance_squared = min_distance * min_distance;
  RNLength max_distance_squared = max_distance * max_distance;

  // Search nodes recursively
  FindAll(root, bbox, 
    query_point, Position(query_point), 
    min_distance_squared, max_distance_squared, 
    IsCompatible, compatible_data, 
    points);

  // Return number of points
  return points.NEntries();
}



template <class PtrType>
int R3Kdtree<PtrType>::
FindAll(PtrType query_point, RNScalar min_distance, RNScalar max_distance, RNArray<PtrType>& points) const
{
  // Find all within some distance
  return FindAll(query_point, min_distance, max_distance, NULL, NULL, points);
}




template <class PtrType>
int R3Kdtree<PtrType>::
FindAll(const R3Point& query_position, RNScalar min_distance, RNScalar max_distance, RNArray<PtrType>& points) const
{
  // Check root
  if (!root) return 0;

  // Use squared distances for efficiency
  RNLength min_distance_squared = min_distance * min_distance;
  RNLength max_distance_squared = max_distance * max_distance;

  // Search nodes recursively
  FindAll(root, bbox, 
    NULL, query_position, 
    min_distance_squared, max_distance_squared, 
    NULL, NULL, 
    points);

  // Return number of points
  return points.NEntries();
}



////////////////////////////////////////////////////////////////////////
// Internal tree creation functions
////////////////////////////////////////////////////////////////////////

template <class PtrType>
int R3Kdtree<PtrType>::
PartitionPoints(PtrType *points, int npoints, RNDimension dim, int imin, int imax)
{
  // Check range
  assert(imin <= imax);
  assert(imin >= 0);
  assert(imin < npoints);
  assert(imax >= 0);
  assert(imax < npoints);
  if (imin == imax) return imin;

  // Choose a coordinate at random to split upon
  int irand = (int) (imin + RNRandomScalar() * (imax - imin + 1));
  if (irand < imin) irand = imin;
  if (irand > imax) irand = imax;
  RNCoord split_coord = Position(points[irand])[dim];

  // Swap values at irand and imax
  PtrType swap = points[irand];
  points[irand] = points[imax];
  points[imax] = swap;

  // Partition points according to coordinate
  int split_index = imin;
  int middle_index = (imin + imax) / 2;
  for (int i = imin; i < imax; i++) {
    assert(split_index <= i);
    RNCoord coord = Position(points[i])[dim];
    if (coord < split_coord) {
      PtrType swap = points[split_index];
      points[split_index] = points[i];
      points[i] = swap;
      split_index++;
    }
    else if (coord == split_coord) {
      if (split_index < middle_index) {
        PtrType swap = points[split_index];
        points[split_index] = points[i];
        points[i] = swap;
        split_index++;
      }
    }
  }

  // Swap values at split_index and imax
  swap = points[split_index];
  points[split_index] = points[imax];
  points[imax] = swap;

  // Now split_index has value split_coord
  // All values to the left of split_index have values < split_coord
  // All values to the right of split_index have values >= split_coord

  // Recurse until we find the median
  if (split_index == imin) {
    if (imin >= npoints/2) return imin;
    else return PartitionPoints(points, npoints, dim, imin+1, imax);
  }
  else if (split_index == imax) {
    if (imax <= npoints/2) return imax;
    else return PartitionPoints(points, npoints, dim, imin, imax-1);
  }
  else {
    if (split_index == npoints/2) return split_index;
    else if (split_index > npoints/2) return PartitionPoints(points, npoints, dim, imin, split_index-1);
    else return PartitionPoints(points, npoints, dim, split_index+1, imax);
  }
}



template <class PtrType>
void R3Kdtree<PtrType>::
InsertPoints(R3KdtreeNode<PtrType> *node, const R3Box& node_box, PtrType *points, int npoints) 
{
  // Make sure node is an empty leaf
  assert(node);
  assert(node->children[0] == NULL);
  assert(node->children[1] == NULL);
  assert(node->npoints == 0);

  // Check number of points
  if (npoints <= R3kdtree_max_points_per_node) {
    // Insert new points into leaf node and return
    for (int i = 0; i < npoints; i++) {
      node->points[node->npoints++] = points[i];
    }
  }
  else {
    // Find dimension to split along
    node->split_dimension = node_box.LongestAxis();
    
    // Partition points according to coordinates in split_dimension
    int split_index = PartitionPoints(points, npoints, node->split_dimension, 0, npoints-1);
    assert((split_index >= 0) && (split_index < npoints));

    // Determine split coordinate
    node->split_coordinate = Position(points[split_index])[node->split_dimension];
    assert(node->split_coordinate >= node_box[RN_LO][node->split_dimension]);
    assert(node->split_coordinate <= node_box[RN_HI][node->split_dimension]);

    // Construct children node boxes
    R3Box node0_box(node_box);
    R3Box node1_box(node_box);
    node0_box[RN_HI][node->split_dimension] = node->split_coordinate;
    node1_box[RN_LO][node->split_dimension] = node->split_coordinate;

    // Create children
    node->children[0] = new R3KdtreeNode<PtrType>(node);
    node->children[1] = new R3KdtreeNode<PtrType>(node);

    // Insert points into children
    InsertPoints(node->children[0], node0_box, points, split_index);
    InsertPoints(node->children[1], node1_box, &points[split_index], npoints - split_index);

    // Increment number of nodes
    nnodes += 2;
  }
}



template <class PtrType>
void R3Kdtree<PtrType>::
Outline(R3KdtreeNode<PtrType> *node, const R3Box& node_box) const
{
  // Draw kdtree nodes recursively
  if (node->children[0]) {
    assert(node->children[1]);
    assert(node->split_coordinate >= node_box[RN_LO][node->split_dimension]);
    assert(node->split_coordinate <= node_box[RN_HI][node->split_dimension]);
    R3Box child0_box(node_box);
    R3Box child1_box(node_box);
    child0_box[RN_HI][node->split_dimension] = node->split_coordinate;
    child1_box[RN_LO][node->split_dimension] = node->split_coordinate;
    Outline(node->children[0], child0_box);
    Outline(node->children[1], child1_box);
  }
  else {
    node_box.Outline();
  }
}



template <class PtrType>
void R3Kdtree<PtrType>::
Outline(void) const
{
  // Draw kdtree nodes recursively
  if (!root) return;
  Outline(root, bbox);
}



template <class PtrType>
int R3Kdtree<PtrType>::
PrintBalance(R3KdtreeNode<PtrType> *node, int depth) const
{
  // Check node
  if (!node) return 0;

  // Initialize number of decendents
  int ndecendents0 = 0;
  int ndecendents1 = 0;

  // Process interior node
  if (node->children[0] && node->children[1]) {
    // Print balance of children
    ndecendents0 = PrintBalance(node->children[0], depth+1);
    ndecendents1 = PrintBalance(node->children[1], depth+1);

    // Print balance of this node
    printf("%d", depth);
    for (int i = 0; i <= depth; i++) printf("  ");
    printf("I %d %d %g\n", ndecendents0, ndecendents1, (double) ndecendents0 / (double) ndecendents1);
  }
  else {
    printf("%d", depth);
    for (int i = 0; i <= depth; i++) printf("  ");
    printf("L %d\n", node->npoints);
  }

  // Return number of nodes rooted in this subtree
  return 1 + ndecendents0 + ndecendents1;
}



template <class PtrType>
void R3Kdtree<PtrType>::
PrintDebugInfo(void) const
{
  // Check root
  if (!root) return;
  PrintBalance(root, 0);
}


#endif


