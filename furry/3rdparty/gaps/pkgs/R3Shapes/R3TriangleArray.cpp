/* Source file for the R3 triangle class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Class type definitions */

RN_CLASS_TYPE_DEFINITIONS(R3TriangleArray);



/* Public functions */

int 
R3InitTriangleArray()
{
    /* Return success */
    return TRUE;
}



void 
R3StopTriangleArray()
{
}



R3TriangleArray::
R3TriangleArray(void)
    : bbox(R3null_box)
{
}



R3TriangleArray::
R3TriangleArray(const R3TriangleArray& array)
  : vertices(array.vertices),
    triangles(array.triangles),
    bbox(array.bbox)
{
}



R3TriangleArray::
R3TriangleArray(const RNArray<R3TriangleVertex *>& vertices, const RNArray<R3Triangle *>& triangles)
  : vertices(vertices),
    triangles(triangles),
    bbox(R3null_box)
{
    // Update bounding box
    Update();
}



const RNBoolean R3TriangleArray::
IsPoint (void) const
{
  // Check if all vertices are at same position
  if (vertices.NEntries() > 1) {
    const R3Point& p0 = vertices[0]->Position();
    for (int i = 1; i < vertices.NEntries(); i++) 
      if (!R3Contains(p0, vertices[i]->Position())) return FALSE;
  }
  return TRUE;
}



const RNBoolean R3TriangleArray::
IsLinear (void) const
{
    // Assume not ???
    return FALSE;
}



const RNBoolean R3TriangleArray::
IsPlanar(void) const
{
  // Check if all triangles are on same plane
  if (triangles.NEntries() > 1) {
    const R3Plane& p0 = triangles[0]->Plane();
    for (int i = 1; i < triangles.NEntries(); i++) 
      if (!R3Contains(p0, triangles[i]->Plane())) return FALSE;
  }
  return TRUE;
}



const RNBoolean R3TriangleArray::
IsConvex(void) const
{
    // If planar - may still be concave ???
    return IsPlanar();
}



const RNInterval R3TriangleArray::
NFacets(void) const
{
    // Return number of trianglets (triangles)
    return RNInterval(triangles.NEntries(), triangles.NEntries());
}



const RNLength R3TriangleArray::
Length (void) const
{
    // Return cumulative perimeter of triangles
    RNLength length = 0.0;
    for (int i = 0; i < triangles.NEntries(); i++)
      length += triangles[i]->Length();
    return length;
}



const RNArea R3TriangleArray::
Area(void) const
{
    // Return cumulative area of triangles
    RNArea area = 0.0;
    for (int i = 0; i < triangles.NEntries(); i++)
      area += triangles[i]->Area();
    return area;
}



const R3Point R3TriangleArray::
Centroid(void) const
{
    // Return centroid of vertices
    R3Point centroid = R3zero_point;
    for (int i = 0; i < vertices.NEntries(); i++) 
      centroid += vertices[i]->Position();
    centroid /= (RNScalar) vertices.NEntries();
    return centroid;
}



const R3Shape& R3TriangleArray::
BShape(void) const
{
    // Return bounding shape
    return bbox;
}



const R3Box R3TriangleArray::
BBox(void) const
{
    // Return bounding box
    return bbox;
}



const R3Sphere R3TriangleArray::
BSphere(void) const
{
    // Return bounding sphere
    return bbox.BSphere();
}



void R3TriangleArray::
Flip (void)
{
    // Flip the triangles
    for (int i = 0; i < triangles.NEntries(); i++)
      triangles[i]->Flip();
}



void R3TriangleArray::
Mirror (const R3Plane& plane)
{
    // Mirror vertices
    for (int i = 0; i < vertices.NEntries(); i++) 
      vertices[i]->Mirror(plane);

    // Update the triangles
    for (int i = 0; i < triangles.NEntries(); i++)
      triangles[i]->Update();

    // Update the bounding box
    Update();
}



void R3TriangleArray::
Transform (const R3Transformation& transformation) 
{
    // Transform vertices
    for (int i = 0; i < vertices.NEntries(); i++) 
      vertices[i]->Transform(transformation);

    // Update the triangles
    for (int i = 0; i < triangles.NEntries(); i++)
      triangles[i]->Update();

    // Update the bounding box
    Update();
}



void R3TriangleArray::
MoveVertex(R3TriangleVertex *vertex, const R3Point& position)
{
    // Move vertex
    vertex->SetPosition(position);

    // Update triangles
    for (int i = 0; i < triangles.NEntries(); i++) {
      R3Triangle *triangle = triangles.Kth(i);
      for (int j = 0; j < 3; j++) {
        if (triangle->Vertex(i) == vertex) {
          triangle->Update();
          break;
        }
      }
    }

    // Update
    Update();
}



void R3TriangleArray::
Update(void)
{
    // Recompute bounding box
    bbox = R3null_box;
    for (int i = 0; i < vertices.NEntries(); i++) {
      R3TriangleVertex *v = vertices.Kth(i);
      bbox.Union(v->Position());
    }
}




