/* Include file for the Catmull Rom spline class */



/* Class definition */

class R3CatmullRomSpline {
public:
  // Constructor functions
  R3CatmullRomSpline(void);
  R3CatmullRomSpline(const R3CatmullRomSpline& curve);
  R3CatmullRomSpline(const RNArray<R3Point *>& points, RNScalar tao = 0.5);
  R3CatmullRomSpline(const R3Point *points, int npoints, RNScalar tao = 0.5);

  // Spline properties
  const int NVertices(void) const;
  RNScalar StartParameter(void) const;
  RNScalar EndParameter(void) const;
  const R3Point& StartPosition(void) const;
  const R3Point& EndPosition(void) const;
  const R3Box& BBox(void) const;
  int Order(void) const;

  // Vertex properties
  R3Point VertexPosition(int k) const;
  R3Vector VertexDirection(int k) const;
  R3Vector VertexDerivative(int k) const;
  void *VertexData(int k) const;

  // Point properties
  int VertexIndex(RNScalar u) const;
  R3Point PointPosition(RNScalar u) const;
  R3Vector PointDirection(RNScalar u) const;
  R3Vector PointDerivative(RNScalar u) const;
  
  // Manipulation functions
  void SetVertexPosition(int k, const R3Point& position);
  void SetVertexData(int k, void *data);

  // Draw functions
  void Draw(RNScalar sample_spacing = 0.1) const;

public:
  // Utility functions
  RNScalar BlendingWeight(RNScalar t, int k) const;
  RNScalar BlendingWeightDerivative(RNScalar t, int k) const;
  R3Point PhantomVertexPosition(int i, int j) const;

private:
  // Upkeep functions
  virtual void UpdateBBox(void);

private:
  R3Point *vertex_positions;
  void **vertex_datas;
  int nvertices;
  RNScalar tao;
  R3Box bbox;
};



/* Inline functions */

inline const int R3CatmullRomSpline::
NVertices(void) const
{
  // Return number of vertices
  return nvertices;
}



inline RNScalar R3CatmullRomSpline::
StartParameter(void) const
{
  // Return parameter at start of spline
  return 0;
}



inline RNScalar R3CatmullRomSpline::
EndParameter(void) const
{
  // Return parameter at end of spline
  return nvertices-1;
}



inline const R3Point& R3CatmullRomSpline::
StartPosition(void) const
{
  // Return position at start of spline
  if (nvertices == 0) return R3zero_point;
  return vertex_positions[0];
}



inline const R3Point& R3CatmullRomSpline::
EndPosition(void) const
{
  // Return position at end of spline
  if (nvertices == 0) return R3zero_point;
  return vertex_positions[nvertices-1];
}



inline const R3Box& R3CatmullRomSpline::
BBox(void) const
{
  // Return bounding box
  if (bbox[0][0] == FLT_MAX) 
    ((R3CatmullRomSpline *) this)->UpdateBBox();
  return bbox;
}



inline int R3CatmullRomSpline::
Order(void) const
{
  // Return order of polynomial
  return 3;
}



// Vertex properties
inline R3Point R3CatmullRomSpline::
VertexPosition(int i) const
{
  // Return position of ith vertex
  assert((i >= 0) && (i < nvertices));
  return vertex_positions[i];
}



inline R3Vector R3CatmullRomSpline::
VertexDerivative(int i) const
{
  // Return derivative of curve at ith vertex
  if (nvertices < 2) return R3zero_vector;
  if (i == nvertices-1) return Order() * (VertexPosition(nvertices-1) - PhantomVertexPosition(nvertices-2, 1));
  else return Order() * (PhantomVertexPosition(i, 0) - VertexPosition(i));
}



inline R3Vector R3CatmullRomSpline::
VertexDirection(int i) const
{
  // Return tangent direction of curve at ith vertex
  R3Vector direction = VertexDerivative(i);
  direction.Normalize();
  return direction;
}



inline void *R3CatmullRomSpline::
VertexData(int i) const
{
  // Return user data associated with vertex
  if (!vertex_datas) return NULL;
  else return vertex_datas[i];
}



inline int R3CatmullRomSpline::
VertexIndex(RNScalar u) const
{
  // Return index at start of segment containing u
  if (u < 0) return 0;
  else if (u >= nvertices-2) return nvertices-2;
  else return (int) u;
}



inline R3Vector R3CatmullRomSpline::
PointDirection(RNScalar u) const
{
  // Return tangent direction of curve at parametric value u
  R3Vector direction = PointDerivative(u);
  direction.Normalize();
  return direction;
}



