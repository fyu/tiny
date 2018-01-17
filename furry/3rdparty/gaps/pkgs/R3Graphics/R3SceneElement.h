/* Include file for the R3 scene element class */



/* Class definitions */

class R3SceneElement {
public:
  // Constructor functions
  R3SceneElement(R3Material *material = NULL);
  ~R3SceneElement(void);

  // Properties
  const R3Box& BBox(void) const;
  R3Point Centroid(void) const;

  // Access functions
  R3Material *Material(void) const;
  int NShapes(void) const;
  R3Shape *Shape(int k) const;
  R3SceneNode *Node(void) const;

  // Manipulation function
  void SetMaterial(R3Material *material);
  void InsertShape(R3Shape *shape);
  void RemoveShape(R3Shape *shape);

  // Ray intersection functions
  RNBoolean Intersects(const R3Ray& ray, R3Shape **hit_shape = NULL,
    R3Point *hit_point = NULL, R3Vector *hit_normal = NULL, RNScalar *hit_t = NULL) const;

  // Draw functions
  void Draw(const R3DrawFlags draw_flags = R3_DEFAULT_DRAW_FLAGS) const;

public:
  // Internal update functions
  void InvalidateBBox(void);
  void UpdateBBox(void);

private:
  friend class R3SceneNode;
  R3SceneNode *node;
  R3Material *material;
  RNArray<R3Shape *> shapes;
  unsigned int opengl_id;
  R3Box bbox;
};



/* Inline functions */

inline R3Point R3SceneElement::
Centroid(void) const
{
  // Return centroid
  return BBox().Centroid();
}



inline R3Material *R3SceneElement::
Material(void) const
{
  // Return material
  return material;
}



inline int R3SceneElement::
NShapes(void) const
{
  // Return number of shapes
  return shapes.NEntries();
}



inline R3Shape *R3SceneElement::
Shape(int k) const
{
  // Return shape
  return shapes[k];
}



inline const R3Box& R3SceneElement::
BBox(void) const
{
  // Return bounding box
  if (bbox[0][0] == FLT_MAX) 
    ((R3SceneElement *) this)->UpdateBBox();
  return bbox;
}



inline R3SceneNode *R3SceneElement::
Node(void) const
{
  // Return node
  return node;
}



