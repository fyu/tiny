/* Include file for the R3 scene node class */



/* Initialization functions */

int R3InitSceneNode();
void R3StopSceneNode();



/* Class definition */

class R3SceneNode {
public:
  // Constructor functions
  R3SceneNode(R3Scene *scene);
  ~R3SceneNode(void);

  // Access functions
  R3Scene *Scene(void) const;
  int SceneIndex(void) const;
  R3SceneNode *Parent(void) const;
  int ParentIndex(void) const;
  int NChildren(void) const;
  R3SceneNode *Child(int k) const;
  int NElements(void) const;
  R3SceneElement *Element(int k) const;
  const R3Affine& Transformation(void) const;
  const R3Box& BBox(void) const;
  R3Point Centroid(void) const;
  const char *Name(void) const;
  void *Data(void) const;

  // Manipulation functions
  void InsertChild(R3SceneNode *node);
  void RemoveChild(R3SceneNode *node);
  void InsertElement(R3SceneElement *element);
  void RemoveElement(R3SceneElement *element);
  void SetTransformation(const R3Affine& transformation);
  void SetName(const char *name);
  void SetData(void *data);

  // Ray intersection functions
  RNBoolean Intersects(const R3Ray& ray, 
    R3SceneNode **hit_node = NULL, R3SceneElement **hit_element = NULL, R3Shape **hit_shape = NULL,
    R3Point *hit_point = NULL, R3Vector *hit_normal = NULL, RNScalar *hit_t = NULL) const;

  // Draw functions
  void Draw(const R3DrawFlags draw_flags = R3_DEFAULT_DRAW_FLAGS) const;

public:
  // Internal update functions
  void InvalidateBBox(void);
  void UpdateBBox(void);

private:
  friend class R3Scene;
  R3Scene *scene;
  int scene_index;
  R3SceneNode *parent;
  int parent_index;
  RNArray<R3SceneNode *> children;
  RNArray<R3SceneElement *> elements;
  R3Affine transformation;
  R3Box bbox;
  char *name;
  void *data;
};



/* Inline functions */

inline R3Scene *R3SceneNode::
Scene(void) const
{
  // Return scene 
  return scene;
}



inline int R3SceneNode::
SceneIndex(void) const
{
  // Return index of node in scene (can be used with scene->Node(xxx))
  return scene_index;
}



inline R3SceneNode *R3SceneNode::
Parent(void) const
{
  // Return parent node
  return parent;
}



inline int R3SceneNode::
ParentIndex(void) const
{
  // Return index of node in parent (can be used with parent->Child(xxx))
  return parent_index;
}



inline int R3SceneNode::
NChildren(void) const
{
  // Return number of children nodes
  return children.NEntries();
}



inline R3SceneNode *R3SceneNode::
Child(int k) const
{
  // Return kth child node
  return children.Kth(k);
}



inline int R3SceneNode::
NElements(void) const
{
  // Return number of elements
  return elements.NEntries();
}



inline R3SceneElement *R3SceneNode::
Element(int k) const
{
  // Return kth element
  return elements.Kth(k);
}



inline const R3Affine& R3SceneNode::
Transformation(void) const
{
  // Return transformation
  return transformation;
}



inline const R3Box& R3SceneNode::
BBox(void) const
{
  // Return bounding box
  if (bbox[0][0] == FLT_MAX) 
    ((R3SceneNode *) this)->UpdateBBox();
  return bbox;
}



inline R3Point R3SceneNode::
Centroid(void) const
{
  // Return centroid of bounding box
  return BBox().Centroid();
}



inline const char *R3SceneNode::
Name(void) const
{
  // Return name
  return name;
}



inline void *R3SceneNode::
Data(void) const
{
  // Return user-defined data
  return data;
}

