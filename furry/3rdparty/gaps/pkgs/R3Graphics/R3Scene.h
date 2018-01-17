/* Include file for the R3 scene class */



/* Initialization functions */

int R3InitScene();
void R3StopScene();



/* Class definition */

class R3Scene {
public:
  // Constructor functions
  R3Scene(void);
  ~R3Scene(void);

  // Property functions
  const R3Box& BBox(void) const;
  R3Point Centroid(void) const;

  // Access functions
  int NNodes(void) const;
  R3SceneNode *Node(int k) const;
  R3SceneNode *Root(void) const;
  int NLights(void) const;
  R3Light *Light(int k) const;
  const R3Camera& Camera(void) const;
  const R2Viewport& Viewport(void) const;
  const R3Viewer& Viewer(void) const;
  const RNRgb& Ambient(void) const;
  const RNRgb& Background(void) const;

  // Manipulation functions
  void InsertNode(R3SceneNode *node);
  void RemoveNode(R3SceneNode *node);
  void InsertLight(R3Light *light);
  void RemoveLight(R3Light *light);
  void SetCamera(const R3Camera& viewer);
  void SetViewport(const R2Viewport& viewport);
  void SetViewer(const R3Viewer& viewer);
  void SetAmbient(const RNRgb& ambient);
  void SetBackground(const RNRgb& background);

  // Ray intersection functions
  RNBoolean Intersects(const R3Ray& ray, 
    R3SceneNode **hit_node = NULL, R3SceneElement **hit_element = NULL, R3Shape **hit_shape = NULL,
    R3Point *hit_point = NULL, R3Vector *hit_normal = NULL, RNScalar *hit_t = NULL) const;

  // I/O functions
  int ReadFile(const char *filename);
  int ReadObjFile(const char *filename);
  int ReadMeshFile(const char *filename);
  int ReadPrincetonFile(const char *filename);
  int ReadParseFile(const char *filename);
  int ReadSupportHierarchyFile(const char *filename);
  int ReadGrammarHierarchyFile(const char *filename);
  int ReadRectangleFile(const char *filename);
  int WriteFile(const char *filename) const;
  int WriteObjFile(const char *filename) const;
  int WritePrincetonFile(const char *filename) const;
  int WriteParseFile(const char *filename) const;
  int WriteSupportHierarchyFile(const char *filename) const;
  int WriteGrammarHierarchyFile(const char *filename) const;

  // Draw functions
  void Draw(const R3DrawFlags draw_flags = R3_DEFAULT_DRAW_FLAGS,
    RNBoolean set_camera = TRUE, RNBoolean set_lights = TRUE) const;

private:
  R3SceneNode *root;
  RNArray<R3SceneNode *> nodes;
  RNArray<R3Light *> lights;
  R3Viewer viewer;
  RNRgb ambient;
  RNRgb background;
};



/* Inline functions */

inline const R3Box& R3Scene::
BBox(void) const
{
  // Return bounding box of root node
  return root->BBox();
}



inline R3Point R3Scene::
Centroid(void) const
{
  // Return centroid of root node
  return root->Centroid();
}



inline R3SceneNode *R3Scene::
Root(void) const
{
  // Return root node
  return root;
}



inline int R3Scene::
NNodes(void) const
{
  // Return number of nodes
  return nodes.NEntries();
}



inline R3SceneNode *R3Scene::
Node(int k) const
{
  // Return kth node
  return nodes.Kth(k);
}



inline int R3Scene::
NLights(void) const
{
  // Return number of lights
  return lights.NEntries();
}



inline R3Light *R3Scene::
Light(int k) const
{
  // Return kth light
  return lights.Kth(k);
}



inline const R3Camera& R3Scene::
Camera(void) const
{
  // Return camera
  return viewer.Camera();
}


inline const R2Viewport& R3Scene::
Viewport(void) const
{
  // Return viewport
  return viewer.Viewport();
}


inline const R3Viewer& R3Scene::
Viewer(void) const
{
  // Return viewer
  return viewer;
}


inline const RNRgb& R3Scene::
Ambient(void) const
{
  // Return ambient light color
  return ambient;
}



inline const RNRgb& R3Scene::
Background(void) const
{
  // Return background color
  return background;
}



inline void R3Scene::
SetAmbient(const RNRgb& ambient) 
{
  // Set ambient light color
  this->ambient = ambient;
}



inline void R3Scene::
SetBackground(const RNRgb& background) 
{
  // Set background color
  this->background = background;
}






