/* Include file for the R3 surfel viewer class */



////////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
////////////////////////////////////////////////////////////////////////

class R3SurfelViewer {
public:
  //////////////////////////////////////////
  //// CONSTRUCTOR/DESTRUCTOR FUNCTIONS ////
  //////////////////////////////////////////

  // Constructor function
  R3SurfelViewer(R3SurfelScene *scene);

  // Destructor function
  virtual ~R3SurfelViewer(void);


  ////////////////////////////////////
  //// UI EVENT HANDLER FUNCTIONS ////
  ////////////////////////////////////

  // Call these at beginning and end of execution
  virtual void Initialize(void);
  virtual void Terminate(void);

  // Call these to handle user input events
  virtual int Redraw(void);
  virtual int Resize(int width, int height);
  virtual int MouseMotion(int x, int y);
  virtual int MouseButton(int x, int y, int button, int state, int shift, int ctrl, int alt);
  virtual int Keyboard(int x, int y, int key, int shift, int ctrl, int alt);


  ////////////////////////////////
  //// SCENE ACCESS FUNCTIONS ////
  ////////////////////////////////

  // Scene properties
  R3SurfelScene *Scene(void) const;


  ///////////////////////////////////
  //// PROPERTY ACCESS FUNCTIONS ////
  ///////////////////////////////////

  // Viewing properties
  const R3Camera& Camera(void) const;
  const R2Viewport& Viewport(void) const;
  const R3Point& CenterPoint(void) const;
  RNScalar SurfelSize(void) const;

  // Visibility properties (0=off, 1=on)
  int SurfelVisibility(void) const;
  int BackgroundVisibility(void) const;
  int ObjectLabelVisibility(void) const;
  int ObjectNameVisibility(void) const;
  int NodeBBoxVisibility(void) const;
  int BlockBBoxVisibility(void) const;
  int CenterPointVisibility(void) const;

  // Color properties
  int SurfelColorScheme(void) const;
  const RNRgb& BackgroundColor(void) const;
  const RNRgb& ObjectNameColor(void) const;
  const RNRgb& NodeBBoxColor(void) const;
  const RNRgb& BlockBBoxColor(void) const;
  const RNRgb& CenterPointColor(void) const;

  // Viewing parameters
  RNScalar TargetResolution(void) const;
  RNScalar FocusRadius(void) const;

  // Frame rate statistics
  RNScalar FrameRate(void) const;
  RNScalar FrameTime(void) const;
  RNScalar CurrentTime(void) const;


  /////////////////////////////////////////
  //// PROPERTY MANIPULATION FUNCTIONS ////
  /////////////////////////////////////////

  // Viewing property manipulation
  void ZoomCamera(RNScalar scale = 10);
  void SetCamera(const R3Camera& camera);
  void SetViewport(const R2Viewport& viewport);
  void SetCenterPoint(const R3Point& point);
  void SetSurfelSize(RNScalar npixels);

  // Visibility manipulation (0=off, 1=on, -1=toggle)
  void SetSurfelVisibility(int visibility);
  void SetBackgroundVisibility(int visibility);
  void SetObjectLabelVisibility(int visibility);
  void SetObjectNameVisibility(int visibility);
  void SetNodeBBoxVisibility(int visibility);
  void SetBlockBBoxVisibility(int visibility);
  void SetCenterPointVisibility(int visibility);

  // Color manipulation
  void SetSurfelColorScheme(int scheme);
  void SetBackgroundColor(const RNRgb& color);
  void SetObjectNameColor(const RNRgb& color);
  void SetNodeBBoxColor(const RNRgb& color);
  void SetBlockBBoxColor(const RNRgb& color);
  void SetCenterPointColor(const RNRgb& color);

  // Viewing parameters
  void SetTargetResolution(RNScalar resolution);
  void SetFocusRadius(RNScalar radius);

  // Image input/output
  int WriteImage(const char *filename);


////////////////////////////////////////////////////////////////////////
// INTERNAL STUFF BELOW HERE
////////////////////////////////////////////////////////////////////////

  // Object editing 
  int SplitLeafNodes(R3SurfelNode *source_node, const R3SurfelConstraint& constraint, 
    RNArray<R3SurfelNode *> *nodesA = NULL, RNArray<R3SurfelNode *> *nodesB = NULL);
  int SplitObject(R3SurfelObject *object, const R3SurfelConstraint& constraint,
    R3SurfelObject **objectA = NULL, R3SurfelObject **objectB = NULL);

  // Working set management
  virtual void UpdateWorkingSet(void);
  virtual void EmptyWorkingSet(void);

protected:
  // Tree manipulation
  void SetScene(R3SurfelScene *scene);

  // Viewing utility functions
  void RotateWorld(RNScalar factor, const R3Point& origin, int, int, int dx, int dy);

  // Pick utility functions
  R3SurfelNode *PickNode(int xcursor, int ycursor, 
    R3Point *hit_position = NULL, R3SurfelBlock **block = NULL, const R3Surfel **surfel = NULL,
    RNBoolean exclude_nonobjects = FALSE, RNBoolean exclude_aerial = FALSE);

  // Draw functions
  void DrawObject(R3SurfelObject *object, RNFlags flags = R3_SURFEL_DEFAULT_DRAW_FLAGS) const;

  // Management functions
  void ReadCoarsestBlocks(RNScalar max_complexity);
  void ReleaseCoarsestBlocks(RNScalar max_complexity);


protected:
  // Tree properties
  R3SurfelScene *scene;

  // Node working set
  R3SurfelNodeSet resident_nodes;

  // Viewing properties
  R3Viewer viewer;
  R3Point center_point;
  RNScalar surfel_size;

  // Visibility properties
  int surfel_visibility;
  int background_visibility;
  int object_label_visibility;
  int object_name_visibility;
  int node_bbox_visibility;
  int block_bbox_visibility;
  int center_point_visibility;

  // Color properties
  int surfel_color_scheme;
  RNRgb background_color;
  RNRgb object_name_color;
  RNRgb node_bbox_color;
  RNRgb block_bbox_color;
  RNRgb center_point_color;

  // Working set parameters
  RNBoolean adapt_working_set_automatically;
  RNScalar target_resolution;
  RNScalar focus_radius;

  // UI state
  int window_height;
  int window_width;
  int mouse_button[3];
  int mouse_position[2];
  int mouse_down_position[2];
  int mouse_drag_distance_squared;
  int shift_down;
  int ctrl_down;
  int alt_down;

  // Timing
  RNTime start_timer;
  RNTime frame_timer;
  RNScalar frame_time;

  // Image capture
  char *image_name;
};



////////////////////////////////////////////////////////////////////////
// KEY CONSTANT DEFINITIONS
////////////////////////////////////////////////////////////////////////

#define R3_SURFEL_VIEWER_ESC_KEY      27
#define R3_SURFEL_VIEWER_SPACE_KEY    32
#define R3_SURFEL_VIEWER_DEL_KEY     127

#define R3_SURFEL_VIEWER_LEFT_KEY   1024
#define R3_SURFEL_VIEWER_RIGHT_KEY  1025
#define R3_SURFEL_VIEWER_DOWN_KEY   1026
#define R3_SURFEL_VIEWER_UP_KEY     1027
#define R3_SURFEL_VIEWER_HOME_KEY   1028
#define R3_SURFEL_VIEWER_END_KEY    1029
#define R3_SURFEL_VIEWER_INSERT_KEY 1030

#define R3_SURFEL_VIEWER_F1_KEY     1064
#define R3_SURFEL_VIEWER_F2_KEY     1065
#define R3_SURFEL_VIEWER_F3_KEY     1066
#define R3_SURFEL_VIEWER_F4_KEY     1067
#define R3_SURFEL_VIEWER_F5_KEY     1068
#define R3_SURFEL_VIEWER_F6_KEY     1069
#define R3_SURFEL_VIEWER_F7_KEY     1070
#define R3_SURFEL_VIEWER_F8_KEY     1071
#define R3_SURFEL_VIEWER_F9_KEY     1072
#define R3_SURFEL_VIEWER_F10_KEY    1073
#define R3_SURFEL_VIEWER_F11_KEY    1074
#define R3_SURFEL_VIEWER_F12_KEY    1075



////////////////////////////////////////////////////////////////////////
// COLOR SCHEME CONSTANT DEFINITIONS
////////////////////////////////////////////////////////////////////////

enum {
  R3_SURFEL_VIEWER_COLOR_BY_RGB,
  R3_SURFEL_VIEWER_COLOR_BY_HEIGHT,
  R3_SURFEL_VIEWER_COLOR_BY_OBJECT,
  R3_SURFEL_VIEWER_COLOR_BY_NODE,
  R3_SURFEL_VIEWER_COLOR_BY_BLOCK,
  R3_SURFEL_VIEWER_COLOR_BY_CURRENT_LABEL,
  R3_SURFEL_VIEWER_COLOR_BY_GROUND_TRUTH_LABEL,
  R3_SURFEL_VIEWER_COLOR_BY_CONFIDENCE,
  R3_SURFEL_VIEWER_NUM_COLOR_SCHEMES
};



////////////////////////////////////////////////////////////////////////
// INLINE FUNCTIONS
////////////////////////////////////////////////////////////////////////

inline const R3Camera& R3SurfelViewer::
Camera(void) const
{
  // Return camera
  return viewer.Camera();
}



inline const R2Viewport& R3SurfelViewer::
Viewport(void) const
{
  // Return viewport
  return viewer.Viewport();
}



inline const R3Point& R3SurfelViewer::
CenterPoint(void) const
{
  // Return center point
  return center_point;
}



inline RNScalar R3SurfelViewer::
SurfelSize(void) const
{
  // Return surfel size
  return surfel_size;
}



inline int R3SurfelViewer::
SurfelVisibility(void) const
{
  // Return surfel visibililty
  return surfel_visibility;
}



inline int R3SurfelViewer::
BackgroundVisibility(void) const
{
  // Return background visibililty
  return background_visibility;
}



inline int R3SurfelViewer::
ObjectLabelVisibility(void) const
{
  // Return object label visibililty
  return object_label_visibility;
}



inline int R3SurfelViewer::
ObjectNameVisibility(void) const
{
  // Return object name visibililty
  return object_name_visibility;
}



inline int R3SurfelViewer::
NodeBBoxVisibility(void) const
{
  // Return node bbox visibililty
  return node_bbox_visibility;
}



inline int R3SurfelViewer::
BlockBBoxVisibility(void) const
{
  // Return block bbox visibililty
  return block_bbox_visibility;
}



inline int R3SurfelViewer::
CenterPointVisibility(void) const
{
  // Return center point visibililty
  return center_point_visibility;
}



inline int R3SurfelViewer::
SurfelColorScheme(void) const
{
  // Return color scheme for drawing surfels
  return surfel_color_scheme;
}



inline const RNRgb& R3SurfelViewer::
BackgroundColor(void) const
{
  // Return background color
  return background_color;
}



inline const RNRgb& R3SurfelViewer::
ObjectNameColor(void) const
{
  // Return object name color
  return object_name_color;
}



inline const RNRgb& R3SurfelViewer::
NodeBBoxColor(void) const
{
  // Return node bbox color
  return node_bbox_color;
}



inline const RNRgb& R3SurfelViewer::
BlockBBoxColor(void) const
{
  // Return block bbox color
  return block_bbox_color;
}



inline const RNRgb& R3SurfelViewer::
CenterPointColor(void) const
{
  // Return center point color
  return center_point_color;
}



inline RNScalar R3SurfelViewer::
TargetResolution(void) const
{
  // Return target resolution
  return target_resolution;
}



inline RNScalar R3SurfelViewer::
FocusRadius(void) const
{
  // Return focus radius
  return focus_radius;
}



inline RNScalar R3SurfelViewer::
FrameRate(void) const
{
  // Return frame rate
  if (frame_time == 0) return 0;
  else return 1.0 / frame_time;
}



inline RNScalar R3SurfelViewer::
FrameTime(void) const
{
  // Return number of seconds per frame refresh
  return frame_time;
}



inline RNScalar R3SurfelViewer::
CurrentTime(void) const
{
  // Return number of seconds since start up
  return start_timer.Elapsed();
}



inline R3SurfelScene *R3SurfelViewer::
Scene(void) const
{
  // Return scene
  return scene;
}



inline void R3SurfelViewer::
SetCamera(const R3Camera& camera)
{
  // Set camera
  viewer.SetCamera(camera);
}



inline void R3SurfelViewer::
SetViewport(const R2Viewport& viewport)
{
  // Set viewport
  viewer.SetViewport(viewport);
}



inline void R3SurfelViewer::
SetCenterPoint(const R3Point& point)
{
  // Set center point
  center_point = point;

  // Update working set
  UpdateWorkingSet();
}



inline void R3SurfelViewer::
SetSurfelSize(RNScalar surfel_size)
{
  // Set surfel size
  this->surfel_size = surfel_size;
}



inline void R3SurfelViewer::
SetSurfelVisibility(int visibility)
{
  // Set surfel visibililty
  if (visibility == -1) surfel_visibility = 1 - surfel_visibility;
  else if (visibility == 0) surfel_visibility = 0;
  else surfel_visibility = 1;
}



inline void R3SurfelViewer::
SetBackgroundVisibility(int visibility)
{
  // Set background visibililty
  if (visibility == -1) background_visibility = 1 - background_visibility;
  else if (visibility == 0) background_visibility = 0;
  else background_visibility = 1;
}



inline void R3SurfelViewer::
SetObjectLabelVisibility(int visibility)
{
  // Set object label visibililty
  if (visibility == -1) object_label_visibility = 1 - object_label_visibility;
  else if (visibility == 0) object_label_visibility = 0;
  else object_label_visibility = 1;
}



inline void R3SurfelViewer::
SetObjectNameVisibility(int visibility)
{
  // Set object name visibililty
  if (visibility == -1) object_name_visibility = 1 - object_name_visibility;
  else if (visibility == 0) object_name_visibility = 0;
  else object_name_visibility = 1;
}



inline void R3SurfelViewer::
SetNodeBBoxVisibility(int visibility)
{
  // Set node bbox visibililty
  if (visibility == -1) node_bbox_visibility = 1 - node_bbox_visibility;
  else if (visibility == 0) node_bbox_visibility = 0;
  else node_bbox_visibility = 1;
}



inline void R3SurfelViewer::
SetBlockBBoxVisibility(int visibility)
{
  // Set block bbox visibililty
  if (visibility == -1) block_bbox_visibility = 1 - block_bbox_visibility;
  else if (visibility == 0) block_bbox_visibility = 0;
  else block_bbox_visibility = 1;
}



inline void R3SurfelViewer::
SetCenterPointVisibility(int visibility)
{
  // Set center point visibililty
  if (visibility == -1) center_point_visibility = 1 - center_point_visibility;
  else if (visibility == 0) center_point_visibility = 0;
  else center_point_visibility = 1;
}



inline void R3SurfelViewer::
SetSurfelColorScheme(int scheme)
{
  // Set color scheme for drawing surfels
  surfel_color_scheme = scheme;
}



inline void R3SurfelViewer::
SetBackgroundColor(const RNRgb& color)
{
  // Set background color
  background_color = color;
}



inline void R3SurfelViewer::
SetObjectNameColor(const RNRgb& color)
{
  // Set object name color
  object_name_color = color;
}



inline void R3SurfelViewer::
SetNodeBBoxColor(const RNRgb& color)
{
  // Set node bbox color
  node_bbox_color = color;
}



inline void R3SurfelViewer::
SetBlockBBoxColor(const RNRgb& color)
{
  // Set block bbox color
  block_bbox_color = color;
}



inline void R3SurfelViewer::
SetCenterPointColor(const RNRgb& color)
{
  // Set center point color
  center_point_color = color;
}



inline void R3SurfelViewer::
SetTargetResolution(RNScalar target_resolution)
{
  // Set target resolution
  this->target_resolution = target_resolution;

  // Update working set
  UpdateWorkingSet();
}



inline void R3SurfelViewer::
SetFocusRadius(RNScalar focus_radius)
{
  // Set focus radius
  this->focus_radius = focus_radius;

  // Update working set
  UpdateWorkingSet();
}



