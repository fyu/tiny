/* Source file for the surfel scene viewer class */



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "R3Graphics/R3Graphics.h"
#include "R3Surfels/R3Surfels.h"
#include "R3SurfelViewer.h"



////////////////////////////////////////////////////////////////////////
// Surfel viewer constructor/destructor
////////////////////////////////////////////////////////////////////////

R3SurfelViewer::
R3SurfelViewer(R3SurfelScene *scene)
  : scene(NULL),
    resident_nodes(),
    viewer(),
    center_point(0,0,0),
    surfel_size(2),
    surfel_visibility(1),
    background_visibility(1),
    object_label_visibility(0),
    object_name_visibility(0),
    node_bbox_visibility(0),
    block_bbox_visibility(0),
    center_point_visibility(0),
    surfel_color_scheme(R3_SURFEL_VIEWER_COLOR_BY_RGB),
    background_color(0,0,0),
    node_bbox_color(0,0,1),
    block_bbox_color(0,1,0),
    center_point_color(1,0,0),
    adapt_working_set_automatically(0),
    target_resolution(10),
    focus_radius(0),
    window_height(0),
    window_width(0),
    shift_down(0),
    ctrl_down(0),
    alt_down(0),
    start_timer(),
    frame_timer(),
    frame_time(-1),
    image_name(NULL)
{
  // Initialize mouse button state
  mouse_button[0] = 0;
  mouse_button[1] = 0;
  mouse_button[2] = 0;

  // Initialize mouse positions
  mouse_position[0] = 0;
  mouse_position[1] = 0;

  // Initialize mouse positions
  mouse_down_position[0] = 0;
  mouse_down_position[1] = 0;

  // Initialize mouse drag distance
  mouse_drag_distance_squared = 0;

  // Initialize timers
  start_timer.Read();
  frame_timer.Read();

  // Set the scene
  if (scene) SetScene(scene);
}



R3SurfelViewer::
~R3SurfelViewer(void)
{
}



////////////////////////////////////////////////////////////////////////
// Text drawing utility functions
////////////////////////////////////////////////////////////////////////

#include "fglut/fglut.h"

static void 
DrawText(const R2Point& p, const char *s, void *font = GLUT_BITMAP_HELVETICA_12)
{
  // Draw text string s and position p
  glRasterPos2d(p[0], p[1]);
  while (*s) glutBitmapCharacter(font, *(s++));
}
  


static void 
DrawText(const R3Point& p, const char *s, void *font = GLUT_BITMAP_HELVETICA_12)
{
  // Draw text string s and position p
  glRasterPos3d(p[0], p[1], p[2]);
  while (*s) glutBitmapCharacter(font, *(s++));
}
  

////////////////////////////////////////////////////////////////////////
// Coloring utility functions
////////////////////////////////////////////////////////////////////////

void
LoadColor(int k)
{
  // Make array of colors
  const int ncolors = 18;
  const RNRgb colors[ncolors] = {
    RNRgb(1, 0, 0), RNRgb(0, 1, 0), RNRgb(0, 0, 1), 
    RNRgb(1, 1, 0), RNRgb(0, 1, 1), RNRgb(1, 0, 1), 
    RNRgb(1, 0.5, 0), RNRgb(0, 1, 0.5), RNRgb(0.5, 0, 1), 
    RNRgb(0.5, 1, 0), RNRgb(0, 0.5, 1), RNRgb(1, 0, 0.5), 
    RNRgb(0.5, 0, 0), RNRgb(0, 0.5, 0), RNRgb(0, 0, 0.5), 
    RNRgb(0.5, 0.5, 0), RNRgb(0, 0.5, 0.5), RNRgb(0.5, 0, 0.5) 
  };

  // Return color
  if (k == 0) RNLoadRgb(colors[0]);
  else RNLoadRgb(colors[1 + (k % (ncolors-1))]);
}



void
LoadColor(double value)
{
  // Compute rgb
  GLdouble r, g, b;
  if (value < 0) {
    r = 0;
    g = 0;
    b = 1;
  }
  else if (value < 0.1) {
    value *= 10;
    r = 0;
    g = value;
    b = 1;
  }
  else if (value < 0.5) {
    value = (value - 0.1) * 2.5;
    r = 0;
    g = 1;
    b = 1 - value;
  }
  else if (value < 0.9) {
    value = (value - 0.5) * 2.5;
    r = value;
    g = 1;
    b = 0;
  }
  else if (value < 1) {
    value = (value - 0.9) * 10;
    r = 1;
    g = 1 - value;
    b = 0;
  }
  else {
    r = 1;
    g = 0;
    b = 0;
  }

  // Load rgb
  glColor3d(r, g, b);
}



////////////////////////////////////////////////////////////////////////
// UI event handler functions
////////////////////////////////////////////////////////////////////////

int R3SurfelViewer::
Redraw(void)
{
  // Check scene
  if (!scene) return 0;

  // Set viewing transformation
  viewer.Camera().Load();

  // Clear window 
  glClearColor(background_color[0], background_color[1], background_color[2], 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Set draw modes
  glDisable(GL_LIGHTING);
  glPointSize(surfel_size);
  glLineWidth(1);

  // Draw surfels
  if (surfel_visibility) {
#if 0
    // Draw labeled objects
    if (surfel_color_scheme == R3_SURFEL_VIEWER_COLOR_BY_CURRENT_LABEL) {
      glPushMatrix();
      R3Vector offset = -0.02 * viewer.Camera().Towards();
      glTranslated(offset[0], offset[1], offset[2]);
      glPointSize(2 * surfel_size);
      for (int i = 0; i < scene->NObjects(); i++) {
        R3SurfelObject *object = scene->Object(i);
        R3SurfelLabel *label = object->HumanLabel();
        if (!label) label = object->PredictedLabel();
        if (!label) continue;
        RNLoadRgb(label->Color());
        DrawObject(object, 0);
      }
      glPointSize(surfel_size);
      glPopMatrix();
    }
#endif

    // Draw resident nodes
    if (surfel_color_scheme == R3_SURFEL_VIEWER_COLOR_BY_CURRENT_LABEL) {
      for (int i = 0; i < resident_nodes.NNodes(); i++) {
        R3SurfelNode *node = resident_nodes.Node(i);
        R3SurfelObject *object = node->Object();
        if (!object) continue;
        R3SurfelLabel *label = object->CurrentLabel();
        int label_index = (label) ? label->SceneIndex() : 0;
        LoadColor(label_index);
        RNLength distance = R3Distance(viewer.Camera().Origin(), node->BBox());
        RNScalar point_size = (distance > 0) ? 10.0 * surfel_size / distance : surfel_size;
        glPointSize(point_size);
        node->Draw(0); 
      }
    }
    else if (surfel_color_scheme == R3_SURFEL_VIEWER_COLOR_BY_GROUND_TRUTH_LABEL) {
      for (int i = 0; i < resident_nodes.NNodes(); i++) {
        R3SurfelNode *node = resident_nodes.Node(i);
        R3SurfelObject *object = node->Object();
        if (!object) continue;
        R3SurfelLabel *label = object->GroundTruthLabel();
        int label_index = (label) ? label->SceneIndex() : 0;
        LoadColor(label_index);
        RNLength distance = R3Distance(viewer.Camera().Origin(), node->BBox());
        RNScalar point_size = (distance > 0) ? 10.0 * surfel_size / distance : surfel_size;
        glPointSize(point_size);
        node->Draw(0); 
      }
    }
    else if (surfel_color_scheme == R3_SURFEL_VIEWER_COLOR_BY_OBJECT) {
      // Draw with colors based on nodes
      for (int i = 0; i < resident_nodes.NNodes(); i++) {
        R3SurfelNode *node = resident_nodes.Node(i);
        R3SurfelObject *object = node->Object();
        int object_index = (object) ? object->SceneIndex() : 0;
        LoadColor(object_index);
        RNLength distance = R3Distance(viewer.Camera().Origin(), node->BBox());
        RNScalar point_size = (distance > 0) ? 10.0 * surfel_size / distance : surfel_size;
        glPointSize(point_size);
        node->Draw(0); 
      }
    }
    else if (surfel_color_scheme == R3_SURFEL_VIEWER_COLOR_BY_NODE) {
      // Draw with colors based on nodes
      for (int i = 0; i < resident_nodes.NNodes(); i++) {
        R3SurfelNode *node = resident_nodes.Node(i);
        LoadColor(i);
        RNLength distance = R3Distance(viewer.Camera().Origin(), node->BBox());
        RNScalar point_size = (distance > 0) ? 10.0 * surfel_size / distance : surfel_size;
        glPointSize(point_size);
        node->Draw(0); 
      }
    }
    else if (surfel_color_scheme == R3_SURFEL_VIEWER_COLOR_BY_BLOCK) {
      // Draw with colors based on blocks
      int count = 0;
      for (int i = 0; i < resident_nodes.NNodes(); i++) {
        R3SurfelNode *node = resident_nodes.Node(i);
        for (int j = 0; j < node->NBlocks(); j++) {
          R3SurfelBlock *block = node->Block(j);
          LoadColor(count++);
          RNLength distance = R3Distance(viewer.Camera().Origin(), block->BBox());
          RNScalar point_size = (distance > 0) ? 10.0 * surfel_size / distance : surfel_size;
          glPointSize(point_size);
          block->Draw(0);
        }
      }
    }
    else if (surfel_color_scheme == R3_SURFEL_VIEWER_COLOR_BY_HEIGHT) {
      // Draw with colors based on heights
      for (int i = 0; i < resident_nodes.NNodes(); i++) {
        R3SurfelNode *node = resident_nodes.Node(i);
        RNLength distance = R3Distance(viewer.Camera().Origin(), node->BBox());
        RNScalar point_size = (distance > 0) ? 10.0 * surfel_size / distance : surfel_size;
        glPointSize(point_size);
        glBegin(GL_POINTS);
        for (int j = 0; j < node->NBlocks(); j++) {
          R3SurfelBlock *block = node->Block(j);
          const R3Point& block_origin = block->Origin();
          for (int k = 0; k < block->NSurfels(); k++) {
            const R3Surfel *surfel = block->Surfel(k);
            double x = surfel->X() + block_origin.X();
            double y = surfel->Y() + block_origin.Y();
            double z = surfel->Z() + block_origin.Z();
            double value = 0.1 * (z - center_point.Z());
            LoadColor(value);
            glVertex3d(x, y, z);
          }
        }
        glEnd();
      }
    }
    else {
      // Draw with regular surfel colors
      for (int i = 0; i < resident_nodes.NNodes(); i++) {
        R3SurfelNode *node = resident_nodes.Node(i);
        RNLength distance = R3Distance(viewer.Camera().Origin(), node->BBox());
        RNScalar point_size = (distance > 0) ? 10.0 * surfel_size / distance : surfel_size;
        glPointSize(point_size);
        node->Draw(); 
      }
    }
  }

  // Reset the point size
  glPointSize(1);

  // Draw object labels
  if (object_label_visibility) {
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, viewer.Viewport().Width(), 0, viewer.Viewport().Height());
    glDisable(GL_DEPTH_TEST);
    glDepthMask(FALSE);
    for (int i = 0; i < scene->NObjects(); i++) {
      R3SurfelObject *object = scene->Object(i);
      for (int i = 0; i < object->NLabelAssignments(); i++) {
        R3SurfelLabelAssignment *assignment = object->LabelAssignment(i);
        if (assignment->Originator() == R3_SURFEL_LABEL_ASSIGNMENT_GROUND_TRUTH_ORIGINATOR) continue;
        RNBoolean confirmed = (assignment->Originator() == R3_SURFEL_LABEL_ASSIGNMENT_HUMAN_ORIGINATOR) ? 1 : 0;
        R3SurfelLabel *label = assignment->Label();
        R3Point position = object->Centroid();
        position[2] = object->BBox().ZMax() + 1;
        R2Point p = viewer.ViewportPoint(position);
        void *font = (confirmed) ? GLUT_BITMAP_HELVETICA_18 : GLUT_BITMAP_HELVETICA_12;
        int width = glutBitmapLength(font, (const unsigned char *) label->Name());
        p[0] -= width / 2;
        RNLoadRgb(label->Color());
        DrawText(p, label->Name(), font);
        break;
      }
    }
    glDepthMask(TRUE);
    glEnable(GL_DEPTH_TEST);
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
  }

  // Draw object names
  if (object_name_visibility) {
    RNLoadRgb(object_name_color);
    for (int i = 0; i < scene->NObjects(); i++) {
      R3SurfelObject *object = scene->Object(i);
      if (object->NParts() > 0) continue;
      if (!object->Name()) continue;
      R3Point position = object->Centroid();
      position[2] = object->BBox().ZMax() + 1;
      DrawText(position, object->Name(), GLUT_BITMAP_HELVETICA_12);
    }
  }

  // Draw node bounding boxes
  if (node_bbox_visibility) {
    RNLoadRgb(node_bbox_color);
    for (int i = 0; i < resident_nodes.NNodes(); i++) {
      R3SurfelNode *node = resident_nodes.Node(i);
      node->BBox().Outline();
      if (node->NParts() > 0) glColor3d(0, 1, 0);
      else glColor3d(1, 0, 0);
    }
  }

  // Draw block bounding boxes
  if (block_bbox_visibility) {
    RNLoadRgb(block_bbox_color);
    for (int i = 0; i < resident_nodes.NNodes(); i++) {
      R3SurfelNode *node = resident_nodes.Node(i);
      for (int j = 0; j < node->NBlocks(); j++) {
        R3SurfelBlock *block = node->Block(j);
        block->BBox().Outline();
      }
    }
  }

  // Draw center point
  if (center_point_visibility) {
    glEnable(GL_LIGHTING);
    GLfloat color[4] = { center_point_color[0], center_point_color[1], center_point_color[2], 1 }; 
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, color);
    R3Sphere(center_point, 1).Draw();
    glDisable(GL_LIGHTING);
  }

  // Capture image and exit
  if (image_name) {
    R2Image image(viewer.Viewport().Width(), viewer.Viewport().Height(), 3);
    image.Capture();
    image.Write(image_name);
    free(image_name);
    image_name = NULL;
  }

  // Update the frame time
  if (frame_time < 0) frame_time = 0.05;
  else frame_time = frame_timer.Elapsed();
  frame_timer.Read();

  // Adapt working set
  if (adapt_working_set_automatically) {
    // Adjust target resolution based on frame time
    if (frame_time > 0.05) {
      SetTargetResolution(0.9 * TargetResolution());
    }
    else if (frame_time < 0.025) {
      SetTargetResolution(1.1 * TargetResolution());
    }

    // Make gross estimate of visible radius
    RNLength camera_height = viewer.Camera().Origin().Z();
    RNLength visible_radius = camera_height;

    // Adjust focus radius 
    SetFocusRadius(visible_radius);

    // Adjust surfel size based on visible radius
    if ((TargetResolution() > 0) && (visible_radius > 0)) {
      RNLength window_width = viewer.Viewport().Width();
      RNLength npixels = window_width / (visible_radius * TargetResolution());
      SetSurfelSize(npixels);
    }
  }

  // Return whether need redraw
  return 0;
}    



int R3SurfelViewer::
Resize(int w, int h)
{
  // Resize window
  glViewport(0, 0, w, h);

  // Resize viewer viewport
  viewer.ResizeViewport(0, 0, w, h);

  // Remember window size
  window_width = w;
  window_height = h;

  // Return whether need redraw
  return 1;
}



int R3SurfelViewer::
MouseMotion(int x, int y)
{
  // Initialize
  int redraw = 0;

  // Compute mouse movement
  int dx = x - mouse_position[0];
  int dy = y - mouse_position[1];

  // Set viewing center point
  R3Point viewing_center_point = center_point;
  const R3Camera& camera = viewer.Camera();
  R3Plane camera_plane(camera.Origin(), camera.Towards());
  RNScalar signed_distance = R3SignedDistance(camera_plane, viewing_center_point);
  if (signed_distance < 0) viewing_center_point -= (signed_distance - 1) * camera.Towards();
  
  // World in hand navigation 
  if (shift_down && mouse_button[0]) viewer.ScaleWorld(2.0, viewing_center_point, x, y, dx, dy);
  else if (ctrl_down && mouse_button[0]) viewer.TranslateWorld(2.0, viewing_center_point, x, y, dx, dy);
  else if (mouse_button[0]) RotateWorld(1.0, viewing_center_point, x, y, dx, dy);
  else if (mouse_button[1]) viewer.ScaleWorld(2.0, viewing_center_point, x, y, dx, dy);
  else if (mouse_button[2]) viewer.TranslateWorld(2.0, viewing_center_point, x, y, dx, dy);
  if (mouse_button[0] || mouse_button[1] || mouse_button[2]) redraw = 1;

  // Remember mouse position 
  mouse_position[0] = x;
  mouse_position[1] = y;

  // Update mouse drag movement
  mouse_drag_distance_squared += dx*dx + dy*dy;

  // Return whether need redraw
  return redraw;
}



int R3SurfelViewer::
MouseButton(int x, int y, int button, int state, int shift, int ctrl, int alt)
{
  // Initialize
  int redraw = 0;

  // Process mouse button event
  if (state == 1) {
    // Button is going down
    mouse_drag_distance_squared = 0;

    // Remember mouse down position 
    mouse_down_position[0] = x;
    mouse_down_position[1] = y;
  }
  else {
    // Button is going up
  }

  // Remember mouse position 
  mouse_position[0] = x;
  mouse_position[1] = y;

  // Remember button state 
  mouse_button[button] = state;

  // Remember modifiers 
  shift_down = shift;
  ctrl_down = ctrl;
  alt_down = alt;

  // Return whether need redraw
  return redraw;
}



int R3SurfelViewer::
Keyboard(int x, int y, int key, int shift, int ctrl, int alt)
{
  // Initialize redraw status
  int redraw = 1;

  // Process debugging commands
  if (alt) {
    // Process debugging commands
    switch (key) {
    case 'B':
      SetBlockBBoxVisibility(-1);
      break;

    case 'b':
      SetNodeBBoxVisibility(-1);
      break;

    case 'C':
    case 'c':
      surfel_color_scheme = (surfel_color_scheme + 1) % R3_SURFEL_VIEWER_NUM_COLOR_SCHEMES;
      break;
      
    case 'L':
    case 'l':
      SetObjectLabelVisibility(-1);
      break;
      
    case 'N':
    case 'n':
      SetObjectNameVisibility(-1);
      break;
      
    case 'O':
    case 'o':
      SetCenterPointVisibility(-1);
      break;
      
    case 'P':
    case 'p':
      SetSurfelVisibility(-1);
      break;

    case 'W':
    case 'w': {
      R3Point pick_position;
      R3SurfelNode *node = PickNode(x, y, &pick_position);
      if (node) SetCenterPoint(pick_position);
      break; }

    case 'X':
    case 'x':
      SetBackgroundVisibility(-1);
      break;
      
    case 'Q': 
    case 'q': {
      R3Point pick_position(0,0,0);
      R3SurfelNode *node = PickNode(x, y, &pick_position);
      if (node) {
        SetCenterPoint(pick_position);
        printf("%g %g %g\n", pick_position[0], pick_position[1], pick_position[2]);
        while (node) {
          const char *node_name = (node->Name()) ? node->Name() : "-";
          R3SurfelObject *object = node->Object();
          const char *object_name = (object && (object->Name())) ? object->Name() : "-";
          char object_index[128];
          if (object) sprintf(object_index, "%d", object->SceneIndex());
          else sprintf(object_index, "%s", "-");
          R3SurfelLabel *ground_truth_label = (object) ? object->GroundTruthLabel() : NULL;
          const char *ground_truth_label_name = (ground_truth_label && (ground_truth_label->Name())) ? ground_truth_label->Name() : "-";
          R3SurfelLabel *current_label = (object) ? object->CurrentLabel() : NULL;
          const char *current_label_name = (current_label && (current_label->Name())) ? current_label->Name() : "-";
          printf("  %4d %4d %-30s  :  %-6s %-30s  :  %-30s %-30s\n",  
                 node->TreeLevel(), node->NParts(), node_name, 
                 object_index, object_name, 
                 current_label_name, ground_truth_label_name);
          node = node->Parent();        
        }
      }
      break; }

    default:
      redraw = 0;
      break;
    }
  }
  else if (ctrl) {
    switch(key) {
    default:
      redraw = 0;
      break;
    }
  }
  else {
    // Process other keyboard events
    switch (key) {
    case R3_SURFEL_VIEWER_F1_KEY:
    case R3_SURFEL_VIEWER_F2_KEY:
    case R3_SURFEL_VIEWER_F3_KEY:
    case R3_SURFEL_VIEWER_F4_KEY: 
      ZoomCamera(pow(0.5, 2.0*(key - R3_SURFEL_VIEWER_F1_KEY)));
      break; 

    case R3_SURFEL_VIEWER_UP_KEY:
      SetTargetResolution(1.5 * TargetResolution());
      break;

    case R3_SURFEL_VIEWER_DOWN_KEY:
      SetTargetResolution(0.67 * TargetResolution());
      break;

    case R3_SURFEL_VIEWER_RIGHT_KEY:
      SetFocusRadius(1.1 * FocusRadius());
      break;

    case R3_SURFEL_VIEWER_LEFT_KEY:
      SetFocusRadius(0.9 * FocusRadius());
      break;

    case '-': 
      SetSurfelSize(0.9 * SurfelSize());
      break;

    case '+': 
      SetSurfelSize(1.1 * SurfelSize());
      break;

    default:
      redraw = 0;
      break;
    }
  }

  // Remember mouse position 
  mouse_position[0] = x;
  mouse_position[1] = y;

  // Remember modifiers 
  shift_down = shift;
  ctrl_down = ctrl;
  alt_down = alt;

  // Return whether need redraw
  return redraw;
}



////////////////////////////////////////////////////////////////////////
// COMMANDS
////////////////////////////////////////////////////////////////////////

void R3SurfelViewer::
Initialize(void)
{
  // Initialize lights
  static GLfloat lmodel_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
  static GLfloat light0_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
  glEnable(GL_LIGHT0);
  static GLfloat light1_diffuse[] = { 0.5, 0.5, 0.5, 1.0 };
  glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
  glEnable(GL_LIGHT1);
  glEnable(GL_NORMALIZE);

  // Initialize graphics modes
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);
}



void R3SurfelViewer::
Terminate(void)
{
}



void R3SurfelViewer::
ZoomCamera(RNScalar scale)
{
  // Zoom into center point
  // scale=1 covers the whole scene
  // scale=0 is zoomed in to a 10m radius

  // Set the camera viewpoint
  R3Box bbox = scene->BBox();
  RNLength r = 10 + scale * bbox.DiagonalRadius();
  R3Point eye = center_point - 2 * r * viewer.Camera().Towards();
  viewer.RepositionCamera(eye);
}



int R3SurfelViewer::
WriteImage(const char *filename)
{
  // Check if can write file
  FILE *fp = fopen(filename, "w");
  if (!fp) return 0;
  else fclose(fp);

  // Remember image name -- capture image next redraw
  if (image_name) free(image_name);
  if (!filename) image_name = NULL;
  else image_name = strdup(filename);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Scene manipulation utility functions
////////////////////////////////////////////////////////////////////////

void R3SurfelViewer::
SetScene(R3SurfelScene *scene)
{
  // Remember scene
  this->scene = scene;

  // Set center point
  center_point = scene->Centroid();
  center_point[2] = scene->BBox().ZMin();

  // Set focus radius
  focus_radius = 400;
  if (focus_radius > scene->BBox().DiagonalRadius()) {
    focus_radius = scene->BBox().DiagonalRadius();
  }

  // Set camera and viewport
  R3Box bbox = scene->BBox();
  RNLength r = bbox.DiagonalRadius();
  static const R3Vector up(0, 1, 0);
  static const R3Vector towards(0, 0, -1);
  R3Point eye = scene->Centroid() - towards * (2 * r); 
  R3Camera camera(eye, towards, up, 0.4, 0.4, 0.01, 100000.0);
  R2Viewport viewport(0, 0, window_width, window_height);
  viewer.SetViewport(viewport);
  viewer.SetCamera(camera);

  // Lock coarsest blocks in memory (~500MB)
  // ReadCoarsestBlocks(32 * 1024 * 1024);

  // Update working set
  UpdateWorkingSet();
}



////////////////////////////////////////////////////////////////////////
// Working set utility functions
////////////////////////////////////////////////////////////////////////

void R3SurfelViewer::
ReadCoarsestBlocks(RNScalar max_complexity)
{
  // Just checking
  if (!scene) return;

  // Get convenient variables
  R3SurfelTree *tree = scene->Tree();
  if (!tree) return;

  // Seed breadth first search with root nodes
  RNQueue<R3SurfelNode *> queue;
  for (int i = 0; i < tree->NNodes(); i++) {
    R3SurfelNode *node = tree->Node(i);
    queue.Insert(node);
  }

  // Visit nodes in breadth first search reading blocks
  RNScalar total_complexity = 0;
  while (!queue.IsEmpty()) {
    R3SurfelNode *node = queue.Pop();
    if (total_complexity + node->Complexity() > max_complexity) break;
    total_complexity += node->Complexity();
    node->ReadBlocks();
    for (int i = 0; i < node->NParts(); i++) {
      R3SurfelNode *part = node->Part(i);
      queue.Push(part);
    }
  }
}



void R3SurfelViewer::
ReleaseCoarsestBlocks(RNScalar max_complexity)
{
  // Just checking
  if (!scene) return;

  // Get convenient variables
  R3SurfelTree *tree = scene->Tree();
  if (!tree) return;

  // Seed breadth first search with root nodes
  RNQueue<R3SurfelNode *> queue;
  for (int i = 0; i < tree->NNodes(); i++) {
    R3SurfelNode *node = tree->Node(i);
    queue.Insert(node);
  }

  // Visit nodes in breadth first search reading blocks
  RNScalar total_complexity = 0;
  while (!queue.IsEmpty()) {
    R3SurfelNode *node = queue.Pop();
    if (total_complexity + node->Complexity() > max_complexity) break;
    total_complexity += node->Complexity();
    node->ReleaseBlocks();
    for (int i = 0; i < node->NParts(); i++) {
      R3SurfelNode *part = node->Part(i);
      queue.Push(part);
    }
  }
}



void R3SurfelViewer::
UpdateWorkingSet(void)
{
  // Just checking
  if (!scene) return;

  // Get convenient variables
  R3SurfelTree *tree = scene->Tree();
  if (!tree) return;

  // Find new set of nodes
  R3SurfelNodeSet new_resident_nodes;
  new_resident_nodes.InsertNodes(tree, center_point, focus_radius, -FLT_MAX, FLT_MAX, target_resolution, RN_EPSILON);

  // Read new working set
  new_resident_nodes.ReadBlocks();

  // Release old working set
  resident_nodes.ReleaseBlocks();

  // Now use newnodes 
  resident_nodes = new_resident_nodes;
}



void R3SurfelViewer::
EmptyWorkingSet(void)
{
  // Just checking
  if (!scene) return;

  // Release blocks from resident nodes
  resident_nodes.ReleaseBlocks();

  // Empty resident nodes
  resident_nodes.Empty();
}



void R3SurfelViewer::
RotateWorld(RNScalar factor, const R3Point& origin, int, int, int dx, int dy)
{
  // Rotate world based on mouse (dx)
  if ((dx == 0) && (dy == 0)) return;
  RNLength vx = (RNLength) dx / (RNLength) viewer.Viewport().Width();
  RNLength vy = (RNLength) dy / (RNLength) viewer.Viewport().Height();
  RNAngle theta = -1 * factor * 4.0 * vx;
  viewer.RotateWorld(origin, R3posz_vector, theta);
  RNAngle phi = factor * 4.0 * vy;
  RNAngle max_phi = R3InteriorAngle(viewer.Camera().Towards(), R3posz_vector) - RN_PI_OVER_TWO;
  RNAngle min_phi = -1.0 * R3InteriorAngle(viewer.Camera().Towards(), R3negz_vector);
  if (phi < min_phi) phi = min_phi;
  if (phi > max_phi) phi = max_phi;
  viewer.RotateWorld(origin, viewer.Camera().Right(), phi);
}



R3SurfelNode *R3SurfelViewer::
PickNode(int x, int y, R3Point *picked_position, 
  R3SurfelBlock **picked_block, const R3Surfel **picked_surfel,
  RNBoolean exclude_nonobjects, RNBoolean exclude_aerial) 
{
  // Initialize result
  if (picked_position) *picked_position = R3zero_point;
  if (picked_block) picked_block = NULL;
  if (picked_surfel) picked_surfel = NULL;

  // Check cursor position
  R2Point cursor_position(x,y);
  if (!R2Contains(viewer.Viewport().BBox(), cursor_position)) {
    return NULL;
  }

  // Allocate select buffer
  const int SELECT_BUFFER_SIZE = 1024;
  GLuint select_buffer[SELECT_BUFFER_SIZE];
  GLint select_buffer_hits;

  // Initialize select buffer
  glSelectBuffer(SELECT_BUFFER_SIZE, select_buffer);
  glRenderMode(GL_SELECT);
  glInitNames();
  glPushName(0);

  // Draw surfels with pick names into selection buffer
  GLint viewport[4];
  glViewport(0, 0, Viewport().Width(), Viewport().Height());
  glGetIntegerv(GL_VIEWPORT, viewport);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPickMatrix((GLdouble) x, (GLdouble) y, 16, 16, viewport);
  viewer.Camera().Load(TRUE);
  glMatrixMode(GL_MODELVIEW);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  for (int i = 0; i < resident_nodes.NNodes(); i++) {
    R3SurfelNode *node = resident_nodes.Node(i);
    glLoadName(i + 1);
    for (int j = 0; j < node->NBlocks(); j++) {
      R3SurfelBlock *block = node->Block(j);
      glPushMatrix();
      const R3Point& origin = block->Origin();
      glTranslated(origin[0], origin[1], origin[2]);
      glBegin(GL_POINTS);
      for (int k = 0; k < block->NSurfels(); k++) {
        const R3Surfel *surfel = block->Surfel(k);
        if (exclude_aerial && surfel->IsAerial()) continue;
        glVertex3fv(surfel->Coords());
      }
      glEnd();
      glPopMatrix();
    }
  }
  glFlush();
  select_buffer_hits = glRenderMode(GL_RENDER);

  // Process select buffer to find front-most hit
  GLuint hit = 0;
  GLuint hit_z = 0xFFFFFFFF;
  GLuint *bufp = select_buffer;
  GLuint numnames, z1, z2;
  for (int i = 0; i < select_buffer_hits; i++) {
    numnames = *bufp++;
    z1 = *bufp++;
    z2 = *bufp++;
    while (numnames--) {
      if (z1 < hit_z) {
        hit = *bufp;
        hit_z = z1/2 + z2/2;
      }
      bufp++;
    }
  }

  // Check if hit anything
  if (hit <= 0) return NULL;

  // Find hit node
  hit--; // subtract the one added to avoid zero
  if (hit < 0) return NULL;
  if (hit >= (GLuint) resident_nodes.NNodes()) return NULL;
  R3SurfelNode *hit_node = resident_nodes.Node(hit);

  // Find node part of an object
  R3SurfelNode *picked_node = hit_node;
  if (exclude_nonobjects) {
    // Find node associated with object
    picked_node = NULL;

    // Check if hit node is part of an object
    if (hit_node->Object()) {
      picked_node = hit_node;
    }

    // Check if hit node has ancestor that is part of an object
    if (picked_node == NULL) {
      R3SurfelNode *ancestor = hit_node->Parent();
      while (ancestor) {
        if (ancestor->Object()) { picked_node = ancestor; break; }
        ancestor = ancestor->Parent();
      }
    }
    
    // Check if hit node has descendent that is part of an object
    if (picked_node == NULL) {
      R3Ray ray = viewer.WorldRay(x, y);
      RNScalar t, picked_t = FLT_MAX;
      RNArray<R3SurfelNode *> stack;
      stack.Insert(hit_node);
      while (!stack.IsEmpty()) {
        R3SurfelNode *node = stack.Tail();
        stack.RemoveTail();
        for (int i = 0; i < node->NParts(); i++) {
          stack.Insert(node->Part(i));
        }
        if (node->Object()) {
          if (R3Intersects(ray, node->BBox(), NULL, NULL, &t)) {
            if (t < picked_t) {
              picked_node = node;
              picked_t = t;
            }
          }
        }
      }
    }
  }
    
  // Find hit position
  GLdouble p[3];
  GLdouble modelview_matrix[16];
  GLdouble projection_matrix[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix);
  glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix);
  GLdouble z = (GLdouble) hit_z / (GLdouble) 0xFFFFFFFF;
  gluUnProject(x, y, z, modelview_matrix, projection_matrix, viewport, &(p[0]), &(p[1]), &(p[2]));
  if (picked_position) picked_position->Reset(p[0], p[1], p[2]);

  // Find hit surfel
  if (picked_block || picked_surfel) {
    // Create pointset in vicinity of picked position
    R3Point position(p[0], p[1], p[2]);
    R3SurfelSphereConstraint sphere_constraint(R3Sphere(position, 0.1));
    R3SurfelPointSet *pointset = CreatePointSet(scene, NULL, &sphere_constraint);
    if (pointset) {
      // Find surfel point closest to picked position
      R3SurfelPoint *closest_point = NULL;
      RNLength closest_distance = FLT_MAX;
      for (int i = 0; i < pointset->NPoints(); i++) {
        R3SurfelPoint *point = pointset->Point(i);
        RNLength distance = R3SquaredDistance(point->Position(), position);
        if (distance < closest_distance) {
          closest_distance = distance;
          closest_point = point;
        }
      }

      // Return closest point
      if (closest_point) {
        if (picked_position) *picked_position = closest_point->Position();
        if (picked_block) *picked_block = closest_point->Block();
        if (picked_surfel) *picked_surfel = closest_point->Surfel();
      }

      // Delete point set
      delete pointset;
    }
  }

  // Return picked node
  return picked_node;
}



void R3SurfelViewer::
DrawObject(R3SurfelObject *object, RNFlags draw_flags) const
{
  // Consider every node of object
  if (object->NNodes() > 0) {
    // Draw nodes
    for (int i = 0; i < object->NNodes(); i++) {
      R3SurfelNode *node = object->Node(i);
      if (!node->DrawResidentDescendents(draw_flags)) {
        if (!node->DrawResidentAncestor(draw_flags)) {
          // const char *object_name = (object->Name()) ? object->Name() : "None";
          // fprintf(stderr, "Did not draw object %s\n", object_name);
        }
      }
    }
  }
  else {
    // Draw parts
    for (int i = 0; i < object->NParts(); i++) {
      R3SurfelObject *part = object->Part(i);
      DrawObject(part, draw_flags);
    }
  }
}
  
  
  
////////////////////////////////////////////////////////////////////////
// OBJECT EDITING 
////////////////////////////////////////////////////////////////////////

int R3SurfelViewer::
SplitLeafNodes(R3SurfelNode *source_node, const R3SurfelConstraint& constraint, 
  RNArray<R3SurfelNode *> *nodesA, RNArray<R3SurfelNode *> *nodesB)
{
  // Get convenient variables
  R3SurfelTree *tree = scene->Tree();
  if (!tree) return 0;
  if (!source_node) source_node = tree->RootNode();
  if (!source_node) return 0;

  // Split leaf nodes, WHILE UPDATING NODES IN VIEWER'S RESIDENT SET
  int countA = 0;
  int countB = 0;
  RNArray<R3SurfelNode *> stack;
  stack.Insert(source_node);
  while (!stack.IsEmpty()) {
    R3SurfelNode *node = stack.Tail();
    stack.RemoveTail();
    if (node->NParts() == 0) {
      // Check if node is resident in working set
      int resident_index = resident_nodes.NodeIndex(node);
    
      // Split leaf node
      RNArray<R3SurfelNode *> partsA, partsB;
      if (tree->SplitLeafNodes(node, constraint, &partsA, &partsB)) {
        // Update resident set
        if (resident_index > 0) {
          resident_nodes.RemoveNode(resident_index);
          for (int j = 0; j < partsA.NEntries(); j++) {
            R3SurfelNode *partA = partsA.Kth(j);
            resident_nodes.InsertNode(partA);
          }
          for (int j = 0; j < partsB.NEntries(); j++) {
            R3SurfelNode *partB = partsB.Kth(j);
            resident_nodes.InsertNode(partB);
          }
        }
      }

      // Insert parts into result
      if (nodesA) nodesA->Append(partsA);
      if (nodesB) nodesB->Append(partsB);
      countA += partsA.NEntries();
      countB += partsB.NEntries();
    }
    else {
      for (int i = 0; i < node->NParts(); i++) {
        R3SurfelNode *part = node->Part(i);
        stack.Insert(part);
      }
    } 
  }

  // Return status
  if (countA == 0) return 0;
  if (countB == 0) return 0;
  return 1;
}



int R3SurfelViewer::
SplitObject(R3SurfelObject *object, const R3SurfelConstraint& constraint,
  R3SurfelObject **resultA, R3SurfelObject **resultB)
{
  // Just checking
  assert(object);
  assert(strcmp(object->Name(), "Root"));
  if (object->NNodes() == 0) return 0;

  // Cureate constraint
  R3SurfelMultiConstraint multi_constraint;
  R3SurfelObjectConstraint object_constraint(object);
  multi_constraint.InsertConstraint(&constraint);
  multi_constraint.InsertConstraint(&object_constraint);

  // Create array of nodes
  RNArray<R3SurfelNode *> nodes;
  for (int i = 0; i < object->NNodes(); i++) {
    R3SurfelNode *node = object->Node(i);
    nodes.Insert(node);
  }

  // Split all nodes
  RNArray<R3SurfelNode *> nodesA, nodesB;
  for (int i = 0; i < nodes.NEntries(); i++) {
    R3SurfelNode *node = nodes.Kth(i);
    SplitLeafNodes(node, constraint, &nodesA, &nodesB);
  }

  // Check nodesA
  if (nodesA.NEntries() == 0) {
    if (resultA) *resultA = NULL;
    if (resultB) *resultB = object; 
    return 0;
  }
    
  // Check nodesA
  if (nodesB.NEntries() == 0) {
    if (resultA) *resultA = object;
    if (resultB) *resultB = NULL; 
    return 0;
  }

  // Create new objects
  R3SurfelObject *objectA = new R3SurfelObject();
  R3SurfelObject *objectB = new R3SurfelObject();
  if (!objectA || !objectB) return 0;
      
  // Remove nodes from object
  while (object->NNodes() > 0) {
    R3SurfelNode *node = object->Node(0);
    object->RemoveNode(node);
  }
  
  // Insert nodes into objectA
  for (int j = 0; j < nodesA.NEntries(); j++) {
    R3SurfelNode *nodeA = nodesA.Kth(j);
    objectA->InsertNode(nodeA);
  }
  
  // Insert nodes into objectB
  for (int j = 0; j < nodesB.NEntries(); j++) {
    R3SurfelNode *nodeB = nodesB.Kth(j);
    objectB->InsertNode(nodeB);
  }
      
  // Copy result
  if (resultA) *resultA = objectA;
  if (resultB) *resultB = objectB;
  
  // Return success
  return 1;
}



