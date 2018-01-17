// Source file for the google viewer program



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV/GSV.h"
#include "R3Graphics/R3Graphics.h"
#include "fglut/fglut.h"



////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////

// Program arguments

static char *scene_name = NULL;
static int print_verbose = 0;
static int print_debug = 0;



// Data and viewing variables

GSVScene *scene = NULL;
R3Viewer *viewer = NULL;



// Selection variables

static GSVRun *selected_run = NULL;
static GSVSegment *selected_segment = NULL;
static GSVPanorama *selected_panorama = NULL;
static GSVImage *selected_image = NULL;
// static GSVColumn *selected_column = NULL;
static int selected_run_index = -1;
static int selected_segment_index = -1;
static int selected_panorama_index = -1;
static int selected_image_index = -1;
static int selected_column_index = -1;



// Display variables

static int show_points = 1;
static int show_scanlines = 0;
static int show_meshes = 0;
static int show_projected_points = 0;
static int show_image = 0;
static int show_image_poses = 0;
static int show_column_poses = 0;
static int show_lasers = 0;
static int show_bbox = 0;
static int show_vdtm = 0;
static int show_debug = 0;

static RNLength point_size = 1;
static RNLength image_plane_distance = 20;
static int snap_viewpoint = 1;



// GLUT variables 

static int GLUTwindow = 0;
static int GLUTwindow_height = 1024;
static int GLUTwindow_width = 765;
static int GLUTmouse[2] = { 0, 0 };
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTmodifiers = 0;



////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

static GSVScene *
ReadScene(const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate google scene
  GSVScene *scene = new GSVScene();
  if (!scene) {
    fprintf(stderr, "Unable to allocate scene\n");
    return NULL;
  }

  // Open google scene files
  if (!scene->ReadFile(filename)) {
    delete scene;
    return NULL;
  }

  // Print statistics
  if (print_verbose) {
    printf("Read scene from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Runs = %d\n", scene->NRuns());
    printf("  # Cameras = %d\n", scene->NCameras());
    printf("  # Lasers = %d\n", scene->NLasers());
    printf("  # Segments = %d\n", scene->NSegments());
    printf("  # Scans = %d\n", scene->NScans());
    printf("  # Tapestries = %d\n", scene->NTapestries());
    printf("  # Panoramas = %d\n", scene->NPanoramas());
    printf("  # Images = %d\n", scene->NImages());
    printf("  # Scanlines = %d\n", scene->NScanlines());
    fflush(stdout);
  }

  // Return scene
  return scene;
}



////////////////////////////////////////////////////////////////////////
// Drawing Functions
////////////////////////////////////////////////////////////////////////

static void 
DrawProjectedPoints(void)
{
  // Get selected image
  if (!selected_image) return;

  // Compute scale factor
  RNScalar scale = (double) GLUTwindow_width / (double) selected_image->Width();

  // Set projection matrix
  glMatrixMode(GL_PROJECTION);  
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(0, GLUTwindow_width, 0, GLUTwindow_height);

  // Set model view matrix
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glScaled(scale, scale, scale);

  // Disable depth buffer 
  glDisable(GL_DEPTH_TEST);
  glDepthMask(FALSE);

  // Draw projected points
  glBegin(GL_POINTS);
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        for (int il = 0; il < scan->NScanlines(); il++) {
          GSVScanline *scanline = scan->Scanline(il);
          for (int ik = 0; ik < scanline->NPoints(); ik++) {
            const R3Point& world_position = scanline->PointPosition(ik);
            R2Point distorted_position = selected_image->DistortedPosition(world_position);
            R2LoadPoint(distorted_position);
          }
        }
      }
    }
  }
  glEnd();

  // Enable depth buffer
  glEnable(GL_DEPTH_TEST);
  glDepthMask(TRUE);

  // Reset projection matrix
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  // Reset model view matrix
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}



static void 
DrawLasers(void)
{
  // Colors
  GLfloat camera_color[3][3] = { 
    { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } 
  };

  // Draw camera pose for every image
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        for (int il = 0; il < scan->NScanlines(); il++) {
          GSVScanline *scanline = scan->Scanline(il);
          GSVPose pose = scanline->Pose();
          glColor3fv(camera_color[ia % 3]);
          R3Span(pose.Viewpoint(), pose.Viewpoint() + pose.Towards()).Draw();
          R3Span(pose.Viewpoint(), pose.Viewpoint() + 0.5 * pose.Up()).Draw();
          R3Span(pose.Viewpoint(), pose.Viewpoint() + 0.25 * pose.Right()).Draw();
        }
      }
    }
  }
}



static void 
DrawImagePoses(void)
{
  // Colors
  GLfloat camera_color[9][3] = { 
    { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 },
    { 0, 1, 1 }, { 1, 0, 1 }, { 1, 1, 0 },
    { 1, 0.5, 0 }, { 0, 1, 0.5 }, { 0.5, 0, 1 } 
  };

  // Draw camera pose for every image
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ip = 0; ip < segment->NPanoramas(); ip++) {
        GSVPanorama *panorama = segment->Panorama(ip);
        for (int ii = 0; ii < panorama->NImages(); ii++) {
          GSVImage *image = panorama->Image(ii);
          GSVPose pose = image->Pose();
          glColor3fv(camera_color[ii % 9]);
          R3Sphere(pose.Viewpoint(), 0.1).Draw();
          R3Span(pose.Viewpoint(), pose.Viewpoint() + pose.Towards()).Draw();
          R3Span(pose.Viewpoint(), pose.Viewpoint() + 0.5 * pose.Up()).Draw();
          R3Span(pose.Viewpoint(), pose.Viewpoint() + 0.25 * pose.Right()).Draw();
        }
      }
    }
  }
}



static void 
DrawTexturedImage(void)
{
  // Persistent variables
  static GLuint texture_id = 0;
  static GSVImage *last_selected_image = NULL;

  // Just checking
  if (!selected_run) return;
  if (!selected_segment) return;
  if (!selected_panorama) return;
  if (!selected_image) return;

  //////////////////////////////////////////

  // Load texture
  if (selected_image != last_selected_image) {
    last_selected_image = selected_image;
    if (texture_id > 0) glDeleteTextures(1, &texture_id);
    texture_id = 0;

    // Create texture
    R2Image *texture = selected_image->UndistortedImage();
    if (texture) {
      // Create texture
      glGenTextures(1, &texture_id);
      glBindTexture(GL_TEXTURE_2D, texture_id);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
      gluBuild2DMipmaps(GL_TEXTURE_2D, 3, texture->Width(), texture->Height(), GL_RGB, GL_UNSIGNED_BYTE, texture->Pixels() );
      assert(texture_id != 0);
      delete texture;
    }
  }

  // Draw textured polygon
  if (texture_id > 0) {
    // Disable depth buffer 
    glDisable(GL_DEPTH_TEST);
    glDepthMask(FALSE);

    // Enable texture
    glBindTexture(GL_TEXTURE_2D, texture_id);
    glEnable(GL_TEXTURE_2D);

    // Draw textured polygon
    GSVCamera *camera = selected_image->Camera();
    RNAngle xfov = 0.5 * camera->XFov();
    RNAngle yfov = 0.5 * camera->YFov();
    GSVPose pose = selected_image->Pose(selected_column_index);
    R3Point viewpoint = pose.Viewpoint();
    R3Vector towards = pose.Towards();
    R3Vector up = pose.Up();
    R3Vector right = pose.Right();
    R3Point origin = viewpoint + towards * image_plane_distance;
    R3Vector dx = right * image_plane_distance * tan(xfov);
    R3Vector dy = up * image_plane_distance * tan(yfov);
    R3Point ur = origin + dx + dy;
    R3Point lr = origin + dx - dy;
    R3Point ul = origin - dx + dy;
    R3Point ll = origin - dx - dy;
    glColor3d(1,1,1);
    glBegin(GL_POLYGON);
    glTexCoord2d(0,0);
    R3LoadPoint(ll);
    glTexCoord2d(1,0);
    R3LoadPoint(lr);
    glTexCoord2d(1, 1);
    R3LoadPoint(ur);
    glTexCoord2d(0,1);
    R3LoadPoint(ul);
    glEnd();

    // Disable texture
    glDisable(GL_TEXTURE_2D);

    // Enable depth buffer
    glEnable(GL_DEPTH_TEST);
    glDepthMask(TRUE);
  }
}



static void 
DrawTexturedMesh(void)
{
  // Persistent variables
  static GSVImage *last_selected_image = NULL;
  static GSVScan *last_selected_scan = NULL;
  static R3Mesh *mesh = NULL;
  static GLuint texture_id = 0;
  static R2Point *texture_coords = NULL;

  // Just checking
  if (!selected_run) return;
  if (!selected_segment) return;
  if (!selected_panorama) return;
  if (!selected_image) return;
  if (selected_segment->NScans() == 0) return;
  GSVScan *selected_scan = selected_segment->Scan(0);
  if (selected_scan->NScanlines() == 0) return;

  //////////////////////////////////////////

  // Load mesh
  if (selected_scan != last_selected_scan) {
    last_selected_scan = selected_scan;
    if (mesh)  delete mesh;
    mesh = selected_scan->Mesh();
    if (!mesh) {
      fprintf(stderr, "Unable to read mesh\n");
      mesh = NULL;
      return;
    }
  }

  //////////////////////////////////////////

  // Load texture
  if (selected_image != last_selected_image) {
    last_selected_image = selected_image;
    if (texture_coords) delete [] texture_coords;
    if (texture_id > 0) glDeleteTextures(1, &texture_id);
    texture_coords = NULL;
    texture_id = 0;

    // Create texture
    R2Image *texture = selected_image->DistortedImage();
    if (texture) {
      // Create texture
      glGenTextures(1, &texture_id);
      glBindTexture(GL_TEXTURE_2D, texture_id);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
      gluBuild2DMipmaps(GL_TEXTURE_2D, 3, texture->Width(), texture->Height(), GL_RGB, GL_UNSIGNED_BYTE, texture->Pixels() );
      assert(texture_id != 0);

      // Create texture coordinates
      texture_coords = new R2Point [ mesh->NVertices() ];
      for (int i = 0; i < mesh->NVertices(); i++) {
        R3MeshVertex *vertex = mesh->Vertex(i);
        const R3Point& position = mesh->VertexPosition(vertex);
        R2Point distorted_position = selected_image->DistortedPosition(position);
        texture_coords[i][0] = distorted_position[0] / texture->Width();
        texture_coords[i][1] = distorted_position[1] / texture->Height();
      }

      // Delete texture
      delete texture;
    }
  }

  //////////////////////////////////////////

  // Draw textured mesh
  if (mesh && texture_coords && texture_id) {
    // Setup OpenGL modes
    glDisable(GL_LIGHTING);
    glColor3d(1,1,1);
    glBindTexture(GL_TEXTURE_2D, texture_id);
    glEnable(GL_TEXTURE_2D);
    glEnable(GL_CULL_FACE);
    
    // Draw mesh
    for (int i = 0; i < mesh->NFaces(); i++) {
      R3MeshFace *face = mesh->Face(i);
      R3MeshVertex *v0 = mesh->VertexOnFace(face, 0);
      R3MeshVertex *v1 = mesh->VertexOnFace(face, 1);
      R3MeshVertex *v2 = mesh->VertexOnFace(face, 2);
      R3Point p0 = mesh->VertexPosition(v0);
      R3Point p1 = mesh->VertexPosition(v1);
      R3Point p2 = mesh->VertexPosition(v2);
      R2Point t0 = texture_coords[mesh->VertexID(v0)];
      R2Point t1 = texture_coords[mesh->VertexID(v1)];
      R2Point t2 = texture_coords[mesh->VertexID(v2)];
      glBegin(GL_POLYGON);
      R3LoadTextureCoords(t0);
      R3LoadPoint(p0);
      R3LoadTextureCoords(t1);
      R3LoadPoint(p1);
      R3LoadTextureCoords(t2);
      R3LoadPoint(p2);
      glEnd();
    }
  }

  // Reset after drawing textured mesh
  glDisable(GL_CULL_FACE);
  glDisable(GL_TEXTURE_2D);
}



static void 
DrawDebug(void)
{
  // Ray trace selected image
  if (!selected_image) return;

  static R3MeshSearchTree *trees[3] = { NULL };
  
  int step = 100;
  glDisable(GL_LIGHTING);
  glColor3d(0, 1, 0);
  glPointSize(5);
  glBegin(GL_POINTS);
  for (int i = 0; i < selected_image->Width(); i += step) {
    for (int j = 0; j < selected_image->Height(); j += step) {
      R2Point undistorted_position(i+0.5*step, j+0.5*step);
      R3Ray world_ray = selected_image->RayThroughUndistortedPosition(undistorted_position);
      R3MeshIntersection closest;
      closest.t = FLT_MAX;
      for (int k = 0; k < selected_segment->NScans(); k++) {
        GSVScan *scan = selected_segment->Scan(k);
        R3Mesh *mesh = scan->Mesh();
        if (!mesh) continue;
        if (!trees[k]) trees[k] = new R3MeshSearchTree(mesh);
        R3MeshSearchTree *tree = trees[k];
        R3MeshIntersection intersection;
        tree->FindIntersection(world_ray, intersection);
        if (intersection.type != R3_MESH_NULL_TYPE) {
          if (intersection.t < closest.t) {
            closest = intersection;
          }
        }
      }
      R3LoadPoint(closest.point);
    }
  }
  glEnd();
  glPointSize(1);
}



////////////////////////////////////////////////////////////////////////
// GLUT callback functions
////////////////////////////////////////////////////////////////////////

void GLUTStop(void)
{
  // Destroy window 
  glutDestroyWindow(GLUTwindow);

  // Exit
  exit(0);
}



void GLUTRedraw(void)
{
  // Snap viewpoint
  if (snap_viewpoint) {
    if (selected_image) { 
      GSVCamera *camera = selected_image->Camera();
      if (camera) {
        GSVPose pose = selected_image->Pose(selected_column_index);
        const R3Point& viewpoint = pose.Viewpoint();
        const R3Vector& towards = pose.Towards();
        const R3Vector& up = pose.Up();
        RNAngle xfov = 0.5 * camera->XFov();
        RNAngle yfov = 0.5 * camera->YFov();
        R3Camera c(viewpoint, towards, up, xfov, yfov, 0.001, 1000000);
        viewer->SetCamera(c);
      }
    }
  }
 
  // Set viewing transformation
  viewer->Camera().Load();

  // Clear window 
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Set point size
  glPointSize(point_size);

  // Draw textured image
  if (show_image) {
    DrawTexturedImage();
  }

  // Draw textured mesh
  if (show_vdtm) {
    DrawTexturedMesh();
  }

  // Draw points
  if (show_points) {
    scene->Draw(GSV_DRAW_POINTS_WITH_VIEWPOINT_DISTANCE_COLOR);
    // scene->Draw(GSV_DRAW_POINTS_WITH_LASER_INDEX_COLOR);
  }


  // Draw scanlines
  if (show_scanlines) {
    scene->Draw(GSV_DRAW_SCANLINES_WITH_POINT_INDEX_COLOR);
  }


  // Draw meshes
  if (show_meshes) {
    glEnable(GL_LIGHTING);
    static GLfloat default_material[] = { 0.8, 0.8, 0.8, 1.0 };
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, default_material); 
    scene->Draw(GSV_DRAW_MESHES_WITHOUT_COLOR);
    glDisable(GL_LIGHTING);
  }

  // Draw projected points
  if (show_projected_points) {
    glColor3d(1, 0, 0);
    DrawProjectedPoints();
  }

  // Draw image poses
  if (show_image_poses) {
    // Draw pose of all images
    DrawImagePoses();

    // Draw pose of selected image
    if (selected_image) {
      glColor3d(1, 1, 1);
      glLineWidth(5);
      selected_image->Pose().Draw();
      glLineWidth(1);
    }
  }

  // Draw column poses
  if (show_column_poses) {
    if (selected_image) {
      glColor3d(1, 1, 1);
      glPointSize(5);
      glBegin(GL_POINTS);
      for (int i = 0; i < selected_image->Width(); i++) {
        GSVPose column_pose = selected_image->Pose(i);
        const R3Point& column_viewpoint = column_pose.Viewpoint();
        R3LoadPoint(column_viewpoint);
      }
      glEnd();
      glPointSize(1);
    }
  }

  // Draw lasers
  if (show_lasers) {
    DrawLasers();
  }

  // Draw bounding box
  if (show_bbox) {
    glColor3d(1, 1, 1);
    scene->BBox().Outline();
  }

  // Draw debug stuff
  if (show_debug) {
    DrawDebug();
  }

  // Swap buffers 
  glutSwapBuffers();
}    



void GLUTResize(int w, int h)
{
  // Resize window
  glViewport(0, 0, w, h);

  // Resize viewer viewport
  viewer->ResizeViewport(0, 0, w, h);

  // Remember window size 
  GLUTwindow_width = w;
  GLUTwindow_height = h;

  // Redraw
  glutPostRedisplay();
}



void GLUTMotion(int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;

  // Compute mouse movement
  int dx = x - GLUTmouse[0];
  int dy = y - GLUTmouse[1];
  
  // World in hand navigation 
  R3Point origin = R3zero_point;
  if (selected_image) origin = selected_image->Pose().Viewpoint();
  else if (selected_panorama) origin = selected_panorama->Viewpoint();
  else if (selected_segment) origin = selected_segment->BBox().Centroid();
  else if (selected_run) origin = selected_run->BBox().Centroid();
  else if (scene) origin = scene->BBox().Centroid();
  origin += 0.1 * R3negz_vector;
  if (GLUTbutton[0]) viewer->RotateWorld(1.0, origin, x, y, dx, dy);
  else if (GLUTbutton[1]) viewer->ScaleWorld(1.0, origin, x, y, dx, dy);
  else if (GLUTbutton[2]) viewer->TranslateWorld(1.0, origin, x, y, dx, dy);
  if (GLUTbutton[0] || GLUTbutton[1] || GLUTbutton[2]) glutPostRedisplay();

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;
}



void GLUTMouse(int button, int state, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;
  
  // Process mouse button event

  // Remember button state 
  int b = (button == GLUT_LEFT_BUTTON) ? 0 : ((button == GLUT_MIDDLE_BUTTON) ? 1 : 2);
  GLUTbutton[b] = (state == GLUT_DOWN) ? 1 : 0;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

   // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;

  // Redraw
  glutPostRedisplay();
}



void GLUTSpecial(int key, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;

  // Process keyboard button event 
  switch (key) {
  case GLUT_KEY_UP:
    if (glutGetModifiers() & GLUT_ACTIVE_SHIFT) {
      if (scene) {
        if (selected_run_index < scene->NRuns() - 1) {
          selected_run = scene->Run(++selected_run_index);
          selected_segment = selected_run->Segment(selected_segment_index);
          selected_panorama = selected_segment->Panorama(selected_panorama_index);
          selected_image = selected_panorama->Image(selected_image_index);
          // selected_column = selected_image->Column(selected_column_index);
        } 
      } 
     }
    else {
      if (selected_segment) {
        if (selected_panorama_index < selected_segment->NPanoramas() - 1) {
          selected_panorama = selected_segment->Panorama(++selected_panorama_index);
          selected_image = selected_panorama->Image(selected_image_index);
          // selected_column = selected_image->Column(selected_column_index);
        }  
      }
    }
    break;

  case GLUT_KEY_DOWN:
    if (glutGetModifiers() & GLUT_ACTIVE_SHIFT) {
      if (scene) {
        if (selected_run_index > 0) {
          selected_run = scene->Run(--selected_run_index);
          selected_segment = selected_run->Segment(selected_segment_index);
          selected_panorama = selected_segment->Panorama(selected_panorama_index);
          selected_image = selected_panorama->Image(selected_image_index);
          // selected_column = selected_image->Column(selected_column_index);
        } 
      } 
    }
    if (selected_segment) {
      if (selected_panorama_index > 0) {
        selected_panorama = selected_segment->Panorama(--selected_panorama_index);
        selected_image = selected_panorama->Image(selected_image_index);
        // selected_column = selected_image->Column(selected_column_index);
      }
    }
    break;

  case GLUT_KEY_RIGHT:
    if (glutGetModifiers() & GLUT_ACTIVE_SHIFT) {
      if (selected_run) {
        if (selected_segment_index < selected_run->NSegments() - 1) {
          selected_segment = selected_run->Segment(++selected_segment_index);
          selected_panorama = selected_segment->Panorama(selected_panorama_index);
          selected_image = selected_panorama->Image(selected_image_index);
          // selected_column = selected_image->Column(selected_column_index);
        }  
      }
    }
    else {
      if (selected_panorama) {
        if (selected_image_index < selected_panorama->NImages() - 1) {
          selected_image = selected_panorama->Image(++selected_image_index);
          // selected_column = selected_image->Column(selected_column_index);
        }
      }
    }
    break;

  case GLUT_KEY_LEFT:
    if (glutGetModifiers() & GLUT_ACTIVE_SHIFT) {
      if (selected_run) {
        if (selected_segment_index > 0) {
          selected_segment = selected_run->Segment(--selected_segment_index);
          selected_panorama = selected_segment->Panorama(selected_panorama_index);
          selected_image = selected_panorama->Image(selected_image_index);
          // selected_column = selected_image->Column(selected_column_index);
        }  
      }
    }
    else {
      if (selected_panorama) {
        if (selected_image_index > 0) {
          selected_image = selected_panorama->Image(--selected_image_index);
          // selected_column = selected_image->Column(selected_column_index);
        }
      }
    }
    break;

  case GLUT_KEY_HOME:
    if (selected_image) {
      selected_column_index = 0;
      // selected_column = selected_image->Column(selected_column_index);
      printf("%d\n", selected_column_index);
    }
    break;

  case GLUT_KEY_END:
    if (selected_image) {
      selected_column_index = selected_image->Width()-1;
      // selected_column = selected_image->Column(selected_column_index);
      printf("%d\n", selected_column_index);
    }
    break;

  case GLUT_KEY_PAGE_UP:
    if (selected_image) {
      selected_column_index += 100;
      if (selected_column_index >= selected_image->Width()) selected_column_index = selected_image->Width()-1;
      // selected_column = selected_image->Column(selected_column_index);
      printf("%d\n", selected_column_index);
    }
    break;

  case GLUT_KEY_PAGE_DOWN:
    if (selected_image) {
      selected_column_index -= 100;
      if (selected_column_index < 0) selected_column_index = 0;
      // selected_column = selected_image->Column(selected_column_index);
      printf("%d\n", selected_column_index);
    }
    break;
  }

  // Check for debug command
  if (glutGetModifiers() & GLUT_ACTIVE_ALT) {
    if (key == GLUT_KEY_F1) show_debug = !show_debug;
  }

  // Print debug information
  if (print_debug) {
    RNScalar timestamp = (selected_image) ? selected_image->Timestamp(selected_column_index) : 0;
    R3Point viewpoint = (selected_image) ? selected_image->Pose(selected_column_index).Viewpoint() : R3zero_point;
    R3Quaternion orientation = (selected_image) ? selected_image->Pose(selected_column_index).Orientation() : R3zero_quaternion;
    printf("r=%-4d s=%-4d p=%-4d i=%-4d c=%-4d tm=%-9.6f vp=%-9.6f %-9.6f %-9.6f or=%-9.6f %-9.6f %-9.6f %-9.6f \n", 
      selected_run_index, selected_segment_index, selected_panorama_index, selected_image_index, selected_column_index,
      timestamp, viewpoint[0], viewpoint[1], viewpoint[2], orientation[0], orientation[1], orientation[2], orientation[3]);
  }

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

  // Redraw
  glutPostRedisplay();
}



void GLUTKeyboard(unsigned char key, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;

  // Process keyboard button event 
  switch (key) {
  case '!':
    show_vdtm = !show_vdtm;
    break;

  case 'B':
  case 'b':
    show_bbox = !show_bbox;
    break;

  case 'C':
    show_column_poses = !show_column_poses;
    break;

  case 'c':
    show_image_poses = !show_image_poses;
    break;

  case 'D':
  case 'd':
    show_projected_points = !show_projected_points;
    break;

  case 'I':
  case 'i':
    show_image = !show_image;
    break;

  case 'L':
  case 'l':
    show_lasers = !show_lasers;
    break;

  case 'M':
  case 'm':
    show_meshes = !show_meshes;
    break;

  case 'P':
  case 'p':
    show_points = !show_points;
    break;

  case 'S':
  case 's':
    show_scanlines = !show_scanlines;
    break;

  case 'V':
  case 'v':
    snap_viewpoint = !snap_viewpoint;
   break;

  case '+':
  case '=':
    point_size *= 1.1;
    break;

  case '_':
  case '-':
    point_size *= 0.9;
    break;

  case 27: // ESCAPE
    GLUTStop();
    break;
  }

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

  // Redraw
  glutPostRedisplay();  
}



////////////////////////////////////////////////////////////////////////
// GLUT callback functions
////////////////////////////////////////////////////////////////////////

static void 
InitInterface(int *argc, char **argv)
{
  // Open window 
  glutInit(argc, argv);
  glutInitWindowPosition(100, 100);
  glutInitWindowSize(GLUTwindow_width, GLUTwindow_height);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); // | GLUT_STENCIL
  GLUTwindow = glutCreateWindow("GSV Street Viewer");

  // Initialize background color
  glClearColor(0, 0, 0, 1);

  // Initialize lights
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
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

  // Initialize GLUT callback functions 
  glutDisplayFunc(GLUTRedraw);
  glutReshapeFunc(GLUTResize);
  glutKeyboardFunc(GLUTKeyboard);
  glutSpecialFunc(GLUTSpecial);
  glutMouseFunc(GLUTMouse);
  glutMotionFunc(GLUTMotion);
}



static void 
RunInterface(void)
{
  // Initialize viewer
  R3Box bbox = scene->BBox();
  assert(!bbox.IsEmpty());
  RNLength r = bbox.DiagonalRadius();
  assert((r > 0.0) && RNIsFinite(r));
  R3Point origin = bbox.Centroid() + R3posz_vector * (2.5 * r);
  RNAngle yfov = RN_PI_OVER_FOUR;
  RNScalar window_aspect = (double) GLUTwindow_height / (double) GLUTwindow_width;
  RNAngle xfov = atan(tan(yfov) / window_aspect);
  R3Camera camera(origin, R3Vector(0, 0, -1), R3Vector(0, 1, 0), xfov, yfov, 0.001 * r, 1000.0 * r);
  R2Viewport viewport(0, 0, GLUTwindow_width, GLUTwindow_height);
  viewer = new R3Viewer(camera, viewport);

  // Initialize selection variables
  if (scene) {
    if (scene->NTapestries() > 0) {
      selected_run_index = 0;
      selected_run = scene->Run(0);
      if (selected_run->NSegments() > 0) {
        selected_segment_index = 0;
        selected_segment = selected_run->Segment(0);
        if (selected_segment->NPanoramas() > 0) {
          selected_panorama_index = 0;
          selected_panorama = selected_segment->Panorama(0);
          if (selected_panorama->NImages() > 0) {
            selected_image_index = 0;
            selected_image = selected_panorama->Image(0);
            if (selected_image->Width() > 0) {
              selected_column_index = selected_image->Width() / 2;
              // selected_column = selected_image->Column(0);
            }
          }
        }
      }
    }
  }

  // Run main loop -- never returns
  glutMainLoop();
}



////////////////////////////////////////////////////////////////////////
// Argument Parsing Functions
////////////////////////////////////////////////////////////////////////

static int 
ParseArgs(int argc, char **argv)
{
  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) { print_verbose = 1; }
      else if (!strcmp(*argv, "-debug")) { print_debug = 1; }
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
    else {
      if (!scene_name) scene_name = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check googles name
  if (!scene_name) {
    fprintf(stderr, "Usage: gsvview scenefile [options]\n");
    return FALSE;
  }

  // Return OK status 
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // Parse program arguments
  if (!ParseArgs(argc, argv)) exit(-1);

  // Read scene
  scene = ReadScene(scene_name);
  if (!scene) exit(-1);

  // Initialize GLUT interface
  InitInterface(&argc, argv);

  // Run GLUT interface
  RunInterface();

  // Return success 
  return 0;
}

















