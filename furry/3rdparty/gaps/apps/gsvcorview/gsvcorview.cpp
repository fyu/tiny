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



// Data variables

static GSVScene *scene = NULL;
static GSVImage *images[2] = { NULL, NULL };
static R2Point features[2] = { R2zero_point, R2zero_point };



// Display variables

static RNBoolean show_distortions[2] = { FALSE, FALSE };
static RNBoolean show_images[2] = { TRUE, TRUE };
static RNBoolean show_features[2] = { TRUE, TRUE };
static RNBoolean show_meshes[2] = { FALSE, FALSE };
static RNBoolean show_points[2] = { FALSE, FALSE };



// GLUT variables 

static int GLUTwindow = 0;
static int GLUTwindow_height = 768;
static int GLUTwindow_width = 1148;
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
// Utility Functions
////////////////////////////////////////////////////////////////////////

static R2Point
ReprojectedPosition(GSVScene *scene, 
  GSVImage *from_image, const R2Point& from_position, RNBoolean from_is_distorted, 
  GSVImage *to_image, RNBoolean to_is_distorted)
{
  // Initialize result
  R2Point to_position(0,0);

  // Get ray through from_position
  R3Ray from_ray;
  if (from_is_distorted) from_ray = from_image->RayThroughDistortedPosition(from_position);
  else from_ray = from_image->RayThroughUndistortedPosition(from_position);

  // This is temporary
  GSVPanorama *from_panorama = from_image->Panorama();
  if (!from_panorama) return to_position;
  GSVSegment *from_segment = from_panorama->Segment();
  if (!from_segment) return to_position;
  R3MeshIntersection closest; closest.t = FLT_MAX;
  for (int ia = 0; ia < from_segment->NScans(); ia++) {
    GSVScan *scan = from_segment->Scan(ia);
    GSVMesh *mesh = scan->Mesh();
    if (!mesh) continue;
    R3MeshIntersection intersection;
    if (mesh->Intersection(from_ray, &intersection)) {
      if (intersection.t < closest.t) {
        closest = intersection;
      }
    }
  }

  // Project to image
  if (to_is_distorted) to_position = to_image->DistortedPosition(closest.point);
  else to_position = to_image->UndistortedPosition(closest.point);

  // Return mapped position in to_image
  return to_position;
}



////////////////////////////////////////////////////////////////////////
// Drawing Functions
////////////////////////////////////////////////////////////////////////

static void 
Draw3D(int k)
{
  // Get convenient variables
  GSVImage *image = images[k];
  if (!image) return;
  GSVCamera *camera = image->Camera();
  if (!camera) return;
  int width = image->Width();
  int height = image->Height();
  if (width * height == 0) return;
  double scale = (double) GLUTwindow_height / (double) height;
  if (scale == 0) return;
  GSVPose pose = image->Pose();
  const R3Point& viewpoint = pose.Viewpoint();
  const R3Vector& towards = pose.Towards();
  const R3Vector& up = pose.Up();
  R3Triad triad(towards, up);
  R3CoordSystem cs = R3CoordSystem(viewpoint, triad);
  RNAngle xfov = 0.5 * camera->XFov();
  RNAngle yfov = 0.5 * camera->YFov();

  // Set viewport
  int viewport_width = (int) (scale * width + 0.5);
  glViewport(k*viewport_width, 0, viewport_width, GLUTwindow_height);

  // Set camera
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluPerspective(2.0 * RN_RAD2DEG(yfov), tan(xfov) / tan(yfov), 0.001, 1000000);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  cs.InverseMatrix().Load();

  // Draw points
  if (show_points[k]) {
    glDisable(GL_LIGHTING);
    glColor3d(0.8, 0.8, 0.8);
    scene->Draw(GSV_DRAW_POINTS_WITH_VIEWPOINT_DISTANCE_COLOR);
  }

  // Draw mesh
  if (show_meshes[k]) {
    glEnable(GL_LIGHTING);
    static GLfloat default_material[] = { 0.8, 0.8, 0.8, 1.0 };
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, default_material); 
    scene->Draw(GSV_DRAW_MESHES_WITHOUT_COLOR);
    glDisable(GL_LIGHTING);
  }

  // Reset window
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  // Reset viewport
  glViewport(0, 0, GLUTwindow_width, GLUTwindow_height);
}



static void 
Draw2D(int k)
{
  // Loaded texture information
  static GSVImage *texture_images[2] = { NULL, NULL };
  static RNBoolean texture_distortions[2] = { FALSE, FALSE };
  static GLuint texture_ids[2] = { 0, 0 };

  // Get convenient variables
  GSVImage *image = images[k];
  if (!image) return;
  int width = image->Width();
  int height = image->Height();
  if (width * height == 0) return;
  double scale = (double) GLUTwindow_height / (double) height;
  if (scale == 0) return;

  // Set viewport
  int viewport_width = (int) (scale * width + 0.5);
  glViewport(k*viewport_width, 0, viewport_width, GLUTwindow_height);

  // Set window
  glMatrixMode(GL_PROJECTION);  
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(0, width-1, 0, height-1);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  
  // Disable depth buffer 
  glDisable(GL_DEPTH_TEST);
  glDepthMask(FALSE);

  // Draw image
  if (show_images[k] && images[k]) {
    // Load texture (if not same as last time)
    if ((texture_ids[k] == 0) || (texture_images[k] != images[k]) || (texture_distortions[k] != show_distortions[k])) {
      if (texture_ids[k] > 0) glDeleteTextures(1, &texture_ids[k]);
      texture_images[k] = NULL;
      texture_ids[k] = 0;
      R2Image *texture;
      if (show_distortions[k]) texture = images[k]->DistortedImage();
      else texture = images[k]->UndistortedImage();
      if (!texture) return;
      texture_images[k] = images[k];
      texture_distortions[k] = show_distortions[k];
      glGenTextures(1, &texture_ids[k]);
      glBindTexture(GL_TEXTURE_2D, texture_ids[k]);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
      gluBuild2DMipmaps(GL_TEXTURE_2D, 3, texture->Width(), texture->Height(), GL_RGB, GL_UNSIGNED_BYTE, texture->Pixels() );
      delete texture;
    }

    if (texture_ids[k]) {
      // Disable lighting
      glDisable(GL_LIGHTING);
      glColor3d(1, 1, 1);
      
      // Enable texture
      glBindTexture(GL_TEXTURE_2D, texture_ids[k]);
      glEnable(GL_TEXTURE_2D);
      
      // Draw image (as textured polygon)
      glBegin(GL_POLYGON);
      glTexCoord2d(0,0);
      glVertex2i(0,0);
      glTexCoord2d(1,0);
      glVertex2i(width-1, 0);
      glTexCoord2d(1, 1);
      glVertex2i(width-1, height-1);
      glTexCoord2d(0,1);
      glVertex2i(0, height-1);
      glEnd();

      // Disable texture
      glDisable(GL_TEXTURE_2D);
    }
  }

  // Draw feature
  if (show_features[k]) {
    glDisable(GL_LIGHTING);
    glColor3d(1, 1, 0);
    glPointSize(10);
    glBegin(GL_POINTS);
    glVertex2d(features[k].X(), features[k].Y());
    glEnd();
    glPointSize(1);
  }

  // Enable depth buffer
  glEnable(GL_DEPTH_TEST);
  glDepthMask(TRUE);

  // Reset window
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  // Reset viewport
  glViewport(0, 0, GLUTwindow_width, GLUTwindow_height);
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
  // Clear window 
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Draw 2D stuff
  Draw2D(0);
  Draw2D(1);

  // Draw 3D stuff
  Draw3D(0);
  Draw3D(1);

  // Swap buffers 
  glutSwapBuffers();
}    



void GLUTResize(int w, int h)
{
  // Resize window
  glViewport(0, 0, w, h);

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
  // int dx = x - GLUTmouse[0];
  // int dy = y - GLUTmouse[1];
  
  // Process event
  // Nothing for now

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;
}



void GLUTMouse(int button, int state, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;
  
  // Process mouse button event
  if (state == 0) {
    if (button == 0) {
      if (images[0] && images[1]) {
        int k = (x > GLUTwindow_width/2) ? 1 : 0; 
        double scale = (double) images[0]->Height() / (double) GLUTwindow_height;
        double px = scale * x;
        double py = scale * y;
        if (k == 1) px -= 0.5 * scale * GLUTwindow_width;
        features[k].Reset(px, py);
        features[1-k] = ReprojectedPosition(scene, 
          images[k], features[k], show_distortions[k], 
          images[1-k], show_distortions[1-k]);
        glutPostRedisplay();
      }
    }
  }

  // Remember button state 
  int b = (button == GLUT_LEFT_BUTTON) ? 0 : ((button == GLUT_MIDDLE_BUTTON) ? 1 : 2);
  GLUTbutton[b] = (state == GLUT_DOWN) ? 1 : 0;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;
}



void GLUTSpecial(int key, int x, int y)
{
  // Determine image index
  int k = (x < GLUTwindow_width/2) ? 0 : 1; 

  // Invert y coordinate
  y = GLUTwindow_height - y;

  // Process keyboard button event 
  switch (key) {
  case GLUT_KEY_UP:
  case GLUT_KEY_DOWN:
    if (images[k]) {
      GSVTapestry *tapestry = images[k]->Tapestry();
      if (tapestry) {
        int index = images[k]->TapestryIndex();
        if ((key == GLUT_KEY_UP) && (index < tapestry->NImages()-1)) {
          images[k] = tapestry->Image(index+1);
          features[k] = ReprojectedPosition(scene, 
            images[1-k], features[1-k], show_distortions[1-k], 
            images[k], show_distortions[k]);
        }
        else if ((key == GLUT_KEY_DOWN) && (index > 0)) {
          images[k] = tapestry->Image(index-1);
          features[k] = ReprojectedPosition(scene, 
            images[1-k], features[1-k], show_distortions[1-k], 
            images[k], show_distortions[k]);
        }
      }
      break;
    }

  case GLUT_KEY_LEFT:
  case GLUT_KEY_RIGHT:
    if (images[k]) {
      GSVPanorama *panorama = images[k]->Panorama();
      if (panorama) {
        int index = images[k]->PanoramaIndex();
        if ((key == GLUT_KEY_RIGHT) && (index < panorama->NImages()-1)) {
          images[k] = panorama->Image(index+1);
          features[k] = ReprojectedPosition(scene, 
            images[1-k], features[1-k], show_distortions[1-k], 
            images[k], show_distortions[k]);
        }
        else if ((key == GLUT_KEY_LEFT) && (index > 0)) {
          images[k] = panorama->Image(index-1);
          features[k] = ReprojectedPosition(scene, 
            images[1-k], features[1-k], show_distortions[1-k], 
            images[k], show_distortions[k]);
        }
      }
    }
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



void GLUTKeyboard(unsigned char key, int x, int y)
{
  // Determine image index
  int k = (x < GLUTwindow_width/2) ? 0 : 1; 

  // Invert y coordinate
  y = GLUTwindow_height - y;

  // Process keyboard button event 
  switch (key) {
  case 'D':
  case 'd':
    show_distortions[k] = !show_distortions[k];
    features[k] = ReprojectedPosition(scene, 
      images[1-k], features[1-k], show_distortions[1-k], 
      images[k], show_distortions[k]);
    break;

  case 'F':
  case 'f':
    show_features[k] = !show_features[k];
    break;

  case 'I':
  case 'i':
    show_images[k] = !show_images[k];
    break;

  case 'M':
  case 'm':
    show_meshes[k] = !show_meshes[k];
    break;

  case 'P':
  case 'p':
    show_points[k] = !show_points[k];
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
  // Initialize images 
  if (!scene) return;
  if (scene->NRuns() == 0) return;
  GSVRun *run = scene->Run(0);
  if (run->NSegments() == 0) return;
  GSVSegment *segment = run->Segment(0);
  if (segment->NPanoramas() < 2) return;
  GSVPanorama *panorama0 = segment->Panorama(0);
  GSVPanorama *panorama1 = segment->Panorama(1);
  if (panorama0->NImages() == 0) return;
  if (panorama1->NImages() == 0) return;
  images[0] = panorama0->Image(0);
  images[1] = panorama1->Image(0);

  // Initialize features
  features[0].Reset(images[0]->Width()/2, images[0]->Height()/2);
  features[1].Reset(images[1]->Width()/2, images[1]->Height()/2);

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

















