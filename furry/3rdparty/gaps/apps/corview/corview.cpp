// Source file for the mesh viewer program



// Include files 

#include "R3Graphics/R3Graphics.h"
#include "fglut/fglut.h"



// Program variables

static char *mesh1_name = NULL;
static char *mesh2_name = NULL;
static char *points1_name = NULL;
static char *points2_name = NULL;
static char *correspondence_name = NULL;
static char *map_name = NULL;
static char *image_name = NULL;
static R3Vector initial_camera_towards(-0.57735, -0.57735, -0.57735);
static R3Vector initial_camera_up(-0.57735, 0.57735, 0.5773);
static R3Point initial_camera_origin(0,0,0);
static RNBoolean initial_camera = FALSE;
static int print_verbose = 0;



// Type definitions

struct SparseVertexCorrespondence {
  // vertices[0] and vertices[1] have same number of entries
  // For every i, vertices[0][i] and vertices[1][i] are corresponding
  // vertices in mesh[0] and mesh[1], respectively
  R3Mesh *mesh[2];
  RNArray<R3MeshVertex *> vertices[2];
  SparseVertexCorrespondence(R3Mesh *mesh0, R3Mesh *mesh1) { 
    mesh[0] = mesh0; 
    mesh[1] = mesh1; 
  };
};

struct DenseVertexCorrespondence {
  // vertices has mesh[0]->NVertices() entries
  // each vertices[i] is pointer to vertex in mesh[1]
  // corresponding to mesh[0]->Vertex(i) in mesh[0]
  // if vertices[i] is NULL, then mesh[0]->Vertex(i) is outlier
  R3Mesh *mesh[2];
  R3MeshVertex **vertices; 
  DenseVertexCorrespondence(R3Mesh *mesh0, R3Mesh *mesh1) { 
    mesh[0] = mesh0; 
    mesh[1] = mesh1; 
    vertices = new R3MeshVertex * [ mesh0->NVertices() ];
    for (int i = 0; i < mesh0->NVertices(); i++) vertices[i] = NULL;
  };
};



// Application variables

static R3Mesh *meshes[2] = { NULL, NULL };
static RNArray<R3MeshVertex *> *points[2] = { NULL, NULL };
static SparseVertexCorrespondence *correspondence = NULL;
static DenseVertexCorrespondence *map = NULL;
static R3Viewer *viewers[2] = { NULL, NULL };



// GLUT variables 

static int GLUTwindow = 0;
static int GLUTwindow_height = 500;
static int GLUTwindow_width = 1000;
static int GLUTmouse[2] = { 0, 0 };
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTmodifiers = 0;
static int GLUTinteraction_window = 0;



// Display variables

static int show_faces = 1;
static int show_edges = 0;
static int show_vertices = 0;
static int show_vertex_names = 0;
static int show_correspondence = 1;
static int show_map = 1;
static int show_points = 0;



// Colors

static int max_colors = 24;
static GLfloat colors[24][4] = {
  {1,0,0,1}, {0,1,0,1}, {0,0,1,1}, {1,0,1,1}, {0,1,1,1}, {1,1,0,1}, 
  {1,.3,.7,1}, {1,.7,.3,1}, {.7,1,.3}, {.3,1,.7,1}, {.7,.3,1,1}, {.3,.7,1,1}, 
  {1,.5,.5,1}, {.5,1,.5,1}, {.5,.5,1,1}, {1,.5,1,1}, {.5,1,1,1}, {1,1,.5,1}, 
  {.5,0,0,1}, {0,.5,0,1}, {0,0,.5,1}, {.5,0,.5,1}, {0,.5,.5,1}, {.5,.5,0,1} 
};



////////////////////////////////////////////////////////////////////////



void GLUTDrawText(const R3Point& p, const char *s)
{
  // Draw text string s and position p
  glRasterPos3d(p[0], p[1], p[2]);
  while (*s) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *(s++));
}
  


void GLUTStop(void)
{
  // Destroy window 
  glutDestroyWindow(GLUTwindow);

  // Exit
  exit(0);
}



void DrawMesh(int m)
{
  // Get mesh
  R3Mesh *mesh = meshes[m];

  // Set viewing transformation
  viewers[m]->Camera().Load();

  // Draw faces
  if (show_faces) {
    glEnable(GL_LIGHTING);
    glPolygonOffset(1, 1);
    glEnable(GL_POLYGON_OFFSET_FILL);
    static GLfloat default_material[] = { 0.8, 0.8, 0.8, 1.0 };
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, default_material); 
    mesh->DrawFaces();
    glDisable(GL_POLYGON_OFFSET_FILL);
  }

  // Draw edges
  if (show_edges) {
    glDisable(GL_LIGHTING);
    glColor3f(1.0, 0.0, 0.0);
    mesh->DrawEdges();
  }

  // Draw vertex names
  if (show_vertex_names) {
    glDisable(GL_LIGHTING);
    glColor3f(0.5, 0.3, 0.1);
    for (int i = 0; i < mesh->NVertices(); i++) {
      R3MeshVertex *vertex = mesh->Vertex(i);
      char buffer[256];
      sprintf(buffer, "%d", mesh->VertexID(vertex));
      GLUTDrawText(mesh->VertexPosition(vertex), buffer);
    }
  }

  // Draw points
  if (show_points) {
    // Draw points
    glEnable(GL_LIGHTING);
    RNLength radius = 0.01 * mesh->BBox().LongestAxisLength();
    for (int i = 0; i < points[m]->NEntries(); i++) {
      R3MeshVertex *vertex = points[m]->Kth(i);
      glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, colors[i%max_colors]);
      R3Sphere(mesh->VertexPosition(vertex), radius).Draw();
    }
  }

  // Draw correspondence
  if (correspondence && show_correspondence) {
    glEnable(GL_LIGHTING);
    RNLength radius = 0.02 * mesh->BBox().LongestAxisLength();
    for (int i = 0; i < correspondence->vertices[m].NEntries(); i++) {
      glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, colors[i%max_colors]);
      R3MeshVertex *vertex = correspondence->vertices[m].Kth(i);
      R3Sphere(mesh->VertexPosition(vertex), radius).Draw();
    }
    glEnd();
  }

  // Draw map
  if (map && show_map) {
    glDisable(GL_LIGHTING);
    glPointSize(5);
    glBegin(GL_POINTS);
    for (int i = 0; i < meshes[0]->NVertices(); i++) {
      int index = (m == 0) ? i : meshes[1]->VertexID(map->vertices[i]);
      R3MeshVertex *vertex = mesh->Vertex(index);
      R3Point position = mesh->VertexPosition(vertex);
      double value = (double) i / (double) meshes[0]->NVertices();
      RNRgb c(0, 0, 0);
      if (value < 0.5) {
        c[0] = 1 - 2 * value;
        c[1] = 2 * value;
      }
      else {
        c[1] = 1 - 2 * (value - 0.5);
        c[2] = 2 * (value - 0.5);
      }
      RNLoadRgb(c);
      R3LoadPoint(position);
    }
    glEnd();
  }
}    



void GLUTRedraw()
{
  // Clear window 
  glClearColor(0.8, 0.8, 1.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Set lights
  static GLfloat light0_position[] = { 3.0, 4.0, 5.0, 0.0 };
  glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
  static GLfloat light1_position[] = { -3.0, -2.0, -3.0, 0.0 };
  glLightfv(GL_LIGHT1, GL_POSITION, light1_position);

  // Draw first mesh
  glViewport(0, 0, GLUTwindow_width/2, GLUTwindow_height);  
  DrawMesh(0);

  // Draw second mesh
  glViewport(GLUTwindow_width/2, 0, GLUTwindow_width/2, GLUTwindow_height);  
  DrawMesh(1);

  // Capture image and exit
  if (image_name) {
    R2Image image(GLUTwindow_width, GLUTwindow_height, 3);
    image.Capture();
    image.Write(image_name);
    GLUTStop();
  }

  // Swap buffers 
  glutSwapBuffers();
}    



void GLUTResize(int w, int h)
{
  // Resize window
  glViewport(0, 0, w, h);

  // Resize viewer viewport
  viewers[0]->ResizeViewport(0, 0, w/2, h);
  viewers[1]->ResizeViewport(0, 0, w/2, h);

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

  // Determine windows
  int w[2] = { 1, 1 };  
  if (GLUTmodifiers & GLUT_ACTIVE_CTRL) {
    w[1 - GLUTinteraction_window] = 0;
  }

  // World in hand navigation 
  R3Point origin = meshes[0]->BBox().Centroid();
  if (GLUTbutton[0]) {
    if (w[0]) viewers[0]->RotateWorld(1.0, origin, x, y, dx, dy);
    if (w[1]) viewers[1]->RotateWorld(1.0, origin, x, y, dx, dy);
  }
  else if (GLUTbutton[1]) {
    if (w[0]) viewers[0]->ScaleWorld(1.0, origin, x, y, dx, dy);
    if (w[1]) viewers[1]->ScaleWorld(1.0, origin, x, y, dx, dy);
  }
  else if (GLUTbutton[2]) {
    if (w[0]) viewers[0]->TranslateWorld(1.0, origin, x, y, dx, dy);
    if (w[1]) viewers[1]->TranslateWorld(1.0, origin, x, y, dx, dy);
  }

  // Check if need redisplay
  if (GLUTbutton[0] || GLUTbutton[1] || GLUTbutton[2]) 
    glutPostRedisplay();

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

  // Remember interaction window (window that mouse was clicked down)
  GLUTinteraction_window = (x < GLUTwindow_width/2) ? 0 : 1;

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
  case 'C':
  case 'c':
    show_correspondence = !show_correspondence;
    break;

  case 'E':
  case 'e':
    show_edges = !show_edges;
    break;

  case 'F':
  case 'f':
    show_faces = !show_faces;
    break;

  case 'M':
  case 'm':
    show_map = !show_map;
    break;

  case 'P':
  case 'p':
    show_points = !show_points;
    break;

  case 'V':
    show_vertex_names = !show_vertex_names;
    break;

  case 'v':
    show_vertices = !show_vertices;
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




#if 0

void GLUTIdle(void)
{
  // Set current window
  if ( glutGetWindow() != GLUTwindow ) 
    glutSetWindow(GLUTwindow);  

  // Redraw
  glutPostRedisplay();
}

#endif



void GLUTInit(int *argc, char **argv)
{
  // Open window 
  glutInit(argc, argv);
  glutInitWindowPosition(100, 100);
  glutInitWindowSize(GLUTwindow_width, GLUTwindow_height);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); // | GLUT_STENCIL
  GLUTwindow = glutCreateWindow("OpenGL Viewer");

  // Initialize background color 
  glClearColor(1.0, 1.0, 1.0, 1.0);

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
  glEnable(GL_LIGHTING); 

  // Initialize graphics modes  
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);

  // Initialize GLUT callback functions 
  glutDisplayFunc(GLUTRedraw);
  glutReshapeFunc(GLUTResize);
  glutKeyboardFunc(GLUTKeyboard);
  glutSpecialFunc(GLUTSpecial);
  glutMouseFunc(GLUTMouse);
  glutMotionFunc(GLUTMotion);
}



void GLUTMainLoop(void)
{
  // Run main loop -- never returns 
  glutMainLoop();
}


 
////////////////////////////////////////////////////////////////////////



static R3Mesh *
ReadMesh(char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate mesh
  R3Mesh *mesh = new R3Mesh();
  assert(mesh);

  // Read mesh from file
  if (!mesh->ReadFile(filename)) {
    fprintf(stderr, "Unable to read mesh %s\n", filename);
    return NULL;
  }

  // Print statistics
  if (print_verbose) {
    printf("Read mesh from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Faces = %d\n", mesh->NFaces());
    printf("  # Edges = %d\n", mesh->NEdges());
    printf("  # Vertices = %d\n", mesh->NVertices());
    fflush(stdout);
  }

  // Return success
  return mesh;
}



static RNArray<R3MeshVertex *> *
ReadPoints(R3Mesh *mesh, char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Parse filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .ply)\n", filename);
    return 0;
  }

  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open vertex file %s\n", filename);
    return NULL;
  }

  // Allocate array
  RNArray<R3MeshVertex *> *vertices = new RNArray<R3MeshVertex *>();
  if (!vertices) {
    fprintf(stderr, "Unable to allocate vertices for %s\n", filename);
    return NULL;
  }

  // Check filename extension
  if (!strcmp(extension, ".vts")) {
    // Read vertex IDs from file
    char buffer[2048];
    while (fgets(buffer, 2048, fp)) {
      int vertex_index;
      if (sscanf(buffer, "%d", &vertex_index) == 1) {
        // Check vertex
        if ((vertex_index < 0) || (vertex_index >= mesh->NVertices())) {
          fprintf(stderr, "Invalid vertex index %d in %s\n", vertex_index, filename);
          return NULL;
        }

        // Find vertex in mesh
        R3MeshVertex *vertex = mesh->Vertex(vertex_index);
        if (!vertex) {
          fprintf(stderr, "Unable to find vertex %d in %s\n", vertex_index, filename);
          return NULL;
        }

        // Insert vertex
        vertices->Insert(vertex);
      }
    }
  }
  else if (!strcmp(extension, ".pid")) {
    // Open points file
    FILE *fp = fopen(filename, "r");
    if (!fp) {
      fprintf(stderr, "Unable to open points file: %s\n", filename);
      return NULL;
    }

    // Read vertex indicies
    int vertex_index;
    while (fscanf(fp, "%d", &vertex_index) == (unsigned int) 1) {
      if ((vertex_index >= 0) && (vertex_index < mesh->NVertices())) {
        R3MeshVertex *vertex = mesh->Vertex(vertex_index);
        vertices->Insert(vertex);
      }
    }

    // Close points file
    fclose(fp);
  }
  else if (!strcmp(extension, ".xyz")) {
    // Open points file
    FILE *fp = fopen(filename, "r");
    if (!fp) {
      fprintf(stderr, "Unable to open points file: %s\n", filename);
      return NULL;
    }

    // Read points
    double x, y, z;
    while (fscanf(fp, "%lf%lf%lf", &x, &y, &z) == (unsigned int) 3) {
      // Find closest vertex
      R3mesh_mark++;
      R3Point input_position(x, y, z);
      R3MeshVertex *closest_vertex = NULL;
      RNScalar closest_squared_distance = FLT_MAX;
      for (int i = 0; i < mesh->NVertices(); i++) {
        R3MeshVertex *vertex = mesh->Vertex(i);
        const R3Point& vertex_position = mesh->VertexPosition(vertex);
        RNScalar squared_distance = R3SquaredDistance(vertex_position, input_position);
        if (squared_distance < closest_squared_distance) {
          closest_squared_distance = squared_distance;
          closest_vertex = vertex;
        }
      }

      // Insert vertex
      if (closest_vertex && (mesh->VertexMark(closest_vertex) != R3mesh_mark)) {
        mesh->SetVertexMark(closest_vertex, R3mesh_mark);
        vertices->Insert(closest_vertex);
      }
    }
  }
  else {
    fprintf(stderr, "Unrecognized point file extension: %s\n", extension);
    return NULL;
  }

  // Close file
  fclose(fp);

  // Print statistics
  if (print_verbose) {
    printf("Read vertices from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Vertices = %d\n", vertices->NEntries());
    fflush(stdout);
  }

  // Return array of vertices
  return vertices;
}



static SparseVertexCorrespondence *
ReadSparseVertexCorrespondence(R3Mesh *mesh0, R3Mesh *mesh1, char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Parse filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .cor)\n", filename);
    return 0;
  }

  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open correspondence file %s\n", filename);
    return NULL;
  }

  // Create correspondence
  SparseVertexCorrespondence *correspondence = new SparseVertexCorrespondence(mesh0, mesh1);
  if (!correspondence) {
    fprintf(stderr, "Unable to allocate correspondence for %s\n", filename); 
    return NULL; 
  }

  // Check filename extension
  if (!strcmp(extension, ".cor")) {
    // Read correspondences
    int id0, id1;
    while (fscanf(fp, "%d%d", &id0, &id1) == (unsigned int) 2) {
      R3MeshVertex *vertex0 = mesh0->Vertex(id0);
      R3MeshVertex *vertex1 = mesh1->Vertex(id1);
      correspondence->vertices[0].Insert(vertex0);
      correspondence->vertices[1].Insert(vertex1);
    }
  }
  else if (!strcmp(extension, ".vk")) {
    // Read header
    double fdummy;
    char format[256], dummy[256];
    int nmeshes, ncorrespondences, idummy;
    if (fscanf(fp, "%s%s%s%s%s%d%s%s%s%d%s%d%d%d%d%d", dummy, format,  dummy, dummy,  dummy, &nmeshes,  
               dummy, dummy,  dummy, &idummy,   dummy, &idummy, &idummy, &idummy, &idummy, &ncorrespondences) != (unsigned int) 16) {
      fprintf(stderr, "Unable to read %s\n", filename);
      return NULL;
    }

    // Read correspondences    
    int id0, id1;
    int count = 0;
    while (fscanf(fp, "%d%lg%d%lg", &id0, &fdummy, &id1, &fdummy) == (unsigned int) 4) {
      R3MeshVertex *vertex0 = mesh0->Vertex(id0);
      R3MeshVertex *vertex1 = mesh1->Vertex(id1);
      correspondence->vertices[0].Insert(vertex0);
      correspondence->vertices[1].Insert(vertex1);
      count++;
    }

    // Check number of correspondences
    if (count != ncorrespondences) {
      fprintf(stderr, "Mismatching number of correspondences in %s\n", filename);
      return NULL;
    }
  }
  else {
    fprintf(stderr, "Unrecognized correspondence file extension: %s\n", extension);
    return NULL;
  }

  // Close file
  fclose(fp);

  // Print statistics
  if (print_verbose) {
    printf("Read sparse correspondences from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Correspondences = %d\n", correspondence->vertices[0].NEntries());
    fflush(stdout);
  }

  // Return correspondence
  return correspondence;
}



static DenseVertexCorrespondence *
ReadDenseVertexCorrespondence(R3Mesh *mesh0, R3Mesh *mesh1, char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Parse filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .cor)\n", filename);
    return 0;
  }

  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open vertex file %s\n", filename);
    return NULL;
  }

  // Create correspondence
  DenseVertexCorrespondence *correspondence = new DenseVertexCorrespondence(mesh0, mesh1);
  if (!correspondence) {
    fprintf(stderr, "Unable to allocate correspondence for %s\n", filename); 
    return NULL; 
  }

  // Check filename extension
  if (!strcmp(extension, ".map")) {
    // Read correspondences
    for (int i = 0; i < mesh0->NVertices(); i++) {
      int id1;
      if (fscanf(fp, "%d", &id1) != (unsigned int) 1) {
        fprintf(stderr, "Unable to read correspondence %d from %s\n", i, filename);
        return NULL;
      }
      R3MeshVertex *vertex1 = mesh1->Vertex(id1);
      correspondence->vertices[i] = vertex1;
    }
  }
  else {
    fprintf(stderr, "Unrecognized correspondence file extension: %s\n", extension);
    return NULL;
  }

  // Close file
  fclose(fp);

  // Print statistics
  if (print_verbose) {
    printf("Read dense correspondences from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Correspondences = %d\n", mesh0->NVertices());
    fflush(stdout);
  }

  // Return correspondence
  return correspondence;
}



static R3Viewer *
CreateBirdsEyeViewer(R3Mesh *mesh)
{
  // Setup camera view looking down the Z axis
  R3Box bbox = mesh->BBox();
  assert(!bbox.IsEmpty());
  RNLength r = bbox.DiagonalRadius();
  assert((r > 0.0) && RNIsFinite(r));
  if (!initial_camera) initial_camera_origin = bbox.Centroid() - initial_camera_towards * (2.5 * r);
  R3Camera camera(initial_camera_origin, initial_camera_towards, initial_camera_up, 0.4, 0.4, 0.01 * r, 100.0 * r);
  R2Viewport viewport(0, 0, GLUTwindow_width, GLUTwindow_height);
  return new R3Viewer(camera, viewport);
}



////////////////////////////////////////////////////////////////////////



static int 
ParseArgs(int argc, char **argv)
{
  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) { 
        print_verbose = 1; 
      }
      else if (!strcmp(*argv, "-image")) { 
        argc--; argv++; image_name = *argv; 
      }
      else if (!strcmp(*argv, "-correspondence")) { 
        argc--; argv++; correspondence_name = *argv; 
      }
      else if (!strcmp(*argv, "-map")) { 
        argc--; argv++; map_name = *argv; 
      }
      else if (!strcmp(*argv, "-points")) { 
        argc--; argv++; points1_name = *argv; 
        argc--; argv++; points2_name = *argv; 
      }
      else { 
        fprintf(stderr, "Invalid program argument: %s", *argv); 
        exit(1); 
      }
      argv++; argc--;
    }
    else {
      if (!mesh1_name) mesh1_name = *argv;
      else if (!mesh2_name) mesh2_name = *argv;
      else if (!correspondence_name) correspondence_name = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check mesh filename
  if (!mesh1_name || !mesh2_name || !correspondence_name) {
    fprintf(stderr, "Usage: corview mesh1 mesh2 correspondence\n");
    return 0;
  }

  // Resolve map vs. correspondence
  if (strstr(correspondence_name, ".map") && !map_name) {
    map_name = correspondence_name;
    correspondence_name = NULL;
  }

  // Return OK status 
  return 1;
}



int main(int argc, char **argv)
{
  // Parse program arguments
  if (!ParseArgs(argc, argv)) exit(-1);

  // Initialize GLUT
  GLUTInit(&argc, argv);

  // Read first mesh
  meshes[0] = ReadMesh(mesh1_name);
  if (!meshes[0]) exit(-1);

  // Read second mesh
  meshes[1] = ReadMesh(mesh2_name);
  if (!meshes[1]) exit(-1);

  // Read first set of points
  if (points1_name) {
    points[0] = ReadPoints(meshes[0], points1_name);
    if (!points[0]) exit(-1);
  }

  // Read second set of points
  if (points2_name) {
    points[1] = ReadPoints(meshes[1], points2_name);
    if (!points[1]) exit(-1);
  }

  // Read sparse correspondence
  if (correspondence_name) {
    correspondence = ReadSparseVertexCorrespondence(meshes[0], meshes[1], correspondence_name);
    if (!correspondence) exit(-1);
  }

  // Read dense correspondence (map)
  if (map_name) {
    map = ReadDenseVertexCorrespondence(meshes[0], meshes[1], map_name);
    if (!map) exit(-1);
  }

  // Create first viewer
  viewers[0] = CreateBirdsEyeViewer(meshes[0]);
  if (!viewers[0]) exit(-1);

  // Create second viewer
  viewers[1] = CreateBirdsEyeViewer(meshes[1]);
  if (!viewers[1]) exit(-1);

  // Run GLUT interface
  GLUTMainLoop();

  // Return success 
  return 0;
}

















