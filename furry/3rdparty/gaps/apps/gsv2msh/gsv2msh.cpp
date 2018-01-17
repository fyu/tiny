// Source file for the google info program



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV/GSV.h"
#include "process.h"



////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////

// Program arguments

static const char *input_scene_name = NULL;
static const char *output_mesh_directory = NULL;
static RNLength min_viewpoint_movement = 0.05;
static RNLength max_depth_discontinuity = 0.5;
static int fill_holes = 0;
static int split_edges = 0;
static int collapse_edges = 0;
static int smooth_vertices = 0;
static int poisson_reconstruction = 0;
int print_verbose = 0;
int print_debug = 0;



////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

static GSVScene *
ReadScene(const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Check if should read points (if .gsv, will be read as needed)
  RNBoolean read_points = (strstr(filename, ".gsv")) ? FALSE : TRUE;

  // Allocate google scene
  GSVScene *scene = new GSVScene();
  if (!scene) {
    fprintf(stderr, "Unable to allocate scene\n");
    return NULL;
  }

  // Open google scene files
  if (!scene->ReadFile(filename, read_points)) {
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



static int
WriteMesh(GSVScan *scan, GSVMesh *mesh, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Write mesh
  if (!mesh->WritePlyFile(scan, filename, TRUE)) return 0;

  // Print debug statistics
  if (print_debug) {
    printf("  Wrote mesh %s ...\n", filename);
    printf("    Time = %.2f seconds\n", start_time.Elapsed());
    printf("    # Faces = %d\n", mesh->NFaces());
    printf("    # Edges = %d\n", mesh->NEdges());
    printf("    # Vertices = %d\n", mesh->NVertices());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
WriteMeshProperties(GSVScan *scan, GSVMesh *mesh, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Write mesh
  if (!mesh->WriteARFFFile(scan, filename)) return 0;

  // Print debug statistics
  if (print_debug) {
    printf("  Wrote mesh properties %s ...\n", filename);
    printf("    Time = %.2f seconds\n", start_time.Elapsed());
    printf("    # Vertices = %d\n", mesh->NVertices());
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Mesh creation 
////////////////////////////////////////////////////////////////////////

static GSVMesh *
CreateMesh(GSVScan *scan, 
  RNLength min_viewpoint_movement,
  RNLength max_depth_discontinuity)
{
  // NOTES: 
  // Should space samples based on distances between vertices, not car movement
  // Should handle fold-overs due to car moving backwards or rotating inwards XXX***XXX

  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Just checking
  if (!scan) return NULL;
  if (scan->NScanlines() == 0) return NULL;

  // Read scan points
  if (!scan->ReadPoints()) {
    fprintf(stderr, "Unable to read points\n");
    return 0;
  }

  // Create mesh
  GSVMesh *mesh = new GSVMesh();
  if (!mesh) {
    fprintf(stderr, "Unable to allocate mesh\n");
    return NULL;
  }

  // Load scan
  if (!mesh->LoadScan(scan, min_viewpoint_movement, max_depth_discontinuity)) {
    fprintf(stderr, "Unable to load scan into mesh\n");
    return NULL;
  }

  // Release scan points
  if (!scan->ReleasePoints()) {
    fprintf(stderr, "Unable to release points\n");
    return 0;
  }

  // Print debug statistics
  if (print_debug) {
    printf("  Created mesh ...\n");
    printf("    Time = %.2f seconds\n", start_time.Elapsed());
    printf("    # Faces = %d\n", mesh->NFaces());
    printf("    # Edges = %d\n", mesh->NEdges());
    printf("    # Vertices = %d\n", mesh->NVertices());
    fflush(stdout);
  }

  // Return mesh
  return mesh;
}



////////////////////////////////////////////////////////////////////////
// Poisson Surface Reconstruction
////////////////////////////////////////////////////////////////////////

static int
WritePoissonPoints(R3Mesh *mesh, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open point file %s\n", filename);
    return 0;
  }

  // Write points
  float buffer[6];
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(i);
    const R3Point& position = mesh->VertexPosition(vertex);
    const R3Vector& normal = mesh->VertexNormal(vertex);

    // Fill buffer
    buffer[0] = position[0];
    buffer[1] = position[1];
    buffer[2] = position[2];
    buffer[3] = -normal[0];
    buffer[4] = -normal[1];
    buffer[5] = -normal[2];

    // Write buffer
    if (fwrite(buffer, sizeof(float), 6, fp) != (unsigned int) 6) {
      fprintf(stderr, "Unable to write point %d to %s\n", i, filename);
      return 0;
    }
  }

  // Close file
  fclose(fp);

  // Print debug statistics
  if (print_debug) {
    printf("  Wrote points to %s ...\n", filename);
    printf("    Time = %.2f seconds\n", start_time.Elapsed());
    printf("    # Points = %d\n", mesh->NVertices());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
DeleteOutlierVertices(R3Mesh *mesh, RNLength min_edge_length_threshold = 1)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int count = 0;

  // Allocate temporary data
  RNLength *min_edge_lengths = new RNLength [ mesh->NVertices() ];

  // Create array of vertices 
  RNArray<R3MeshVertex *> vertices;
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(i);
    vertices.Insert(vertex);
    min_edge_lengths[i] = FLT_MAX;
    for (int j = 0; j < mesh->VertexValence(vertex); j++) {
      R3MeshEdge *edge = mesh->EdgeOnVertex(vertex, j);
      RNLength edge_length = mesh->EdgeLength(edge);
      if (edge_length < min_edge_lengths[i]) {
        min_edge_lengths[i] = edge_length;
      }
    }
  }

  // Delete vertices with minimum edge lengths over threshold
  for (int i = 0; i < vertices.NEntries(); i++) {
    R3MeshVertex *vertex = vertices.Kth(i);
    if (min_edge_lengths[i] < min_edge_length_threshold) continue;
    mesh->DeleteVertex(vertex);
    count++;
  }

  // Delete temporary data
  delete [] min_edge_lengths;

  // Print debug statistics
  if (print_debug) {
    printf("  Deleted outlier vertices ...\n");
    printf("    Time = %.2f seconds\n", start_time.Elapsed());
    printf("    # Deleted = %d\n", count);
    fflush(stdout);
  }

  // Return success
  return 1;
}



R3Mesh *
CreatePoissonMesh(GSVScan *scan, 
  RNLength min_viewpoint_movement,
  RNLength max_depth_discontinuity)
{
  // Create initial mesh
  R3Mesh *mesh1 = CreateMesh(scan, min_viewpoint_movement, max_depth_discontinuity);
  if (!mesh1) return NULL;
  
  // Write points
  char points_name[4096];
  sprintf(points_name, "foo.bnpts");
  if (!WritePoissonPoints(mesh1, points_name)) { 
    delete mesh1; 
    return NULL; 
  }

  // Delete initial mesh
  delete mesh1;

  // Run Misha's code
  char mesh_name[4096];
  sprintf(mesh_name, "foo.ply");
  char command[4096];
  sprintf(command, "c:/Funk/downloads/PoissonRecon/PoissonRecon.64.exe --in %s --out %s --verbose --depth 8", points_name, mesh_name);
  system(command);

  // Read reconstructed mesh
  R3Mesh *mesh2 = new R3Mesh();
  if (!mesh2->ReadFile(mesh_name)) {
    delete mesh2;
    return NULL;
  }

  // Flip faces
  for (int i = 0; i < mesh2->NFaces(); i++) {
    mesh2->FlipFace(mesh2->Face(i));
  }

  // Delete outlier vertices
  if (!DeleteOutlierVertices(mesh2)) {
    delete mesh2;
    return NULL;
  }

  // Return mesh
  return mesh2;
}



static int
WriteMeshes(GSVScene *scene, 
  const char *output_mesh_directory,
  RNLength min_viewpoint_movement,
  RNLength max_depth_discontinuity)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int count = 0;
  if (print_verbose) {
    printf("Creating meshes ...\n");
    fflush(stdout);
  }

  // Write mesh for each scan
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (scan->NScanlines() < 2) continue;

        // Create mesh
        GSVMesh *mesh = CreateMesh(scan, min_viewpoint_movement, max_depth_discontinuity);
        if (!mesh) return 0;

        // Write mesh
        char mesh_name[4096];
        sprintf(mesh_name, "%s/%s/%02d_%02d_DA_Mesh.ply", output_mesh_directory, run->Name(), is, ia);
        if (!WriteMesh(scan, mesh, mesh_name)) return 0;

        // Write mesh properties
        char properties_name[4096];
        sprintf(properties_name, "%s/%s/%02d_%02d_DA_Mesh.arff", output_mesh_directory, run->Name(), is, ia);
        // if (!WriteMeshProperties(scan, mesh, properties_name)) return 0;

        // Process mesh with basic functions
        if (fill_holes && !FillHoles(mesh)) continue;
        for (RNLength edge_length = 2; edge_length > 0.2; edge_length *= 0.5) {
          if (split_edges && !SplitEdges(mesh, 2*edge_length)) continue;
          if (collapse_edges && !CollapseEdges(mesh, edge_length)) continue;
          if (smooth_vertices && !SmoothVertices(mesh, edge_length)) continue;
        }
 
        // Write processed mesh
        sprintf(mesh_name, "%s/%s/%02d_%02d_DA_ProcessedMesh.ply", output_mesh_directory, run->Name(), is, ia);
        if (!WriteMesh(scan, mesh, mesh_name)) return 0;

        // Delete mesh
        delete mesh;

        // Poisson reconstruction
        if (poisson_reconstruction) {
          // Create Poisson mesh
          R3Mesh *poisson_mesh = CreatePoissonMesh(scan, min_viewpoint_movement, max_depth_discontinuity);
          if (!poisson_mesh) return 0; 
          
          // Write Poisson mesh
          char mesh_name[4096];
          sprintf(mesh_name, "%s/%s/%02d_%02d_DA_PoissonMesh.ply", output_mesh_directory, run->Name(), is, ia);
          if (!poisson_mesh->WriteFile(mesh_name)) return 0;

          // Delete Poisson mesh
          delete poisson_mesh;
        }

        // Update statistics
        count++;
      }
    }
  }

  // Print statistics
  if (print_verbose) {
    printf("  Finished creating meshes\n");
    printf("  Overall Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Meshes = %d\n", count);
    fflush(stdout);
  }

  // Return success
  return 1;
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
      if (!strcmp(*argv, "-v")) print_verbose = 1;
      else if (!strcmp(*argv, "-debug")) print_debug = 1;
      else if (!strcmp(*argv, "-fill_holes")) fill_holes = 1;
      else if (!strcmp(*argv, "-split_edges")) split_edges = 1;
      else if (!strcmp(*argv, "-collapse_edges")) collapse_edges = 1;
      else if (!strcmp(*argv, "-smooth_vertices")) smooth_vertices = 1;
      else if (!strcmp(*argv, "-poisson_reconstruction")) poisson_reconstruction = 1;
      else if (!strcmp(*argv, "-min_viewpoint_movement")) { argc--; argv++; min_viewpoint_movement = atof(*argv); }
      else if (!strcmp(*argv, "-max_depth_discontinuity")) { argc--; argv++; max_depth_discontinuity = atof(*argv); }
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
    else {
      if (!input_scene_name) input_scene_name = *argv;
      else if (!output_mesh_directory) output_mesh_directory = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check scene name
  if (!input_scene_name || !output_mesh_directory) {
    fprintf(stderr, "Usage: gsv2msh input_scene output_directory [options]\n");
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
  GSVScene *scene = ReadScene(input_scene_name);
  if (!scene) exit(-1);

  // Create output directories
  char cmd[1024];
  sprintf(cmd, "mkdir -p %s", output_mesh_directory);  
  system(cmd);
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    sprintf(cmd, "mkdir -p %s/%s", output_mesh_directory, run->Name());
    system(cmd);
  }
    
  // Create mesh
  if (!WriteMeshes(scene, output_mesh_directory, min_viewpoint_movement, max_depth_discontinuity)) exit(-1);

  // Return success 
  return 0;
}

















