// Source file for the google street view pose optimization program 



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV/GSV.h"
#include "R3Graphics/R3Graphics.h"
#include "R3CatmullRomSpline.h"
#include "fglut/fglut.h"
#include "poseoptimization.h"
// #include "DrawText.cpp"



////////////////////////////////////////////////////////////////////////
// Program arguments
////////////////////////////////////////////////////////////////////////

// Input/output

static char *input_scene_name = NULL;
static char *output_scene_name = NULL;
static char *input_correspondences_name = NULL;
static char *output_correspondences_name = NULL;
static char *input_transformations_name = NULL;
static char *output_transformations_name = NULL;
static char *input_inertias_name = NULL;
static char *output_inertias_name = NULL;
static char *input_options_name = NULL;
static char *output_options_name = NULL;
static char *output_world_distance_error_plot_name = NULL;
static char *output_image_distance_error_plot_name = NULL;
static int print_verbose = 0;
static int print_debug = 0;
static int interactive = 0;


// Feature creation

static int create_scan_pole_features = 0;
static int create_scan_curb_features = 0;
static int create_scan_ridge_and_valley_features = 0;
static int create_scan_edge_features = 0;
static int create_scan_plane_features = 0;
static int create_image_corner_features = 0;
static int create_image_sift_features = 0;
static int create_image_line_features = 0;
static int create_features = 0;


// Feature update

static int update_scan_features = 0;
static int update_pixel_features = 0;
static int update_features = 0;


// Pair creation

static int create_pairs = 0;
static int create_all_pairs = 0;
static int create_feasible_pairs = 0;
static int create_close_pairs = 0;
static int create_coplanar_pairs = 0;
static int create_pixel_pixel_plane_pairs = 0;
static int remove_pairs = 0;
static int remove_all_pairs = 0;
static int remove_infeasible_pairs = 0;



// Correspondence creation

static int create_correspondences = 0;
static int create_all_correspondences = 0;
static int create_feasible_correspondences = 0;
static int create_close_correspondences = 0;
static int create_coplanar_correspondences = 0;
static int create_transitive_closure_correspondences = 0;
static int create_pixel_pixel_plane_correspondences = 0;
static int create_icp_correspondences = 0;
static int remove_correspondences = 0;
static int remove_all_correspondences = 0;
static int remove_infeasible_correspondences = 0;



// Cluster creation

static int create_clusters = 0;
static int create_all_clusters = 0;
static int create_close_clusters = 0;
static int create_icp_clusters = 0;
static int create_coplanar_clusters = 0;
static int create_merged_clusters = 0;
static int create_hierarchical_clusters = 0;
static int remove_clusters = 0;
static int remove_all_clusters = 0;



// Transformation creation

static int update_transformations = 0;
static int update_all_transformations = 0;
static int update_laser_translations = 0;
static int update_laser_rotations = 0;
static int update_camera_translations = 0;
static int update_camera_rotations = 0;
static int update_path_translations = 0;
static int update_path_rotations = 0;
static int update_cluster_positions = 0;
static int update_cluster_rotations = 0;
static int update_pixel_depths = 0;



// Optimization parameters

static int include_correspondence_equations = -1;
static int include_cluster_equations = -1;
static int solver = RN_CERES_SOLVER;
static RNLength path_vertex_spacing = 2;
static RNLength optimization_radius = RN_UNKNOWN;
static R2Point optimization_center(RN_UNKNOWN, RN_UNKNOWN);
static RNScalar rigidity = 0.5;



// Data management
static int load_meshes = 0;
static int load_segments = 0;
static int load_planar_grids = 0;
static int load_all_points = 0;
static const char *segment_grid_name = "LabelMajorStructure";



////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////

// Data variables

static GSVScene *scene = NULL;
static GSVPoseOptimization *optimization = NULL;
static GSVFeatureCorrespondence *last_correspondence = NULL;
static GSVFeature *last_feature = NULL;
static int num_initial_pairs = 0;
static int num_initial_correspondences = 0;



// Display variables

static RNRgb background_color(0, 0, 0);
static int show_scans[2] = { 1, 1 };
static int show_images[2] = { 0, 0 };
static int show_meshes[2] = { 0, 0 };
static int show_cameras[2] = { 0, 0 };
static int show_lasers[2] = { 1, 1 };
static int show_paths[2] = { 0, 0 };
static int show_bboxes[2] = { 0, 0 };
static int show_segments[2] = { 0, 0 };
static int show_planar_grids[2] = { 0, 0 };
static int show_features[2] = { 0, 0 };
static int show_pairs[2] = { 0, 0 };
static int show_correspondences[2] = { 1, 1 };
static int show_clusters[2] = { 0, 0 };
static int show_optimization[2] = { 1, 1 };
static R3Viewer *viewer[2] = { NULL, NULL };
static GSVImage *selected_image[2] = { NULL, NULL };
static GSVScanline *selected_scanline[2] = { NULL, NULL };
static GSVFeature *selected_feature[2] = { NULL, NULL };
static R3Point center_point[2] = { R3Point(0,0,0), R3Point(0,0,0) };
static int viewpoint_scheme[2] = { 0, 0 };
static int color_scheme[2] = { 1, 1 };
static float point_size[2] = { 1, 1 };
static int timestamp_radius = 10;
static double image_plane_distance = 20;
static int split_screen = 0;



// Useful constants

static int GSV_IMAGE_HEIGHT = 2592;
static int GSV_IMAGE_WIDTH = 1936;



// GLUT variables 

static int GLUTwindow = 0;
static int GLUTwindow_height = 1024;
static int GLUTwindow_width = 2 * GLUTwindow_height * GSV_IMAGE_WIDTH / GSV_IMAGE_HEIGHT;
static int GLUTmouse[2] = { 0, 0 };
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTmodifiers = 0;
static int GLUTside = 0;



// Color schemes

enum {
  NO_COLOR_SCHEME,
  SEGMENT_COLOR_SCHEME,
  VIEWPOINT_DEPTH_COLOR_SCHEME,
  HEIGHT_COLOR_SCHEME,
  PICK_COLOR_SCHEME
};



////////////////////////////////////////////////////////////////////////
// I/O functions
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

  // Read scene 
  if (!scene->ReadFile(filename, load_all_points)) {
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
WriteScene(GSVScene *scene, GSVPoseOptimization *optimization, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Read all points (for now)
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        scan->ReadPoints();
      }
    }
  }

  // Apply pose transformations
  if (optimization) {
    if (!optimization->ApplyOptimizedTransformationsToScene()) return 0;
  }

  // Write scene
  if (!scene->WriteFile(filename)) return 0;

  // Release all points (for now)
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        scan->ReleasePoints();
      }
    }
  }

  // Print statistics
  if (print_verbose) {
    printf("Wrote scene to %s ...\n", filename);
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

  // Return success
  return 1;
}



static int
ReadCorrespondences(GSVPoseOptimization *optimization, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Read correspondences
  if (!optimization->ReadFile(filename)) return 0;

  // Remember number of initial pairs and correspondences
  num_initial_pairs = optimization->NPairs();
  num_initial_correspondences = optimization->NCorrespondences();

  // Print statistics
  if (print_verbose) {
    printf("Read %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Features = %d\n", optimization->NFeatures());
    printf("  # Pairs = %d\n", optimization->NPairs());
    printf("  # Correspondences = %d\n", optimization->NCorrespondences());
    printf("  # Clusters = %d\n", optimization->NClusters());
    printf("  World RMSD = %g\n", optimization->WorldRMSD());
    printf("  Image RMSD = %g\n", optimization->ImageRMSD());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
WriteCorrespondences(GSVPoseOptimization *optimization, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Determine if should apply transformations when write
  RNBoolean apply_pose_transformations = (output_scene_name) ? TRUE : FALSE;

  // Write correspondences
  if (!optimization->WriteFile(filename, apply_pose_transformations)) return 0;

  // Print statistics
  if (print_verbose) {
    printf("Wrote %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Features = %d\n", optimization->NFeatures());
    printf("  # Pairs = %d\n", optimization->NPairs());
    printf("  # Correspondences = %d\n", optimization->NCorrespondences());
    printf("  # Clusters = %d\n", optimization->NClusters());
    printf("  World RMSD = %g\n", optimization->WorldRMSD());
    printf("  Image RMSD = %g\n", optimization->ImageRMSD());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
ReadTransformations(GSVPoseOptimization *optimization, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Read transformations
  if (!optimization->ReadTransformationsFile(filename)) return 0;

  // Print statistics
  if (print_verbose) {
    printf("Read transformations from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  World RMSD = %g\n", optimization->WorldRMSD());
    printf("  Image RMSD = %g\n", optimization->ImageRMSD());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
WriteTransformations(GSVPoseOptimization *optimization, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Write transformations
  if (!optimization->WriteTransformationsFile(filename)) return 0;

  // Print statistics
  if (print_verbose) {
    printf("Wrote transformations to %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  World RMSD = %g\n", optimization->WorldRMSD());
    printf("  Image RMSD = %g\n", optimization->ImageRMSD());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
ReadInertias(GSVPoseOptimization *optimization, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Read inertias
  if (!optimization->ReadInertiaFile(filename)) return 0;

  // Print statistics
  if (print_verbose) {
    printf("Read inertias from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
WriteInertias(GSVPoseOptimization *optimization, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Write inertias
  if (!optimization->WriteInertiaFile(filename)) return 0;

  // Print statistics
  if (print_verbose) {
    printf("Wrote inertias to %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
UpdateInertias(GSVPoseOptimization *optimization)
{
  // Set inertias based on optimization center and radius
  if ((optimization_radius > 0) && (optimization_center != R2unknown_point)) {
    // Compute useful variables for computing proximity weights
    RNScalar sigma = optimization_radius / 3.0;
    RNScalar exp_factor = -0.5 / (sigma * sigma);
    RNLength radius_squared  = optimization_radius * optimization_radius;
    
    // Set inertias for path variables
    for (int i = 0; i < optimization->vertices.NEntries(); i++) {
      GSVPathVertex *vertex = optimization->vertices.Kth(i);
      R3Point vertex_position = vertex->pose.Viewpoint();
      RNLength dx = optimization_center.X() - vertex_position.X();
      RNLength dy = optimization_center.Y() - vertex_position.Y();
      RNLength dd = dx*dx + dy*dy;
      RNScalar inertia = (dd < radius_squared) ? 1.0 - exp(dd * exp_factor) : 1.0;
      for (int j = 0; j < 6; j++) {
        vertex->inertia[j] = inertia;
      }
    }
    
    // Set inertias for cluster variables
    for (int i = 0; i < optimization->clusters.NEntries(); i++) {
      GSVFeatureCluster *cluster = optimization->clusters.Kth(i);
      R3Point cluster_position = cluster->scan_position;
      RNLength dx = optimization_center.X() - cluster_position.X();
      RNLength dy = optimization_center.Y() - cluster_position.Y();
      RNLength dd = dx*dx + dy*dy;
      RNScalar inertia = (dd < radius_squared) ? 1.0 - exp(dd * exp_factor) : 1.0;
      for (int j = 0; j < 6; j++) {
        cluster->inertia[j] = inertia;
      }
    }
  }
    
  // Fix inertias of axis aligned vectors
  for (int i = 0; i < optimization->NClusters(); i++) {
    GSVFeatureCluster *cluster = optimization->Cluster(i);
    if (cluster->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) {
      cluster->inertia[optimization->RX] = 1.0;
      cluster->inertia[optimization->RY] = 1.0;
      cluster->inertia[optimization->RZ] = 1.0;
    }
    else if (cluster->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) {
      cluster->inertia[optimization->TZ] = 1.0;
      if (RNIsEqual(fabs(cluster->scan_direction.X()), 1.0) ||
          RNIsEqual(fabs(cluster->scan_direction.Y()), 1.0) || 
          RNIsEqual(fabs(cluster->scan_direction.Z()), 1.0)) {
        cluster->inertia[optimization->RX] = 1.0;
        cluster->inertia[optimization->RY] = 1.0;
        cluster->inertia[optimization->RZ] = 1.0;
      }
    }
    else if (cluster->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) {
      cluster->inertia[optimization->TY] = 1.0;
      cluster->inertia[optimization->TZ] = 1.0;
      if (RNIsEqual(fabs(cluster->scan_normal.X()), 1.0) ||
          RNIsEqual(fabs(cluster->scan_normal.Y()), 1.0) || 
          RNIsEqual(fabs(cluster->scan_normal.Z()), 1.0)) {
        cluster->inertia[optimization->RX] = 1.0;
        cluster->inertia[optimization->RY] = 1.0;
        cluster->inertia[optimization->RZ] = 1.0;
      }
    }
  }

  // Return success
  return 1;
}



static int
ReadOptions(GSVPoseOptimization *optimization, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open options file %s\n", filename);
    return 0;
  }

  // Parse file
  char buffer[4096];
  while (fgets(buffer, 4096, fp)) {
    char *bufferp = buffer;
    if (*bufferp == '#') continue;
    if (*bufferp == '\0') continue;
    char *keyword = strtok(bufferp, " \t\n=");
    if (!keyword) continue;
    char *value = strtok(NULL, " \t\n=");
    if (!value) continue;

    // Check keyword
    if (!strcmp(keyword, "#")) continue;
    
    // Feature creation options
    else if (!strcmp(keyword, "create_scan_pole_features"))
      create_features |= create_scan_pole_features = atoi(value);
    else if (!strcmp(keyword, "create_scan_curb_features"))
      create_features |= create_scan_curb_features = atoi(value);
    else if (!strcmp(keyword, "create_scan_ridge_and_valley_features"))
      create_features |= create_scan_ridge_and_valley_features = atoi(value);
    else if (!strcmp(keyword, "create_scan_edge_features"))
      create_features |= create_scan_edge_features = atoi(value);
    else if (!strcmp(keyword, "create_scan_plane_features"))
      create_features |= create_scan_plane_features = atoi(value);
    else if (!strcmp(keyword, "create_image_corner_features"))
      create_features |= create_image_corner_features = atoi(value);
    else if (!strcmp(keyword, "create_image_sift_features"))
      create_features |= create_image_sift_features = atoi(value);
    else if (!strcmp(keyword, "create_image_line_features"))
      create_features |= create_image_line_features = atoi(value);

    // Pair creation selections
    else if (!strcmp(keyword, "create_all_pairs"))
      create_pairs |= create_all_pairs = atoi(value);
    else if (!strcmp(keyword, "create_feasible_pairs"))
      create_pairs |= create_feasible_pairs = atoi(value);
    else if (!strcmp(keyword, "create_close_pairs"))
      create_pairs |= create_close_pairs = atoi(value);
    else if (!strcmp(keyword, "create_coplanar_pairs"))
      create_pairs |= create_coplanar_pairs = atoi(value);
    else if (!strcmp(keyword, "create_pixel_pixel_plane_pairs"))
      create_pairs |= create_pixel_pixel_plane_pairs = atoi(value);

    // Pair removal selections
    else if (!strcmp(keyword, "remove_all_pairs"))
      remove_pairs |= remove_all_pairs = atoi(value);
    else if (!strcmp(keyword, "remove_infeasible_pairs"))
      remove_pairs |= remove_infeasible_pairs = atoi(value);

    // Correspondence creation selections
    else if (!strcmp(keyword, "create_all_correspondences"))
      create_correspondences |= create_all_correspondences = atoi(value);
    else if (!strcmp(keyword, "create_feasible_correspondences"))
      create_correspondences |= create_feasible_correspondences = atoi(value);
    else if (!strcmp(keyword, "create_close_correspondences"))
      create_correspondences |= create_close_correspondences = atoi(value);
    else if (!strcmp(keyword, "create_coplanar_correspondences"))
      create_correspondences |= create_coplanar_correspondences = atoi(value);
    else if (!strcmp(keyword, "create_transitive_closure_correspondences"))
      create_correspondences |= create_transitive_closure_correspondences = atoi(value);
    else if (!strcmp(keyword, "create_pixel_pixel_plane_correspondences"))
      create_correspondences |= create_pixel_pixel_plane_correspondences = atoi(value);
    else if (!strcmp(keyword, "create_icp_correspondences"))
      create_correspondences |= create_icp_correspondences = atoi(value);

    // Correspondence removal selections
    else if (!strcmp(keyword, "remove_all_correspondences"))
      remove_correspondences |= remove_all_correspondences = atoi(value);
    else if (!strcmp(keyword, "remove_infeasible_correspondences"))
      remove_correspondences |= remove_infeasible_correspondences = atoi(value);

    // Cluster creation selections
    else if (!strcmp(keyword, "create_all_clusters"))
      create_clusters |= create_all_clusters = atoi(value);
    else if (!strcmp(keyword, "create_close_clusters"))
      create_clusters |= create_close_clusters = atoi(value);
    else if (!strcmp(keyword, "create_merged_clusters"))
      create_clusters |= create_merged_clusters = atoi(value);
    else if (!strcmp(keyword, "create_hierarchical_clusters"))
      create_clusters |= create_hierarchical_clusters = atoi(value);
    else if (!strcmp(keyword, "create_coplanar_clusters"))
      create_clusters |= create_coplanar_clusters = atoi(value);
    else if (!strcmp(keyword, "create_icp_clusters"))
      create_clusters |= create_icp_clusters = atoi(value);

    // Cluster removal selections
    else if (!strcmp(keyword, "remove_all_clusters"))
      remove_clusters |= remove_all_clusters = atoi(value);

    // Feature update options 
    else if (!strcmp(keyword, "update_scan_features"))
      update_features |= update_scan_features = atoi(value);
    else if (!strcmp(keyword, "update_pixel_features"))
      update_features |= update_pixel_features = atoi(value);

    // Optimization options
    else if (!strcmp(keyword, "update_laser_translations"))
      update_transformations = update_laser_translations = atoi(value);
    else if (!strcmp(keyword, "update_laser_rotations"))
      update_transformations = update_laser_rotations = atoi(value);
    else if (!strcmp(keyword, "update_camera_translations"))
      update_transformations = update_camera_translations = atoi(value);
    else if (!strcmp(keyword, "update_camera_rotations"))
      update_transformations = update_camera_rotations = atoi(value);
    else if (!strcmp(keyword, "update_path_translations"))
      update_transformations = update_path_translations = atoi(value);
    else if (!strcmp(keyword, "update_path_rotations"))
      update_transformations = update_path_rotations = atoi(value);
    else if (!strcmp(keyword, "update_cluster_positions"))
      update_transformations = update_cluster_positions = atoi(value);
    else if (!strcmp(keyword, "update_cluster_rotations"))
      update_transformations = update_cluster_rotations = atoi(value);
    else if (!strcmp(keyword, "update_pixel_depths"))
      update_transformations = update_pixel_depths = atoi(value);
    else if (!strcmp(keyword, "include_correspondence_equations"))
      include_correspondence_equations = atoi(value);
    else if (!strcmp(keyword, "include_cluster_equations"))
      include_cluster_equations = atoi(value);

    // Optimization parameters
    else if (!strcmp(keyword, "optimization_radius")) 
      optimization_radius = atof(value);
    else if (!strcmp(keyword, "optimization_center_x")) 
      optimization_center[0] = atof(value);
    else if (!strcmp(keyword, "optimization_center_y")) 
      optimization_center[1] = atof(value);
    else if (!strcmp(keyword, "rigidity"))
      rigidity = atof(value);
    else if (!strcmp(keyword, "solver"))
      solver = atoi(value);
    
    // Other options parsed by poseoptimization
    else if (!optimization->ParseOption(keyword, value)) {
      fprintf(stderr, "Unable to parse option %s from file %s\n", keyword, filename);
      return 0;
    }
  }

  // Close file
  fclose(fp);

  // Print statistics
  if (print_verbose) {
    printf("Read options from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
WriteOptions(GSVPoseOptimization *optimization, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open options file %s\n", filename);
    return 0;
  }

  // Feature creation options
  fprintf(fp, "# Feature creation options\n");
  fprintf(fp, "create_scan_pole_features = %d\n", create_scan_pole_features);
  fprintf(fp, "create_scan_curb_features = %d\n", create_scan_curb_features);
  fprintf(fp, "create_scan_ridge_and_valley_features = %d\n", create_scan_ridge_and_valley_features);
  fprintf(fp, "create_scan_edge_features = %d\n", create_scan_edge_features);
  fprintf(fp, "create_scan_plane_features = %d\n", create_scan_plane_features);
  fprintf(fp, "create_image_corner_features = %d\n", create_image_corner_features);
  fprintf(fp, "create_image_sift_features = %d\n", create_image_sift_features);
  fprintf(fp, "create_image_line_features = %d\n", create_image_line_features);
  fprintf(fp, "\n");

  // Pair creation/removal selections
  fprintf(fp, "# Pair creation options\n");
  fprintf(fp, "create_all_pairs = %d\n", create_all_pairs);
  fprintf(fp, "create_feasible_pairs = %d\n", create_feasible_pairs);
  fprintf(fp, "create_close_pairs = %d\n", create_close_pairs);
  fprintf(fp, "create_coplanar_pairs = %d\n", create_coplanar_pairs);
  fprintf(fp, "create_pixel_pixel_plane_pairs = %d\n", create_pixel_pixel_plane_pairs);
  fprintf(fp, "remove_all_pairs = %d\n", remove_all_pairs);
  fprintf(fp, "remove_infeasible_pairs = %d\n", remove_infeasible_pairs);
  fprintf(fp, "\n");

  // Correspondence creation/removal selections
  fprintf(fp, "# Correspondence creation options\n");
  fprintf(fp, "create_all_correspondences = %d\n", create_all_correspondences);
  fprintf(fp, "create_feasible_correspondences = %d\n", create_feasible_correspondences);
  fprintf(fp, "create_close_correspondences = %d\n", create_close_correspondences);
  fprintf(fp, "create_coplanar_correspondences = %d\n", create_coplanar_correspondences);
  fprintf(fp, "create_transitive_closure_correspondences = %d\n", create_transitive_closure_correspondences);
  fprintf(fp, "create_pixel_pixel_plane_correspondences = %d\n", create_pixel_pixel_plane_correspondences);
  fprintf(fp, "create_icp_correspondences = %d\n", create_icp_correspondences);
  fprintf(fp, "remove_all_correspondences = %d\n", remove_all_correspondences);
  fprintf(fp, "remove_infeasible_correspondences = %d\n", remove_infeasible_correspondences);
  fprintf(fp, "\n");

  // Cluster creation/removal selections
  fprintf(fp, "# Cluster creation options\n");
  fprintf(fp, "create_all_clusters = %d\n", create_all_clusters);
  fprintf(fp, "create_close_clusters = %d\n", create_close_clusters);
  fprintf(fp, "create_merged_clusters = %d\n", create_merged_clusters);
  fprintf(fp, "create_hierarchical_clusters = %d\n", create_hierarchical_clusters);
  fprintf(fp, "create_coplanar_clusters = %d\n", create_coplanar_clusters);
  fprintf(fp, "create_icp_clusters = %d\n", create_icp_clusters);
  fprintf(fp, "remove_all_clusters = %d\n", remove_all_clusters);
  fprintf(fp, "\n");

  // Feature update options 
  fprintf(fp, "# Feature update options\n");
  fprintf(fp, "update_scan_features = %d\n", update_scan_features);
  fprintf(fp, "update_pixel_features = %d\n", update_pixel_features);
  fprintf(fp, "\n");

  // Path creation parameters
  fprintf(fp, "# Path creation parameters\n");
  fprintf(fp, "max_path_vertex_spacing = %g\n", optimization->max_path_vertex_spacing);
  fprintf(fp, "\n");

  // Optimization options
  fprintf(fp, "# Optimization parameters\n");
  fprintf(fp, "update_laser_translations = %d\n", update_laser_translations);
  fprintf(fp, "update_laser_rotations = %d\n", update_laser_rotations);
  fprintf(fp, "update_camera_translations = %d\n", update_camera_translations);
  fprintf(fp, "update_camera_rotations = %d\n", update_camera_rotations);
  fprintf(fp, "update_path_translations = %d\n", update_path_translations);
  fprintf(fp, "update_path_rotations = %d\n", update_path_rotations);
  fprintf(fp, "update_cluster_positions = %d\n", update_cluster_positions);
  fprintf(fp, "update_cluster_rotations = %d\n", update_cluster_rotations);
  fprintf(fp, "update_pixel_depths = %d\n", update_pixel_depths);
  fprintf(fp, "include_correspondence_equations = %d\n", include_correspondence_equations);
  fprintf(fp, "include_cluster_equations = %d\n", include_cluster_equations);

  // Optimization parameters
  fprintf(fp, "optimization_radius = %g\n", optimization_radius);
  fprintf(fp, "optimization_center_x = %g\n", optimization_center.X());
  fprintf(fp, "optimization_center_y = %g\n", optimization_center.Y());
  fprintf(fp, "rigidity = %g\n", rigidity);
  fprintf(fp, "solver = %d\n", solver);
  fprintf(fp, "\n");
  
  // Write poseoptimization options
  optimization->WriteOptions(fp);

  // Close file
  fclose(fp);

  // Print statistics
  if (print_verbose) {
    printf("Wrote options to %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
WriteErrorPlot(GSVPoseOptimization *optimization, const char *filename = NULL,  
  RNBoolean include_scan_scan_correspondences = TRUE, 
  RNBoolean include_scan_image_correspondences = TRUE, 
  RNBoolean include_image_image_correspondences = TRUE, 
  RNBoolean world_distance = TRUE)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Set parameters
  double bin_spacing = (world_distance) ? 0.1 : 1.0;
  RNLength max_error = (world_distance) ? 4 : 100;

  // Update clusters and pixel depths
  if (!optimization->Solve(FALSE, FALSE,  FALSE, FALSE,  FALSE, FALSE,  TRUE, TRUE, TRUE)) return 0;

  // Create histogram
  int count = 0;
  int nbins = (int) (max_error / bin_spacing);
  double *bins = new double [nbins];
  for (int i = 0; i < nbins; i++) bins[i] = 0;

  // Add correspondence errors
  for (int i = 0; i < optimization->NCorrespondences(); i++) {
    GSVFeatureCorrespondence *correspondence = optimization->Correspondence(i);
    for (int j = 0; j < 2; j++) {
      GSVFeature *feature0 = correspondence->Feature(j);
      GSVFeature *feature1 = correspondence->Feature(1-j);
      if (feature0->image && feature1->image && !include_image_image_correspondences) continue;
      if (!feature0->image && feature1->image && !include_scan_image_correspondences) continue;
      if (feature0->image && !feature1->image && !include_scan_image_correspondences) continue;
      if (!feature0->image && !feature1->image && !include_scan_scan_correspondences) continue;
      RNLength dd = 0;
      if (world_distance) dd = optimization->SquaredWorldDistance(feature0, feature1);
      else dd = optimization->SquaredImageDistance(feature0, feature1);
      if (dd == RN_INFINITY) continue;
      RNLength error = sqrt(dd);
      double bin = error / bin_spacing;
      int bin0 = (int) bin;
      int bin1 = bin0 + 1;
      if (bin0 >= nbins) bin0 = nbins - 1;
      if (bin1 >= nbins) bin1 = nbins - 1;
      double t = bin - bin0;
      bins[bin0] += 1-t;
      bins[bin1] += t;
      count++;
    }
  }

  // Add cluster errors
  for (int i = 0; i < optimization->NClusters(); i++) {
    GSVFeatureCluster *cluster = optimization->Cluster(i);
    for (int j = 0; j < cluster->NFeatures(); j++) {
      GSVFeature *feature = cluster->Feature(j);
      if (feature->image && !include_image_image_correspondences) continue;
      if (feature->image && !include_scan_image_correspondences) continue;
      RNLength dd = 0;
      if (world_distance) dd = optimization->SquaredWorldDistance(cluster, feature);
      else dd = optimization->SquaredImageDistance(cluster, feature);
      if (dd == RN_INFINITY) continue;
      RNLength error = sqrt(dd);
      double bin = error / bin_spacing;
      int bin0 = (int) bin;
      int bin1 = bin0 + 1;
      if (bin0 >= nbins) bin0 = nbins - 1;
      if (bin1 >= nbins) bin1 = nbins - 1;
      double t = bin - bin0;
      bins[bin0] += 1-t;
      bins[bin1] += t;
      count++;
    }
    if (world_distance) {
      for (int j = 0; j < cluster->NClusters(); j++) {
        GSVFeatureCluster *child = cluster->Cluster(j);
        if (!include_scan_scan_correspondences) continue;
        RNLength dd = 0;
        dd = optimization->SquaredWorldDistance(cluster, child);
        if (dd == RN_INFINITY) continue;
        RNLength error = sqrt(dd);
        double bin = error / bin_spacing;
        int bin0 = (int) bin;
        int bin1 = bin0 + 1;
        if (bin0 >= nbins) bin0 = nbins - 1;
        if (bin1 >= nbins) bin1 = nbins - 1;
        double t = bin - bin0;
        bins[bin0] += 1-t;
        bins[bin1] += t;
        count++;
      }
    }
  }

  // Open file
  FILE *fp = stdout;
  if (filename) fp = fopen(filename, "w");
  if (!fp) { 
    fprintf(stderr, "Unable to open error plot %s\n", filename); 
    return 0; 
  }

  // Write CDF error plot
  double sum = 0;
  if (count == 0) count = 1;
  for (int i = 0; i < nbins; i++) {
    sum += bins[i];
    fprintf(fp, "%g %g\n", (i+0.5) * bin_spacing, sum / count);
  }

  // Close file
  if (filename) fclose(fp);

  // Delete histogram 
  delete [] bins;

  // Print statistics
  if (print_verbose) {
    const char *name = (world_distance) ? "world" : "image";
    printf("Wrote %s distance error plot to %s ...\n", name, filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Features = %d\n", optimization->NFeatures());
    printf("  # Pairs = %d\n", optimization->NPairs());
    printf("  # Correspondences = %d\n", optimization->NCorrespondences());
    printf("  # Clusters = %d\n", optimization->NClusters());
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Optimization functions
////////////////////////////////////////////////////////////////////////

static int
Solve(GSVPoseOptimization *optimization, 
  RNBoolean laser_translations = -1, RNBoolean laser_rotations = -1, 
  RNBoolean camera_translations = -1, RNBoolean camera_rotations = -1, 
  RNBoolean path_translations = -1, RNBoolean path_rotations = -1, 
  RNBoolean cluster_positions = -1, RNBoolean cluster_rotations = -1,
  RNBoolean pixel_depths = -1, 
  RNBoolean correspondence_equations = -1, RNBoolean cluster_equations = -1)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Update stuff by default
  if (pixel_depths < 0) pixel_depths = 1;
  if (cluster_positions < 0) cluster_positions = cluster_equations;
  if (cluster_rotations < 0) cluster_rotations = cluster_equations;
  if (cluster_equations < 0) cluster_equations = cluster_positions;

  // Update optimization
  if (!optimization->Solve(
    (laser_translations >= 0) ? laser_translations : update_laser_translations, 
    (laser_rotations >= 0) ? laser_rotations : update_laser_rotations,
    (camera_translations >= 0) ? camera_translations : update_camera_translations, 
    (camera_rotations >= 0) ? camera_rotations : update_camera_rotations,
    (path_translations >= 0) ? path_translations : update_path_translations, 
    (path_rotations >= 0) ? path_rotations : update_path_rotations,
    (cluster_positions >= 0) ? cluster_positions : update_cluster_positions,
    (cluster_rotations >= 0) ? cluster_rotations : update_cluster_rotations,
    (pixel_depths >= 0) ? pixel_depths : update_pixel_depths,
    (correspondence_equations >= 0) ? correspondence_equations : include_correspondence_equations, 
    (cluster_equations >= 0) ? cluster_equations : include_cluster_equations, 
    rigidity, solver)) {
    fprintf(stderr, "Unable to update optimization\n");
    return 0;
  }

  // Print statistics
  if (print_verbose) {
    printf("Optimized pose transformations ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Features = %d\n", optimization->NFeatures());
    printf("  # Pairs = %d\n", optimization->NPairs());
    printf("  # Correspondences = %d\n", optimization->NCorrespondences());
    printf("  # Clusters = %d\n", optimization->NClusters());
    printf("  World RMSD = %g\n", optimization->WorldRMSD());
    printf("  Image RMSD = %g\n", optimization->ImageRMSD());
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Feature creation functions
////////////////////////////////////////////////////////////////////////

static int
CreateFeatures(GSVPoseOptimization *optimization)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int saved_nfeatures = optimization->NFeatures();
  if (print_verbose) {
    printf("Creating features ...\n");
    fflush(stdout);
  }

  // Create scan pole features
  if (create_scan_pole_features) {
    if (!optimization->CreateScanPoleFeatures()) return 0;
  }
    
  // Create scan curb features
  if (create_scan_curb_features) {
    if (!optimization->CreateScanCurbFeatures()) return 0;
  }
    
  // Create scan edge features
  if (create_scan_ridge_and_valley_features) {
    if (!optimization->CreateScanRidgeAndValleyFeatures()) return 0;
  }
    
  // Create scan edge features
  if (create_scan_edge_features) {
    if (!optimization->CreateScanEdgeFeatures()) return 0;
  }
    
  // Create scan plane features
  if (create_scan_plane_features) {
    if (!optimization->CreateScanPlaneFeatures()) return 0;
  }
    
  // Create image corner features
  if (create_image_corner_features) {
    if (!optimization->CreateImageCornerFeatures()) return 0;
  }

  // Create image sift features
  if (create_image_sift_features) {
    if (!optimization->CreateImageSiftFeatures()) return 0;
  }

  // Create image line features
  if (create_image_line_features) {
    if (!optimization->CreateImageLineFeatures()) return 0;
  }

  // Print statistics
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # New Features = %d\n", optimization->NFeatures() - saved_nfeatures);
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Feature update
////////////////////////////////////////////////////////////////////////

static int 
UpdateScanFeatures(GSVPoseOptimization *optimization)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int count = 0;

  // Snap scan features to scan point
  for (int i = 0; i < optimization->NFeatures(); i++) {
    GSVFeature *feature = optimization->Feature(i);
    if ((feature->feature_type != GSV_SCAN_POINT_FEATURE_TYPE) &&
        (feature->feature_type != GSV_SCAN_PLANE_FEATURE_TYPE) &&
        (feature->feature_type != GSV_SCAN_LINE_FEATURE_TYPE)) continue;
    GSVScanline *scanline = feature->scanline;
    if (!scanline) continue;
    int point_index = feature->scan_point_index;
    if (point_index < 0) continue;
    GSVScan *scan = scanline->Scan();
    scan->ReadPoints();
    feature->scan_position = scanline->PointPosition(point_index);
    scan->ReleasePoints();
    count++;
  }
  
  // Print statistics
  if (print_verbose) {
    printf("Updated scan features ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Scan Features = %d\n", count);
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int 
UpdatePixelFeatures(GSVPoseOptimization *optimization)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_verbose) {
    printf("Updating pixel features ...\n");
    fflush(stdout);
  }

  // Solve for cluster and pixel locations
  int saved_ncorrespondences = optimization->NCorrespondences();
  optimization->CreateAllCorrespondences();
  optimization->Solve(FALSE, FALSE,  FALSE, FALSE,  FALSE, FALSE,  TRUE, TRUE, TRUE);
  optimization->TruncateCorrespondences(saved_ncorrespondences);

#if 0
  // Delete features that are not in bbox
  const R3Box& bbox = R3null_box;
  int nfeatures1 = optimization->NFeatures();
  if (!bbox.IsEmpty()) {
    RNArray<GSVFeature *> features = optimization->features;
    for (int i = 0; i < features.NEntries(); i++) {
      GSVFeature *feature = features.Kth(i);
      if (feature->feature_type != GSV_IMAGE_POINT_FEATURE_TYPE) continue;
      const R3Point& viewpoint = feature->image->Pose().Viewpoint();
      if (!R3Contains(bbox, viewpoint)) {
        optimization->RemoveFeature(feature); 
      }
    }
  }
  
  // Delete features that are not in center of image
  int nfeatures2 = optimization->NFeatures();
  RNScalar center_size = 0.8;
  if (center_size < 1.0) {
    RNScalar lo = 0.5 * (1 - center_size);
    RNScalar hi = 1 - lo;
    RNArray<GSVFeature *> features = optimization->features;
    for (int i = 0; i < features.NEntries(); i++) {
      GSVFeature *feature = features.Kth(i);
      if (feature->feature_type != GSV_IMAGE_POINT_FEATURE_TYPE) continue;
      RNScalar x = feature->image_position.X() / feature->image->Width();
      RNScalar y = feature->image_position.X() / feature->image->Width();
      if ((x < lo) || (x > hi) || (y < lo) || (y > hi)) {
        optimization->RemoveFeature(feature); 
      }
    }
  }
  
  // Update positions of pixel features
  if (!optimization->UpdatePixelPositionsFromMeshIntersections()) exit(-1);

  // Delete features with no mesh intersection
  int nfeatures3 = optimization->NFeatures();
  if (1) {
    RNArray<GSVFeature *> features = optimization->features;
    for (int i = 0; i < features.NEntries(); i++) {
      GSVFeature *feature = features.Kth(i);
      if (feature->feature_type != GSV_IMAGE_POINT_FEATURE_TYPE) continue;
      if (feature->image_t > 0.0) continue;
      optimization->RemoveFeature(feature);
    }
  }

  // Print statistics
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Features = %d\n", optimization->NFeatures());
    printf("  # Features outside box = %d\n", nfeatures1 - nfeatures2);
    printf("  # Features near border = %d\n", nfeatures2 - nfeatures3);
    printf("  # Features not on mesh = %d\n", nfeatures3 - optimization->NFeatures());
    fflush(stdout);
  }
#endif

  // Print statistics
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Features = %d\n", optimization->NFeatures());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
UpdateFeatures(GSVPoseOptimization *optimization)
{
  // Update scan feature positions based on scan point positions
  if (update_scan_features) {
    if (!UpdateScanFeatures(optimization)) exit(-1);
  }

  // Update pixel feature positions based on ray intersections
  if (update_pixel_features) {
    if (!UpdatePixelFeatures(optimization)) exit(-1);
  }

  // Return success
  return 1;
}


////////////////////////////////////////////////////////////////////////
// Pair creation functions
////////////////////////////////////////////////////////////////////////

static int
CreatePairs(GSVPoseOptimization *optimization)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int saved_npairs = optimization->NPairs();
  if (print_verbose) {
    printf("Creating pairs ...\n");
    fflush(stdout);
  }

  // Create coplanar pairs
  if (create_all_pairs) {
    optimization->CreateAllPairs();
  }
  else if (create_feasible_pairs) {
    optimization->CreateFeasiblePairs();
  }
  else if (create_close_pairs) {
    optimization->CreateClosePairs();
  }

  // Create pixel pixel plane pairs
  if (create_pixel_pixel_plane_pairs) {
    if (!optimization->CreatePixelPixelPlanePairs()) return 0;
  }

  // Create coplanar pairs
  if (create_coplanar_pairs) {
    optimization->CreateCoplanarPairs();
  }

  // Print statistics
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Features = %d\n", optimization->NFeatures());
    printf("  # Pairs = %d\n", optimization->NPairs());
    printf("  # Correspondences = %d\n", optimization->NCorrespondences());
    printf("  # Clusters = %d\n", optimization->NClusters());
    printf("  # New Pairs = %d\n", optimization->NPairs() - saved_npairs);
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
RemovePairs(GSVPoseOptimization *optimization)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int count = 0;
  if (print_verbose) {
    printf("Removing pairs ...\n");
    fflush(stdout);
  }

  // Check which pairs to remove
  if (remove_all_pairs) {
    // Remove all pairs
    count = optimization->NPairs();
    optimization->TruncatePairs(0);
  }
  else if (remove_infeasible_pairs) {
    // Copy array of pairs
    RNArray<GSVFeaturePair *> original_pairs;
    for (int i = 0; i < optimization->NPairs(); i++) {
      GSVFeaturePair *pair = optimization->Pair(i);
      original_pairs.Insert(pair);
    }

    // Remove infeasible pairs
    for (int i = 0; i < original_pairs.NEntries(); i++) {
      GSVFeaturePair *pair = original_pairs.Kth(i);
      GSVFeature *feature0 = pair->Feature(0);
      GSVFeature *feature1 = pair->Feature(1);
      if (optimization->PairScore(feature0, feature1) < 0) {
        optimization->RemovePair(pair);
        count++;
      }
    }
  }

  // Print statistics
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Features = %d\n", optimization->NFeatures());
    printf("  # Pairs = %d\n", optimization->NPairs());
    printf("  # Correspondences = %d\n", optimization->NCorrespondences());
    printf("  # Clusters = %d\n", optimization->NClusters());
    printf("  # Pairs Removed = %d\n", count);
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Correspondence creation functions
////////////////////////////////////////////////////////////////////////

static int
CreateCorrespondences(GSVPoseOptimization *optimization)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int saved_ncorrespondences = optimization->NCorrespondences();
  if (print_verbose) {
    printf("Creating correspondences ...\n");
    fflush(stdout);
  }

  // Create correspondences from pairs
  if (create_all_correspondences) {
    // Create all correspondences
    optimization->CreateAllCorrespondences();
  }
  else if (create_feasible_correspondences) {
    // Create feasible correspondences
    optimization->CreateFeasibleCorrespondences();
  }
  else if (create_close_correspondences) {
    // Create close correspondences
    if (!optimization->CreateCloseCorrespondences()) return 0;
  }
  else if (create_icp_correspondences) {
    if (!optimization->CreateICPCorrespondences()) return 0;
  }

  // Create coplanar correspondences
  if (create_coplanar_correspondences) {
    optimization->CreateCoplanarCorrespondences();
  }

  // Create pixel pixel plane correspondences
  if (create_pixel_pixel_plane_correspondences) {
    if (!optimization->CreatePixelPixelPlaneCorrespondences()) return 0;
  }

  // Create transitive closure correspondences
  if (create_transitive_closure_correspondences) {
    if (!optimization->CreateTransitiveClosureCorrespondences()) return 0;
  }

  // Print statistics
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Features = %d\n", optimization->NFeatures());
    printf("  # Pairs = %d\n", optimization->NPairs());
    printf("  # Correspondences = %d\n", optimization->NCorrespondences());
    printf("  # Clusters = %d\n", optimization->NClusters());
    printf("  # New Correspondences = %d\n", optimization->NCorrespondences() - saved_ncorrespondences);
    printf("  World RMSD = %g\n", optimization->WorldRMSD());
    printf("  Image RMSD = %g\n", optimization->ImageRMSD());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
RemoveCorrespondences(GSVPoseOptimization *optimization)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int count = 0;
  if (print_verbose) {
    printf("Removing correspondences ...\n");
    fflush(stdout);
  }

  // Check which correspondences to remove
  if (remove_all_correspondences) {
    // Remove all correspondences
    count = optimization->NCorrespondences();
    optimization->TruncateCorrespondences(0);
  }
  else if (remove_infeasible_correspondences) {
    // Update correspondences
    optimization->UpdatePixelPositionsFromRayIntersections();
    
    // Copy array of correspondences
    RNArray<GSVFeatureCorrespondence *> original_correspondences;
    for (int i = 0; i < optimization->NCorrespondences(); i++) {
      GSVFeatureCorrespondence *correspondence = optimization->Correspondence(i);
      original_correspondences.Insert(correspondence);
    }

    // Remove infeasible correspondences
    for (int i = 0; i < original_correspondences.NEntries(); i++) {
      GSVFeatureCorrespondence *correspondence = original_correspondences.Kth(i);
      GSVFeature *feature0 = correspondence->Feature(0);
      GSVFeature *feature1 = correspondence->Feature(1);
      if (optimization->CorrespondenceScore(feature0, feature1) < 0) {
        optimization->RemoveCorrespondence(correspondence);
        count++;
      }
    }
  }

  // Print statistics
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Features = %d\n", optimization->NFeatures());
    printf("  # Pairs = %d\n", optimization->NPairs());
    printf("  # Correspondences = %d\n", optimization->NCorrespondences());
    printf("  # Clusters = %d\n", optimization->NClusters());
    printf("  # Correspondences Removed = %d\n", count);
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Cluster creation functions
////////////////////////////////////////////////////////////////////////

static int
CreateClusters(GSVPoseOptimization *optimization)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int saved_nclusters = optimization->NClusters();
  if (print_verbose) {
    printf("Creating clusters ...\n");
    fflush(stdout);
  }

  // Create clusters from pairs
  if (create_all_clusters) {
    optimization->CreateAllClusters();
  }
  else if (create_close_clusters) {
    optimization->CreateCloseClusters();
  }
  else if (create_icp_clusters) {
    optimization->CreateICPClusters();
  }

  // Create coplanar clusters
  if (create_coplanar_clusters) {
    optimization->CreateCoplanarClusters();
  }

  // Create merged clusters
  if (create_merged_clusters) {
    if (!optimization->CreateMergedClusters()) return 0;
  }

  // Create hierarchical clusters
  if (create_hierarchical_clusters) {
    if (!optimization->CreateHierarchicalClusters()) return 0;
  }

  // Print statistics
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Features = %d\n", optimization->NFeatures());
    printf("  # Pairs = %d\n", optimization->NPairs());
    printf("  # Correspondences = %d\n", optimization->NCorrespondences());
    printf("  # Clusters = %d\n", optimization->NClusters());
    printf("  # New Clusters = %d\n", optimization->NClusters() - saved_nclusters);
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
RemoveClusters(GSVPoseOptimization *optimization)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int count = 0;
  if (print_verbose) {
    printf("Removing clusters ...\n");
    fflush(stdout);
  }

  // Check which clusters to remove
  if (remove_all_clusters) {
    // Remove all clusters
    count = optimization->NClusters();
    optimization->TruncateClusters(0);
  }

  // Print statistics
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Features = %d\n", optimization->NFeatures());
    printf("  # Pairs = %d\n", optimization->NPairs());
    printf("  # Correspondences = %d\n", optimization->NCorrespondences());
    printf("  # Clusters = %d\n", optimization->NClusters());
    printf("  # Clusters Removed = %d\n", count);
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Transformation creation functions
////////////////////////////////////////////////////////////////////////

static int
UpdateTransformations(GSVPoseOptimization *optimization)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_verbose) {
    printf("Updating transformations ...\n");
    fflush(stdout);
  }

#if 1
  // Update pixel depths
  if (update_all_transformations || update_pixel_depths || update_cluster_positions || update_cluster_rotations) {
    optimization->UpdatePixelPositionsFromRayIntersections();
  }

  // Update cluster positions
  if (update_all_transformations || update_cluster_positions || update_cluster_rotations) {
    optimization->UpdateClusterPositionsFromFeaturePositions(); 
  }

  // Print errors at start
  if (!Solve(optimization, FALSE, FALSE,  FALSE, FALSE,  FALSE, FALSE,  FALSE, FALSE, FALSE)) return 0;
#else
  // Update pixel depths and cluster positions
  if (update_all_transformations || update_pixel_depths || update_cluster_positions || update_cluster_rotations) {
    if (!Solve(optimization, FALSE, FALSE,  FALSE, FALSE,  FALSE, FALSE,  TRUE, TRUE, TRUE)) return 0;
  }
#endif

  // Update paths
  if (update_all_transformations || update_path_translations) {
    if (!Solve(optimization, FALSE, FALSE,  FALSE, FALSE,  TRUE, FALSE,  
      update_all_transformations || update_cluster_positions,
      update_all_transformations || update_cluster_rotations,
      update_all_transformations || update_pixel_depths)) return 0;
  }
  if (update_all_transformations || update_path_rotations) {
    if (!Solve(optimization, FALSE, FALSE,  FALSE, FALSE,  update_path_translations, TRUE,  
      update_all_transformations || update_cluster_positions,
      update_all_transformations || update_cluster_rotations,
      update_all_transformations || update_pixel_depths)) return 0;
  }
    
  if (update_all_transformations || update_laser_translations) {
    if (!Solve(optimization, TRUE, FALSE,  FALSE, FALSE,  FALSE, FALSE,  
      update_all_transformations || update_cluster_positions,
      update_all_transformations || update_cluster_rotations,
      update_all_transformations || update_pixel_depths)) return 0;
  }
  if (update_all_transformations || update_laser_rotations) {
    if (!Solve(optimization, update_laser_translations, TRUE,  FALSE, FALSE,  FALSE, FALSE,  
      update_all_transformations || update_cluster_positions,
      update_all_transformations || update_cluster_rotations,
      update_all_transformations || update_pixel_depths)) return 0;
  }

  if (update_all_transformations || update_camera_translations) {
    if (!Solve(optimization, FALSE, FALSE,  TRUE, FALSE,  FALSE, FALSE,   
      update_all_transformations || update_cluster_positions,
      update_all_transformations || update_cluster_rotations,
      update_all_transformations || update_pixel_depths)) return 0;
  }
  if (update_all_transformations || update_camera_rotations) {
    if (!Solve(optimization, FALSE, FALSE,  update_camera_translations, TRUE,  FALSE, FALSE,   
      update_all_transformations || update_cluster_positions,
      update_all_transformations || update_cluster_rotations,
      update_all_transformations || update_pixel_depths)) return 0;
  }
  if (update_all_transformations || update_camera_rotations) {
    if (!Solve(optimization, FALSE, FALSE,  update_camera_translations, TRUE,  FALSE, FALSE,  
      update_all_transformations || update_cluster_positions,
      update_all_transformations || update_cluster_rotations,
      update_all_transformations || update_pixel_depths)) return 0;
  }

  // Print statistics
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Features = %d\n", optimization->NFeatures());
    printf("  # Pairs = %d\n", optimization->NPairs());
    printf("  # Correspondences = %d\n", optimization->NCorrespondences());
    printf("  # Clusters = %d\n", optimization->NClusters());
    printf("  World RMSD = %g\n", optimization->WorldRMSD());
    printf("  Image RMSD = %g\n", optimization->ImageRMSD());
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Scanline data functions
////////////////////////////////////////////////////////////////////////

static int
ReadScanData(GSVScan *scan, RNScalar timestamp) 
{
  // Check scan 
  if (!scan) return 1;
  ScanData *scan_data = (ScanData *) scan->Data();
  if (scan_data) { scan_data->timestamp = timestamp; return 1; }
  if (scan->NScanlines() == 0) return 1;
   
  // Get convenient variables
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  GSVRun *run = segment->Run();
  if (!run) return 0;
  int segment_index = segment->RunIndex();
  static const char *image_directory = "gsv_data/laser_images";
  char image_name[4096];

  // Read points
  if (!scan->ReadPoints()) {
    fprintf(stderr, "Unable to read points for %s %02d %02d\n", run->Name(), segment_index, scan_index);
    return 0;
  }

  // Fill scan data
  scan_data = new ScanData();
  scan->SetData(scan_data);
  scan_data->scan = scan;
  scan_data->mesh = NULL; 
  scan_data->timestamp = timestamp;
  scan_data->index = scan->SceneIndex();

  // Read mesh
  if (load_meshes) {
    scan_data->mesh = scan->Mesh();
  }

  // Read segmentation grids
  if (load_segments) {
    sprintf(image_name, "%s/%s/%02d_%02d_DA_Scanline.grd", image_directory, run->Name(), segment_index, scan_index);
    scan_data->scanline_grid.Read(image_name);
    sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionX.grd", image_directory, run->Name(), segment_index, scan_index);
    scan_data->position_x_grid.Read(image_name);
    sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionY.grd", image_directory, run->Name(), segment_index, scan_index);
    scan_data->position_y_grid.Read(image_name);
    sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionZ.grd", image_directory, run->Name(), segment_index, scan_index);
    scan_data->position_z_grid.Read(image_name);
    sprintf(image_name, "%s/%s/%02d_%02d_DA_Timestamp.grd", image_directory, run->Name(), segment_index, scan_index);
    scan_data->timestamp_grid.Read(image_name);
    sprintf(image_name, "%s/%s/%02d_%02d_DA_%s.grd", image_directory, run->Name(), segment_index, scan_index, segment_grid_name);
    scan_data->segment_id_grid.Read(image_name); 
  }
  
  // Read planar grids
  if (load_planar_grids) {
    for (int j = 0; j < 100; j++) {
      sprintf(image_name, "%s/%s/%02d_%02d_%06d_ST_Displacement.grd", image_directory, run->Name(), segment_index, scan_index, j);
      if (!RNFileExists(image_name)) break;
      R3PlanarGrid *planar_grid = new R3PlanarGrid();
      if (!planar_grid->ReadFile(image_name)) { delete planar_grid; continue; }
      if (planar_grid->Cardinality() < 16 * 1024) { delete planar_grid; continue; }
      if (fabs(planar_grid->Plane().Normal().Z()) > 0.1) { delete planar_grid; continue; }
      scan_data->planar_grids.Insert(planar_grid);
    }
  }

  // Return success
  return 1;
}



static int
ReleaseScanData(GSVScan *scan)
{
  // Check scan
  if (!scan) return 1;
  ScanData *scan_data = (ScanData *) scan->Data();
  if (!scan_data) return 1;

  // Release points
  if (!scan->ReleasePoints()) {
    GSVSegment *segment = scan->Segment();
    GSVRun *run = segment->Run();
    fprintf(stderr, "Unable to release points for %s %02d %02d\n", run->Name(), segment->RunIndex(), scan->SegmentIndex());
    return 0;
  }

  // Delete mesh
  if (scan_data->mesh) delete scan_data->mesh;

  // Delete planar grids
  for (int j = 0; j < scan_data->planar_grids.NEntries(); j++) {
    R3PlanarGrid *planar_grid = scan_data->planar_grids.Kth(j);
    delete planar_grid;
  }      

  // Delete scan data
  delete scan_data;
  scan->SetData(NULL);

  // Return success
  return 1;
}



static int
ToggleScanData(GSVSegment *segment, RNScalar timestamp) 
{
  // Check segment
  if (!segment) return 1;

  // Fill scan data
  for (int i = 0; i < segment->NScans(); i++) {
    GSVScan *scan = segment->Scan(i);
    if (scan->Data()) ReleaseScanData(scan);
    else ReadScanData(scan, timestamp);
  }      

  // Return success
  return 1;
}



static int
SelectScanData(GSVSegment *selected_segment, RNScalar timestamp) 
{
  // Fill scan data for segment, and empty all others
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (segment == selected_segment) ReadScanData(scan, timestamp);
        else ReleaseScanData(scan);
      }
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Color functions
////////////////////////////////////////////////////////////////////////

static void
LoadSegmentColor(GSVSegment *segment, int color_scheme)
{
  // Check segment
  if (!segment) return;

  // Segment colors
  const int max_segment_colors = 12;
  const RNRgb segment_colors[max_segment_colors] = {
    RNRgb(1.0, 0.0, 0.0), RNRgb(0.0, 1.0, 0.0), RNRgb(0.0, 0.0, 1.0), 
    RNRgb(0.8, 0.8, 0.0), RNRgb(0.0, 0.8, 0.8), RNRgb(0.8, 0.0, 0.8),
    RNRgb(0.8, 0.5, 0.2), RNRgb(0.8, 0.2, 0.5), RNRgb(0.2, 0.8, 0.5),
    RNRgb(0.5, 0.8, 0.2), RNRgb(0.2, 0.5, 0.8), RNRgb(0.5, 0.2, 0.8),

  };

  // Check color scheme
  if (color_scheme == SEGMENT_COLOR_SCHEME) {
    // Find run index and segment index
    GSVRun *run = segment->Run();
    if (!run) glColor3d(0.0, 0.0, 0.0);
    int segment_index = segment->RunIndex();
    int run_index = run->SceneIndex();

    // Load segment color
    RNLoadRgb(segment_colors[(7*run_index + segment_index) % max_segment_colors]);
  }
  else if (color_scheme == PICK_COLOR_SCHEME) {
    int segment_index = segment->SceneIndex();
    unsigned char rgba[4];
    rgba[0] = (segment_index >> 16) & 0xFF;
    rgba[1] = (segment_index >> 8) & 0xFF;
    rgba[2] = segment_index & 0xFF;
    rgba[3] = 0xFE;
    glColor4ubv(rgba);
  }
}



static void
LoadScanlineColor(GSVScanline *scanline, int side, int color_scheme)
{
  // Check scanline
  if (!scanline) return;

  // Check color scheme
  if (color_scheme == SEGMENT_COLOR_SCHEME) {
    GSVScan *scan = scanline->Scan();
    if (!scan) glColor3d(0.0, 0.0, 0.0);
    GSVSegment *segment = scan->Segment();
    LoadSegmentColor(segment, color_scheme);
  }
  else if (color_scheme == PICK_COLOR_SCHEME) {
    int scanline_index = 0;
    ScanlineData *scanline_data = (ScanlineData *) scanline->Data();
    if (scanline_data) scanline_index = scanline_data->index;
    else scanline_index = scanline->SceneIndex();
    unsigned char rgba[4];
    rgba[0] = (scanline_index >> 16) & 0xFF;
    rgba[1] = (scanline_index >> 8) & 0xFF;
    rgba[2] = scanline_index & 0xFF;
    rgba[3] = 0xFD;
    glColor4ubv(rgba);
  }
}



static void
LoadImageColor(GSVImage *image, int side, int color_scheme)
{
  // Check image
  if (!image) return;

  // Check color scheme
  if (color_scheme == SEGMENT_COLOR_SCHEME) {
    if (image == selected_image[side]) { glColor3d(1.0, 1.0, 1.0); return; }
    GSVPanorama *panorama = image->Panorama();
    if (!panorama) { glColor3d(0.0, 0.0, 0.0); return; }
    GSVSegment *segment = panorama->Segment();
    LoadSegmentColor(segment, color_scheme);
  }
  else if (color_scheme == PICK_COLOR_SCHEME) {
    int image_index = 0;
    ImageData *image_data = (ImageData *) image->Data();
    if (image_data) image_index = image_data->index;
    else image_index = image->SceneIndex();
    unsigned char rgba[4];
    rgba[0] = (image_index >> 16) & 0xFF;
    rgba[1] = (image_index >> 8) & 0xFF;
    rgba[2] = image_index & 0xFF;
    rgba[3] = 0xFC;
    glColor4ubv(rgba);
  }
}



static void
LoadFeatureColor(GSVFeature *feature, int side, int color_scheme, int feature_index)
{
  // Check feature
  if (!feature) return;

  // Check color scheme
  if (color_scheme == SEGMENT_COLOR_SCHEME) {
    if (feature == selected_feature[side]) { glColor3d(1.0, 1.0, 1.0); return; }
    GSVScanline *scanline = feature->scanline;
    if (!scanline) { glColor3d(0.0, 0.0, 0.0); return; }
    LoadScanlineColor(scanline, side, color_scheme);
  }
  else if (color_scheme == PICK_COLOR_SCHEME) {
    unsigned char rgba[4];
    rgba[0] = (feature_index >> 16) & 0xFF;
    rgba[1] = (feature_index >> 8) & 0xFF;
    rgba[2] = feature_index & 0xFF;
    rgba[3] = 0xFB;
    glColor4ubv(rgba);
  }
}



static void
LoadClusterColor(GSVFeatureCluster *cluster, int side, int color_scheme, int cluster_index)
{
  // Check cluster
  if (!cluster) return;

  // Check color scheme
  if (color_scheme == SEGMENT_COLOR_SCHEME) {
    glColor3d(0.5, 0.5, 0.5);
  }
  else if (color_scheme == PICK_COLOR_SCHEME) {
    unsigned char rgba[4];
    rgba[0] = (cluster_index >> 16) & 0xFF;
    rgba[1] = (cluster_index >> 8) & 0xFF;
    rgba[2] = cluster_index & 0xFF;
    rgba[3] = 0xFA;
    glColor4ubv(rgba);
  }
}



static void
LoadHeatMapColor(RNScalar value)
{
  // Compute color
  GLdouble r, g, b;
  if (value < 0.2) {
    value *= 4;
    r = 1;
    g = value;
    b = 0;
  }
  else if (value < 0.4) {
    value = (value - 0.25) * 4;
    r = 1 - value;
    g = 1;
    b = 0;
  }
  else if (value < 0.6) {
    value = (value - 0.5) * 4;
    r = 0;
    g = 1;
    b = value;
  }
  else if (value < 0.8) {
    value = (value - 0.75) * 4;
    r = 0;
    g = 1 - value;
    b = 1;
  }
  else {
    value = (value - 0.75) * 4;
    r = 0;
    g = 0;
    b = 1;
  }

  // Load color
  RNLoadRgb(r, g, b);
}



////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////

static R3Vector
ComputeNormal(GSVScanline *scanline, int point_index)
{
  // Get convenient variables
  GSVScan *scan = scanline->Scan();
  int scanline_index = scanline->ScanIndex();
  R3Point p = scanline->PointPosition(point_index);

  // Find point on previous row
  R3Point pA = (point_index > 0) ? scanline->PointPosition(point_index-1) : p;

  // Find point on next row
  R3Point pB = (point_index < scanline->NPoints()-1) ? scanline->PointPosition(point_index+1) : p;

  // Find point on previous scanline
  R3Point p0 = p;
  GSVScanline *scanline0 = NULL;
  if (scanline_index > 0) {
    RNScalar closest_dd = FLT_MAX;
    scanline0 = scan->Scanline(scanline_index-1);
    for (int i = 0; i < scanline0->NPoints(); i++) {
      const R3Point& q = scanline0->PointPosition(i);
      RNLength dd = R3SquaredDistance(q, p);
      if (dd < closest_dd) { p0 = q; closest_dd = dd; }
    }
  }

  // Find point on next scanline
  R3Point p1 = p;
  GSVScanline *scanline1 = NULL;
  if (scanline_index < scan->NScanlines()-1) {
    RNScalar closest_dd = FLT_MAX;
    scanline1 = scan->Scanline(scanline_index-1);
    for (int i = 0; i < scanline1->NPoints(); i++) {
      const R3Point& q = scanline1->PointPosition(i);
      RNLength dd = R3SquaredDistance(q, p);
      if (dd < closest_dd) { p0 = q; closest_dd = dd; }
    }
  }

  // Compute normal with cross product
  R3Vector v01 = p1 - p0;
  v01.Normalize();
  R3Vector vAB = pB - pA;
  vAB.Normalize();
  R3Vector n = v01 % vAB;
  n.Normalize();

  // Flip to face laser
  R3Vector v = p - scanline->Pose().Viewpoint();
  if (n.Dot(v) > 0) n.Flip();

  // Return normal
  return n;
}



////////////////////////////////////////////////////////////////////////
// Draw functions
////////////////////////////////////////////////////////////////////////

static void
DrawFeature(GSVFeature *feature, int side, RNBoolean transform)
{
  // Set transformation
  static R3Affine transformation = R3identity_affine;
  if (transform) {
    transformation = optimization->OptimizedTransformation(feature);
    transformation.Push();
  }

  // Draw feature
  feature->Draw();

  // Add stuff for selected feature
  if (feature == selected_feature[side]) {
    if ((feature->image) && (viewpoint_scheme[side] != 1)) {
      // Draw line to viewpoint
      glLineWidth(1);
      glColor3d(1.0, 1.0, 1.0);
      glBegin(GL_LINES);
      R3LoadPoint(feature->image->Pose().Viewpoint());
      R3LoadPoint(feature->scan_position);
      glEnd();
    }
  }

  // Reset transformation
  if (transform) {
    transformation.Pop();
  }
}



static void
DrawCluster(GSVFeatureCluster *cluster, int side, RNBoolean transform)
{
  // Set transformation
  static R3Affine transformation = R3identity_affine;
  if (transform) {
    transformation = optimization->OptimizedTransformation(cluster);
    transformation.Push();
  }

  // Draw cluster
  cluster->Draw();

  // Reset transformation
  if (transform) {
    transformation.Pop();
  }
}



static void 
DrawImage(GSVImage *image, int side, RNBoolean transform)
{
  // Persistent variables
  static GLuint texture_id[2] = { 0, 0 };
  static GSVImage *last_image[2] = { NULL, NULL };

  // Load texture
  if (image != last_image[side]) {
    last_image[side] = image;
    if (texture_id[side] > 0) glDeleteTextures(1, &texture_id[side]);
    texture_id[side] = 0;

    // Create texture
    R2Image *texture = image->UndistortedImage();
    if (texture) {
      // Create texture
      glGenTextures(1, &texture_id[side]);
      glBindTexture(GL_TEXTURE_2D, texture_id[side]);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
      gluBuild2DMipmaps(GL_TEXTURE_2D, 3, texture->Width(), texture->Height(), GL_RGB, GL_UNSIGNED_BYTE, texture->Pixels() );
      assert(texture_id[side] != 0);
      delete texture;
    }
  }

  // Draw textured polygon
  if (texture_id[side] > 0) {
    // Set transformation
    static R3Affine transformation = R3identity_affine;
    if (transform) {
      R2Point center(image->Width()/2, image->Height()/2);
      transformation = optimization->OptimizedTransformation(image, center);
      transformation.Push();
    }

    // Disable depth buffer 
    // glDisable(GL_DEPTH_TEST);
    // glDepthMask(FALSE);

    // Enable texture
    glBindTexture(GL_TEXTURE_2D, texture_id[side]);
    glEnable(GL_TEXTURE_2D);

    // Get camera parameters
    GSVCamera *camera = image->Camera();
    RNAngle xfov = 0.5 * camera->XFov();
    RNAngle yfov = 0.5 * camera->YFov();
    GSVPose pose = image->Pose();
    R3Point viewpoint = pose.Viewpoint();
    R3Vector towards = pose.Towards();
    R3Vector up = pose.Up();
    R3Vector right = pose.Right();

    // Draw textured polygon
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
    // glEnable(GL_DEPTH_TEST);
    // glDepthMask(TRUE);

    // Reset transformation
    if (transform) {
      transformation.Pop();
    }
  }
}



static void 
DrawFeatures(int side, RNBoolean transform, int color_scheme)
{
  // Draw features
  for (int i = 0; i < optimization->NFeatures(); i++) {
    GSVFeature *feature = optimization->Feature(i);

    // Check if should draw feature
    if (feature->image && selected_image[side] && (feature->image != selected_image[side])) continue;

    // Draw feature
    LoadFeatureColor(feature, side, color_scheme, feature->index);
    DrawFeature(feature, side, transform);
  }
}



static void 
DrawPairs(int side, RNBoolean transform, int color_scheme)
{
  // Draw pairs
  for (int i = 0; i < optimization->NPairs(); i++) {
    GSVFeaturePair *pair = optimization->Pair(i);
    GSVFeature *feature0 = pair->Feature(0);
    GSVFeature *feature1 = pair->Feature(1);

    // Check if should draw pair
    if (selected_image[side]) {
      if (feature0->image && (feature0->image != selected_image[side]) && 
          feature1->image && (feature1->image != selected_image[side])) continue;
    }
    if (selected_image[1-side]) {
      if (feature0->image && (feature0->image != selected_image[1-side]) && 
          feature1->image && (feature1->image != selected_image[1-side])) continue;
    }

    // Draw features
    LoadFeatureColor(feature0, side, color_scheme, feature0->index);
    DrawFeature(feature0, side, transform);
    LoadFeatureColor(feature1, side, color_scheme, feature1->index);
    DrawFeature(feature1, side, transform);

    // Draw line between features
    glBegin(GL_LINES);
    LoadFeatureColor(feature0, side, color_scheme, 0);
    R3LoadPoint(optimization->WorldPosition(feature0, transform));
    LoadFeatureColor(feature1, side, color_scheme, 0);
    R3LoadPoint(optimization->WorldPosition(feature1, transform));
    glEnd();
  }
}


static void 
DrawCorrespondences(int side, RNBoolean transform, int color_scheme)
{
  // Draw correspondences
  for (int i = 0; i < optimization->NCorrespondences(); i++) {
    GSVFeatureCorrespondence *correspondence = optimization->Correspondence(i);
    GSVFeature *feature0 = correspondence->Feature(0);
    GSVFeature *feature1 = correspondence->Feature(1);

    // Check if should draw correspondence
    if (selected_image[side]) {
      if (feature0->image && (feature0->image != selected_image[side]) && 
          feature1->image && (feature1->image != selected_image[side])) continue;
    }
    if (selected_image[1-side]) {
      if (feature0->image && (feature0->image != selected_image[1-side]) && 
          feature1->image && (feature1->image != selected_image[1-side])) continue;
    }

    // Draw features
    LoadFeatureColor(feature0, side, color_scheme, feature0->index);
    DrawFeature(feature0, side, transform);
    LoadFeatureColor(feature1, side, color_scheme, feature1->index);
    DrawFeature(feature1, side, transform);

    // Draw line between features
    glBegin(GL_LINES);
    LoadFeatureColor(feature0, side, color_scheme, 0);
    R3LoadPoint(optimization->WorldPosition(feature0, transform));
    LoadFeatureColor(feature1, side, color_scheme, 0);
    R3LoadPoint(optimization->WorldPosition(feature1, transform));
    glEnd();

#if 0
    // Draw text with distance between features
    char buffer[1024];
    R3Point p0 = optimization->WorldPosition(feature0, transform);
    R3Point p1 = optimization->WorldPosition(feature1, transform);
    R3Point p = 0.5 * (p0 + p1);
    RNScalar dd = optimization->SquaredWorldDistance(feature0, feature1, transform);
    RNScalar d = sqrt(dd);
    RNScalar r = (d > 1.0) ? d - 1.0 : 0;
    RNScalar g = ((d > 0.5) && (d < 1.5)) ? d - 0.5 : 0;
    RNScalar b = ((d > 0.0) && (d < 1.0)) ? 0.5 + 0.5 * d : 0.5;
    glColor3d(r, g, b);
    sprintf(buffer, "%.0f", 10 * d);
    DrawText(p, buffer);
#endif
  }
}



static void 
DrawClusters(int side, RNBoolean transform, int color_scheme)
{
  // Draw clusters
  for (int i = 0; i < optimization->NClusters(); i++) {
    GSVFeatureCluster *cluster = optimization->Cluster(i);
    const R3Point& cluster_position = optimization->WorldPosition(cluster, transform);

    // Draw cluster
    LoadClusterColor(cluster, side, color_scheme, cluster->index);
    DrawCluster(cluster, side, transform);
    
    // Draw lines to features
    for (int j = 0; j < cluster->NFeatures(); j++) {
      GSVFeature *feature = cluster->Feature(j);
      const R3Point& feature_position = optimization->WorldPosition(feature, transform);
      LoadFeatureColor(feature, side, color_scheme, feature->index);
      DrawFeature(feature, side, transform);
      if (cluster->scan_position != R3unknown_point) {
        glBegin(GL_LINES);
        R3LoadPoint(cluster_position);
        R3LoadPoint(feature_position);
        glEnd();
      }
    }

    // Draw lines to children clusters
    glLineWidth(5);
    for (int j = 0; j < cluster->NClusters(); j++) {
      GSVFeatureCluster *child = cluster->Cluster(j);
      const R3Point& child_position = optimization->WorldPosition(child, transform);
      LoadClusterColor(child, side, color_scheme, child->index);
      if (child->scan_position != R3unknown_point) {
        glBegin(GL_LINES);
        R3LoadPoint(cluster_position);
        R3LoadPoint(child_position);
        glEnd();
      }
    }
    glLineWidth(1);
  }
}



static void 
DrawScans(int side, RNBoolean transform, int color_scheme)
{
  // Draw scans
  glBegin(GL_POINTS);
  if (color_scheme == SEGMENT_COLOR_SCHEME) {
    for (int ir = 0; ir < scene->NRuns(); ir++) {
      GSVRun *run = scene->Run(ir);
      for (int is = 0; is < run->NSegments(); is++) {
        GSVSegment *segment = run->Segment(is);
        LoadSegmentColor(segment, color_scheme);
        for (int ia = 0; ia < segment->NScans(); ia++) {
          if (ia == 1) continue;
          GSVScan *scan = segment->Scan(ia);
          if (scan->NScanlines() == 0) continue;
          ScanData *scan_data = (ScanData *) scan->Data();
          if (!scan_data) continue;
          for (int ie = 0; ie < scan->NScanlines(); ie++) {
            GSVScanline *scanline = scan->Scanline(ie);
            double timestamp_delta = (scan_data->timestamp >= 0) ? fabs(scanline->Timestamp() - scan_data->timestamp) : 0;
            if (timestamp_delta > timestamp_radius) continue;
            const R3Affine& scanline_transformation = optimization->OptimizedTransformation(scanline);
            for (int ik = 0; ik < scanline->NPoints(); ik++) {
              R3Point position = scanline->PointPosition(ik);
              if (transform) position.Transform(scanline_transformation);
              R3LoadPoint(position);
            }
          }
        }
      }
    }
  }
  else if (color_scheme == VIEWPOINT_DEPTH_COLOR_SCHEME) {
    for (int ir = 0; ir < scene->NRuns(); ir++) {
      GSVRun *run = scene->Run(ir);
      for (int is = 0; is < run->NSegments(); is++) {
        GSVSegment *segment = run->Segment(is);
        for (int ia = 0; ia < segment->NScans(); ia++) {
          if (ia == 1) continue;
          GSVScan *scan = segment->Scan(ia);
          if (scan->NScanlines() == 0) continue;
          ScanData *scan_data = (ScanData *) scan->Data();
          if (!scan_data) continue;
          for (int ie = 0; ie < scan->NScanlines(); ie++) {
            GSVScanline *scanline = scan->Scanline(ie);
            double timestamp_delta = (scan_data->timestamp >= 0) ? fabs(scanline->Timestamp() - scan_data->timestamp) : 0;
            if (timestamp_delta > timestamp_radius) continue;
            const R3Affine& scanline_transformation = optimization->OptimizedTransformation(scanline);
            R3Point viewpoint = scanline->Pose().Viewpoint();
            R3Vector towards = scanline->Pose().Towards();
            if (transform) viewpoint.Transform(scanline_transformation);
            if (transform) towards.Transform(scanline_transformation);
            for (int ik = 0; ik < scanline->NPoints(); ik++) {
              R3Point position = scanline->PointPosition(ik);
              if (transform) position.Transform(scanline_transformation);
              R3Vector vector = position - viewpoint;
              RNScalar depth = vector.Dot(towards);
              RNScalar value = depth / image_plane_distance;
              LoadHeatMapColor(value);
              R3LoadPoint(position);
            }
          }
        }
      }
    }
  }
  else if (color_scheme == HEIGHT_COLOR_SCHEME) {
    for (int ir = 0; ir < scene->NRuns(); ir++) {
      GSVRun *run = scene->Run(ir);
      for (int is = 0; is < run->NSegments(); is++) {
        GSVSegment *segment = run->Segment(is);
        for (int ia = 0; ia < segment->NScans(); ia++) {
          if (ia == 1) continue;
          GSVScan *scan = segment->Scan(ia);
          if (scan->NScanlines() == 0) continue;
          ScanData *scan_data = (ScanData *) scan->Data();
          if (!scan_data) continue;
          for (int ie = 0; ie < scan->NScanlines(); ie++) {
            GSVScanline *scanline = scan->Scanline(ie);
            double timestamp_delta = (scan_data->timestamp >= 0) ? fabs(scanline->Timestamp() - scan_data->timestamp) : 0;
            if (timestamp_delta > timestamp_radius) continue;
            const R3Affine& scanline_transformation = optimization->OptimizedTransformation(scanline);
            RNCoord ground_z = scanline->EstimatedGroundZ();
            for (int ik = 0; ik < scanline->NPoints(); ik++) {
              R3Point position = scanline->PointPosition(ik);
              RNScalar height = position.Z() - ground_z;
              if (transform) position.Transform(scanline_transformation);
              RNScalar value = 2 * height / image_plane_distance;
              LoadHeatMapColor(value);
              R3LoadPoint(position);
            }
          }
        }
      }
    }
  }
  else if (color_scheme == PICK_COLOR_SCHEME) {
    for (int ir = 0; ir < scene->NRuns(); ir++) {
      GSVRun *run = scene->Run(ir);
      for (int is = 0; is < run->NSegments(); is++) {
        GSVSegment *segment = run->Segment(is);
        for (int ia = 0; ia < segment->NScans(); ia++) {
          if (ia == 1) continue;
          GSVScan *scan = segment->Scan(ia);
          if (scan->NScanlines() == 0) continue;
          ScanData *scan_data = (ScanData *) scan->Data();
          if (!scan_data) continue;
          for (int ie = 0; ie < scan->NScanlines(); ie++) {
            GSVScanline *scanline = scan->Scanline(ie);
            double timestamp_delta = (scan_data->timestamp >= 0) ? fabs(scanline->Timestamp() - scan_data->timestamp) : 0;
            if (timestamp_delta > timestamp_radius) continue;
            const R3Affine& scanline_transformation = optimization->OptimizedTransformation(scanline);
            ScanlineData *scanline_data = (ScanlineData *) scanline->Data();
            if (!scanline_data) continue;
            int scanline_index = scanline_data->index;
            unsigned char rgba[4];
            rgba[0] = (scanline_index >> 16) & 0xFF;
            rgba[1] = (scanline_index >> 8) & 0xFF;
            rgba[2] = scanline_index & 0xFF;
            for (int ik = 0; ik < scanline->NPoints(); ik++) {
              R3Point position = scanline->PointPosition(ik);
              if (transform) position.Transform(scanline_transformation);
              rgba[3] = (ik + 1) & 0xFF;
              glColor4ubv(rgba);
              R3LoadPoint(position);
            }
          }
        }
      }
    }
  }
  glEnd();
}



static void 
DrawMeshes(int side, RNBoolean transform, int color_scheme)
{
  // Colors
  const int max_mesh_colors = 12;
  const RNRgb mesh_colors[max_mesh_colors] = {
    RNRgb(1.0, 0.0, 0.0), RNRgb(0.0, 1.0, 0.0), RNRgb(0.0, 0.0, 1.0), 
    RNRgb(0.8, 0.8, 0.0), RNRgb(0.0, 0.8, 0.8), RNRgb(0.8, 0.0, 0.8),
    RNRgb(0.8, 0.5, 0.2), RNRgb(0.8, 0.2, 0.5), RNRgb(0.2, 0.8, 0.5),
    RNRgb(0.5, 0.8, 0.2), RNRgb(0.2, 0.5, 0.8), RNRgb(0.5, 0.2, 0.8),

  };

  // Set OpenGL modes
  glEnable(GL_LIGHTING);
  static GLfloat material[] = { 0.8, 0.8, 0.8, 1.0 };

  // Draw meshes
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      RNRgb color = mesh_colors[(7*ir + is) % max_mesh_colors];
      material[0] = color[0]; material[1] = color[1]; material[2] = color[2];
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, material); 
      for (int ia = 0; ia < segment->NScans(); ia++) {
        if (ia == 1) continue;
        GSVScan *scan = segment->Scan(ia);
        ScanData *scan_data = (ScanData *) scan->Data();
        if (!scan_data) continue;
        if (scan->NScanlines() == 0) continue;
        if (!scan_data->mesh) continue;
        scan_data->mesh->Draw();
      }
    }
  }

  // Reset OpenGL modes
  glDisable(GL_LIGHTING);
}



static void 
DrawLasers(int side, RNBoolean transform, int color_scheme)
{
  // Scanline drawing parameters
  static int scanline_step = 0;
  if (scanline_step == 0) {
    const int max_scanlines = 32 * 1024;
    scanline_step = 1 + scene->NScanlines() / max_scanlines;
  }

  // Draw pose for every laser
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        for (int ie = 0; ie < scan->NScanlines(); ie += scanline_step) {
          GSVScanline *scanline = scan->Scanline(ie);
          GSVPose pose = scanline->Pose();
          R3Point position = pose.Viewpoint();
          R3Vector towards = pose.Towards();
          R3Vector up = pose.Up();
          if (transform) {
            R3Affine transformation = optimization->OptimizedTransformation(scanline);
            position.Transform(transformation);
            towards.Transform(transformation);
            up.Transform(transformation);
          }
          if (scanline == selected_scanline[side]) glColor3d(1.0, 1.0, 1.0);
          else LoadScanlineColor(scanline, side, color_scheme);
          R3Span(position, position + towards).Draw();
          R3Span(position, position + 0.5*up).Draw();
        }
      }
    }
  }
}



static void 
DrawCameras(int side, RNBoolean transform, int color_scheme)
{
  // Draw camera pose for every image
  int image_index = 0;
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ip = 0; ip < segment->NPanoramas(); ip++) {
        GSVPanorama *panorama = segment->Panorama(ip);
        for (int ii = 0; ii < panorama->NImages(); ii++) {
          GSVImage *image = panorama->Image(ii);
          GSVPose pose = image->Pose();
          R3Point position = pose.Viewpoint();
          R3Vector towards = pose.Towards();
          R3Vector up = pose.Up();
          if (transform) {
            R2Point center(image->Width()/2, image->Height()/2);
            R3Affine transformation = optimization->OptimizedTransformation(image, center);
            position.Transform(transformation);
            towards.Transform(transformation);
            up.Transform(transformation);
          }
          LoadImageColor(image, side, color_scheme);
          R3Span(position, position + towards).Draw();
          R3Span(position, position + 0.5 * up).Draw();
          image_index++;
        }
      }
    }
  }
}



static void 
DrawPaths(int side, RNBoolean transform, int color_scheme)
{
  // Draw path for every segment
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      SegmentData *segment_data = (SegmentData *) segment->Data();
      if (!segment_data) continue;
      GSVPath *path = segment_data->path;
      LoadSegmentColor(segment, color_scheme);

      // Draw vertices
      glPointSize(3);
      glBegin(GL_POINTS);
      for (int i = 0; i < path->NVertices(); i++) {
        R3Point point = path->VertexPosition(i);
        if (transform) point += path->VertexTranslation(i);
        R3LoadPoint(point);
      }
      glEnd();
      glPointSize(1);

      // Draw spline curve
      glBegin(GL_LINE_STRIP);
      for (RNScalar u = 0; u <= path->VertexParameter(path->NVertices()-1); u += 0.1) {
        R3Point point = path->spline->PointPosition(u);
        if (transform) point += path->Translation(u);
        R3LoadPoint(point);
      }
      glEnd();
    }
  }
}



static void 
DrawSegments(int side, RNBoolean transform, int color_scheme)
{
  // Colors
  const int max_segment_colors = 12;
  const RNRgb segment_colors[max_segment_colors] = {
    RNRgb(1.0, 0.0, 0.0), RNRgb(0.0, 1.0, 0.0), RNRgb(0.0, 0.0, 1.0), 
    RNRgb(0.8, 0.8, 0.0), RNRgb(0.0, 0.8, 0.8), RNRgb(0.8, 0.0, 0.8),
    RNRgb(0.5, 0.8, 0.2), RNRgb(0.2, 0.5, 0.8), RNRgb(0.5, 0.2, 0.8),
    RNRgb(0.8, 0.5, 0.2), RNRgb(0.8, 0.2, 0.5), RNRgb(0.2, 0.8, 0.5),
  };

  // Check if loaded segment data
  if (!load_segments) return;

  // Draw segments
  glDisable(GL_LIGHTING);
  glBegin(GL_POINTS);
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        ScanData *scan_data = (ScanData *) scan->Data();
        if (!scan_data) continue;
        if (scan->NScanlines() == 0) continue;
        if (ia == 1) continue;
        if (scan_data->position_x_grid.YResolution() == 0) continue;
        for (int i = 0; i < scan_data->position_x_grid.XResolution(); i++) {
          RNScalar timestamp = scan_data->timestamp_grid.GridValue(i, 0);
          if (timestamp == R2_GRID_UNKNOWN_VALUE) continue;
          double timestamp_delta = (scan_data->timestamp >= 0) ? fabs(timestamp - scan_data->timestamp) : 0;
          if (timestamp_delta > timestamp_radius) continue;
          for (int j = 0; j < scan_data->position_x_grid.YResolution(); j++) {
            RNScalar x = scan_data->position_x_grid.GridValue(i, j);
            if (x == R2_GRID_UNKNOWN_VALUE) continue;
            RNScalar y = scan_data->position_y_grid.GridValue(i, j);
            if (y == R2_GRID_UNKNOWN_VALUE) continue;
            RNScalar z = scan_data->position_z_grid.GridValue(i, j);
            if (z == R2_GRID_UNKNOWN_VALUE) continue;
            RNScalar scanline_index_value = scan_data->scanline_grid.GridValue(i, j);
            if (scanline_index_value == R2_GRID_UNKNOWN_VALUE) continue;
            RNScalar segment_index_value = scan_data->segment_id_grid.GridValue(i, j);
            if (segment_index_value == R2_GRID_UNKNOWN_VALUE) continue;
            int scanline_index = (int) (scanline_index_value + 0.5);
            int segment_index = (int) (segment_index_value + 0.5);
            GSVScanline *scanline = scan->Scanline(scanline_index);
            R3Point position(x, y, z);
            if (transform) {
              const R3Affine& scanline_transformation = optimization->OptimizedTransformation(scanline);
              position.Transform(scanline_transformation);
            }
            RNLoadRgb(segment_colors[segment_index % max_segment_colors]);
            R3LoadPoint(position);
          }
        }
      }
    }
  }
  glEnd();
}



static void 
DrawPlanarGrids(int side, RNBoolean transform, int color_scheme)
{
  // Draw planar grids
  glDisable(GL_LIGHTING);
  glBegin(GL_TRIANGLES);
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        ScanData *scan_data = (ScanData *) scan->Data();
        if (!scan_data) continue;
        if (scan->NScanlines() == 0) continue;
        if (ia == 1) continue;
        for (int g = 0; g < scan_data->planar_grids.NEntries(); g++) {
          R3PlanarGrid *planar_grid = scan_data->planar_grids.Kth(g);
          R3Vector normal = planar_grid->Plane().Normal();
          for (int i = 1; i < planar_grid->XResolution(); i++) {
            for (int j = 1; j < planar_grid->YResolution(); j++) {
              RNScalar d00 = planar_grid->GridValue(i-1, j-1);
              RNScalar d11 = planar_grid->GridValue(i, j);
              RNScalar d01 = planar_grid->GridValue(i-1, j);
              RNScalar d10 = planar_grid->GridValue(i, j-1);
              if ((d00 != R2_GRID_UNKNOWN_VALUE) && (d11 != R2_GRID_UNKNOWN_VALUE)) {
                R3Point p00 = planar_grid->WorldPosition(i-1, j-1) + d00 * normal;
                R3Point p11 = planar_grid->WorldPosition(i, j) + d11 * normal;
                RNScalar c00 = 0.5 * (1.0 + d00);
                RNScalar c11 = 0.5 * (1.0 + d11);
                if (d01 != R2_GRID_UNKNOWN_VALUE) {
                  R3Point p01 = planar_grid->WorldPosition(i-1, j) + d01 * normal;
                  RNScalar c01 = 0.5 * (1.0 + d01);
                  glColor3d(c00, c00, c00);
                  R3LoadPoint(p00);
                  glColor3d(c01, c01, c01);
                  R3LoadPoint(p01);
                  glColor3d(c11, c11, c11);
                  R3LoadPoint(p11);
                }
                if (d10 != R2_GRID_UNKNOWN_VALUE) {
                  R3Point p10 = planar_grid->WorldPosition(i, j-1) + d10 * normal;
                  RNScalar c10 = 0.5 * (1.0 + d10);
                  glColor3d(c00, c00, c00);
                  R3LoadPoint(p00);
                  glColor3d(c11, c11, c11);
                  R3LoadPoint(p11);
                  glColor3d(c10, c10, c10);
                  R3LoadPoint(p10);
                }
              }
              else if ((d01 != R2_GRID_UNKNOWN_VALUE) && (d10 != R2_GRID_UNKNOWN_VALUE)) {
                R3Point p01 = planar_grid->WorldPosition(i-1, j) + d01 * normal;
                R3Point p10 = planar_grid->WorldPosition(i, j-1) + d10 * normal;
                RNScalar c01 = 0.5 * (1.0 + d01);
                RNScalar c10 = 0.5 * (1.0 + d10);
                if (d00 != R2_GRID_UNKNOWN_VALUE) {
                  R3Point p00 = planar_grid->WorldPosition(i-1, j-1) + d00 * normal;
                  RNScalar c00 = 0.5 * (1.0 + d00);
                  glColor3d(c00, c00, c00);
                  R3LoadPoint(p00);
                  glColor3d(c01, c01, c01);
                  R3LoadPoint(p01);
                  glColor3d(c10, c10, c10);
                  R3LoadPoint(p10);
                }
                if (d11 != R2_GRID_UNKNOWN_VALUE) {
                  R3Point p11 = planar_grid->WorldPosition(i, j) + d11 * normal;
                  RNScalar c11 = 0.5 * (1.0 + d11);
                  glColor3d(c01, c01, c01);
                  R3LoadPoint(p01);
                  glColor3d(c11, c11, c11);
                  R3LoadPoint(p11);
                  glColor3d(c10, c10, c10);
                  R3LoadPoint(p10);
                }
              }
            }
          }
        }
      }
    }
  }
  glEnd();
}



static void 
DrawImages(int side, RNBoolean transform, int color_scheme)
{
  // Check selected image
  if (!selected_image[side]) return;

  // Draw selected image
  DrawImage(selected_image[side], side, transform);
}



static void 
DrawBBoxes(int side, RNBoolean transform, int color_scheme)
{
  // Draw bbox for every segment
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      LoadSegmentColor(segment, color_scheme);
      segment->BBox().Outline();
    }
  }
}



////////////////////////////////////////////////////////////////////////
// Pick functions
////////////////////////////////////////////////////////////////////////

static int 
Pick(GSVPoseOptimization *optimization, int x, int y, 
  RNBoolean pick_features = TRUE, RNBoolean pick_cameras = TRUE, RNBoolean pick_lasers = TRUE, RNBoolean pick_scans = TRUE,
  GSVFeature **picked_feature = NULL, GSVImage **picked_image = NULL, GSVScanline **picked_scanline = NULL, int *picked_point_index = NULL, 
  R3Point *picked_position = NULL)
{
  // How close the cursor has to be to a point (in pixels)
  int pick_tolerance = 10;
  int pick_clusters = FALSE; 

  // Clear window 
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Set viewing transformation
  int side = split_screen * ((x < GLUTwindow_width/2) ? 0 : 1);
  int viewport_width = GLUTwindow_width/(split_screen+1);
  glViewport(side*viewport_width, 0, viewport_width, GLUTwindow_height);
  viewer[side]->Camera().Load();

  // Set OpenGL stuff
  glLineWidth(2 * pick_tolerance);
  glPointSize(2 * pick_tolerance);    

  // Draw everything
  if (pick_clusters && show_clusters[side]) DrawClusters(side, show_optimization[side], PICK_COLOR_SCHEME);
  if (pick_features && show_features[side]) DrawFeatures(side, show_optimization[side], PICK_COLOR_SCHEME);
  if (pick_features && show_pairs[side]) DrawPairs(side, show_optimization[side], PICK_COLOR_SCHEME);
  if (pick_features && show_correspondences[side]) DrawCorrespondences(side, show_optimization[side], PICK_COLOR_SCHEME);
  if (pick_cameras && show_cameras[side]) DrawCameras(side, show_optimization[side], PICK_COLOR_SCHEME);
  if (pick_lasers && show_lasers[side]) DrawLasers(side, show_optimization[side], PICK_COLOR_SCHEME);
  if (pick_scans && show_scans[side]) DrawScans(side, show_optimization[side], PICK_COLOR_SCHEME);

  // Reset OpenGL stuff
  glPointSize(1);
  glLineWidth(1);
  glFinish();

  // Read color buffer at cursor position
  unsigned char rgba[4];
  glReadPixels(x, y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, rgba);
  if (rgba[3] == 0) return 0;

  // Return scanline
  int r = rgba[0] & 0xFF;
  int g = rgba[1] & 0xFF;
  int b = rgba[2] & 0xFF;
  int a = rgba[3] & 0xFF;

  // Determine pick type
  int pick_type = 0;
  GSVScanline *scanline = NULL;
  if (a == 0xFA) {
    // Picked cluster
    int cluster_index = (r << 16) | (g << 8) | b;
    if (cluster_index < 0) return 0;
    if (cluster_index >= optimization->NClusters()) return 0;
    // GSVFeatureCluster *cluster = optimization->Cluster(cluster_index);
    if (picked_feature) *picked_feature = NULL;
    if (picked_image) *picked_image = NULL;
    if (picked_scanline) *picked_scanline = NULL;
    if (picked_point_index) *picked_point_index = -1;
    pick_type = 5;
  }
  else if (a == 0xFB) {
    // Picked feature
    int feature_index = (r << 16) | (g << 8) | b;
    if (feature_index < 0) return 0;
    if (feature_index >= optimization->NFeatures()) return 0;
    GSVFeature *feature = optimization->Feature(feature_index);
    if (picked_feature) *picked_feature = feature;
    if (picked_image) *picked_image = NULL;
    if (picked_scanline) *picked_scanline = feature->scanline;
    if (picked_point_index) *picked_point_index = feature->scan_point_index;
    scanline = feature->scanline;
    pick_type = 1;
  }
  else if (a == 0xFC) {
    // Picked camera
    int image_index = (r << 16) | (g << 8) | b;
    if (image_index < 0) return 0;
    if (image_index > optimization->images.NEntries()) return 0;
    GSVImage *image = optimization->images.Kth(image_index);
    if (picked_feature) *picked_feature = NULL;
    if (picked_image) *picked_image = image;
    if (picked_scanline) *picked_scanline = NULL;
    if (picked_point_index) *picked_point_index = -1;
    pick_type = 2;
  }
  else if (a == 0xFD) {
    // Picked laser
    int scanline_index = (r << 16) | (g << 8) | b;
    if (scanline_index < 0) return 0;
    if (scanline_index > optimization->scanlines.NEntries()) return 0;
    scanline = optimization->scanlines.Kth(scanline_index);
    if (picked_feature) *picked_feature = NULL;
    if (picked_image) *picked_image = NULL;
    if (picked_scanline) *picked_scanline = scanline;
    if (picked_point_index) *picked_point_index = -1;
    pick_type = 3;
  }
  else {
    // Picked scan point
    pick_type = 1;
    int scanline_index = (r << 16) | (g << 8) | b;
    if (scanline_index < 0) return 0;
    if (scanline_index > optimization->scanlines.NEntries()) return 0;
    scanline = optimization->scanlines.Kth(scanline_index);
    if (picked_feature) *picked_feature = NULL;
    if (picked_image) *picked_image = NULL;
    if (picked_scanline) *picked_scanline = scanline;
    if (picked_point_index) *picked_point_index = a - 1;
    pick_type = 4;
  }

  // Return position
  if (picked_position) {
    // Find hit position
    GLfloat depth;
    GLdouble p[3];
    GLint viewport[4];
    GLdouble modelview_matrix[16];
    GLdouble projection_matrix[16];
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix);
    glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix);
    glReadPixels(x, y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
    gluUnProject(x, y, depth, modelview_matrix, projection_matrix, viewport, &(p[0]), &(p[1]), &(p[2]));
    R3Point position(p[0], p[1], p[2]);
    *picked_position = position;
  }
  
  // Return pick type
  return pick_type;
}



////////////////////////////////////////////////////////////////////////
// GLUT callback functions
////////////////////////////////////////////////////////////////////////

void AtExit(void)
{
  // Write correspondences
  if (output_correspondences_name) {
    WriteCorrespondences(optimization, output_correspondences_name);
  }

  // Write transformations
  if (output_transformations_name) {
    if (!WriteTransformations(optimization, output_transformations_name)) exit(-1);
  }

  // Write inertias
  if (output_inertias_name) {
    if (!WriteInertias(optimization, output_inertias_name)) exit(-1);
  }

  // Write options
  if (output_options_name) {
    if (!WriteOptions(optimization, output_options_name)) exit(-1);
  }

  // Write world distance error plot
  if (output_world_distance_error_plot_name) {
    WriteErrorPlot(optimization, output_world_distance_error_plot_name, TRUE, TRUE, TRUE, TRUE);
  }

  // Write image distance error plot
  if (output_image_distance_error_plot_name) {
    WriteErrorPlot(optimization, output_image_distance_error_plot_name, TRUE, TRUE, TRUE, FALSE);
  }

  // Write scene
  if (output_scene_name) {
    WriteScene(scene, optimization, output_scene_name);
  }
}



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

  // Draw both sides of screen
  for (int side = 0; side <= split_screen; side++) {
    // Update viewer to match selected camera
    if (viewpoint_scheme[side] && selected_image[side]) {
      GSVCamera *camera = selected_image[side]->Camera();
      if (camera) {
        RNAngle xfov = 0.5 * camera->XFov();
        RNAngle yfov = 0.5 * camera->YFov();
        GSVPose pose = selected_image[side]->Pose();
        R3Point viewpoint = pose.Viewpoint();
        R3Vector towards = pose.Towards();
        R3Vector up = pose.Up();
        if (show_optimization[side]) {
          GSVImage *image = selected_image[side];
          R2Point center(image->Width()/2, image->Height()/2);
          R3Affine transformation = optimization->OptimizedTransformation(image, center);
          viewpoint.Transform(transformation);
          towards.Transform(transformation);
          up.Transform(transformation);
        }
        center_point[side] = viewpoint + towards;
        R3Camera c(viewpoint, towards, up, xfov, yfov, 0.001, 1000000);
        viewer[side]->SetCamera(c);
      }
    }

    // Set viewport
    if (!split_screen) glViewport(0, 0, GLUTwindow_width, GLUTwindow_height);    
    else glViewport(side*GLUTwindow_width/2, 0, GLUTwindow_width/2, GLUTwindow_height);
    
    // Set viewing transformation
    viewer[side]->Camera().Load();
    glPointSize(point_size[side]);

    // Draw scans
    if (show_scans[side]) {
      DrawScans(side, show_optimization[side], color_scheme[side]);
    }

    // Draw clusters
    if (show_clusters[side]) {
      DrawClusters(side, show_optimization[side], color_scheme[side]);
    }

    // Draw correspondences
    if (show_correspondences[side]) {
      DrawCorrespondences(side, show_optimization[side], color_scheme[side]);
    }

   // Draw pairs
    if (show_pairs[side]) {
      DrawPairs(side, show_optimization[side], color_scheme[side]);
    }

    // Draw features
    if (show_features[side]) {
      DrawFeatures(side, show_optimization[side], color_scheme[side]);
    }

    // Draw selected feature
    if (selected_feature[side]) {
      RNLoadRgb(1.0, 1.0, 1.0);
      DrawFeature(selected_feature[side], side, show_optimization[side]);
    }
        
    // Draw meshes
    if (show_meshes[side]) {
      DrawMeshes(side, show_optimization[side], color_scheme[side]);
    }
    
    // Draw paths
    if (show_paths[side]) {
      DrawPaths(side, show_optimization[side], color_scheme[side]);
    }
    
    // Draw cameras
    if (show_cameras[side]) {
      DrawCameras(side, show_optimization[side], color_scheme[side]);
    }
    
    // Draw lasers
    if (show_lasers[side]) {
      DrawLasers(side, show_optimization[side], color_scheme[side]);
  }
    
    // Draw images
    if (show_images[side]) {
      DrawImages(side, show_optimization[side], color_scheme[side]);
    }
    
    // Draw segments
    if (show_segments[side]) {
      DrawSegments(side, show_optimization[side], color_scheme[side]);
    }
    
    // Draw planar grids
    if (show_planar_grids[side]) {
      DrawPlanarGrids(side, show_optimization[side], color_scheme[side]);
    }
    
    // Draw bounding box
    if (show_bboxes[side]) {
      DrawBBoxes(side, show_optimization[side], color_scheme[side]);
    }
  }

  // Swap buffers 
  glutSwapBuffers();
}    



void GLUTResize(int w, int h)
{
  // Resize window
  glViewport(0, 0, w, h);

  // Resize viewer viewport
  if (split_screen) {
    viewer[0]->ResizeViewport(0, 0, w/2, h);
    viewer[1]->ResizeViewport(w/2, 0, w/2, h);
  }
  else {
    viewer[0]->ResizeViewport(0, 0, w, h);
    viewer[1]->ResizeViewport(0, 0, w, h);
  }

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
  if (GLUTbutton[0]) viewer[GLUTside]->RotateWorld(1.0, center_point[GLUTside], x, y, dx, dy);
  else if (GLUTbutton[1]) viewer[GLUTside]->ScaleWorld(1.0, center_point[GLUTside], x, y, dx, dy);
  else if (GLUTbutton[2]) viewer[GLUTside]->TranslateWorld(1.0, center_point[GLUTside], x, y, dx, dy);
  if (GLUTbutton[0] || GLUTbutton[1] || GLUTbutton[2]) glutPostRedisplay();

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;
}



void GLUTMouse(int button, int state, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;
  
  // Compute side
  int side = split_screen * ((x < GLUTwindow_width/2) ? 0 : 1);

  // Check for drag
  static int last_mouse_down[2] = { -999, -999 };
  RNBoolean drag = (fabs(x - last_mouse_down[0]) > 10) | (fabs(y - last_mouse_down[1]) > 10);
  last_mouse_down[0] = (state == GLUT_DOWN) ? x : -999;
  last_mouse_down[1] = (state == GLUT_DOWN) ? y : -999;

  // Process mouse button event
  if (!drag && (state == GLUT_UP)) {
    if (button == GLUT_LEFT_BUTTON) {
      // Check for double click
      static RNBoolean double_click = FALSE;
      static RNTime last_mouse_down_time;
      double_click = !double_click && (last_mouse_down_time.Elapsed() < 1);
      last_mouse_down_time.Read();

      if (double_click) {
        // Toggle display of selected segment
        R3Point position;
        GSVImage *image = NULL;
        GSVScanline *scanline = NULL;
        if (Pick(optimization, x, y, FALSE, show_cameras[side], show_lasers[side], show_scans[side], NULL, &image, &scanline, NULL, &position)) {
          center_point[side] = position;
          if (scanline) {
            RNScalar timestamp = scanline->Timestamp();
            GSVScan *scan = scanline->Scan();
            GSVSegment *segment = scan->Segment();
            ToggleScanData(segment, timestamp);
            glutPostRedisplay();
            printf("Scanline %20s %02d %02d %06d %12.6f\n", 
              segment->Run()->Name(), segment->RunIndex(), scan->SegmentIndex(), scanline->ScanIndex(), timestamp);
          }
          else if (image) {
            RNScalar timestamp = image->Timestamp();
            GSVPanorama *panorama = image->Panorama();
            GSVSegment *segment = panorama->Segment();
            ToggleScanData(segment, timestamp);
            glutPostRedisplay();
            printf("Image %20s %02d %06d %02d %12.6f\n",  
              segment->Run()->Name(), segment->RunIndex(), panorama->SegmentIndex(), image->PanoramaIndex(), timestamp);
          }
        }
      }
      else {
        // Set center point
        R3Point position;
        if (Pick(optimization, x, y, show_features[side], show_cameras[side], show_lasers[side], show_scans[side], NULL, NULL, NULL, NULL, &position)) {
          center_point[side] = position;
          glutPostRedisplay();
        }
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

  // Remember side 
  GLUTside = side;
}



void GLUTKeyboard(unsigned char key, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;

  // Compute side
  int side = split_screen * ((x < GLUTwindow_width/2) ? 0 : 1);

  // Process keyboard button event 
  switch (key) {
  case 'A':
  case 'a': 
  case 'N':
  case 'n': 
    if (1) {
      // Add scan feature
      R3Point position;
      int point_index = -1;
      GSVScanline *scanline = NULL;
      if (Pick(optimization, x, y, FALSE, FALSE, FALSE, show_scans[side], NULL, NULL, &scanline, &point_index, &position)) {
        if (point_index >= 0) {
          center_point[side] = position;
          scanline->Scan()->ReadPoints();
          position = scanline->PointPosition(point_index);
          R3Vector normal = ComputeNormal(scanline, point_index);
          scanline->Scan()->ReleasePoints();
          if (show_optimization[side]) position.InverseTransform(optimization->OptimizedTransformation(scanline));
          int feature_type = ((key == 'a') || (key == 'A')) ? GSV_SCAN_POINT_FEATURE_TYPE : GSV_SCAN_PLANE_FEATURE_TYPE;
          GSVFeature *feature = new GSVFeature(feature_type, scanline, point_index, position, R3zero_vector, normal, 1.0, 1.0);
          optimization->InsertFeature(feature);
          last_feature = feature;
          last_correspondence = NULL;
          if ((key == 'A') || (key == 'N')) {
            if (split_screen && selected_feature[1-side]) {
              GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(selected_feature[1-side], feature, 1);
              optimization->InsertCorrespondence(correspondence);
              last_correspondence = correspondence;
            }
            else if (selected_feature[side]) {
              GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(selected_feature[side], feature, 1);
              optimization->InsertCorrespondence(correspondence);
              last_correspondence = correspondence;
            }
          }
          selected_feature[side] = feature;
        }
      }
    }
    break; 

  case 'z': 
    if (1) {
      // Add scan lline feature
      R3Point position;
      int point_index = -1;
      GSVScanline *scanline = NULL;
      if (Pick(optimization, x, y, FALSE, FALSE, FALSE, show_scans[side], NULL, NULL, &scanline, &point_index, &position)) {
        if (point_index >= 0) {
          center_point[side] = position;
          if (show_optimization[side]) position.InverseTransform(optimization->OptimizedTransformation(scanline));
          int feature_type = GSV_SCAN_LINE_FEATURE_TYPE;
          GSVFeature *feature = new GSVFeature(feature_type, scanline, point_index, position, R3posz_vector, R3zero_vector, 1.0, 1.0);
          optimization->InsertFeature(feature);
          last_feature = feature;
          last_correspondence = NULL;
          selected_feature[side] = feature;
        }
      }
    }
    break; 

  case 'X':
  case 'x': 
    if (selected_image[side]) {
      // Add image feature
      int sx = (side == 0) ? x : x - (GLUTwindow_width/2);
      double ix = (double) (sx * selected_image[side]->Width()) / (GLUTwindow_width/2);
      double iy = (double) (y * selected_image[side]->Height()) / GLUTwindow_height;
      GSVFeature *feature = new GSVFeature(GSV_IMAGE_POINT_FEATURE_TYPE, 
        selected_image[side], R2Point(ix, iy), R2zero_vector, 1.0, 1.0);
      optimization->InsertFeature(feature);
      last_feature = feature;
      last_correspondence = NULL;
      if (key == 'X') {
        if (split_screen && selected_feature[1-side]) {
          GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(selected_feature[1-side], feature, 1);
          optimization->InsertCorrespondence(correspondence);
          last_correspondence = correspondence;
        }
        else if (selected_feature[side]) {
          GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(selected_feature[side], feature, 1);
          optimization->InsertCorrespondence(correspondence);
          last_correspondence = correspondence;
        }
      }
      selected_feature[side] = feature;
    }
    break; 

  case 'B':
  case 'b':
    show_bboxes[side] = !show_bboxes[side];
    break;

  case 'C':
  case 'c':
    show_cameras[side] = !show_cameras[side];
    break;

  case 'D':
  case 'd':
    if (color_scheme[side] == SEGMENT_COLOR_SCHEME) color_scheme[side] = VIEWPOINT_DEPTH_COLOR_SCHEME;
    else if (color_scheme[side] == VIEWPOINT_DEPTH_COLOR_SCHEME) color_scheme[side] = HEIGHT_COLOR_SCHEME;
    else color_scheme[side] = SEGMENT_COLOR_SCHEME;
    break;

  case 'F':
  case 'f':
    show_features[side] = !show_features[side];
    break;

  case 'G':
  case 'g':
    show_clusters[side] = !show_clusters[side];
    break;

  case 'I':
  case 'i':
    show_images[side] = !show_images[side];
    break;

  case 'L':
  case 'l':
    show_lasers[side] = !show_lasers[side];
    break;

  case 'M':
  case 'm':
    show_correspondences[side] = !show_correspondences[side];
    break;

  case 'O':
  case 'o':
    show_optimization[side] = !show_optimization[side];
    if (show_optimization[side]) printf("Optimized\n");
    else printf("Not optimized\n");
    break;

  case 'P':
  case 'p':
    show_pairs[side] = !show_pairs[side];
    break;

  case 'r':
    rigidity = 0.75 * rigidity;
    printf("Rigidity = %g\n", rigidity);
    break;

  case 'R':
    rigidity = 1.5 * rigidity;
    printf("Rigidity = %g\n", rigidity);
    break;

  case 'S':
  case 's':
    show_scans[side] = !show_scans[side];
    break;

  case 't':
    optimization->max_pair_world_distance *= 0.75;
    optimization->max_correspondence_world_distance *= 0.75;
    printf("Max distance = %g %g\n", 
      optimization->max_pair_world_distance, 
      optimization->max_correspondence_world_distance);
    break;

  case 'T':
    optimization->max_pair_world_distance *= 1.5;
    optimization->max_correspondence_world_distance *= 1.5;
    printf("Max distance = %g %g\n", 
      optimization->max_pair_world_distance, 
      optimization->max_correspondence_world_distance);
    break;

  case 'V':
  case 'v':
    viewpoint_scheme[side] = !viewpoint_scheme[side];
    break;

  case 'W':
  case 'w':
    split_screen = 1 - split_screen;
    GLUTResize(GLUTwindow_width, GLUTwindow_height);
    break;

  case '!':
    for (int i = 0; i < optimization->LASER_NV; i++) optimization->laser_inertia_weights[i] *= 0.1;
    for (int i = 0; i < optimization->CAMERA_NV; i++) optimization->camera_inertia_weights[i] *= 0.1;
    printf("%g\n", optimization->laser_inertia_weights[0]);
    break;

  case '@':
    for (int i = 0; i < optimization->LASER_NV; i++) optimization->laser_inertia_weights[i] *= 10;
    for (int i = 0; i < optimization->CAMERA_NV; i++) optimization->camera_inertia_weights[i] *= 10;
    printf("%g\n", optimization->laser_inertia_weights[0]);
    break;

  case '#':
    show_paths[side] = !show_paths[side];
    break;

  case '$': {
    R3Point position;
    GSVImage *image = NULL;
    GSVScanline *scanline = NULL;
    if (Pick(optimization, x, y, FALSE, show_cameras[side], show_lasers[side], show_scans[side], NULL, &image, &scanline, NULL, &position)) {
      if (scanline) {
        RNScalar timestamp = scanline->Timestamp();
        GSVScan *scan = scanline->Scan();
        GSVSegment *segment = scan->Segment();
        GSVRun *run = segment->Run();
        printf("%30s %02d %12g %12.6f %12.6f %12.6f\n", run->Name(), segment->RunIndex(), timestamp, position.X(), position.Y(), position.Z());
      }
      else if (image) {
        RNScalar timestamp = image->Timestamp();
        GSVPanorama *panorama = image->Panorama();
        GSVSegment *segment = panorama->Segment();
        GSVRun *run = segment->Run();
        printf("%30s %02d %12g %12.6f %12.6f %12.6f\n", run->Name(), segment->RunIndex(), timestamp, position.X(), position.Y(), position.Z());
      }
    }
    break; }

  case '%':
    solver = (solver + 1) % RN_NUM_SOLVERS;
    printf("Solver = %d\n", solver);
    break;

  case '^':
    show_segments[side] = !show_segments[side];
    break;

  case '&':
    show_meshes[side] = !show_meshes[side];
    break;

  case '*':
    show_planar_grids[side] = !show_planar_grids[side];
    break;

  case '(':
    optimization->expression_type = 0;
    break;

  case ')':
    optimization->expression_type = 1;
    break;

  case '+':
   point_size[side] += 1;
   break;

  case '_':
    point_size[side] -= 1;
    break;

  case '=':
    if (image_plane_distance < 1000) image_plane_distance *= 1.1;
    break;

  case '-':
    if (image_plane_distance > 0.01) image_plane_distance *= 0.9;
    break;

  case '/':
    if (optimization->max_correspondence_shape_context_descriptor_distance == 0)
      optimization->max_correspondence_shape_context_descriptor_distance = 0.01;
    else optimization->max_correspondence_shape_context_descriptor_distance *= 2;
    printf("SD = %g\n", optimization->max_correspondence_shape_context_descriptor_distance);
    break;

  case '\\':
    if (optimization->max_correspondence_shape_context_descriptor_distance <= 0.01)
      optimization->max_correspondence_shape_context_descriptor_distance = 0;
    else optimization->max_correspondence_shape_context_descriptor_distance /= 2;
    printf("SD = %g\n", optimization->max_correspondence_shape_context_descriptor_distance);
    break;

  case 1:  // ctrl-A
    // Update pose translations and rotations
    Solve(optimization, FALSE, FALSE,  FALSE, FALSE,  TRUE, TRUE);
    break; 

  case 2:  // ctrl-B
    // Update all translations (including camera and laser offsets)
    Solve(optimization, TRUE, FALSE,  TRUE, FALSE,  TRUE, FALSE);
    break;

  case 3:  // ctrl-C
    // Update camera translations only
    Solve(optimization, FALSE, FALSE,  TRUE, FALSE,  FALSE, FALSE);
    break;

  case 4:  // ctrl-D
    // Update pixel depths only
    Solve(optimization, FALSE, FALSE,  FALSE, FALSE,  FALSE, FALSE);
    break;

  case 5: // ctrl-E
    // Clear pose transformations
    // optimization->ClearOptimizationVariables();
    optimization->UpdatePixelPositionsFromRayIntersections();
    optimization->UpdateClusterPositionsFromFeaturePositions();
    Solve(optimization, FALSE, FALSE,  FALSE, FALSE,  FALSE, FALSE,  FALSE, FALSE,  FALSE, FALSE, FALSE);
    break;

  case 6:  // ctrl-F
    // Update path translations and rotations
    Solve(optimization, FALSE, FALSE,  FALSE, FALSE,  TRUE, TRUE);
    break; 

  case 7:  // ctrl-G
    // Empty all correspondences
    if (TRUE) {
      while (optimization->NCorrespondences() > 0) {
        GSVFeatureCorrespondence *correspondence = optimization->Correspondence(optimization->NCorrespondences()-1);
        optimization->RemoveCorrespondence(correspondence);
      }
    }
    break; 

  case 12:  // ctrl-L
    // Update laser only
    Solve(optimization, TRUE, FALSE,  TRUE, FALSE,  FALSE, TRUE);
    break; 

  case 18: // ctrl-R
    // Update pose rotations
    Solve(optimization, FALSE, FALSE,  FALSE, FALSE,  FALSE, TRUE);
    break; 

  case 19:  // ctrl-S
    if (output_correspondences_name) {
      WriteCorrespondences(optimization, output_correspondences_name);
    }
    break;

  case 20: // ctrl-T
    // Update pose translations
    Solve(optimization, FALSE, FALSE,  FALSE, FALSE,  TRUE, FALSE);
    break; 

  case 26:  // ctrl-Z
    if (last_correspondence) {
      // Remove last correspondence
      printf("HEREA\n");
      optimization->RemoveCorrespondence(last_correspondence);
      delete last_correspondence;
      last_correspondence = NULL;
    }
    if (last_feature) {
      // Remove last feature
      printf("HEREB\n");
      if (last_feature == selected_feature[0]) selected_feature[0] = NULL;
      if (last_feature == selected_feature[1]) selected_feature[1] = NULL;
      optimization->RemoveFeature(last_feature);
      delete last_feature;
      last_feature = NULL;
    }
    break;

  case 17: // ctrl-Q
    // Quit
    GLUTStop();
    break;

  case ' ':
    if (1) {
      // Select segment and unselect all others
      R3Point position;
      GSVImage *image = NULL;
      GSVScanline *scanline = NULL;
      if (Pick(optimization, x, y, FALSE, show_cameras[side], show_lasers[side], show_scans[side], NULL, &image, &scanline, NULL, &position)) {
        center_point[side] = position;
        GSVSegment *segment = NULL;
        RNScalar timestamp = -1;
        if (scanline) {
          selected_scanline[side] = scanline;
          segment = scanline->Scan()->Segment();
          timestamp = scanline->Timestamp();
        }
        else if (image) {
          selected_image[side] = image;
          segment = image->Tapestry()->Segment();
          timestamp = image->Timestamp();
        }
        if (segment) {
          if (glutGetModifiers() & GLUT_ACTIVE_SHIFT) ToggleScanData(segment, timestamp);
          else SelectScanData(segment, timestamp);
        }
      }
    }
    break;

  case 27: // ESC
    // Select no segments
    selected_feature[0] = NULL;
    selected_feature[1] = NULL;
    selected_scanline[0] = NULL;
    selected_scanline[1] = NULL;
    selected_image[0] = NULL;
    selected_image[1] = NULL;
    SelectScanData(NULL, -1);
    break;
  }

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

  // Remember side 
  GLUTside = side;

  // Redraw
  glutPostRedisplay();  
}



void GLUTSpecial(int key, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;

  // Compute side
  int side = split_screen * ((x < GLUTwindow_width/2) ? 0 : 1);

  // Process keyboard button event 
  switch (key) {
  case GLUT_KEY_F1:
    // Select feature
    if (show_features[side] || show_correspondences[side] || show_clusters[side]) {
      R3Point position;
      selected_feature[side] = NULL;
      if (Pick(optimization, x, y, TRUE, FALSE, FALSE, FALSE, &selected_feature[side], NULL, NULL, NULL, &position)) {
        if (selected_feature[side]) {
          center_point[side] = position;
          if ((selected_feature[side]->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) ||
              (selected_feature[side]->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) ||
              (selected_feature[side]->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) {
            GSVScanline *scanline = selected_feature[side]->scanline;
            GSVScan *scan = (scanline) ? scanline->Scan() : NULL;
            GSVSegment *segment = (scan) ? scan->Segment() : NULL;
            GSVRun *run = (segment) ? segment->Run() : NULL;
            int scanline_index = (scanline) ? scanline->ScanIndex() : -1;
            int scan_index = (scan) ? scan->SegmentIndex() : -1;
            int segment_index = (segment) ? segment->RunIndex() : -1;
            const char *run_name = (run) ? run->Name() : "None";
            R3Point scan_position = optimization->WorldPosition(selected_feature[side], FALSE); 
            R3Vector scan_direction = optimization->WorldNormal(selected_feature[side], FALSE); 
            RNScalar timestamp = scanline->Timestamp();
            GSVDescriptor *descriptor = &selected_feature[side]->descriptor;
            printf("%2d   %6d    %30s %2d %2d %6d %3d     %8.6f     %12.6f %12.6f %12.6f    %8.6f %8.6f %8.6f    %8.6f \n", 
                   selected_feature[side]->feature_type, selected_feature[side]->index, 
                   run_name, segment_index, scan_index, scanline_index, 
                   selected_feature[side]->scan_point_index, selected_feature[side]->score,
                   scan_position.X(), scan_position.Y(), scan_position.Z(), 
                   scan_direction.X(), scan_direction.Y(), scan_direction.Z(), timestamp);
            for (int i = 0; i < descriptor->nvalues; i++) 
              printf("%6.3f ", descriptor->values[i]);
            printf("\n");
          }
          else if ((selected_feature[side]->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) ||
                   (selected_feature[side]->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) {
            GSVImage *image = selected_feature[side]->image;
            GSVPanorama *panorama = (image) ? image->Panorama() : NULL;
            GSVSegment *segment = (panorama) ? panorama->Segment() : NULL;
            GSVRun *run = (segment) ? segment->Run() : NULL;
            int image_index = (image) ? image->PanoramaIndex() : -1;
            int panorama_index = (panorama) ? panorama->SegmentIndex() : -1;
            int segment_index = (segment) ? segment->RunIndex() : -1;
            const char *run_name = (run) ? run->Name() : "None";
            RNScalar timestamp = image->Timestamp();
            GSVDescriptor *descriptor = &selected_feature[side]->descriptor;
            printf("%2d   %6d    %30s %2d %6d %2d    %8.6f    %12.6f %12.6f    %8.6f %8.6f    %8.6f    %8.6f \n", 
                   selected_feature[side]->feature_type, selected_feature[side]->index, 
                   run_name, segment_index,  panorama_index, image_index, selected_feature[side]->score,
                   selected_feature[side]->image_position.X(), selected_feature[side]->image_position.Y(), 
                   selected_feature[side]->image_direction.X(), selected_feature[side]->image_direction.Y(),
                   selected_feature[side]->image_t, timestamp);
            for (int i = 0; i < descriptor->nvalues; i++) 
              printf("%6.3f ", descriptor->values[i]);
            printf("\n");
          }
        }
      }
    }
    break; 

  case GLUT_KEY_F2:
    // Create correspondence to selected feature 
    if (show_features[side] || show_correspondences[side] || show_clusters[side]) {
      R3Point position;
      GSVFeature *feature = NULL;
      if (split_screen && selected_feature[1-side]) {
        if (Pick(optimization, x, y, TRUE, FALSE, FALSE, FALSE, &feature, NULL, NULL, NULL, &position)) {
          if (feature != selected_feature[1-side]) {
            GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(selected_feature[1-side], feature, 1);
            optimization->InsertCorrespondence(correspondence);
            center_point[side] = position;
            last_correspondence = correspondence;
            last_feature = NULL;
          }
        }
      }
      else if (selected_feature[side]) {
        if (Pick(optimization, x, y, TRUE, FALSE, FALSE, FALSE, &feature, NULL, NULL, NULL, &position)) {
          if (feature != selected_feature[1-side]) {
            GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(selected_feature[side], feature, 1);
            optimization->InsertCorrespondence(correspondence);
            center_point[side] = position;
            last_correspondence = correspondence;
            last_feature = NULL;
          }
        }
      }
    }
    break; 

  case GLUT_KEY_F3:
    // Delete all correspondences with selected feature
    if (show_features[side] || show_correspondences[side] || show_clusters[side]) {
      R3Point position;
      GSVFeature *feature = NULL;
      if (Pick(optimization, x, y, TRUE, FALSE, FALSE, FALSE, &feature, NULL, NULL, NULL, &position)) {
        if (feature) {
          center_point[side] = position;
          RNArray<GSVFeatureCorrespondence *> correspondences;
          for (int i = 0; i < optimization->NCorrespondences(); i++)
            correspondences.Insert(optimization->Correspondence(i));
          for (int i = 0; i < correspondences.NEntries(); i++) {
            GSVFeatureCorrespondence *correspondence = correspondences.Kth(i);
            if ((correspondence->Feature(0) == feature) || (correspondence->Feature(1) == feature)) {
              optimization->RemoveCorrespondence(correspondence);
              delete correspondence;
            }
          }
        }
      }
    }
    break; 

  case GLUT_KEY_F4: {
    optimization->CreateCloseClusters();
    // int saved_ncorrespondences = optimization->NCorrespondences();
    // optimization->CreatePixelPixelPlaneCorrespondences();
    // printf("Created %d correspondences\n", optimization->NCorrespondences() - saved_ncorrespondences);
    break; }

  case GLUT_KEY_F5: {
    int saved_ncorrespondences = optimization->NCorrespondences();
    optimization->CreateTransitiveClosureCorrespondences();
    printf("Created %d correspondences\n", optimization->NCorrespondences() - saved_ncorrespondences);
    break; }

  case GLUT_KEY_F6:
    if (show_features[side] || show_correspondences[side] || show_clusters[side]) {
      R3Point position;
      selected_feature[side] = NULL;
      if (Pick(optimization, x, y, TRUE, FALSE, FALSE, FALSE, &selected_feature[side], NULL, NULL, NULL, &position)) {
        if (selected_feature[side]) {
          for (int i = 0; i < optimization->NCorrespondences(); i++) {
            GSVFeatureCorrespondence *correspondence = optimization->Correspondence(i);
            GSVFeature *feature0 = correspondence->Feature(0);
            GSVFeature *feature1 = correspondence->Feature(1);
            if ((feature0 != selected_feature[side]) && (feature1 != selected_feature[side])) continue;
            GSVDescriptor& descriptor0 = feature0->descriptor;
            GSVDescriptor& descriptor1 = feature1->descriptor;
            printf("%d %d : %g %g %g : %g\n", feature0->index, feature1->index, 
                   sqrt(optimization->SquaredWorldDistance(feature0, feature1)), 
                   sqrt(optimization->SquaredImageDistance(feature0, feature1)),
                   sqrt(descriptor0.SquaredDistance(descriptor1)),
                   optimization->CorrespondenceScore(feature0, feature1));
          }
        }
      }
    }
    break;

   case GLUT_KEY_F7: {
     optimization->TruncateClusters(0);
     optimization->CreateAllClusters();
     printf("Created %d clusters\n", optimization->NClusters());
     break; }

   case GLUT_KEY_F9: {
    static int correspondence_index = 0;
    if (correspondence_index < optimization->NCorrespondences()) {
      split_screen = 1;
      GSVFeatureCorrespondence *correspondence = optimization->Correspondence(correspondence_index);
      GSVFeature *feature0 = correspondence->Feature(0);
      GSVFeature *feature1 = correspondence->Feature(1);
      selected_feature[0] = feature0;
      selected_feature[1] = feature1;
      selected_image[0] = feature0->image;
      selected_image[1] = feature1->image;
      viewpoint_scheme[0] = 1;
      viewpoint_scheme[1] = 1;
      GLUTResize(GLUTwindow_width, GLUTwindow_height);
      correspondence_index = (correspondence_index < optimization->NCorrespondences()-1) ? correspondence_index + 1 : 0;
    }
    break; }

  case GLUT_KEY_F10:
    optimization->TruncateCorrespondences(num_initial_correspondences);
    optimization->CreateCloseCorrespondences();
    printf("%g %d %d\n", optimization->max_correspondence_world_distance, 
       optimization->NPairs(), optimization->NCorrespondences());
    break;

  case GLUT_KEY_F11: {
    optimization->TruncateCorrespondences(num_initial_correspondences);
    optimization->CreateCloseCorrespondences();
    optimization->Solve(FALSE, FALSE,  FALSE, FALSE,  TRUE, TRUE);
    printf("%g %d %d\n", optimization->max_correspondence_world_distance, 
       optimization->NPairs(), optimization->NCorrespondences());
    break; }

  case GLUT_KEY_F12: {
    static RNScalar max_distance = 20;
    rigidity = 0.1 * max_distance;
    optimization->max_pair_world_distance = max_distance;
    optimization->max_correspondence_world_distance = max_distance;
    optimization->TruncateCorrespondences(num_initial_correspondences);
    optimization->CreateCloseCorrespondences();
    optimization->Solve(FALSE, FALSE,  FALSE, FALSE,  TRUE, FALSE);
    if (max_distance < 2) optimization->Solve(FALSE, FALSE,  FALSE, FALSE,  TRUE, TRUE);
    printf("%g %g %d\n", optimization->max_correspondence_world_distance, rigidity, optimization->NCorrespondences());
    max_distance *= 0.5;
    break; }

  case GLUT_KEY_DOWN:
  case GLUT_KEY_UP:
    if (selected_image[side]) {
      GSVTapestry *tapestry = selected_image[side]->Tapestry();
      int index = selected_image[side]->TapestryIndex();
      if ((key == GLUT_KEY_DOWN) && (index > 0)) 
        selected_image[side] = tapestry->Image(index-1);
      else if ((key == GLUT_KEY_UP) && (index < tapestry->NImages()-1)) 
        selected_image[side] = tapestry->Image(index+1);
    }
    break;

  case GLUT_KEY_RIGHT:
  case GLUT_KEY_LEFT:
    if (selected_image[side]) {
      GSVPanorama *panorama = selected_image[side]->Panorama();
      int index = selected_image[side]->PanoramaIndex();
      if (key == GLUT_KEY_LEFT) index = (index > 0) ? index-1 : 7;
      if (key == GLUT_KEY_RIGHT) index = (index < 7) ? index+1 : 0;
      selected_image[side] = panorama->Image(index);
    }
    break;
  }

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

  // Remember side 
  GLUTside = side;

  // Redraw
  glutPostRedisplay();
}



////////////////////////////////////////////////////////////////////////
// Top-Level functions
////////////////////////////////////////////////////////////////////////

static void 
InitInterface(int *argc, char **argv)
{
  // Open window 
  glutInit(argc, argv);
  glutInitWindowPosition(100, 100);
  glutInitWindowSize(GLUTwindow_width, GLUTwindow_height);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_ALPHA); 
  GLUTwindow = glutCreateWindow("GSV Street Viewer");

  // Initialize background color
  glClearColor(background_color[0], background_color[1], background_color[2], 0);
  // glClearColor(1, 1, 1, 0);

  // Initialize lights
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  static GLfloat lmodel_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
  glEnable(GL_NORMALIZE);

  // Initialize headlight
  static GLfloat light0_diffuse[] = { 0.5, 0.5, 0.5, 1.0 };
  static GLfloat light0_specular[] = { 0.1, 0.1, 0.1, 1.0 };
  static GLfloat light0_position[] = { 0.0, 0.0, 1.0, 0.0 };
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
  glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
  glEnable(GL_LIGHT0);

  // Initialize side light 
  static GLfloat light1_diffuse[] = { 0.5, 0.5, 0.5, 1.0 };
  static GLfloat light1_specular[] = { 0.5, 0.5, 0.5, 1.0 };
  static GLfloat light1_position[] = { 0.5, -0.5, 0.866, 0 };
  glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
  glLightfv(GL_LIGHT1, GL_SPECULAR, light1_specular);
  glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
  glEnable(GL_LIGHT1);

  // Initialize materials 
  static GLfloat ambient_material[] = { 0.1, 0.1, 0.1, 1.0 };
  static GLfloat diffuse_material[] = { 0.5, 0.5, 0.5, 1.0 };
  static GLfloat specular_material[] = { 1, 1, 1, 1 };
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient_material);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse_material);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular_material);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 10);

  // Initialize graphics modes  
  glEnable(GL_DEPTH_TEST);

  // Initialize GLUT callback functions 
  glutDisplayFunc(GLUTRedraw);
  glutReshapeFunc(GLUTResize);
  glutKeyboardFunc(GLUTKeyboard);
  glutSpecialFunc(GLUTSpecial);
  glutMouseFunc(GLUTMouse);
  glutMotionFunc(GLUTMotion);
  atexit(AtExit);
}



static void 
RunInterface(void)
{
  // Initialize center points
  center_point[0] = scene->BBox().Centroid();
  center_point[1] = scene->BBox().Centroid();

  // Initialize viewers
  for (int side = 0; side < 2; side++) {
    R3Box bbox = scene->BBox();
    assert(!bbox.IsEmpty());
    RNLength r = bbox.DiagonalRadius();
    assert((r > 0.0) && RNIsFinite(r));
    R3Point origin = center_point[side] + R3posz_vector * (2.5 * r);
    RNAngle yfov = RN_PI_OVER_FOUR;
    RNScalar window_aspect = (double) GLUTwindow_height / (double) GLUTwindow_width;
    RNAngle xfov = atan(tan(yfov) / window_aspect);
    R3Camera camera(origin, R3Vector(0, 0, -1), R3Vector(0, 1, 0), xfov, yfov, 0.001 * r, 1000.0 * r);
    R2Viewport viewport(0, 0, GLUTwindow_width, GLUTwindow_height);
    viewer[side] = new R3Viewer(camera, viewport);
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
      else if (!strcmp(*argv, "-interactive")) { interactive = 1; }
      else if (!strcmp(*argv, "-input_scene")) { argc--; argv++; input_scene_name = *argv; }
      else if (!strcmp(*argv, "-output_scene")) { argc--; argv++; output_scene_name = *argv; }
      else if (!strcmp(*argv, "-input_features")) { argc--; argv++; input_correspondences_name = *argv; }
      else if (!strcmp(*argv, "-output_features")) { argc--; argv++; output_correspondences_name = *argv; }
      else if (!strcmp(*argv, "-input_pairs")) { argc--; argv++; input_correspondences_name = *argv; }
      else if (!strcmp(*argv, "-output_pairs")) { argc--; argv++; output_correspondences_name = *argv; }
      else if (!strcmp(*argv, "-input_correspondences")) { argc--; argv++; input_correspondences_name = *argv; }
      else if (!strcmp(*argv, "-output_correspondences")) { argc--; argv++; output_correspondences_name = *argv; }
      else if (!strcmp(*argv, "-input_transformations")) { argc--; argv++; input_transformations_name = *argv; }
      else if (!strcmp(*argv, "-output_transformations")) { argc--; argv++; output_transformations_name = *argv; }
      else if (!strcmp(*argv, "-input_inertias")) { argc--; argv++; input_inertias_name = *argv; }
      else if (!strcmp(*argv, "-output_inertias")) { argc--; argv++; output_inertias_name = *argv; }
      else if (!strcmp(*argv, "-input_options")) { argc--; argv++; input_options_name = *argv; }
      else if (!strcmp(*argv, "-output_options")) { argc--; argv++; output_options_name = *argv; }
      else if (!strcmp(*argv, "-output_world_distance_error_plot")) { argc--; argv++; output_world_distance_error_plot_name = *argv; }
      else if (!strcmp(*argv, "-output_image_distance_error_plot")) { argc--; argv++; output_image_distance_error_plot_name = *argv; }
      else if (!strcmp(*argv, "-create_scan_pole_features")) { create_features = create_scan_pole_features = 1; }
      else if (!strcmp(*argv, "-create_scan_curb_features")) { create_features = create_scan_curb_features = 1; }
      else if (!strcmp(*argv, "-create_scan_ridge_and_valley_features")) { create_features = create_scan_ridge_and_valley_features = 1; }      
      else if (!strcmp(*argv, "-create_scan_edge_features")) { create_features = create_scan_edge_features = 1; }      
      else if (!strcmp(*argv, "-create_scan_plane_features")) { create_features = create_scan_plane_features = 1; }
      else if (!strcmp(*argv, "-create_image_corner_features")) { create_features = create_image_corner_features = 1; }
      else if (!strcmp(*argv, "-create_image_sift_features")) { create_features = create_image_sift_features = 1; }
      else if (!strcmp(*argv, "-create_image_line_features")) { create_features = create_image_line_features = 1; }
      else if (!strcmp(*argv, "-create_scan_features")) { create_features = create_scan_pole_features = create_scan_plane_features = 1; }
      else if (!strcmp(*argv, "-create_image_features")) { create_features = create_image_corner_features = 1; }
      else if (!strcmp(*argv, "-create_features")) { create_features = create_scan_pole_features = create_scan_plane_features = create_image_corner_features = 1; }
      else if (!strcmp(*argv, "-create_pairs")) { create_pairs = create_close_pairs = 1; }
      else if (!strcmp(*argv, "-create_all_pairs")) { create_pairs = create_all_pairs = 1; }
      else if (!strcmp(*argv, "-create_feasible_pairs")) { create_pairs = create_feasible_pairs = 1; }
      else if (!strcmp(*argv, "-create_close_pairs")) { create_pairs = create_close_pairs = 1; }
      else if (!strcmp(*argv, "-create_coplanar_pairs")) { create_pairs = create_coplanar_pairs = 1; }
      else if (!strcmp(*argv, "-create_pixel_pixel_plane_pairs")) { create_pairs = create_pixel_pixel_plane_pairs = 1; }
      else if (!strcmp(*argv, "-create_correspondences")) { create_correspondences = create_icp_correspondences = 1; }
      else if (!strcmp(*argv, "-create_all_correspondences")) { create_correspondences = create_all_correspondences = 1; }
      else if (!strcmp(*argv, "-create_feasible_correspondences")) { create_correspondences = create_feasible_correspondences = 1; }
      else if (!strcmp(*argv, "-create_icp_correspondences")) { create_correspondences = create_icp_correspondences = 1; }
      else if (!strcmp(*argv, "-create_close_correspondences")) { create_correspondences = create_close_correspondences = 1; }
      else if (!strcmp(*argv, "-create_coplanar_correspondences")) { create_correspondences = create_coplanar_correspondences = 1; }
      else if (!strcmp(*argv, "-create_pixel_pixel_plane_correspondences")) { create_correspondences = create_pixel_pixel_plane_correspondences = 1; }
      else if (!strcmp(*argv, "-create_transitive_closure_correspondences")) { create_correspondences = create_transitive_closure_correspondences = 1; }
      else if (!strcmp(*argv, "-create_all_clusters")) { create_clusters = create_all_clusters = 1; }
      else if (!strcmp(*argv, "-create_close_clusters")) { create_clusters = create_close_clusters = 1; }
      else if (!strcmp(*argv, "-create_merged_clusters")) { create_clusters = create_merged_clusters = 1; }
      else if (!strcmp(*argv, "-create_hierarchical_clusters")) { create_clusters = create_hierarchical_clusters = 1; }
      else if (!strcmp(*argv, "-create_coplanar_clusters")) { create_clusters = create_coplanar_clusters = 1; }
      else if (!strcmp(*argv, "-create_icp_clusters")) { create_clusters = create_icp_clusters = 1; }
      else if (!strcmp(*argv, "-create_transformations")) { update_transformations = update_all_transformations = 1; }
      else if (!strcmp(*argv, "-update_transformations")) { update_transformations = update_all_transformations = 1; }
      else if (!strcmp(*argv, "-update_path_transformations")) { update_transformations = update_path_translations = update_path_rotations = update_pixel_depths = 1; }
      else if (!strcmp(*argv, "-update_path_translations")) { update_transformations = update_path_translations = update_pixel_depths = 1; }
      else if (!strcmp(*argv, "-update_path_rotations")) { update_transformations = update_path_rotations = update_pixel_depths = 1; }
      else if (!strcmp(*argv, "-update_cluster_positions")) { update_transformations = update_cluster_positions = include_cluster_equations = 1; }
      else if (!strcmp(*argv, "-update_cluster_rotations")) { update_transformations = update_cluster_rotations = include_cluster_equations = 1; }
      else if (!strcmp(*argv, "-update_pixel_depths")) { update_transformations = update_pixel_depths = 1; }
      else if (!strcmp(*argv, "-update_scan_features")) { update_features = update_scan_features = 1; }
      else if (!strcmp(*argv, "-update_pixel_features")) { update_features = update_pixel_features = 1; }
      else if (!strcmp(*argv, "-remove_pairs")) { remove_pairs = 1; }
      else if (!strcmp(*argv, "-remove_all_pairs")) { remove_pairs = remove_all_pairs = 1; }
      else if (!strcmp(*argv, "-remove_infeasible_pairs")) { remove_pairs = remove_infeasible_pairs = 1; }
      else if (!strcmp(*argv, "-remove_correspondences")) { remove_correspondences = remove_all_correspondences = 1; }
      else if (!strcmp(*argv, "-remove_all_correspondences")) { remove_correspondences = remove_all_correspondences = 1; }
      else if (!strcmp(*argv, "-remove_all_clusters")) { remove_clusters = remove_all_clusters = 1; }
      else if (!strcmp(*argv, "-remove_infeasible_correspondences")) { remove_correspondences = remove_infeasible_correspondences = 1; }
      else if (!strcmp(*argv, "-include_cluster_equations")) { update_cluster_positions = include_cluster_equations = 1; }
      else if (!strcmp(*argv, "-exclude_cluster_equations")) { update_cluster_positions = include_cluster_equations = 0; }
      else if (!strcmp(*argv, "-include_correspondence_equations")) { include_correspondence_equations = 1; }
      else if (!strcmp(*argv, "-exclude_correspondence_equations")) { include_correspondence_equations = 0; }
      else if (!strcmp(*argv, "-load_meshes")) { load_meshes = 1; }
      else if (!strcmp(*argv, "-load_planar_grids")) { load_planar_grids = 1; }
      else if (!strcmp(*argv, "-load_all_points")) { load_all_points = 1; }
      else if (!strcmp(*argv, "-load_segments")) { load_segments = 1; }
      else if (!strcmp(*argv, "-segment_grid_name")) { argc--; argv++; segment_grid_name = *argv; }
      else if (!strcmp(*argv, "-timestamp_radius")) { argc--; argv++; timestamp_radius = atof(*argv); }
      else if (!strcmp(*argv, "-path_vertex_spacing")) { argc--; argv++; path_vertex_spacing = atof(*argv); }
      else if (!strcmp(*argv, "-rigidity")) { argc--; argv++; rigidity = atof(*argv); }
      else if (!strcmp(*argv, "-solver")) { argc--; argv++; solver = atoi(*argv); }
      else if (!strcmp(*argv, "-background")) { 
        argc--; argv++; background_color[0] = atof(*argv); 
        argc--; argv++; background_color[1] = atof(*argv); 
        argc--; argv++; background_color[2] = atof(*argv); 
      }
      else { 
        fprintf(stderr, "Invalid program argument: %s\n", *argv); 
        exit(1); 
      }
      argv++; argc--;
    }
    else {
      if (!input_scene_name) input_scene_name = *argv;
      else if (!output_scene_name) output_scene_name = *argv;
      else if (!input_correspondences_name) input_correspondences_name = *argv;
      else if (!output_correspondences_name) output_correspondences_name = *argv;
      else { fprintf(stderr, "Invalid program argument: %s\n", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check scene name
  if (!input_scene_name) {
    fprintf(stderr, "Usage: gsvalign input_scenefile [output_scenefile] [-interactive] [options]\n");
    return FALSE;
  }

  // Update point loading parameters
  if (!strstr(input_scene_name, ".gsv")) {
    load_all_points = TRUE;
  }

  // Make interactive if no output
  if (!output_scene_name && !output_correspondences_name && !output_transformations_name && !output_inertias_name && 
      !output_world_distance_error_plot_name && !output_image_distance_error_plot_name) {
    interactive = TRUE;
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
  scene = ReadScene(input_scene_name);
  if (!scene) exit(-1);

  // Allocate optimization data structures
  optimization = new GSVPoseOptimization(scene, path_vertex_spacing);
  if (!optimization) exit(-1);

  // Read options
  if (input_options_name) {
    if (!ReadOptions(optimization, input_options_name)) exit(-1);
  }

  // Set parameters
  if (include_correspondence_equations < 0) include_correspondence_equations = TRUE;
  if (include_cluster_equations < 0) include_cluster_equations = TRUE;
  if (update_cluster_positions) include_cluster_equations = TRUE;

  // Read correspondences
  if (input_correspondences_name) {
    if (!ReadCorrespondences(optimization, input_correspondences_name)) exit(-1);
  }

  // Read transformations
  if (input_transformations_name) {
    if (!ReadTransformations(optimization, input_transformations_name)) exit(-1);
  }

  // Read inertias
  UpdateInertias(optimization);
  if (input_inertias_name) {
    if (!ReadInertias(optimization, input_inertias_name)) exit(-1);
  }

  // Remove pairs
  if (remove_pairs) {
    if (!RemovePairs(optimization)) exit(-1);
  }

  // Remove correspondences
  if (remove_correspondences) {
    if (!RemoveCorrespondences(optimization)) exit(-1);
  }

  // Remove clusters
  if (remove_clusters) {
    if (!RemoveClusters(optimization)) exit(-1);
  }

  // Create features
  if (create_features) {
    if (!CreateFeatures(optimization)) exit(-1);
  }

  // Update feature positions based on scan point positions
  if (update_features) {
    if (!UpdateFeatures(optimization)) exit(-1);
  }

  // Create pairs
  if (create_pairs) {
    if (!CreatePairs(optimization)) exit(-1);
  }

  // Create correspondences
  if (create_correspondences) {
    if (!CreateCorrespondences(optimization)) exit(-1);
  }

  // Create clusters
  if (create_clusters) {
    if (!CreateClusters(optimization)) exit(-1);
  }

  // Update transformations
  if (update_transformations) {
    if (!UpdateTransformations(optimization)) exit(-1);
  }

  // Check if interactive
  if (interactive) {
    // Initialize GLUT interface
    InitInterface(&argc, argv);
    
    // Run GLUT interface
    RunInterface();
  }
  else {
    // Write correspondences
    if (output_correspondences_name) {
      if (!WriteCorrespondences(optimization, output_correspondences_name)) exit(-1);
    }

    // Write transformations
    if (output_transformations_name) {
      if (!WriteTransformations(optimization, output_transformations_name)) exit(-1);
    }

    // Write inertias
    if (output_inertias_name) {
      if (!WriteInertias(optimization, output_inertias_name)) exit(-1);
    }

    // Write options
    if (output_options_name) {
      if (!WriteOptions(optimization, output_options_name)) exit(-1);
    }

    // Write world distance error plot
    if (output_world_distance_error_plot_name) {
      if (!WriteErrorPlot(optimization, output_world_distance_error_plot_name, TRUE, TRUE, TRUE, TRUE)) exit(-1);
    }

    // Write image distance error plot
    if (output_image_distance_error_plot_name) {
      if (!WriteErrorPlot(optimization, output_image_distance_error_plot_name, TRUE, TRUE, TRUE, FALSE)) exit(-1);
    }

    // Write scene
    if (output_scene_name) {
      if (!WriteScene(scene, optimization, output_scene_name)) exit(-1);
    }
  }

  // Return success 
  return 0;
}


