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

static const char *input_scene_name = NULL;
static const char *output_image_directory = NULL;
static int capture_XY_images = 0;
static int capture_SA_images = 0;
static int capture_DA_images = 0;
static int capture_DH_images = 0;
static int capture_UV_images = 0;
static int capture_opengl_images = 0;
static int include_base_images = 0;
static int include_normal_images = 0;
static int include_curvature_images = 0;
static int include_boundary_images = 0;
static int include_pca_images = 0;
static int include_color_images = 0;
static int include_depth_images = 0;
static int include_hull_images = 0;
static int include_label_images = 0;
static int print_verbose = 0;
static int print_debug = 0;


// XY, DA, DH parameters

static int scan_ground_index = 25;
static RNLength SA_viewpoint_spacing = 0.01;
static RNLength XY_image_spacing = 0.1;
static RNLength DA_image_spacing = 0.1;
static RNLength DH_image_spacing = 0.1;
static RNLength DH_max_height = 30;
static RNLength DH_max_hole_size = 1;

// Depth and color parameters

static double min_frustum_depth = 1;
static double max_frustum_depth = 1000;
static double max_viewpoint_distance = 100;
static double max_timestamp_difference = 100;



////////////////////////////////////////////////////////////////////////
// Grid names
////////////////////////////////////////////////////////////////////////

// Grid names
const char *base_image_names [] = {
  "Scanline", "PointIndex", "Timestamp",
  "ViewpointAngle", "ViewpointDistance", "ViewpointDepth",
  "Height", "GroundZ", 
  "PositionX", "PositionY", "PositionZ",
  "ViewpointX", "ViewpointY", "ViewpointZ",
  "TowardsX", "TowardsY", "TowardsZ",
  "UpX", "UpY", "UpZ"
};

const char *xy_image_names [] = {
  "Density", "ViewpointDistance", 
  "MinZ", "MaxZ", "GroundZ"
};

int xy_image_combination_methods [] = {
  2, 1, 1, 2, 0
};

const char *normal_image_names [] = {
  "NormalX", "NormalY", "NormalZ", "NdotV"
};

const char *curvature_image_names [] = {
  "HorizontalCurvature", "VerticalCurvature"
};

const char *boundary_image_names [] = {
  "BoundaryType", 
};

const char *pca_image_names [] = {
  "Density",
  "NearestNeighbor1", "NearestNeighbor2", "NearestNeighbor4", "NearestNeighbor8", "NearestNeighbor16",
  "PrincipleAxis1X", "PrincipleAxis1Y", "PrincipleAxis1Z",
  "PrincipleAxis2X", "PrincipleAxis2Y", "PrincipleAxis2Z",
  "PrincipleAxis3X", "PrincipleAxis3Y", "PrincipleAxis3Z",
  "Lambda1", "Lambda2", "Lambda3",
  "Lambda21", "Lambda31", "Lambda32"
};

const char *color_image_names [] = {
  "ColorRed", "ColorGreen", "ColorBlue"
};

int num_base_image_names = sizeof(base_image_names) / sizeof(const char *);
int num_xy_image_names = sizeof(xy_image_names) / sizeof(const char *);
int num_normal_image_names = sizeof(normal_image_names) / sizeof(const char *);
int num_curvature_image_names = sizeof(curvature_image_names) / sizeof(const char *);
int num_boundary_image_names = sizeof(boundary_image_names) / sizeof(const char *);
int num_pca_image_names = sizeof(pca_image_names) / sizeof(const char *);
int num_color_image_names = sizeof(color_image_names) / sizeof(const char *);



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

  // Read scene 
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



////////////////////////////////////////////////////////////////////////
// Utility types and functions
////////////////////////////////////////////////////////////////////////

static void
InterpolateMissingColumns(R2Grid& grid)
{
  // Initialize previous column with values
  int ix0 = -1;

  // Consider every column
  for (int ix1 = 0; ix1 < grid.XResolution(); ix1++) {
    // Check if column has values
    RNBoolean found_value = FALSE;
    for (int iy = 0; iy < grid.YResolution(); iy++) {
      if (grid.GridValue(ix1, iy) != R2_GRID_UNKNOWN_VALUE) { 
        found_value = TRUE; 
        break; 
      }
    }

    // Skip column, if has no values
    if (!found_value) continue;

    // Iterpolate values in skipped columns
    if ((ix0 >= 0) && (ix0 < ix1-1)) {
      for (int iy = 0; iy < grid.YResolution(); iy++) {
        RNScalar value0 = grid.GridValue(ix0, iy);
        if (value0 == R2_GRID_UNKNOWN_VALUE) continue;
        RNScalar value1 = grid.GridValue(ix1, iy);
        if (value1 == R2_GRID_UNKNOWN_VALUE) continue;
        for (int ix = ix0+1; ix < ix1; ix++) {
          RNScalar t = (double) (ix - ix0) / (double) (ix1 - ix0);
          RNScalar value = (1-t)*value0 + t*value1;
          grid.SetGridValue(ix, iy, value);
        }
      }
    }

    // Remember last column with values
    ix0 = ix1;
  }
}



static int
CombineImages(R2Grid& result, const RNArray<R2Grid *>& inputs, int combination_method = 0)
{
  // Check number of inputs
  if (inputs.IsEmpty()) return 0;

  // Initialize result
  if (combination_method == 0) result.Clear(0);
  else if (combination_method == 1) result.Clear(FLT_MAX);
  else if (combination_method == 2) result.Clear(-FLT_MAX);
  else result.Clear(0); 

  // Initialize count
  R2Grid count(result.XResolution(), result.YResolution());
    
  // Process grids
  for (int i = 0; i < inputs.NEntries(); i++) {
    R2Grid *grid = inputs.Kth(i);
    for (int iy = 0; iy < result.YResolution(); iy++) {
      for (int ix = 0; ix < result.XResolution(); ix++) {
        R2Point world_position = result.WorldPosition(ix, iy);
        RNScalar value = grid->WorldValue(world_position);
        if (value == R2_GRID_UNKNOWN_VALUE) continue;
        if (value == 0) continue; // Not right, should check if outside bounding box
        count.AddGridValue(ix, iy, 1.0);
        if (combination_method == 0) {
          result.AddGridValue(ix, iy, value);
        }
        else if (combination_method == 1) {
          RNScalar old_value = result.GridValue(ix, iy);
          if ((old_value != R2_GRID_UNKNOWN_VALUE) && (value > old_value)) continue;
          result.SetGridValue(ix, iy, value);
        }
        else if (combination_method == 2) {
          RNScalar old_value = result.GridValue(ix, iy);
          if ((old_value != R2_GRID_UNKNOWN_VALUE) && (value < old_value)) continue;
          result.SetGridValue(ix, iy, value);
        }
      }
    }
  }

  // Substitute unknowns
  if (combination_method == 0) { result.Substitute(0, R2_GRID_UNKNOWN_VALUE); result.Divide(count); }
  else if (combination_method == 1) result.Substitute(FLT_MAX, R2_GRID_UNKNOWN_VALUE);
  else if (combination_method == 2) result.Substitute(-FLT_MAX, R2_GRID_UNKNOWN_VALUE);

  // Return success
  return 1;
}    
    


////////////////////////////////////////////////////////////////////////
// Parameterization independent grid functions
////////////////////////////////////////////////////////////////////////

static int
WriteNormalImages(GSVScan *scan, const char *output_image_directory, const char *image_parameterization)
{
  // Get convenient variables
  if (!scan) return 0;
  if (scan->NScanlines() == 0) return 1;
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;

  // Skip if already done
  char image_name[4096];
  sprintf(image_name, "%s/%s/%02d_%02d_%s_NormalX.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (RNFileExists(image_name)) return 1;
  
  // Print message
  if (print_debug) {
    printf("    Creating normal images ...\n");
    fflush(stdout);
  }

  // Read position grids
  R2Grid positionX_grid, positionY_grid, positionZ_grid;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionX.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!positionX_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionY.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!positionY_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionZ.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!positionZ_grid.Read(image_name)) return 0;

  // Read viewpoint grids
  R2Grid viewpointX_grid, viewpointY_grid, viewpointZ_grid;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_ViewpointX.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!viewpointX_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_ViewpointY.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!viewpointY_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_ViewpointZ.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!viewpointZ_grid.Read(image_name)) return 0;

  // Create normal grids
  R2Grid normalX_grid(positionX_grid);   normalX_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid normalY_grid(normalX_grid);
  R2Grid normalZ_grid(normalX_grid);
  R2Grid NdotV_grid(normalX_grid);

  // Fill normal grids
  int xres = normalX_grid.XResolution();
  int yres = normalX_grid.YResolution();
  for (int i = 0; i < xres; i++) {
    int i0 = (i > 0) ? i-1 : 0;
    int i1 = (i < xres-1) ? i+1 : xres-1;
    for (int j = 0; j < yres; j++) {
      int j0 = (j > 0) ? j-1 : 0;
      int j1 = (j < yres-1) ? j+1 : yres-1;

      // Check position values
      if (positionX_grid.GridValue(i, j) == R2_GRID_UNKNOWN_VALUE) continue;
      if (positionY_grid.GridValue(i, j) == R2_GRID_UNKNOWN_VALUE) continue;
      if (positionZ_grid.GridValue(i, j) == R2_GRID_UNKNOWN_VALUE) continue;

      // Get 3D vector between horizontal neighbors
      RNScalar xA0 = positionX_grid.GridValue(i0, j);
      if (xA0 == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar yA0 = positionY_grid.GridValue(i0, j);
      if (yA0 == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar zA0 = positionZ_grid.GridValue(i0, j);
      if (zA0 == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar xA1 = positionX_grid.GridValue(i1, j);
      if (xA1 == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar yA1 = positionY_grid.GridValue(i1, j);
      if (yA1 == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar zA1 = positionZ_grid.GridValue(i1, j);
      if (zA1 == R2_GRID_UNKNOWN_VALUE) continue;
      R3Vector vA(xA1 - xA0, yA1 - yA0, zA1 - zA0);
      vA.Normalize();

      // Get 3D vector between vertical neighbors
      RNScalar xB0 = positionX_grid.GridValue(i, j0);
      if (xB0 == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar yB0 = positionY_grid.GridValue(i, j0);
      if (yB0 == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar zB0 = positionZ_grid.GridValue(i, j0);
      if (zB0 == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar xB1 = positionX_grid.GridValue(i, j1);
      if (xB1 == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar yB1 = positionY_grid.GridValue(i, j1);
      if (yB1 == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar zB1 = positionZ_grid.GridValue(i, j1);
      if (zB1 == R2_GRID_UNKNOWN_VALUE) continue;
      R3Vector vB(xB1 - xB0, yB1 - yB0, zB1 - zB0);
      vB.Normalize();

      // Compute view vector
      RNScalar dx = positionX_grid.GridValue(i, j) - viewpointX_grid.GridValue(i, j);
      RNScalar dy = positionY_grid.GridValue(i, j) - viewpointY_grid.GridValue(i, j);
      RNScalar dz = positionZ_grid.GridValue(i, j) - viewpointZ_grid.GridValue(i, j);
      R3Vector view_vector(dx, dy, dz); 
      view_vector.Normalize();

      // Compute normal
      R3Vector normal = vA % vB;
      normal.Normalize();

      // Compute NdotV
      RNScalar NdotV = normal.Dot(-view_vector);
      if (NdotV < 0) { NdotV = -NdotV; normal.Flip(); }

      // Update normal and NdotV grids
      normalX_grid.SetGridValue(i, j, normal.X());
      normalY_grid.SetGridValue(i, j, normal.Y());
      normalZ_grid.SetGridValue(i, j, normal.Z());
      NdotV_grid.SetGridValue(i, j, NdotV);
    }
  }

  // Write grids
  sprintf(image_name, "%s/%s/%02d_%02d_%s_NormalX.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  normalX_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_NormalY.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  normalY_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_NormalZ.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  normalZ_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_NdotV.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  NdotV_grid.Write(image_name);

  // Return success
  return 1;
}



static int
WriteCurvatureImages(GSVScan *scan, const char *output_image_directory, const char *image_parameterization)
{
  // Get convenient variables
  if (!scan) return 0;
  if (scan->NScanlines() == 0) return 1;
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;

  // Skip if already done
  char image_name[4096];
  sprintf(image_name, "%s/%s/%02d_%02d_%s_HorizontalCurvature.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (RNFileExists(image_name)) return 1;
  
  // Print message
  if (print_debug) {
    printf("    Creating curvature images ...\n");
    fflush(stdout);
  }

  // Read position grids
  R2Grid positionX_grid, positionY_grid, positionZ_grid;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionX.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!positionX_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionY.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!positionY_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionZ.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!positionZ_grid.Read(image_name)) return 0;
  positionX_grid.Blur(1);
  positionY_grid.Blur(1);
  positionZ_grid.Blur(1);

  // Read normal grids
  R2Grid normalX_grid, normalY_grid, normalZ_grid;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_NormalX.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!normalX_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_NormalY.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!normalY_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_NormalZ.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!normalZ_grid.Read(image_name)) return 0;
  normalX_grid.Blur(1);
  normalY_grid.Blur(1);
  normalZ_grid.Blur(1);

  // Create curvature grids
  R2Grid horizontal_curvature_grid(positionX_grid);   
  horizontal_curvature_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid vertical_curvature_grid(horizontal_curvature_grid);

  // Compute whether left-right directions should be flipped 
  int dir = (scan->SegmentIndex() != 2) ? 1 : -1;

  // Fill horizontal curvature grid
  for (int i = 1; i < horizontal_curvature_grid.XResolution()-1; i++) {
    for (int j = 0; j < horizontal_curvature_grid.YResolution(); j++) {
      // Get position values
      RNScalar positionX = positionX_grid.GridValue(i, j);
      if (positionX == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar positionY = positionY_grid.GridValue(i, j);
      if (positionY == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar positionZ = positionZ_grid.GridValue(i, j);
      if (positionZ == R2_GRID_UNKNOWN_VALUE) continue;
      R3Point position(positionX, positionY, positionZ);

      // Get normal values
      RNScalar normalX = normalX_grid.GridValue(i, j);
      if (normalX == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar normalY = normalY_grid.GridValue(i, j);
      if (normalY == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar normalZ = normalZ_grid.GridValue(i, j);
      if (normalZ == R2_GRID_UNKNOWN_VALUE) continue;
      R3Vector normal(normalX, normalY, normalZ);

      // Get neighbor positions
      RNScalar positionAX = positionX_grid.GridValue(i-dir, j);
      if (positionAX == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar positionAY = positionY_grid.GridValue(i-dir, j);
      if (positionAY == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar positionAZ = positionZ_grid.GridValue(i-dir, j);
      if (positionAZ == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar positionBX = positionX_grid.GridValue(i+dir, j);
      if (positionBX == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar positionBY = positionY_grid.GridValue(i+dir, j);
      if (positionBY == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar positionBZ = positionZ_grid.GridValue(i+dir, j);
      if (positionBZ == R2_GRID_UNKNOWN_VALUE) continue;
      R3Point positionA(positionAX, positionAY, positionAZ);
      R3Point positionB(positionBX, positionBY, positionBZ);
      R3Vector vA = positionA - position;
      R3Vector vB = positionB - position;
      RNLength lenA = vA.Length();
      RNLength lenB = vB.Length();
      if (RNIsZero(lenA) || RNIsZero(lenB)) continue;
      vA /= lenA; vB /= lenB;
      if (lenA < 1) lenA = 1;
      if (lenB < 1) lenB = 1;
      RNScalar angle = R3InteriorAngle(vA, vB);
      R3Vector dir_vector = vB + vA;
      RNScalar dot = dir_vector.Dot(normal);
      RNScalar sign = (dot > 0) ? -1 : 1;
      RNScalar curvature = sign * (RN_PI - angle) / (lenA + lenB);
      horizontal_curvature_grid.SetGridValue(i, j, curvature);
    }
  }

  // Fill vertical curvature grid
  for (int i = 0; i < horizontal_curvature_grid.XResolution(); i++) {
    for (int j = 1; j < horizontal_curvature_grid.YResolution()-1; j++) {
      // Get position values
      RNScalar positionX = positionX_grid.GridValue(i, j);
      if (positionX == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar positionY = positionY_grid.GridValue(i, j);
      if (positionY == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar positionZ = positionZ_grid.GridValue(i, j);
      if (positionZ == R2_GRID_UNKNOWN_VALUE) continue;
      R3Point position(positionX, positionY, positionZ);

      // Get normal values
      RNScalar normalX = normalX_grid.GridValue(i, j);
      if (normalX == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar normalY = normalY_grid.GridValue(i, j);
      if (normalY == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar normalZ = normalZ_grid.GridValue(i, j);
      if (normalZ == R2_GRID_UNKNOWN_VALUE) continue;
      R3Vector normal(normalX, normalY, normalZ);

      // Get neighbor positions
      RNScalar positionAX = positionX_grid.GridValue(i, j-1);
      if (positionAX == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar positionAY = positionY_grid.GridValue(i, j-1);
      if (positionAY == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar positionAZ = positionZ_grid.GridValue(i, j-1);
      if (positionAZ == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar positionBX = positionX_grid.GridValue(i, j+1);
      if (positionBX == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar positionBY = positionY_grid.GridValue(i, j+1);
      if (positionBY == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar positionBZ = positionZ_grid.GridValue(i, j+1);
      if (positionBZ == R2_GRID_UNKNOWN_VALUE) continue;
      R3Point positionA(positionAX, positionAY, positionAZ);
      R3Point positionB(positionBX, positionBY, positionBZ);
      R3Vector vA = positionA - position;
      R3Vector vB = positionB - position;
      RNLength lenA = vA.Length();
      RNLength lenB = vB.Length();
      if (RNIsZero(lenA) || RNIsZero(lenB)) continue;
      vA /= lenA; vB /= lenB;
      if (lenA < 1) lenA = 1;
      if (lenB < 1) lenB = 1;
      RNScalar angle = R3InteriorAngle(vA, vB);
      R3Vector dir_vector = vB + vA;
      RNScalar dot = dir_vector.Dot(normal);
      RNScalar sign = (dot > 0) ? -1 : 1;
      RNScalar curvature = sign * (RN_PI - angle) / (lenA + lenB);
      vertical_curvature_grid.SetGridValue(i, j, curvature);
    }
  }

  // Write grids
  sprintf(image_name, "%s/%s/%02d_%02d_%s_HorizontalCurvature.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  horizontal_curvature_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_VerticalCurvature.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  vertical_curvature_grid.Write(image_name);

  // Return success
  return 1;
}



static int
WriteBoundaryImages(GSVScan *scan, const char *output_image_directory, const char *image_parameterization,
  RNLength max_depth_discontinuity = 1.0, RNScalar curvature_threshold = 0.25)
{
  // Get convenient variables
  if (!scan) return 0;
  if (scan->NScanlines() == 0) return 1;
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;

  // Skip if already done
  char image_name[4096];
  sprintf(image_name, "%s/%s/%02d_%02d_%s_BoundaryType.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (RNFileExists(image_name)) return 1;
  
  // Print message
  if (print_debug) {
    printf("    Creating boundary images ...\n");
    fflush(stdout);
  }

  // Compute whether left-right directions should be flipped 
  int dir = (scan->SegmentIndex() != 2) ? 1 : -1;

  // Read viewpoint depth grid
  R2Grid viewpoint_depth_grid;
  R2Grid horizontal_curvature_grid;
  R2Grid vertical_curvature_grid;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_ViewpointDepth.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!viewpoint_depth_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_HorizontalCurvature.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!horizontal_curvature_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_VerticalCurvature.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!vertical_curvature_grid.Read(image_name)) return 0;

  // Create boundary grid
  R2Grid boundary_type_grid(viewpoint_depth_grid.XResolution(), viewpoint_depth_grid.YResolution());   

  // Fill boundary grids
  for (int i = 1; i < boundary_type_grid.XResolution()-1; i++) {
    for (int j = 1; j < boundary_type_grid.YResolution()-1; j++) {
      RNScalar depth = viewpoint_depth_grid.GridValue(i, j);
      if (depth == R2_GRID_UNKNOWN_VALUE) continue;

      // Initialize boundary type 
      int boundary_type = 0;
      RNScalar neighbor_depth;

      // Check left
      neighbor_depth = viewpoint_depth_grid.GridValue(i-dir, j);
      if (neighbor_depth == R2_GRID_UNKNOWN_VALUE) boundary_type += GSV_LEFT_UNKNOWN_BOUNDARY | GSV_LEFT_SILHOUETTE_BOUNDARY;
      else if (neighbor_depth > depth + max_depth_discontinuity) boundary_type += GSV_LEFT_SILHOUETTE_BOUNDARY;
      else if (neighbor_depth < depth - max_depth_discontinuity) boundary_type += GSV_LEFT_SHADOW_BOUNDARY;

      // Check right
      neighbor_depth = viewpoint_depth_grid.GridValue(i+dir, j);
      if (neighbor_depth == R2_GRID_UNKNOWN_VALUE) boundary_type += GSV_RIGHT_UNKNOWN_BOUNDARY | GSV_RIGHT_SILHOUETTE_BOUNDARY;
      else if (neighbor_depth > depth + max_depth_discontinuity) boundary_type += GSV_RIGHT_SILHOUETTE_BOUNDARY;
      else if (neighbor_depth < depth - max_depth_discontinuity) boundary_type += GSV_RIGHT_SHADOW_BOUNDARY;

      // Check down
      neighbor_depth = viewpoint_depth_grid.GridValue(i, j-1);
      if (neighbor_depth == R2_GRID_UNKNOWN_VALUE) boundary_type += GSV_DOWN_UNKNOWN_BOUNDARY | GSV_DOWN_SILHOUETTE_BOUNDARY;
      else if (neighbor_depth > depth + max_depth_discontinuity) boundary_type += GSV_DOWN_SILHOUETTE_BOUNDARY;
      else if (neighbor_depth < depth - max_depth_discontinuity) boundary_type += GSV_DOWN_SHADOW_BOUNDARY;

      // Check up
      neighbor_depth = viewpoint_depth_grid.GridValue(i, j+1);
      if (neighbor_depth == R2_GRID_UNKNOWN_VALUE) boundary_type += GSV_UP_UNKNOWN_BOUNDARY | GSV_UP_SILHOUETTE_BOUNDARY;
      else if (neighbor_depth > depth + max_depth_discontinuity) boundary_type += GSV_UP_SILHOUETTE_BOUNDARY;
      else if (neighbor_depth < depth - max_depth_discontinuity) boundary_type += GSV_UP_SHADOW_BOUNDARY;

      // Check other types
      if (boundary_type == 0) {
        // Check horizontal curvature
        RNScalar hcurv = horizontal_curvature_grid.GridValue(i, j);
        if (hcurv != R2_GRID_UNKNOWN_VALUE) {
          if (hcurv > curvature_threshold) boundary_type += GSV_HORIZONTAL_RIDGE_BOUNDARY;
          else if (hcurv < -curvature_threshold) boundary_type += GSV_HORIZONTAL_VALLEY_BOUNDARY;
        }

        // Check vertical curvature
        RNScalar vcurv = vertical_curvature_grid.GridValue(i, j);
        if (vcurv != R2_GRID_UNKNOWN_VALUE) {
          if (vcurv > curvature_threshold) boundary_type += GSV_VERTICAL_RIDGE_BOUNDARY;
          else if (vcurv < -curvature_threshold) boundary_type += GSV_VERTICAL_VALLEY_BOUNDARY;
        }
      }

      // Set grid value
      boundary_type_grid.SetGridValue(i, j, boundary_type);
    }
  }

  // Write grids
  sprintf(image_name, "%s/%s/%02d_%02d_%s_BoundaryType.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  boundary_type_grid.Write(image_name);

  // Return success
  return 1;
}



static int
WritePCAImages(GSVScan *scan, const char *output_image_directory, const char *image_parameterization,
  int max_image_distance = 16, double max_world_distance = 2.0)
{
  // Get convenient variables
  if (!scan) return 0;
  if (scan->NScanlines() == 0) return 1;
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;
  double max_dd = max_world_distance * max_world_distance;

  // Skip if already done
  char image_name[4096];
  sprintf(image_name, "%s/%s/%02d_%02d_%s_Density.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (RNFileExists(image_name)) return 1;
  
  // Print message
  if (print_debug) {
    printf("    Creating pca images ...\n");
    fflush(stdout);
  }

  // Read position grids
  R2Grid positionX_grid, positionY_grid, positionZ_grid;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionX.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!positionX_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionY.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!positionY_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PositionZ.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!positionZ_grid.Read(image_name)) return 0;

  // Read viewpoint grids
  R2Grid viewpointX_grid, viewpointY_grid, viewpointZ_grid;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_ViewpointX.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!viewpointX_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_ViewpointY.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!viewpointY_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_%s_ViewpointZ.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!viewpointZ_grid.Read(image_name)) return 0;

  // Initialize output grids
  R2Grid density_grid(positionX_grid);       
  density_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid nearest_neighbor_1_grid(density_grid);       
  R2Grid nearest_neighbor_2_grid(density_grid);       
  R2Grid nearest_neighbor_4_grid(density_grid);       
  R2Grid nearest_neighbor_8_grid(density_grid);       
  R2Grid nearest_neighbor_16_grid(density_grid);       
  R2Grid principle_axis_1_x_grid(density_grid);
  R2Grid principle_axis_1_y_grid(density_grid);
  R2Grid principle_axis_1_z_grid(density_grid);
  R2Grid principle_axis_2_x_grid(density_grid);
  R2Grid principle_axis_2_y_grid(density_grid);
  R2Grid principle_axis_2_z_grid(density_grid);
  R2Grid principle_axis_3_x_grid(density_grid);
  R2Grid principle_axis_3_y_grid(density_grid);
  R2Grid principle_axis_3_z_grid(density_grid);
  R2Grid lambda_1_grid(density_grid);
  R2Grid lambda_2_grid(density_grid);
  R2Grid lambda_3_grid(density_grid);

  // Allocate temporary memory for points
  R3Point *points = new R3Point[(2*max_image_distance+1) * (2*max_image_distance+1)];

  // For every point
  for (int ci = 0; ci < positionX_grid.XResolution(); ci++) {
    for (int cj = 0; cj < positionX_grid.YResolution(); cj++) {
      // Get point position
      RNScalar cx = positionX_grid.GridValue(ci, cj);
      if (cx == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar cy = positionY_grid.GridValue(ci, cj);
      if (cy == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar cz = positionZ_grid.GridValue(ci, cj);
      if (cz == R2_GRID_UNKNOWN_VALUE) continue;
      R3Point cp(cx, cy, cz);

      // Find neighbor points
      int npoints = 0;
      for (int di = -max_image_distance; di <= max_image_distance; di++) {
        int i = ci + di;
        if ((i < 0) || (i >= positionX_grid.XResolution())) continue;
        for (int dj = -max_image_distance; dj <= max_image_distance; dj++) {
          int j = cj + dj;
          if ((j < 0) || (j >= positionX_grid.YResolution())) continue;
          RNScalar x = positionX_grid.GridValue(i, j);
          if (x == R2_GRID_UNKNOWN_VALUE) continue;
          RNScalar y = positionY_grid.GridValue(i, j);
          if (y == R2_GRID_UNKNOWN_VALUE) continue;
          RNScalar z = positionZ_grid.GridValue(i, j);
          if (z == R2_GRID_UNKNOWN_VALUE) continue;
          R3Point p(x, y, z);
          RNLength dd = R3SquaredDistance(p, cp);
          if (dd > max_dd) continue;
          points[npoints++] = p;
        }
      }

      // Check if found enough neighbor points
      if (npoints == 0) continue;

      // Sort the 16 closest points
      for (int k1 = 0; k1 < npoints; k1++) {
        RNLength d1 = R3SquaredDistance(points[k1], cp);
        int k2max = k1;
        if (k2max > 17) {
          k2max = 17;
          R3Point swap = points[k2max];
          points[k2max] = points[k1];
          points[k1] = swap;
        }
        for (int k2 = k2max; k2 > 0; k2--) {
          RNLength d2 = R3SquaredDistance(points[k2-1], cp);
          if (d1 > d2) break;
          R3Point swap = points[k2];
          points[k2] = points[k2-1];
          points[k2-1] = swap;
        }
      }
          
      // Compute density
      RNLength r3 = max_world_distance * max_world_distance * max_world_distance;
      RNVolume volume = (4.0/3.0) * RN_PI * r3;
      RNScalar density = (volume > 0) ? npoints / volume : 0;

      // Compute nearest neighbor distances
      RNScalar nearest_neighbor_1 =  (npoints > 1)  ? R3Distance(points[1],  cp) : max_world_distance;
      RNScalar nearest_neighbor_2 =  (npoints > 2)  ? R3Distance(points[2],  cp) : max_world_distance;
      RNScalar nearest_neighbor_4 =  (npoints > 4)  ? R3Distance(points[4],  cp) : max_world_distance;
      RNScalar nearest_neighbor_8 =  (npoints > 8)  ? R3Distance(points[8],  cp) : max_world_distance;
      RNScalar nearest_neighbor_16 = (npoints > 16) ? R3Distance(points[16], cp) : max_world_distance;

      // Compute principle axes
      RNScalar variances[3];
      R3Point centroid = R3Centroid(npoints, points);
      R3Triad triad = R3PrincipleAxes(centroid, npoints, points, NULL, variances);
      
      // Compute view vector
      RNScalar dx = positionX_grid.GridValue(ci, cj) - viewpointX_grid.GridValue(ci, cj);
      RNScalar dy = positionY_grid.GridValue(ci, cj) - viewpointY_grid.GridValue(ci, cj);
      RNScalar dz = positionZ_grid.GridValue(ci, cj) - viewpointZ_grid.GridValue(ci, cj);
      R3Vector view_vector(dx, dy, dz); 
      view_vector.Normalize();

      // Get axes oriented based on view vector
      R3Vector axes[3];
      axes[0] = (triad[0].Dot(view_vector) < 0) ? triad[0] : -triad[0];
      axes[2] = (triad[2].Dot(view_vector) < 0) ? triad[2] : -triad[2];
      axes[1] = triad[2] % triad[1];
      axes[1].Normalize();

      // Fill grids
      density_grid.SetGridValue(ci, cj, density);
      nearest_neighbor_1_grid.SetGridValue(ci, cj, nearest_neighbor_1);
      nearest_neighbor_2_grid.SetGridValue(ci, cj, nearest_neighbor_2);
      nearest_neighbor_4_grid.SetGridValue(ci, cj, nearest_neighbor_4);
      nearest_neighbor_8_grid.SetGridValue(ci, cj, nearest_neighbor_8);
      nearest_neighbor_16_grid.SetGridValue(ci, cj, nearest_neighbor_16);
      principle_axis_1_x_grid.SetGridValue(ci, cj, axes[0][0]);
      principle_axis_1_y_grid.SetGridValue(ci, cj, axes[0][1]);
      principle_axis_1_z_grid.SetGridValue(ci, cj, axes[0][2]);
      principle_axis_2_x_grid.SetGridValue(ci, cj, axes[1][0]);
      principle_axis_2_y_grid.SetGridValue(ci, cj, axes[1][1]);
      principle_axis_2_z_grid.SetGridValue(ci, cj, axes[1][2]);
      principle_axis_3_x_grid.SetGridValue(ci, cj, axes[2][0]);
      principle_axis_3_y_grid.SetGridValue(ci, cj, axes[2][1]);
      principle_axis_3_z_grid.SetGridValue(ci, cj, axes[2][2]);
      lambda_1_grid.SetGridValue(ci, cj, variances[0]);
      lambda_2_grid.SetGridValue(ci, cj, variances[1]);
      lambda_3_grid.SetGridValue(ci, cj, variances[2]);
    }
  }

  // Delete temporary memory
  delete [] points;

  // Write images
  sprintf(image_name, "%s/%s/%02d_%02d_%s_Density.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  density_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_NearestNeighbor1.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  nearest_neighbor_1_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_NearestNeighbor2.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  nearest_neighbor_2_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_NearestNeighbor4.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  nearest_neighbor_4_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_NearestNeighbor8.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  nearest_neighbor_8_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_NearestNeighbor16.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  nearest_neighbor_16_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PrincipleAxis1X.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  principle_axis_1_x_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PrincipleAxis1Y.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  principle_axis_1_y_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PrincipleAxis1Z.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  principle_axis_1_z_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PrincipleAxis2X.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  principle_axis_2_x_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PrincipleAxis2Y.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  principle_axis_2_y_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PrincipleAxis2Z.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  principle_axis_2_z_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PrincipleAxis3X.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  principle_axis_3_x_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PrincipleAxis3Y.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  principle_axis_3_y_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_PrincipleAxis3Z.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  principle_axis_3_z_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_Lambda1.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  lambda_1_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_Lambda2.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  lambda_2_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_Lambda3.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  lambda_3_grid.Write(image_name);

  // Write more images
  R2Grid lambda_21_grid(lambda_2_grid);  lambda_21_grid.Divide(lambda_1_grid);
  R2Grid lambda_31_grid(lambda_3_grid);  lambda_31_grid.Divide(lambda_1_grid);
  R2Grid lambda_32_grid(lambda_3_grid);  lambda_32_grid.Divide(lambda_2_grid);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_Lambda21.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  lambda_21_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_Lambda31.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  lambda_31_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_Lambda32.grd", output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  lambda_32_grid.Write(image_name);

  // Return success
  return 1;
}



static int
WriteColorImages(GSVScan *scan, const char *output_image_directory, const char *image_parameterization)
{
  // Parameters
  RNScalar color_image_zoom_factor = 1;

  // Get convenient variables
  if (!scan) return 1;
  if (scan->NScanlines() == 0) return 1;
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;
  GSVScene *scene = run->Scene();
  const char *cache_directory = scene->CacheDataDirectoryName();
  char image_name[4096];
 
  // Skip if already done
  int camera_index = (scan_index == 0) ? 6 : 2;
  sprintf(image_name, "%s/%s/%02d_%02d_%02d_%s_Color.bmp", 
    output_image_directory, run->Name(), segment_index, scan_index, camera_index, image_parameterization);
  if (RNFileExists(image_name)) return 1;
  
  // Print message
  if (print_debug) {
    printf("    Creating color images ...\n");
    fflush(stdout);
  }

  // Read scanline image
  R2Grid scanline_grid, timestamp_grid;
  sprintf(image_name, "%s/laser_images/%s/%02d_%02d_%s_Scanline.grd", cache_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!scanline_grid.Read(image_name)) return 0;

  // Read position images
  R2Grid position_x_grid, position_y_grid, position_z_grid;
  sprintf(image_name, "%s/laser_images/%s/%02d_%02d_%s_PositionX.grd", cache_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!position_x_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/laser_images/%s/%02d_%02d_%s_PositionY.grd", cache_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!position_y_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/laser_images/%s/%02d_%02d_%s_PositionZ.grd", cache_directory, run->Name(), segment_index, scan_index, image_parameterization);
  if (!position_z_grid.Read(image_name)) return 0;

  // Compute world to grid transformation
  int xres = (int) (color_image_zoom_factor * scanline_grid.XResolution());
  int yres = (int) (color_image_zoom_factor * scanline_grid.YResolution());

  // Create color images for scan
  R2Affine world_to_grid = R2identity_affine;
  world_to_grid.Scale(color_image_zoom_factor);
  world_to_grid.Transform(scanline_grid.WorldToGridTransformation());
  R2Grid scan_red_grid(xres, yres, world_to_grid);         
  scan_red_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid scan_green_grid(scan_red_grid);
  R2Grid scan_blue_grid(scan_red_grid);
  R2Grid scan_timestamp_grid(scan_red_grid);

  // Consider tapestries relevant to scan
  int ntapestries[3] = { 3, 3, 3 };
  int tapestries[3][3] = { { 5, 6, 7 }, { 7, 0, 1 }, { 1, 2, 3 } };
  for (int k = 0; k < ntapestries[scan_index]; k++) {
    int it = tapestries[scan_index][k];
    GSVTapestry *tapestry = segment->Tapestry(it);
    GSVCamera *camera = tapestry->Camera();
    if (!camera) continue;

    // Print message
    if (print_verbose) {
      printf("      Creating color images for tapestry %02d\n", it);
      fflush(stdout);
    }

    // Initialize tapestry grids
    R2Grid tapestry_red_grid(xres, yres, world_to_grid);         
    tapestry_red_grid.Clear(R2_GRID_UNKNOWN_VALUE);
    R2Grid tapestry_green_grid(tapestry_red_grid);
    R2Grid tapestry_blue_grid(tapestry_red_grid);
    R2Grid tapestry_timestamp_grid(tapestry_red_grid);

    // Consider every image in tapestry
    for (int ii = 0; ii < tapestry->NImages(); ii++) {
      GSVImage *image = tapestry->Image(ii);

      // Get image info
      const GSVPose& image_pose = image->Pose();
      RNScalar image_timestamp = image->Timestamp();
      const R3Point& image_viewpoint = image_pose.Viewpoint();
      R3Frustum image_frustum(image_viewpoint, image_pose.Towards(), image_pose.Up(), 
        camera->XFov(), camera->YFov(), min_frustum_depth, max_frustum_depth);

      // Get distorted image
      R2Image *distorted_image = image->DistortedImage();
      if (!distorted_image) continue;

      // Consider every grid pixel
      for (int ix = 0; ix < xres; ix++) {
        for (int iy = 0; iy < yres; iy++) {
          int px = (int) (ix / color_image_zoom_factor + 0.5);
          if (px >= position_x_grid.XResolution()) continue;
          int py = (int) (iy / color_image_zoom_factor + 0.5);
          if (py >= position_x_grid.YResolution()) continue;
          RNScalar x = position_x_grid.GridValue(px, py);
          if (x == R2_GRID_UNKNOWN_VALUE) continue;
          RNScalar y = position_y_grid.GridValue(px, py);
          if (y == R2_GRID_UNKNOWN_VALUE) continue;
          RNScalar z = position_z_grid.GridValue(px, py);
          if (z == R2_GRID_UNKNOWN_VALUE) continue;
          RNScalar scanline_value = scanline_grid.GridValue(px, py);
          if (scanline_value == R2_GRID_UNKNOWN_VALUE) continue;
          int scanline_index = (int) (scanline_value + 0.5);
          GSVScanline *scanline = scan->Scanline(scanline_index);
          R3Point position(x, y, z);

          // Check if scanline timestamp is within range
          RNScalar scanline_timestamp = scanline->Timestamp();
          if (max_timestamp_difference > 0) {
            RNScalar timestamp_difference = fabs(scanline_timestamp - image->Timestamp());
            if (timestamp_difference > max_timestamp_difference) continue;
          }

          // Check if scanline viewpoint is within range
          const GSVPose& scanline_pose = scanline->Pose();
          if (max_viewpoint_distance > 0) {
            if (R3Distance(scanline_pose.Viewpoint(), image_viewpoint) > max_viewpoint_distance) continue;
          }
 
          // Check if point is within image view frustum
          if (!image_frustum.Intersects(position)) continue;
        
          // Retrieve color of point from distorted image
          R2Point image_position = image->DistortedPosition(position);
          if (image_position.X() < 0) continue;
          if (image_position.Y() < 0) continue;
          int image_ix1 = (int) image_position.X();
          int image_iy1 = (int) image_position.Y();
          if ((image_ix1 < 0) || (image_ix1 >= distorted_image->Width())) continue;
          if ((image_iy1 < 0) || (image_iy1 >= distorted_image->Height())) continue;
          int image_ix2 = image_ix1 + 1;
          int image_iy2 = image_iy1 + 1;
          if ((image_ix2 < 0) || (image_ix2 >= distorted_image->Width())) continue;
          if ((image_iy2 < 0) || (image_iy2 >= distorted_image->Height())) continue;
          RNRgb color11 = distorted_image->PixelRGB(image_ix1, image_iy1);
          RNRgb color12 = distorted_image->PixelRGB(image_ix1, image_iy2);
          RNRgb color21 = distorted_image->PixelRGB(image_ix2, image_iy1);
          RNRgb color22 = distorted_image->PixelRGB(image_ix2, image_iy2);
          RNScalar tx = image_position.X() - image_ix1;
          RNScalar ty = image_position.Y() - image_iy1;
          RNRgb colorA = (1-tx)*color11 + tx*color21;
          RNRgb colorB = (1-tx)*color12 + tx*color22;
          RNRgb color = (1-ty)*colorA + ty*colorB; 

          // Update tapestry grid
          RNScalar tapestry_timestamp = tapestry_timestamp_grid.GridValue(ix, iy);
          if ((tapestry_timestamp == R2_GRID_UNKNOWN_VALUE) || 
              (fabs(scanline_timestamp - image_timestamp) < fabs(scanline_timestamp - tapestry_timestamp))) {
            tapestry_red_grid.SetGridValue(ix, iy, color.R());
            tapestry_green_grid.SetGridValue(ix, iy, color.G());
            tapestry_blue_grid.SetGridValue(ix, iy, color.B());
            tapestry_timestamp_grid.SetGridValue(ix, iy, image_timestamp);
          }
              
          // Update scan grid
          RNScalar scan_timestamp = scan_timestamp_grid.GridValue(ix, iy);
          if ((scan_timestamp == R2_GRID_UNKNOWN_VALUE) || 
              (fabs(scanline_timestamp - image_timestamp) < fabs(scanline_timestamp - scan_timestamp))) {
            scan_red_grid.SetGridValue(ix, iy, color.R());
            scan_green_grid.SetGridValue(ix, iy, color.G());
            scan_blue_grid.SetGridValue(ix, iy, color.B());
            scan_timestamp_grid.SetGridValue(ix, iy, image_timestamp);
          }
        }
      }
      
      // Delete distorted image
      delete distorted_image;
    }

    // Create tapestry color image
    R2Image tapestry_rgb_image(xres, yres);
    for (int ix = 0; ix < xres; ix++) {
      for (int iy = 0; iy < yres; iy++) {
        RNScalar r = tapestry_red_grid.GridValue(ix, iy);
        if (r == R2_GRID_UNKNOWN_VALUE) continue;
        RNScalar g = tapestry_green_grid.GridValue(ix, iy);
        if (g == R2_GRID_UNKNOWN_VALUE) continue;
        RNScalar b = tapestry_blue_grid.GridValue(ix, iy);
        if (b == R2_GRID_UNKNOWN_VALUE) continue;
        RNRgb color(r, g, b);
        tapestry_rgb_image.SetPixelRGB(ix, iy, color);
      }
    }

    // Write tapestry grids 
    sprintf(image_name, "%s/%s/%02d_%02d_%02d_%s_Color.bmp", 
      output_image_directory, run->Name(), segment_index, scan_index, it, image_parameterization);
    tapestry_rgb_image.Write(image_name);
    sprintf(image_name, "%s/%s/%02d_%02d_%02d_%s_ColorRed.grd", 
      output_image_directory, run->Name(), segment_index, scan_index, it, image_parameterization);
    tapestry_red_grid.Write(image_name);
    sprintf(image_name, "%s/%s/%02d_%02d_%02d_%s_ColorGreen.grd", 
      output_image_directory, run->Name(), segment_index, scan_index, it, image_parameterization);
    tapestry_green_grid.Write(image_name);
    sprintf(image_name, "%s/%s/%02d_%02d_%02d_%s_ColorBlue.grd", 
      output_image_directory, run->Name(), segment_index, scan_index, it, image_parameterization);
    tapestry_blue_grid.Write(image_name);
    sprintf(image_name, "%s/%s/%02d_%02d_%02d_%s_ColorTimestamp.grd", 
      output_image_directory, run->Name(), segment_index, scan_index, it, image_parameterization);
    tapestry_timestamp_grid.Write(image_name);
  }

  // Create scan color image
  R2Image scan_rgb_image(xres, yres);
  for (int ix = 0; ix < xres; ix++) {
    for (int iy = 0; iy < yres; iy++) {
      RNScalar r = scan_red_grid.GridValue(ix, iy);
      if (r == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar g = scan_green_grid.GridValue(ix, iy);
      if (g == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar b = scan_blue_grid.GridValue(ix, iy);
      if (b == R2_GRID_UNKNOWN_VALUE) continue;
      RNRgb color(r, g, b);
      scan_rgb_image.SetPixelRGB(ix, iy, color);
    }
  }

  // Write scan grids 
  sprintf(image_name, "%s/%s/%02d_%02d_%s_Color.bmp", 
    output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  scan_rgb_image.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_ColorRed.grd", 
    output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  scan_red_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_ColorGreen.grd", 
    output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  scan_green_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_ColorBlue.grd", 
    output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  scan_blue_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_%s_ColorTimestamp.grd", 
    output_image_directory, run->Name(), segment_index, scan_index, image_parameterization);
  scan_timestamp_grid.Write(image_name);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// XY images (overhead)
////////////////////////////////////////////////////////////////////////

static int
WriteXYBaseImages(GSVScan *scan, const char *output_image_directory)
{
  // Get convenient variables
  if (!scan) return 1;
  if (scan->NScanlines() == 0) return 1;
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;

  // Skip if already done
  char image_name[4096];
  sprintf(image_name, "%s/%s/%02d_%02d_XY_Density.grd", output_image_directory, run->Name(), segment_index, scan_index);
  if (RNFileExists(image_name)) return 1;
  
  // Print message
  if (print_debug) {
    printf("    Creating base images ...\n");
    fflush(stdout);
  }

  // Compute grid resolution
  const R3Box& scan_bbox = scan->BBox();
  R2Box grid_bbox(scan_bbox[0][0], scan_bbox[0][1], scan_bbox[1][0], scan_bbox[1][1]);
  int xres = (int) (grid_bbox.XLength() / XY_image_spacing + 0.5);
  int yres = (int) (grid_bbox.YLength() / XY_image_spacing + 0.5);
  if ((xres == 0) || (yres == 0)) return 0;

  // Initialize grids
  R2Grid density_grid(xres, yres, grid_bbox);       density_grid.Clear(0);
  R2Grid distance_grid(xres, yres, grid_bbox);      distance_grid.Clear(0);
  R2Grid minZ_grid(xres, yres, grid_bbox);          minZ_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid maxZ_grid(xres, yres, grid_bbox);          maxZ_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid groundZ_grid(xres, yres, grid_bbox);       groundZ_grid.Clear(R2_GRID_UNKNOWN_VALUE);

  // Read scan points
  if (!scan->ReadPoints()) {
    fprintf(stderr, "Unable to read points for %s %02d %02d\n", run->Name(), segment_index, scan_index);
    return 0;
  }

  // Fill other grids
  for (int ie = 0; ie < scan->NScanlines(); ie++) {
    GSVScanline *scanline = scan->Scanline(ie);
    const GSVPose& pose = scanline->Pose();
    const R3Point& viewpoint = pose.Viewpoint();
    // const R3Vector& towards = pose.Towards();
    // const R3Vector& up = pose.Up();
    // R3Vector right = towards % up;
    RNCoord groundZ = scanline->EstimatedGroundZ();
    if (groundZ == RN_UNKNOWN) continue;
    for (int j = 0; j < scanline->NPoints(); j++) {
      R3Point position = scanline->PointPosition(j);
      RNLength distance = R2Distance(R2Point(position.X(), position.Y()), R2Point(viewpoint.X(), viewpoint.Y()));
      R2Point grid_position = density_grid.GridPosition(R2Point(position.X(), position.Y()));
      int ix = (int) (grid_position.X() + 0.5);
      int iy = (int) (grid_position.Y() + 0.5);
      RNScalar previous_distance = distance_grid.GridValue(ix, iy);
      RNScalar previous_minZ = minZ_grid.GridValue(ix, iy);
      RNScalar previous_maxZ = maxZ_grid.GridValue(ix, iy);
      RNScalar previous_groundZ = groundZ_grid.GridValue(ix, iy);
      density_grid.RasterizeWorldPoint(position[0], position[1], 1);
      if ((previous_distance == R2_GRID_UNKNOWN_VALUE) || (distance < previous_distance))
        distance_grid.SetGridValue(ix, iy, distance);
      if ((previous_minZ == R2_GRID_UNKNOWN_VALUE) || (position.Z() < previous_minZ))
        minZ_grid.SetGridValue(ix, iy, position.Z());
      if ((previous_maxZ == R2_GRID_UNKNOWN_VALUE) || (position.Z() > previous_maxZ))
        maxZ_grid.SetGridValue(ix, iy, position.Z());
      if ((previous_groundZ == R2_GRID_UNKNOWN_VALUE) || (previous_distance == R2_GRID_UNKNOWN_VALUE) || (distance < previous_distance))
        groundZ_grid.SetGridValue(ix, iy, groundZ);
    }
  }

  // Re-Fill distance grid
  distance_grid.Clear(0);
  for (int ie = 0; ie < scan->NScanlines(); ie++) {
    GSVScanline *scanline = scan->Scanline(ie);
    const GSVPose& pose = scanline->Pose();
    const R3Point& viewpoint = pose.Viewpoint();
    R2Point grid_position = density_grid.GridPosition(R2Point(viewpoint.X(), viewpoint.Y()));
    int ix = (int) (grid_position.X() + 0.5);
    int iy = (int) (grid_position.Y() + 0.5);
    distance_grid.SetGridValue(ix, iy, 1.0);
  }

  // Re-compute distance grid
  distance_grid.SquaredDistanceTransform();
  distance_grid.Sqrt();

  // Release scan points
  if (!scan->ReleasePoints()) {
    fprintf(stderr, "Unable to release points for %s %02d %02d\n", run->Name(), segment_index, scan_index);
    return 0;
  }

  // Write grids
  sprintf(image_name, "%s/%s/%02d_%02d_XY_Density.grd", output_image_directory, run->Name(), segment_index, scan_index);
  density_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_XY_ViewpointDistance.grd", output_image_directory, run->Name(), segment_index, scan_index);
  distance_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_XY_MinZ.grd", output_image_directory, run->Name(), segment_index, scan_index);
  minZ_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_XY_MaxZ.grd", output_image_directory, run->Name(), segment_index, scan_index);
  maxZ_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_XY_GroundZ.grd", output_image_directory, run->Name(), segment_index, scan_index);
  groundZ_grid.Write(image_name);

  // Return success
  return 1;
}



static int
WriteXYImages(GSVScan *scan, const char *output_image_directory)
{
  // Get convenient variables
  if (!scan) return 1;
  if (scan->NScanlines() == 0) return 1;
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;

  // Print message
  if (print_verbose) {
    printf("  Creating XY images for %s %02d %02d\n", run->Name(), segment_index, scan_index);
    fflush(stdout);
  }

  // Capture XY images
  if (include_base_images) {
    if (!WriteXYBaseImages(scan, output_image_directory)) return 0;
  }

  // Return success
  return 1;
}


static int
WriteXYImages(GSVSegment *segment, const char *output_image_directory)
{
  // Get convenient variables
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;

  // Print message
  if (print_verbose) {
    printf("  Creating XY images for %s %02d\n", run->Name(), segment_index);
    fflush(stdout);
  }

  // Combine XY images from scans
  if (include_base_images) {
    for (int i = 0; i < num_xy_image_names; i++) {
      // Skip if already done
      char image_name[1024];
      sprintf(image_name, "%s/%s/%02d_XY_%s.grd", output_image_directory, run->Name(), segment_index, xy_image_names[i]);
      if (RNFileExists(image_name)) continue;

      // Read scan grids
      RNArray<R2Grid *> scan_grids;
      for (int ia = 0; ia < segment->NScans(); ia++) {
        if (ia == 1) continue;
        R2Grid *scan_grid = new R2Grid();
        char scan_image_name[1024];
        sprintf(scan_image_name, "%s/%s/%02d_%02d_XY_%s.grd", output_image_directory, run->Name(), segment_index, ia, xy_image_names[i]);
        if (!scan_grid->Read(scan_image_name)) return 0;
        scan_grids.Insert(scan_grid);
      }

      // Compute grid resolution
      const R3Box& segment_bbox = segment->BBox();
      R2Box grid_bbox(segment_bbox[0][0], segment_bbox[0][1], segment_bbox[1][0], segment_bbox[1][1]);
      int xres = (int) (grid_bbox.XLength() / XY_image_spacing + 0.5);
      int yres = (int) (grid_bbox.YLength() / XY_image_spacing + 0.5);
      if ((xres == 0) || (yres == 0)) return 0;

      // Create combined grid
      R2Grid combined_grid(xres, yres, grid_bbox);
      if (!CombineImages(combined_grid, scan_grids, xy_image_combination_methods[i])) {
        return 0;
      }

      // Fill missing values
      if (!strcmp(xy_image_names[i], "GroundZ")) {
        combined_grid.FillHoles(9999999);
      }

      // Delete scan grids
      for (int ia = 0; ia < scan_grids.NEntries(); ia++) {
        delete scan_grids[ia];
      }

      // Write combined grid
      if (!combined_grid.Write(image_name)) return 0;
    }
  }

  // Return success
  return 1;
}

static int
WriteXYImages(GSVRun *run, const char *output_image_directory)
{
  // Print message
  if (print_verbose) {
    printf("  Creating XY images for %s\n", run->Name());
    fflush(stdout);
  }

  // Combine XY images from segments
  if (include_base_images) {
    for (int i = 0; i < num_xy_image_names; i++) {
      // Skip if already done
      char image_name[1024];
      sprintf(image_name, "%s/%s/XY_%s.grd", output_image_directory, run->Name(), xy_image_names[i]);
      if (RNFileExists(image_name)) continue;

      // Read segment grids
      RNArray<R2Grid *> segment_grids;
      for (int is = 0; is < run->NSegments(); is++) {
        R2Grid *segment_grid = new R2Grid();
        char segment_image_name[1024];
        sprintf(segment_image_name, "%s/%s/%02d_XY_%s.grd", output_image_directory, run->Name(), is, xy_image_names[i]);
        if (!segment_grid->Read(segment_image_name)) return 0;
        segment_grids.Insert(segment_grid);
      }

      // Compute grid resolution
      const R3Box& run_bbox = run->BBox();
      R2Box grid_bbox(run_bbox[0][0], run_bbox[0][1], run_bbox[1][0], run_bbox[1][1]);
      int xres = (int) (grid_bbox.XLength() / XY_image_spacing + 0.5);
      int yres = (int) (grid_bbox.YLength() / XY_image_spacing + 0.5);
      if ((xres == 0) || (yres == 0)) return 0;

      // Create combined grid
      R2Grid combined_grid(xres, yres, grid_bbox);
      if (!CombineImages(combined_grid, segment_grids, xy_image_combination_methods[i])) {
        return 0;
      }

      // Delete segment grids
      for (int is = 0; is < segment_grids.NEntries(); is++) {
        delete segment_grids[is];
      }

      // Write combined grid
      if (!combined_grid.Write(image_name)) return 0;
    }
  }

  // Return success
  return 1;
}



static int
WriteXYImages(GSVScene *scene, const char *output_image_directory)
{
  // Print message
  if (print_verbose) {
    printf("  Creating XY images\n");
    fflush(stdout);
  }

  // Combine XY images from segments
  if (include_base_images) {
    for (int i = 0; i < num_xy_image_names; i++) {
      // Skip if already done
      char image_name[1024];
      sprintf(image_name, "%s/XY_%s.grd", output_image_directory, xy_image_names[i]);
      if (RNFileExists(image_name)) continue;

      // Read run grids
      RNArray<R2Grid *> run_grids;
      for (int ir = 0; ir < scene->NRuns(); ir++) {
        GSVRun *run = scene->Run(ir);
        R2Grid *run_grid = new R2Grid();
        char run_image_name[1024];
        sprintf(run_image_name, "%s/%s/XY_%s.grd", output_image_directory, run->Name(), xy_image_names[i]);
        if (!run_grid->Read(run_image_name)) return 0;
        run_grids.Insert(run_grid);
      }

      // Compute grid resolution
      const R3Box& scene_bbox = scene->BBox();
      R2Box grid_bbox(scene_bbox[0][0], scene_bbox[0][1], scene_bbox[1][0], scene_bbox[1][1]);
      int xres = (int) (grid_bbox.XLength() / XY_image_spacing + 0.5);
      int yres = (int) (grid_bbox.YLength() / XY_image_spacing + 0.5);
      if ((xres == 0) || (yres == 0)) return 0;

      // Create combined grid
      R2Grid combined_grid(xres, yres, grid_bbox);
      if (!CombineImages(combined_grid, run_grids, xy_image_combination_methods[i])) {
        return 0;
      }

      // Delete run grids
      for (int ir = 0; ir < run_grids.NEntries(); ir++) {
        delete run_grids[ir];
      }

      // Write combined grid
      if (!combined_grid.Write(image_name)) return 0;
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// SA images (Scanline vs. Angle)
////////////////////////////////////////////////////////////////////////

static int
WriteSABaseImages(GSVScan *scan, const char *output_image_directory)
{
  // Get convenient variables
  if (!scan) return 1;
  if (scan->NScanlines() == 0) return 1;
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;
  
  // Skip if already done
  char image_name[4096];
  sprintf(image_name, "%s/%s/%02d_%02d_SA_Timestamp.grd", output_image_directory, run->Name(), segment_index, scan_index);
  if (RNFileExists(image_name)) return 1;
  
  // Print message
  if (print_debug) {
    printf("    Creating base images ...\n");
    fflush(stdout);
  }

  // Count scanlines
  int nscanlines = 0;
  R3Point prev_viewpoint = scan->Scanline(0)->Pose().Viewpoint();
  for (int ie = 0; ie < scan->NScanlines(); ie++) {
    GSVScanline *scanline = scan->Scanline(ie);
    const GSVPose& pose = scanline->Pose();
    const R3Point& viewpoint = pose.Viewpoint();
    if ((ie > 0) && (R3Distance(viewpoint, prev_viewpoint) < SA_viewpoint_spacing)) continue;
    prev_viewpoint = viewpoint;
    nscanlines++;
  }

  // Compute grid resolution
  int xres = nscanlines;
  int yres = 180;

  // Initialize grids
  R2Box grid_box(0, 0, scan->NScanlines(), 180);
  R2Grid timestamp_grid(xres, yres, grid_box);       timestamp_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid angle_grid(xres, yres, grid_box);           angle_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid depth_grid(xres, yres, grid_box);           depth_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid distance_grid(xres, yres, grid_box);        distance_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid height_grid(xres, yres, grid_box);          height_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid travel_distance_grid(xres, yres, grid_box); travel_distance_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid scanline_grid(xres, yres, grid_box);        scanline_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid pointindex_grid(xres, yres, grid_box);      pointindex_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid groundZ_grid(xres, yres, grid_box);         groundZ_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid positionX_grid(xres, yres, grid_box);       positionX_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid positionY_grid(xres, yres, grid_box);       positionY_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid positionZ_grid(xres, yres, grid_box);       positionZ_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid viewpointX_grid(xres, yres, grid_box);      viewpointX_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid viewpointY_grid(xres, yres, grid_box);      viewpointY_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid viewpointZ_grid(xres, yres, grid_box);      viewpointZ_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid towardsX_grid(xres, yres, grid_box);        towardsX_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid towardsY_grid(xres, yres, grid_box);        towardsY_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid towardsZ_grid(xres, yres, grid_box);        towardsZ_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid upX_grid(xres, yres, grid_box);             upX_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid upY_grid(xres, yres, grid_box);             upY_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid upZ_grid(xres, yres, grid_box);             upZ_grid.Clear(R2_GRID_UNKNOWN_VALUE);

  // Read scan points
  if (!scan->ReadPoints()) {
    fprintf(stderr, "Unable to read points for %s %02d %02d\n", run->Name(), segment_index, scan_index);
    return 0;
  }

  // Fill grids
  nscanlines = 0;
  RNScalar travel_distance = 0;
  prev_viewpoint = scan->Scanline(0)->Pose().Viewpoint();
  for (int ie = 0; ie < scan->NScanlines(); ie++) {
    GSVScanline *scanline = scan->Scanline(ie);
    const GSVPose& pose = scanline->Pose();
    const R3Point& viewpoint = pose.Viewpoint();
    if ((ie > 0) && (R3Distance(viewpoint, prev_viewpoint) < SA_viewpoint_spacing)) continue;
    travel_distance += R3Distance(viewpoint, prev_viewpoint);
    prev_viewpoint = viewpoint;
    int ix = nscanlines++;
    if (ix >= xres) continue;
    RNScalar timestamp = scanline->Timestamp();
    RNCoord groundZ = scanline->EstimatedGroundZ();
    const R3Vector& towards = pose.Towards();
    const R3Vector& up = pose.Up();
    for (int j = 0; j < scanline->NPoints(); j++) {
      R3Point position = scanline->PointPosition(j);
      RNLength height = position.Z() - groundZ;
      R3Vector v = position - viewpoint;
      RNScalar depth = v.Dot(towards);
      RNLength distance = v.Length();
      if (RNIsZero(distance)) continue;
      v /= distance;
      RNScalar dot = v.Dot(towards);
      RNAngle angle = (dot < 1) ? ((dot > -1) ? acos(dot) : RN_PI) : 0;
      if (v.Dot(up) > 0) angle = RN_PI_OVER_TWO + angle;
      else angle = RN_PI_OVER_TWO - angle;
      int iy = (int) (yres * angle / RN_PI + 0.5);
      if ((iy < 0) || (iy >= yres)) continue;

      // Remember values for closest scan point
      RNScalar old_distance_value = distance_grid.GridValue(ix, iy);
      if ((old_distance_value == R2_GRID_UNKNOWN_VALUE) || (distance < old_distance_value)) {
        timestamp_grid.SetGridValue(ix, iy, timestamp);
        angle_grid.SetGridValue(ix, iy, angle);
        depth_grid.SetGridValue(ix, iy, depth);
        distance_grid.SetGridValue(ix, iy, distance);
        height_grid.SetGridValue(ix, iy, height);
        travel_distance_grid.SetGridValue(ix, iy, travel_distance);
        scanline_grid.SetGridValue(ix, iy, ie);
        pointindex_grid.SetGridValue(ix, iy, j);
        groundZ_grid.SetGridValue(ix, iy, groundZ);
        positionX_grid.SetGridValue(ix, iy, position.X());
        positionY_grid.SetGridValue(ix, iy, position.Y());
        positionZ_grid.SetGridValue(ix, iy, position.Z());
        viewpointX_grid.SetGridValue(ix, iy, viewpoint.X());
        viewpointY_grid.SetGridValue(ix, iy, viewpoint.Y());
        viewpointZ_grid.SetGridValue(ix, iy, viewpoint.Z());
        towardsX_grid.SetGridValue(ix, iy, towards.X());
        towardsY_grid.SetGridValue(ix, iy, towards.Y());
        towardsZ_grid.SetGridValue(ix, iy, towards.Z());
        upX_grid.SetGridValue(ix, iy, up.X());
        upY_grid.SetGridValue(ix, iy, up.Y());
        upZ_grid.SetGridValue(ix, iy, up.Z());
      }
    }
  }

  // Release scan points
  if (!scan->ReleasePoints()) {
    fprintf(stderr, "Unable to release points for %s %02d %02d\n", run->Name(), segment_index, scan_index);
    return 0;
  }

  // Create hole grid
  R2Grid hole_grid(distance_grid);
  hole_grid.Clear(0);
  hole_grid.Substitute(R2_GRID_UNKNOWN_VALUE, 1);

  // Write grids
  sprintf(image_name, "%s/%s/%02d_%02d_SA_Timestamp.grd", output_image_directory, run->Name(), segment_index, scan_index);
  timestamp_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_Hole.grd", output_image_directory, run->Name(), segment_index, scan_index);
  hole_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_ViewpointAngle.grd", output_image_directory, run->Name(), segment_index, scan_index);
  angle_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_ViewpointDistance.grd", output_image_directory, run->Name(), segment_index, scan_index);
  distance_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_ViewpointDepth.grd", output_image_directory, run->Name(), segment_index, scan_index);
  depth_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_Height.grd", output_image_directory, run->Name(), segment_index, scan_index);
  height_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_TravelDistance.grd", output_image_directory, run->Name(), segment_index, scan_index);
  travel_distance_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_Scanline.grd", output_image_directory, run->Name(), segment_index, scan_index);
  scanline_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_PointIndex.grd", output_image_directory, run->Name(), segment_index, scan_index);
  pointindex_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_GroundZ.grd", output_image_directory, run->Name(), segment_index, scan_index);
  groundZ_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_PositionX.grd", output_image_directory, run->Name(), segment_index, scan_index);
  positionX_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_PositionY.grd", output_image_directory, run->Name(), segment_index, scan_index);
  positionY_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_PositionZ.grd", output_image_directory, run->Name(), segment_index, scan_index);
  positionZ_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_ViewpointX.grd", output_image_directory, run->Name(), segment_index, scan_index);
  viewpointX_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_ViewpointY.grd", output_image_directory, run->Name(), segment_index, scan_index);
  viewpointY_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_ViewpointZ.grd", output_image_directory, run->Name(), segment_index, scan_index);
  viewpointZ_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_TowardsX.grd", output_image_directory, run->Name(), segment_index, scan_index);
  towardsX_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_TowardsY.grd", output_image_directory, run->Name(), segment_index, scan_index);
  towardsY_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_TowardsZ.grd", output_image_directory, run->Name(), segment_index, scan_index);
  towardsZ_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_UpX.grd", output_image_directory, run->Name(), segment_index, scan_index);
  upX_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_UpY.grd", output_image_directory, run->Name(), segment_index, scan_index);
  upY_grid.Write(image_name);
  sprintf(image_name, "%s/%s/%02d_%02d_SA_UpZ.grd", output_image_directory, run->Name(), segment_index, scan_index);
  upZ_grid.Write(image_name);

  // Return success
  return 1;
}



static int
WriteSAImages(GSVScan *scan, const char *output_image_directory)
{
  // Get convenient variables
  if (!scan) return 1;
  if (scan->NScanlines() == 0) return 1;
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;
  
  // Print message
  if (print_verbose) {
    printf("  Creating SA images for %s %02d %02d\n", run->Name(), segment_index, scan_index);
    fflush(stdout);
  }

  // Capture SA images
  if (include_base_images) {
    if (!WriteSABaseImages(scan, output_image_directory)) return 0;
  }
  if (include_normal_images) {
    if (!WriteNormalImages(scan, output_image_directory, "SA")) return 0;
  }
  if (include_curvature_images) {
    if (!WriteCurvatureImages(scan, output_image_directory, "SA")) return 0;
  }
  if (include_boundary_images) {
    if (!WriteBoundaryImages(scan, output_image_directory, "SA")) return 0;
  }
  if (include_pca_images) {
    if (!WritePCAImages(scan, output_image_directory, "SA")) return 0;
  }
  if (include_color_images) {
    if (!WriteColorImages(scan, output_image_directory, "SA")) return 0;
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// DA images (Distance vs. Angle)
////////////////////////////////////////////////////////////////////////

static int
WriteDAImage(GSVScan *scan, const char *output_image_directory, 
  const R2Grid& sa_travel_distance_grid, const R2Grid& sa_viewpoint_distance_grid, 
  const char *grid_name)
{
  // Get convenient variables
  if (!scan) return 0;
  if (scan->NScanlines() == 0) return 0;
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;

  // Skip if already done
  char da_name[4096];
  sprintf(da_name, "%s/%s/%02d_%02d_DA_%s.grd", output_image_directory, run->Name(), segment_index, scan_index, grid_name);
  if (RNFileExists(da_name)) return 1;
  
  // Read SA image
  R2Grid sa_grid;
  char sa_name[4096];
  sprintf(sa_name, "%s/%s/%02d_%02d_SA_%s.grd", output_image_directory, run->Name(), segment_index, scan_index, grid_name);
  sa_grid.Read(sa_name);

  // Compute grid resolution
  RNScalar total_travel_distance = sa_travel_distance_grid.Maximum();
  if (total_travel_distance == 0) return 0;
  int xres = (int) (total_travel_distance / DA_image_spacing + 0.5);
  if (xres == 0) return 0;
  int yres = 180;

  // Compute world to grid transformation
  R2Affine world_to_grid = R2identity_affine;
  world_to_grid.XScale(1.0 / DA_image_spacing);

  // Create DA image
  R2Grid da_grid(xres, yres, world_to_grid);
  da_grid.Clear(R2_GRID_UNKNOWN_VALUE);

  // Create DA distance image
  R2Grid da_viewpoint_distance_grid(xres, yres, world_to_grid);
  da_viewpoint_distance_grid.Clear(FLT_MAX);

  // Fill DA image
  for (int i = 0; i < sa_grid.XResolution(); i++) {
    for (int iy = 0; iy < sa_grid.YResolution(); iy++) {
      RNScalar value = sa_grid.GridValue(i, iy);
      if (value == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar viewpoint_distance = sa_viewpoint_distance_grid.GridValue(i, iy);
      if (viewpoint_distance == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar travel_distance = sa_travel_distance_grid.GridValue(i, iy);
      if (travel_distance == R2_GRID_UNKNOWN_VALUE) continue;
      int ix = (int) (xres * travel_distance / total_travel_distance + 0.5);
      if ((ix < 0) || (ix >= xres)) continue;
      RNScalar old_viewpoint_distance = da_viewpoint_distance_grid(ix, iy);
      if (viewpoint_distance < old_viewpoint_distance)  {
        da_viewpoint_distance_grid.SetGridValue(ix, iy, viewpoint_distance);
        da_grid.SetGridValue(ix, iy, value);
      }
    }
  }
  
  // Fill holes in DA image
  InterpolateMissingColumns(da_grid);

  // Write DA image
  da_grid.Write(da_name);

  // Return success
  return 1;
}



static int
WriteDAImages(GSVScan *scan, const char *output_image_directory)
{
  // Get convenient variables
  if (!scan) return 1;
  if (scan->NScanlines() == 0) return 1;
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;
  
  // Print message
  if (print_verbose) {
    printf("  Creating DA images for %s %02d %02d\n", run->Name(), segment_index, scan_index);
    fflush(stdout);
  }

  // Read SA scanline and distance images
  char sa_travel_distance_name[4096];
  char sa_viewpoint_distance_name[4096];
  R2Grid sa_travel_distance_grid;
  R2Grid sa_viewpoint_distance_grid;
  sprintf(sa_travel_distance_name, "%s/%s/%02d_%02d_SA_TravelDistance.grd", output_image_directory, run->Name(), segment_index, scan_index);
  sprintf(sa_viewpoint_distance_name, "%s/%s/%02d_%02d_SA_ViewpointDistance.grd", output_image_directory, run->Name(), segment_index, scan_index);
  if (!RNFileExists(sa_travel_distance_name) && (!WriteSAImages(scan, output_image_directory))) return 0;
  sa_travel_distance_grid.Read(sa_travel_distance_name);
  sa_viewpoint_distance_grid.Read(sa_viewpoint_distance_name);

  // Write DA images
  if (include_base_images) {
    for (int i = 0; i < num_base_image_names; i++) {
      if (!WriteDAImage(scan, output_image_directory, sa_travel_distance_grid, sa_viewpoint_distance_grid, base_image_names[i])) return 0;
    }
  }
#if 1
  if (include_normal_images) {
    if (!WriteNormalImages(scan, output_image_directory, "DA")) return 0;
  }
  if (include_curvature_images) {
    if (!WriteCurvatureImages(scan, output_image_directory, "DA")) return 0;
  }
  if (include_boundary_images) {
    if (!WriteBoundaryImages(scan, output_image_directory, "DA")) return 0;
  }
  if (include_pca_images) {
    if (!WritePCAImages(scan, output_image_directory, "DA")) return 0;
  }
#else
  if (include_normal_images) {
    for (int i = 0; i < num_normal_image_names; i++) {
      if (!WriteDAImage(scan, output_image_directory, sa_travel_distance_grid, sa_viewpoint_distance_grid, normal_image_names[i])) return 0;
    }
  }
  if (include_curvature_images) {
    for (int i = 0; i < num_curvature_image_names; i++) {
      if (!WriteDAImage(scan, output_image_directory, sa_travel_distance_grid, sa_viewpoint_distance_grid, curvature_image_names[i])) return 0;
    }
  }
  if (include_boundary_images) {
    for (int i = 0; i < num_boundary_image_names; i++) {
      if (!WriteDAImage(scan, output_image_directory, sa_travel_distance_grid, sa_viewpoint_distance_grid, boundary_image_names[i])) return 0;
    }
  }
  if (include_pca_images) {
    for (int i = 0; i < num_pca_image_names; i++) {
      if (!WriteDAImage(scan, output_image_directory, sa_travel_distance_grid, sa_viewpoint_distance_grid, pca_image_names[i])) return 0;
    }
  }
#endif
  if (include_color_images) {
    for (int i = 0; i < num_color_image_names; i++) {
      if (!WriteDAImage(scan, output_image_directory, sa_travel_distance_grid, sa_viewpoint_distance_grid, color_image_names[i])) return 0;
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// DH images (Distance vs. Height)
////////////////////////////////////////////////////////////////////////

static int
WriteDHImage(GSVScan *scan, const char *output_image_directory, 
  const R2Grid& da_distance_grid, const R2Grid& da_height_grid,
  const char *grid_name)
{
  // Get convenient variables
  if (!scan) return 0;
  if (scan->NScanlines() == 0) return 0;
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;

  // Skip if already done
  char dh_name[4096];
  sprintf(dh_name, "%s/%s/%02d_%02d_DH_%s.grd", output_image_directory, run->Name(), segment_index, scan_index, grid_name);
  if (RNFileExists(dh_name)) return 1;
  
  // Read DA image
  R2Grid da_grid;
  char da_name[4096];
  sprintf(da_name, "%s/%s/%02d_%02d_DA_%s.grd", output_image_directory, run->Name(), segment_index, scan_index, grid_name);
  da_grid.Read(da_name);

  // Compute world to grid transformation
  int yres = (int) (DH_max_height / DH_image_spacing + 0.5);
  R2Affine world_to_grid = R2identity_affine;
  world_to_grid.Scale(1.0 / DH_image_spacing);

  // Fill holes in DA image
  int max_hole_size = (int) (DH_max_hole_size / DH_image_spacing + 0.5);
  da_grid.FillHoles(max_hole_size);

  // Create DH image
  R2Grid dh_grid(da_grid.XResolution(), yres, world_to_grid);
  dh_grid.Clear(R2_GRID_UNKNOWN_VALUE);

  // Create DH distance image
  R2Grid dh_distance_grid(da_grid.XResolution(), yres, world_to_grid);
  dh_distance_grid.Clear(FLT_MAX);

  // Fill DH image
  for (int ix = 0; ix < da_grid.XResolution(); ix++) {
    for (int j = scan_ground_index; j < da_grid.YResolution(); j++) {
      RNScalar height = da_height_grid.GridValue(ix, j);
      if (height == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar distance = da_distance_grid.GridValue(ix, j);
      if (distance == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar value = da_grid.GridValue(ix, j);
      if (value == R2_GRID_UNKNOWN_VALUE) continue;
      int iy = (int) (yres * height / DH_max_height + 0.5);
      if ((iy < 0) || (iy >= yres)) continue;
      RNScalar old_distance = dh_distance_grid(ix, iy);
      if (distance < old_distance)  {
        dh_grid.SetGridValue(ix, iy, value);
        dh_distance_grid.SetGridValue(ix, iy, distance);
      }
    }
  }
  
  // Fill holes in DH image
  dh_grid.FillHoles(max_hole_size);

  // Write DH image
  dh_grid.Write(dh_name);

  // Return success
  return 1;
}



static int
WriteDHImages(GSVScan *scan, const char *output_image_directory)
{
  // Get convenient variables
  if (!scan) return 1;
  if (scan->NScanlines() == 0) return 1;
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;
  
  // Print message
  if (print_verbose) {
    printf("  Creating DH images for %s %02d %02d\n", run->Name(), segment_index, scan_index);
    fflush(stdout);
  }

  // Read DA height and distance images
  char da_height_name[4096];
  char da_distance_name[4096];
  R2Grid da_height_grid;
  R2Grid da_distance_grid;
  sprintf(da_height_name, "%s/%s/%02d_%02d_DA_Height.grd", output_image_directory, run->Name(), segment_index, scan_index);
  sprintf(da_distance_name, "%s/%s/%02d_%02d_DA_ViewpointDistance.grd", output_image_directory, run->Name(), segment_index, scan_index);
  if (!RNFileExists(da_height_name) && (!WriteDAImages(scan, output_image_directory))) return 0;
  da_height_grid.Read(da_height_name);
  da_distance_grid.Read(da_distance_name);

  // Fill holes in DA height and distance images
  int max_hole_size = (int) (DH_max_hole_size / DH_image_spacing + 0.5);
  da_height_grid.FillHoles(max_hole_size);
  da_distance_grid.FillHoles(max_hole_size);

  // Write DH images
  if (include_base_images) {
    for (int i = 0; i < num_base_image_names; i++) {
      if (!WriteDHImage(scan, output_image_directory, da_distance_grid, da_height_grid, base_image_names[i])) return 0;
    }
  }
  if (include_normal_images) {
    for (int i = 0; i < num_normal_image_names; i++) {
      if (!WriteDHImage(scan, output_image_directory, da_distance_grid, da_height_grid, normal_image_names[i])) return 0;
    }
  }
  if (include_curvature_images) {
    for (int i = 0; i < num_curvature_image_names; i++) {
      if (!WriteDHImage(scan, output_image_directory, da_distance_grid, da_height_grid, curvature_image_names[i])) return 0;
    }
  }
  if (include_boundary_images) {
    for (int i = 0; i < num_boundary_image_names; i++) {
      if (!WriteDHImage(scan, output_image_directory, da_distance_grid, da_height_grid, boundary_image_names[i])) return 0;
    }
  }
  if (include_pca_images) {
    for (int i = 0; i < num_pca_image_names; i++) {
      if (!WriteDHImage(scan, output_image_directory, da_distance_grid, da_height_grid, pca_image_names[i])) return 0;
    }
  }
  if (include_color_images) {
    for (int i = 0; i < num_color_image_names; i++) {
      if (!WriteDHImage(scan, output_image_directory, da_distance_grid, da_height_grid, color_image_names[i])) return 0;
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// UV images (undistorted image pixel coordinates)
////////////////////////////////////////////////////////////////////////

static void
RasterizeTriangle(R2Grid& depth_grid, R2Grid& sa_index_grid, 
  const R2Point position[3], RNScalar depth[3], int grid_index[3])
{
  // Sort vertex indices by Y coordinate
  int iv0, iv1, iv2;
  if (position[0].Y() < position[1].Y()) {
    if (position[0].Y() < position[2].Y()) { 
      if (position[1].Y() < position[2].Y()) { iv0 = 0; iv1 = 1; iv2 = 2; }
      else { iv0 = 0; iv1 = 2; iv2 = 1; }
    }
    else { iv0 = 2; iv1 = 0; iv2 = 1; }
  }
  else {
    if (position[1].Y() < position[2].Y()) { 
      if (position[0].Y() < position[2].Y()) { iv0 = 1; iv1 = 0; iv2 = 2; }
      else { iv0 = 1; iv1 = 2; iv2 = 0; }
    }
    else { iv0 = 2; iv1 = 1; iv2 = 0; }
  }

  // Sort vertex coordinates
  double dy = position[iv2].Y() - position[iv0].Y();
  double t1 = (dy > 0) ? (position[iv1].Y() - position[iv0].Y()) / dy : 0;
  double y0 = position[iv0].Y();
  double y1 = position[iv1].Y();
  double y2 = position[iv2].Y();
  double x0 = position[iv0].X();
  double x2 = position[iv2].X();
  double x1a = position[iv1].X();
  double x1b = (1-t1)*x0 + t1*x2;
  double depth0 = depth[iv0];
  double depth2 = depth[iv2];
  double depth1a = depth[iv1];
  double depth1b = (1-t1)*depth0 + t1*depth2;
  double gridindex0 = grid_index[iv0];
  double gridindex1 = grid_index[iv1];
  double gridindex2 = grid_index[iv2];
  if (x1a > x1b) { 
    double swap;
    swap = x1a; x1a = x1b; x1b = swap; 
    swap = depth1a; depth1a = depth1b; depth1b = swap; 
  }

  // Rasterize lower half of triangle
  int iy0 = (int) (y0 + 0.5);
  int iy1 = (int) (y1 + 0.5);
  int nysteps = (iy1 - iy0) + 1;
  double xa_step = (x1a - x0) / nysteps;
  double xb_step = (x1b - x0) / nysteps;
  double deptha_step = (depth1a  - depth0) / nysteps;
  double depthb_step = (depth1b  - depth0) / nysteps;
  double xa = x0;
  double xb = x0;
  double deptha = depth0;
  double depthb = depth0;
  for (int iy = iy0; iy <= iy1; iy++) {
    if ((iy >= 0) && (iy < depth_grid.YResolution())) {
      int ixa = (int) (xa + 0.5);
      int ixb = (int) (xb + 0.5);
      int nxsteps = (ixb - ixa) + 1;
      double depth_step = (depthb - deptha) / nxsteps;
      double depth = deptha;
      double gridindex = (iy - iy0 < iy1 - iy) ? gridindex0 : gridindex1;
      for (int ix = ixa; ix <= ixb; ix++) {
        if ((ix >= 0) && (ix < depth_grid.XResolution())) {
          if (depth < depth_grid.GridValue(ix, iy)) {
            depth_grid.SetGridValue(ix, iy, depth);
            sa_index_grid.SetGridValue(ix, iy, gridindex);
          }
        }
        depth += depth_step;
      }
    }
    xa += xa_step;
    xb += xb_step;
    deptha += deptha_step;
    depthb += depthb_step;
  }

  // Rasterize upper half of triangle
  int iy2 = (int) (y2 + 0.5);
  nysteps = (iy2 - iy1) + 1;
  xa_step = (x2 - x1a) / nysteps;
  xb_step = (x2 - x1b) / nysteps;
  deptha_step = (depth2  - depth1a) / nysteps;
  depthb_step = (depth2  - depth1b) / nysteps;
  xa = x1a;
  xb = x1b;
  deptha = depth1a;
  depthb = depth1b;
  for (int iy = iy1; iy <= iy2; iy++) {
    if ((iy >= 0) && (iy < depth_grid.YResolution())) {
      int ixa = (int) (xa + 0.5);
      int ixb = (int) (xb + 0.5);
      int nxsteps = (ixb - ixa) + 1;
      double depth_step = (depthb - deptha) / nxsteps;
      double depth = deptha;
      double gridindex = (iy2 - iy < iy - iy1) ? gridindex2 : gridindex1;
      for (int ix = ixa; ix <= ixb; ix++) {
        if ((ix >= 0) && (ix < depth_grid.XResolution())) {
          if (depth < depth_grid.GridValue(ix, iy)) {
            depth_grid.SetGridValue(ix, iy, depth);
            sa_index_grid.SetGridValue(ix, iy, gridindex);
          }
        }
        depth += depth_step;
      }
    }
    xa += xa_step;
    xb += xb_step;
    deptha += deptha_step;
    depthb += depthb_step;
  }
}




static int
WriteUVDepthImages(GSVImage *image, const char *output_image_directory, 
  const char *prefix_name = "", RNBoolean fill_holes = FALSE)
{
  // Get convenient variables
  GSVCamera *camera = image->Camera();
  if (!camera) return 0;
  GSVPanorama *panorama = image->Panorama();
  if (!panorama) return 0;
  int image_index = image->PanoramaIndex();
  GSVSegment *segment = panorama->Segment();
  if (!segment) return 0;
  int panorama_index = panorama->SegmentIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;
  int segment_index = segment->RunIndex();
  char image_name[4096];

  // Get filenames 
  char depth_grid_name[4096], sa_index_grid_name[4096];
  sprintf(depth_grid_name, "%s/%s/%02d_%06d_%02d_UV_%sDepth.pfm", 
    output_image_directory, run->Name(), segment_index, panorama_index, image_index, prefix_name);
  sprintf(sa_index_grid_name, "%s/%s/%02d_%06d_%02d_UV_%sSAIndex.pfm", 
    output_image_directory, run->Name(), segment_index, panorama_index, image_index, prefix_name);

  // Skip if already done
  if (RNFileExists(depth_grid_name)) return 1;
  
  // Print message
  if (print_verbose) {
    printf("  Creating UV images for %s %02d %06d %02d\n", run->Name(), segment_index, panorama_index, image_index);
    fflush(stdout);
  }

  // Get undistorted image
  R2Image *undistorted_image = image->UndistortedImage();
  if (!undistorted_image) {
    fprintf(stderr, "Unable to get undistorted image\n");
    return 0;
  }

  // Get more convenient variables
  R2Box image_bbox(0, 0, undistorted_image->Width()-1, undistorted_image->Height()-1);
  const GSVPose& pose = image->Pose();
  R3Frustum frustum(pose.Viewpoint(), pose.Towards(), pose.Up(), 
    camera->XFov(), camera->YFov(), min_frustum_depth, max_frustum_depth);

  // Create grids
  R2Grid sa_index_grid(undistorted_image->Width(), undistorted_image->Height());
  sa_index_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid depth_grid(undistorted_image->Width(), undistorted_image->Height());
  depth_grid.Clear(max_frustum_depth);

  // Draw scans in same segment
  for (int ia = 0; ia < segment->NScans(); ia++) {
    GSVScan *scan = segment->Scan(ia);

    // Check if scan bounding box intersects image view frustum
    const R3Box& scan_bbox = scan->BBox();
    if (!frustum.Intersects(scan_bbox)) continue;

    // Read SA grids
    R2Grid sa_position_x_grid, sa_position_y_grid, sa_position_z_grid;
    sprintf(image_name, "%s/%s/%02d_%02d_SA_PositionX.grd", output_image_directory, run->Name(), segment_index, ia);
    if (!sa_position_x_grid.Read(image_name)) return 0;
    sprintf(image_name, "%s/%s/%02d_%02d_SA_PositionY.grd", output_image_directory, run->Name(), segment_index, ia);
    if (!sa_position_y_grid.Read(image_name)) return 0;
    sprintf(image_name, "%s/%s/%02d_%02d_SA_PositionZ.grd", output_image_directory, run->Name(), segment_index, ia);
    if (!sa_position_z_grid.Read(image_name)) return 0;

    // Fill holes
    if (fill_holes) {
      sa_position_x_grid.FillHoles();
      sa_position_y_grid.FillHoles();
      sa_position_z_grid.FillHoles();
    }

    // Rasterize SA grids into UV depth images
    for (int ix = 1; ix < sa_position_x_grid.XResolution(); ix++) {
      for (int iy = 1; iy < sa_position_x_grid.YResolution(); iy++) {
        // Rasterize 4 triangles per quad
        for (int it = 0; it < 4; it++) {
          // Get vertex data
          int nvertices = 0;
          int idx = 0; int idy = 0;
          R3Point world_position[3];
          R2Point undistorted_position[3];
          R2Point distorted_position[3];
          RNScalar depth[3] = { 0.0, 0.0, 0.0 };
          int grid_index[3] = { -1, -1, -1 };
          for (int iv = 0; iv < 3; iv++) {
            // Get SA index offsets based on triangle number
            if      ((it == 0) && (iv == 0)) { idx = -1; idy = -1; }
            else if ((it == 0) && (iv == 1)) { idx =  0; idy = -1; }
            else if ((it == 0) && (iv == 2)) { idx =  0; idy =  0; }
            else if ((it == 1) && (iv == 0)) { idx = -1; idy = -1; }
            else if ((it == 1) && (iv == 1)) { idx =  0; idy =  0; }
            else if ((it == 1) && (iv == 2)) { idx = -1; idy =  0; }
            else if ((it == 2) && (iv == 0)) { idx =  0; idy = -1; }
            else if ((it == 2) && (iv == 1)) { idx = -1; idy =  0; }
            else if ((it == 2) && (iv == 2)) { idx = -1; idy = -1; }
            else if ((it == 3) && (iv == 0)) { idx =  0; idy = -1; }
            else if ((it == 3) && (iv == 1)) { idx =  0; idy =  0; }
            else if ((it == 3) && (iv == 2)) { idx = -1; idy =  0; }

            // Get world position
            RNScalar x = sa_position_x_grid.GridValue(ix+idx, iy+idy);
            if (x == R2_GRID_UNKNOWN_VALUE) continue;
            RNScalar y = sa_position_y_grid.GridValue(ix+idx, iy+idy);
            if (y == R2_GRID_UNKNOWN_VALUE) continue;
            RNScalar z = sa_position_z_grid.GridValue(ix+idx, iy+idy);
            if (z == R2_GRID_UNKNOWN_VALUE) continue;
            world_position[nvertices].Reset(x, y, z);

            // Get undistorted image position
            undistorted_position[nvertices] = image->UndistortedPosition(world_position[nvertices]);
            if (undistorted_position[nvertices].X() == RN_UNKNOWN) continue;
            if (undistorted_position[nvertices].Y() == RN_UNKNOWN) continue;

            // Get depth
            distorted_position[nvertices] = image->DistortedPosition(undistorted_position[nvertices]);
            if (distorted_position[nvertices].X() == RN_UNKNOWN) continue;
            if (distorted_position[nvertices].Y() == RN_UNKNOWN) continue;
            int column_index = (int) (distorted_position[nvertices].X() + 0.5);
            GSVPose pose = image->Pose(column_index);
            const R3Point& viewpoint = pose.Viewpoint();
            const R3Vector& towards = pose.Towards();
            R3Vector vector = world_position[nvertices] - viewpoint;
            depth[nvertices] = vector.Dot(towards);
            if (depth[nvertices] <= 0) continue;

            // Get scanline index
            sa_position_x_grid.IndicesToIndex(ix+idx, iy+idy, grid_index[nvertices]);
            if (grid_index[nvertices] < 0) continue;
            if (ia > 0) grid_index[nvertices] += 180 * segment->Scan(0)->NScanlines();
            if (ia > 1) grid_index[nvertices] += 180 * segment->Scan(1)->NScanlines();

            // Got everything for this vertex
            nvertices++;
          }

          // Accommodate spans and points
          if (nvertices == 0) continue;
          if (nvertices < 2) {
            undistorted_position[1] = undistorted_position[0];
            depth[1] = depth[0];
            grid_index[1] = grid_index[0];
          }
          if (nvertices < 3) {
            undistorted_position[2] = undistorted_position[1];
            depth[2] = depth[1];
            grid_index[2] = grid_index[1];
          }

          // Rasterize
          RasterizeTriangle(depth_grid, sa_index_grid, undistorted_position, depth, grid_index);
        }
      }
    }
  }

  // Delete undistorted image
  delete undistorted_image;

  // Update depth grid
  depth_grid.Threshold(max_frustum_depth - RN_EPSILON, R2_GRID_KEEP_VALUE, R2_GRID_UNKNOWN_VALUE);

  // Write grids
  depth_grid.Write(depth_grid_name);
  sa_index_grid.Write(sa_index_grid_name);

  // Return success
  return 1;
}



static int
WriteUVImage(GSVImage *image, const char *output_image_directory, const char *grid_name)
{
  // Get convenient variables
  GSVCamera *camera = image->Camera();
  if (!camera) return 0;
  GSVPanorama *panorama = image->Panorama();
  if (!panorama) return 0;
  int image_index = image->PanoramaIndex();
  GSVSegment *segment = panorama->Segment();
  if (!segment) return 0;
  int panorama_index = panorama->SegmentIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;
  int segment_index = segment->RunIndex();
  char image_name[4096];

  // Get filenames 
  char output_grid_name[4096];
  sprintf(output_grid_name, "%s/%s/%02d_%06d_%02d_UV_%s.pfm", 
    output_image_directory, run->Name(), segment_index, panorama_index, image_index, grid_name);

  // Skip if already done
  if (RNFileExists(output_grid_name)) return 1;
  
  // Print message
  if (print_verbose) {
    printf("  Creating UV image for %s %02d %06d %02d %s\n", run->Name(), segment_index, panorama_index, image_index, grid_name);
    fflush(stdout);
  }

  // Read UV grids
  R2Grid sa_grids[3];
  R2Grid uv_index_grid;
  sprintf(image_name, "%s/%s/%02d_%06d_%02d_UV_GridIndex.pfm", 
    output_image_directory, run->Name(), segment_index, panorama_index, image_index);
  if (!uv_index_grid.Read(image_name)) return 0;
  for (int ia = 0; ia < segment->NScans(); ia++) {
    sprintf(image_name, "%s/%s/%02d_%02d_SA_%s.grd", 
      output_image_directory, run->Name(), segment_index, ia, grid_name);
    if (!sa_grids[ia].Read(image_name)) return 0;
  }

  // Create output grid
  R2Grid output_grid(uv_index_grid.XResolution(), uv_index_grid.YResolution());
  output_grid.Clear(R2_GRID_UNKNOWN_VALUE);

  // Write output grid
  for (int iy = 0; iy < uv_index_grid.YResolution(); iy++) {
    for (int ix = 0; ix < uv_index_grid.XResolution(); ix++) {
      RNScalar grid_index_value = uv_index_grid.GridValue(ix, iy);
      if (grid_index_value == R2_GRID_UNKNOWN_VALUE) continue;

      // Determine scan index and grid index
      int scan_index = 0;
      int grid_index = (int) (grid_index_value + 0.5);
      if (grid_index > 180*segment->Scan(0)->NScanlines()) {
        grid_index -= 180*segment->Scan(0)->NScanlines();
        scan_index++;
      }
      if (grid_index > 180*segment->Scan(1)->NScanlines()) {
        grid_index -= 180*segment->Scan(1)->NScanlines();
        scan_index++;
      }
      if (grid_index > 180*segment->Scan(2)->NScanlines()) {
        RNAbort("Invalid grid index\n");
        return 0;
      }

      // Get value
      RNScalar value = sa_grids[scan_index].GridValue(grid_index);
      if (value == R2_GRID_UNKNOWN_VALUE) continue;

      // Write value
      output_grid.SetGridValue(ix, iy, value);
    }
  }

  // Write grid
  output_grid.Write(output_grid_name);

  // Return success
  return 1;
}




static int
WriteUVImages(GSVImage *image, const char *output_image_directory)
{
  // Write depth images 
  if (include_depth_images) {
    if (!WriteUVDepthImages(image, output_image_directory, "", FALSE)) return 0;
  }

  // Write hull images
  if (include_hull_images) {
    if (!WriteUVDepthImages(image, output_image_directory, "Hull_", TRUE)) return 0;
  }

  // Temporary
  if (include_label_images) {
    if (!WriteUVImage(image, output_image_directory, "MajorStructureLabel")) return 0;
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Depth images
////////////////////////////////////////////////////////////////////////



// GLUT state variables

static int run_index = 0;
static int segment_index = 0;
static int panorama_index = 0;
static int image_index = 0;



// GLUT window variables 

static int GLUTwindow = 0;
static int GLUTwindow_height = 1000;
static int GLUTwindow_width = 1000;
static GSVScene *GLUTscene = NULL;




static int
WriteOpenGLImage(GSVScene *scene, GSVImage *image)
{
  // Get convenient variables
  GSVCamera *camera = image->Camera();
  if (!camera) return 0;
  GSVPanorama *panorama = image->Panorama();
  if (!panorama) return 0;
  int image_index = image->PanoramaIndex();
  GSVSegment *segment = panorama->Segment();
  if (!segment) return 0;
  int panorama_index = panorama->SegmentIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;
  int segment_index = segment->RunIndex();

  // Get filenames
  char depth_grid_name[4096], scan_grid_name[4096], scanline_grid_name[4096], point_index_grid_name[4096];
  sprintf(depth_grid_name, "%s/%s/%02d_%02d_%06d_UV_OpenGL_Depth.pfm", 
    output_image_directory, run->Name(), segment_index, panorama_index, image_index);
  sprintf(scan_grid_name, "%s/%s/%02d_%02d_%06d_UV_OpenGL_Scan.pfm", 
    output_image_directory, run->Name(), segment_index, panorama_index, image_index);
  sprintf(scanline_grid_name, "%s/%s/%02d_%02d_%06d_UV_OpenGL_Scanline.pfm", 
    output_image_directory, run->Name(), segment_index, panorama_index, image_index);
  sprintf(point_index_grid_name, "%s/%s/%02d_%02d_%06d_UV_OpenGL_PointIndex.pfm", 
    output_image_directory, run->Name(), segment_index, panorama_index, image_index);

  // Check if files already exist
  if (RNFileExists(depth_grid_name)) return 1;

  // Get distorted image
  R2Image *distorted_image = image->DistortedImage();
  if (!distorted_image) {
    fprintf(stderr, "Unable to get distorted image\n");
    return 0;
  }

  // Get more convenient variables
  R2Box image_bbox(0, 0, distorted_image->Width()-1, distorted_image->Height()-1);
  const GSVPose& pose = image->Pose();
  R3Frustum frustum(pose.Viewpoint(), pose.Towards(), pose.Up(), 
    camera->XFov(), camera->YFov(), min_frustum_depth, max_frustum_depth);

  // Create grids
  R2Grid depth_grid(distorted_image->Width(), distorted_image->Height());
  depth_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  R2Grid scan_grid(depth_grid);
  R2Grid scanline_grid(depth_grid);
  R2Grid point_index_grid(depth_grid);

  // Set viewport
  glViewport(0, 0, GLUTwindow_width, GLUTwindow_height);
  
  // Set projection matrix
  glMatrixMode(GL_PROJECTION);  
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, GLUTwindow_width-1, 0, GLUTwindow_height-1, 0, max_frustum_depth);
      
  // Set model view matrix
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  // Read back transformation matrices
  GLint viewport[16];
  GLdouble projection_matrix[16];
  GLdouble modelview_matrix[16];
  glGetIntegerv(GL_VIEWPORT, viewport);
  glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix);
  glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix);

  // Allocate buffer
  unsigned char *color_buffer = new unsigned char [ 4 * GLUTwindow_width * GLUTwindow_height ];
  GLfloat *depth_buffer = new GLfloat [ GLUTwindow_width * GLUTwindow_height ];

  // Draw the image in pieces that fit within the GLUT window
  for (int ox = 0; ox < distorted_image->Width(); ox += GLUTwindow_width) {
    for (int oy = 0; oy < distorted_image->Height(); oy += GLUTwindow_height) {
      // Clear window 
      glClearColor(0, 0, 0, 0);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      // Draw scans in same segment
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        GSVMesh *mesh = scan->Mesh();
        if (!mesh) continue;

        // Check if mesh bounding box intersects image view frustum
        const R3Box& mesh_bbox = mesh->BBox();
        if (!frustum.Intersects(mesh_bbox)) {
          delete mesh;
          continue;
        }
        
        // Draw vertices
        glBegin(GL_POINTS);
        for (int k = 0; k < mesh->NVertices(); k++) {
          GSVMeshVertex *vertex = mesh->Vertex(k);
          const R3Point& world_position = mesh->VertexPosition(vertex);
          if (!frustum.Intersects(world_position)) continue;
          R2Point distorted_position = image->DistortedPosition(world_position);
          if (!R2Contains(image_bbox, distorted_position)) continue;
          GSVScanline *scanline = mesh->VertexScanline(vertex);
          int scanline_index = scanline->SegmentIndex();
          int point_index = mesh->VertexPointIndex(vertex);
          int column_index = (int) (distorted_position.X() + 0.5);
          GSVPose pose = image->Pose(column_index);
          const R3Point& viewpoint = pose.Viewpoint();
          const R3Vector& towards = pose.Towards();
          R3Vector vector = world_position - viewpoint;
          RNScalar depth = vector.Dot(towards);
          unsigned char r = (point_index+1) & 0XFF;
          unsigned char g = (scanline_index >> 16) & 0XFF;
          unsigned char b = (scanline_index >> 8) & 0XFF;
          unsigned char a = (scanline_index) & 0XFF;
          glColor4ub(r, g, b, a);
          glVertex3d(distorted_position[0] - ox, distorted_position[1] - oy, -depth);
        }
        glEnd();

        // Draw edges
        for (int k = 0; k < mesh->NEdges(); k++) {
          GSVMeshEdge *edge = mesh->Edge(k);
          R3Span edge_span = mesh->EdgeSpan(edge);
          if (!frustum.Intersects(edge_span)) continue;
          glBegin(GL_LINES);
          for (int iv = 0; iv < 2; iv++) {
            GSVMeshVertex *vertex = mesh->VertexOnEdge(edge, iv);
            const R3Point& world_position = mesh->VertexPosition(vertex);
            R2Point distorted_position = image->DistortedPosition(world_position);
            if (!R2Contains(image_bbox, distorted_position)) continue;
            GSVScanline *scanline = mesh->VertexScanline(vertex);
            int scanline_index = scanline->SegmentIndex();
            int point_index = mesh->VertexPointIndex(vertex);
            int column_index = (int) (distorted_position.X() + 0.5);
            GSVPose pose = image->Pose(column_index);
            const R3Point& viewpoint = pose.Viewpoint();
            const R3Vector& towards = pose.Towards();
            R3Vector vector = world_position - viewpoint;
            RNScalar depth = vector.Dot(towards);
            unsigned char r = (point_index+1) & 0XFF;
            unsigned char g = (scanline_index >> 16) & 0XFF;
            unsigned char b = (scanline_index >> 8) & 0XFF;
            unsigned char a = (scanline_index) & 0XFF;
            glColor4ub(r, g, b, a);
            glVertex3d(distorted_position[0] - ox, distorted_position[1] - oy, -depth);
          }
          glEnd();
        }

        // Draw faces
        for (int k = 0; k < mesh->NFaces(); k++) {
          GSVMeshFace *face = mesh->Face(k);
          const R3Box& face_bbox = mesh->FaceBBox(face);
          if (!frustum.Intersects(face_bbox)) continue;
          glBegin(GL_POLYGON);
          for (int iv = 0; iv < 3; iv++) {
            GSVMeshVertex *vertex = mesh->VertexOnFace(face, iv);
            const R3Point& world_position = mesh->VertexPosition(vertex);
            R2Point distorted_position = image->DistortedPosition(world_position);
            if (!R2Contains(image_bbox, distorted_position)) continue;
            GSVScanline *scanline = mesh->VertexScanline(vertex);
            int segment_scanline_index = scanline->SegmentIndex();
            int point_index = mesh->VertexPointIndex(vertex);
            int column_index = (int) (distorted_position.X() + 0.5);
            GSVPose pose = image->Pose(column_index);
            const R3Point& viewpoint = pose.Viewpoint();
            const R3Vector& towards = pose.Towards();
            R3Vector vector = world_position - viewpoint;
            RNScalar depth = vector.Dot(towards);
            unsigned char r = (point_index+1) & 0XFF;
            unsigned char g = (segment_scanline_index >> 16) & 0XFF;
            unsigned char b = (segment_scanline_index >> 8) & 0XFF;
            unsigned char a = (segment_scanline_index) & 0XFF;
            glColor4ub(r, g, b, a);
            glVertex3d(distorted_position[0] - ox, distorted_position[1] - oy, -depth);
          }
          glEnd();
        }

        // Delete mesh
        delete mesh;
      }

      // Read color and depth buffers
      glReadPixels(0, 0, GLUTwindow_width, GLUTwindow_height, GL_RGBA, GL_UNSIGNED_BYTE, color_buffer);
      glReadPixels(0, 0, GLUTwindow_width, GLUTwindow_height, GL_DEPTH_COMPONENT, GL_FLOAT, depth_buffer);

      // Fill depth values into grid
      GLfloat *depthp = depth_buffer;
      unsigned char *colorp = color_buffer;
      for (int iy = 0; iy < GLUTwindow_height; iy++) {
        for (int ix = 0; ix < GLUTwindow_width; ix++) {
          double x, y, z;
          unsigned char r = *(colorp++);
          unsigned char g = *(colorp++);
          unsigned char b = *(colorp++);
          unsigned char a = *(colorp++);
          GLfloat d = *(depthp++);
          if (r == 0) continue;
          if (ix + ox >= depth_grid.XResolution()) continue;
          if (iy + oy >= depth_grid.YResolution()) continue;
          int segment_scanline_index = (g<<16) + (b<<8) + a;
          int point_index = r-1;
          GSVScanline *scanline = segment->Scanline(segment_scanline_index);
          if (!scanline) continue;
          GSVScan *scan = scanline->Scan();
          if (!scan) continue;
          gluUnProject(ix, iy, d, modelview_matrix, projection_matrix, viewport, &x, &y, &z);
          scan_grid.SetGridValue(ix + ox, iy + oy, scan->SegmentIndex());
          scanline_grid.SetGridValue(ix + ox, iy + oy, scanline->ScanIndex());
          point_index_grid.SetGridValue(ix + ox, iy + oy, point_index);
          depth_grid.SetGridValue(ix + ox, iy + oy, -z);
        }
      }
    }
  }

  // Delete distorted image
  delete distorted_image;

  // Delete buffers
  delete [] color_buffer;
  delete [] depth_buffer;

  // Reset projection matrix
  glMatrixMode(GL_PROJECTION);  
  glPopMatrix();

  // Reset modelview matrix
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  // Write grids
  depth_grid.Write(depth_grid_name);
  scanline_grid.Write(scanline_grid_name);
  scan_grid.Write(scan_grid_name);
  point_index_grid.Write(point_index_grid_name);

  // Print message
  if (print_verbose) {
    printf("  Created opengl images for %s %02d %02d %d\n", run->Name(), segment_index, panorama_index, image_index);
    fflush(stdout);
  }

  // Return success
  return 1;
}

void GLUTStop(void)
{
  // Destroy window 
  glutDestroyWindow(GLUTwindow);

  // Exit
  exit(0);
}



void GLUTKeyboard(unsigned char key, int x, int y)
{
  // Process quit (ESC) command
  if (key == 27) GLUTStop();
}



void GLUTResize(int w, int h)
{
  // Disallow resize
  glutReshapeWindow(GLUTwindow_width, GLUTwindow_height);
}



void GLUTIdle(void)
{
  // Redraw
  glutPostRedisplay();  
}



void GLUTRedraw(void)
{
  // Check scene
  if (!GLUTscene) return;

  // Get image
  if (run_index >= GLUTscene->NRuns()) GLUTStop();
  GSVRun *run = GLUTscene->Run(run_index);
  if (segment_index >= run->NSegments()) GLUTStop();
  GSVSegment *segment = run->Segment(segment_index);
  if (panorama_index >= segment->NPanoramas()) GLUTStop();
  GSVPanorama *panorama = segment->Panorama(panorama_index);
  if (image_index >= panorama->NImages()) GLUTStop();
  GSVImage *image = panorama->Image(image_index);

  // Write depth image
  if (!WriteOpenGLImage(GLUTscene, image)) GLUTStop();

  // Update indices
  if (++image_index >= panorama->NImages()) { panorama_index++; image_index = 0; }
  if (panorama_index >= segment->NPanoramas()) { segment_index++; panorama_index = 0; }
  if (segment_index >= run->NSegments()) { run_index++; segment_index = 0; }

  // Swap buffers 
  glutSwapBuffers();
}    



static void 
WriteOpenGLImages(GSVScene *scene, const char *)
{
  // Print message
  if (print_verbose) {
    printf("Creating depth images ...\n");
    fflush(stdout);
  }

  // Set scene
  GLUTscene = scene;

  // Open window 
  int argc = 1;
  char program_name[1024];
  strcpy(program_name, "gsv2img");
  char **argv = new char * [ 1 ];
  argv[0] = program_name;
  glutInit(&argc, argv);
  glutInitWindowPosition(10, 10);
  glutInitWindowSize(GLUTwindow_width, GLUTwindow_height);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  GLUTwindow = glutCreateWindow("gsv2img");

  // Initialize background color
  glClearColor(0, 0, 0, 1);

  // Initialize graphics modes  
  glEnable(GL_DEPTH_TEST);

  // Initialize GLUT callback functions 
  glutDisplayFunc(GLUTRedraw);
  glutReshapeFunc(GLUTResize);
  glutKeyboardFunc(GLUTKeyboard);
  glutIdleFunc(GLUTIdle);

  // Run main loop -- never returns
  glutMainLoop();
}




////////////////////////////////////////////////////////////////////////
// Top-level image computation function
////////////////////////////////////////////////////////////////////////

static int
WriteImages(GSVScene *scene, const char *output_image_directory)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_verbose) {
    printf("Creating images ...\n");
    fflush(stdout);
  }

  // Write scan images
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (capture_SA_images && !WriteSAImages(scan, output_image_directory)) return 0;
        if (capture_DA_images && !WriteDAImages(scan, output_image_directory)) return 0;
        if (capture_DH_images && !WriteDHImages(scan, output_image_directory)) return 0;
        if (capture_XY_images && !WriteXYImages(scan, output_image_directory)) return 0;
      }
      if (capture_XY_images && !WriteXYImages(segment, output_image_directory)) return 0;
    }
    if (capture_XY_images && !WriteXYImages(run, output_image_directory)) return 0;
  }
  if (capture_XY_images && !WriteXYImages(scene, output_image_directory)) return 0;

  // Write image images
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ip = 0; ip < segment->NPanoramas(); ip++) {
        GSVPanorama *panorama = segment->Panorama(ip);
        for (int ii = 0; ii < panorama->NImages(); ii++) {
          GSVImage *image = panorama->Image(ii);
          if (capture_UV_images && !WriteUVImages(image, output_image_directory)) return 0;
        }
      }
    }
  }

  // Write opengl images
  if (capture_opengl_images) {
    // This does not return ???
    WriteOpenGLImages(scene, output_image_directory);
  }

  // Print statistics
  if (print_verbose) {
    printf("  Done in %.2f seconds\n", start_time.Elapsed());
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
  // By default ...
  int default_parameterizations = 1;
  int default_images = 1;
  int capture_everything = 0;

  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) { print_verbose = 1; }
      else if (!strcmp(*argv, "-debug")) { print_debug = 1; }
      else if (!strcmp(*argv, "-XY")) { capture_XY_images = 1; default_parameterizations = 0; }
      else if (!strcmp(*argv, "-SA")) { capture_SA_images = 1; default_parameterizations = 0; }
      else if (!strcmp(*argv, "-DA")) { capture_DA_images = 1; default_parameterizations = 0; }
      else if (!strcmp(*argv, "-DH")) { capture_DH_images = 1; default_parameterizations = 0; }
      else if (!strcmp(*argv, "-UV")) { capture_UV_images = 1; default_parameterizations = 0; }
      else if (!strcmp(*argv, "-opengl")) { capture_opengl_images = 1; default_parameterizations = 0; }
      else if (!strcmp(*argv, "-base")) { include_base_images = 1; default_images = 0; }
      else if (!strcmp(*argv, "-normal")) { include_normal_images = 1; default_images = 0; }
      else if (!strcmp(*argv, "-curvature")) { include_curvature_images = 1; default_images = 0; }
      else if (!strcmp(*argv, "-boundary")) { include_boundary_images = 1; default_images = 0; }
      else if (!strcmp(*argv, "-pca")) { include_pca_images = 1; default_images = 0; }
      else if (!strcmp(*argv, "-color")) { include_color_images = 1; default_images = 0; }
      else if (!strcmp(*argv, "-depth")) { include_depth_images = 1; default_images = 0; }
      else if (!strcmp(*argv, "-hull")) { include_hull_images = 1; default_images = 0; }
      else if (!strcmp(*argv, "-label")) { include_label_images = 1; default_images = 0; }
      else if (!strcmp(*argv, "-all")) { capture_everything = 1; default_parameterizations = default_images = 0; }
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
    else {
      if (!input_scene_name) input_scene_name = *argv;
      else if (!output_image_directory) output_image_directory = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check scene name
  if (!input_scene_name) {
    fprintf(stderr, "Usage: gsv2img input_scene [output_directory] [options]\n");
    return FALSE;
  }

  // Check output image directory
  if (!output_image_directory) {
    output_image_directory = "gsv_data/laser_images";
  }

  // Check image selection
  if (capture_everything) {
    capture_XY_images = 1;
    capture_SA_images = 1;
    capture_DA_images = 1;
    capture_DH_images = 1;
    capture_UV_images = 1;
    // capture_opengl_images = 1;
    include_base_images = 1;
    include_normal_images = 1;
    include_curvature_images = 1;
    include_boundary_images = 1;
    include_pca_images = 1;
    include_color_images = 1;
    include_depth_images = 1;
    include_hull_images = 1;
  }
  if (capture_UV_images) {
    capture_SA_images = 1;
    include_base_images = 1;
  }
  if (default_parameterizations) {
    capture_XY_images = 1;
    capture_SA_images = 1;
    capture_DA_images = 1;
    capture_DH_images = 1;
  }
  if (default_images) {
    include_base_images = 1;
    include_normal_images = 1;
    include_curvature_images = 1;
    include_boundary_images = 1;
    include_pca_images = 1;
  }

  // Check dependencies
  if (include_pca_images) {
    include_base_images = 1;
  }
  if (include_normal_images) {
    include_base_images = 1;
  }
  if (include_curvature_images) {
    include_base_images = 1;
    include_normal_images = 1;
  }
  if (include_boundary_images) {
    include_base_images = 1;
  }
  if (include_color_images) {
    include_base_images = 1;
  }

  // Check verbosity 
  if (print_debug) {
    print_verbose = 1;
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

  // Create output image directories
  char cmd[1024];
  sprintf(cmd, "mkdir -p %s", output_image_directory);  
  system(cmd);
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    sprintf(cmd, "mkdir -p %s/%s", output_image_directory, run->Name());
    system(cmd);
  }
    
  // Write images
  if (!WriteImages(scene, output_image_directory)) exit(-1);

  // Return success 
  return 0;
}

















