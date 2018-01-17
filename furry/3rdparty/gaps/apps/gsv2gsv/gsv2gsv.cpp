// Source file for the google street view processing program



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV/GSV.h"



////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////

// Program arguments

static char *input_scene_name = NULL;
static char *output_scene_name = NULL;
static char *selected_run_name = NULL;
static char *commands_name = NULL;
static int fill_cache = 0;
static int remove_points = 0;
static int print_verbose = 0;
static int print_debug = 0;



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

  // Read google scene files
  if (!scene->ReadFile(filename, !remove_points)) {
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
WriteScene(GSVScene *scene, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Write google scene files
  if (!scene->WriteFile(filename)) return 0;

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



////////////////////////////////////////////////////////////////////////
// Processing Functions
////////////////////////////////////////////////////////////////////////

static GSVScene *
CopyScene(GSVScene *input_scene, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int command_count = 0;

  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }

  // Create output scene
  GSVScene *output_scene = new GSVScene();

  // Set output raw directory name
  char buffer[2048];
  strcpy(buffer, output_scene_name);
  char *bufferp = strrchr(buffer, '/');
  if (bufferp) *bufferp = '\0';
  else strcpy(buffer, ".");
  strcat(buffer, "/raw_data");
  output_scene->SetRawDataDirectoryName(buffer);

  // Set output cache directory name
  strcpy(buffer, output_scene_name);
  bufferp = strrchr(buffer, '/');
  if (bufferp) *bufferp = '\0';
  else strcpy(buffer, ".");
  strcat(buffer, "/gsv_data");
  output_scene->SetCacheDataDirectoryName(buffer);

  // Check directory names
  if ((!strcmp(input_scene->RawDataDirectoryName(), output_scene->RawDataDirectoryName())) ||
      (!strcmp(input_scene->CacheDataDirectoryName(), output_scene->CacheDataDirectoryName()))) {
    fprintf(stderr, "Unable to execute commands for input and output scenes in same directory\n");
    return NULL;
  }

  // Copy elements of input scene to output scene
  char keyword[1024];
  while (fscanf(fp, "%s", keyword) == (unsigned int) 1) {
    command_count++;
    if (!strcmp(keyword, "copy_segment")) {
      // Parse subsegment stuff
      char input_run_name[1024];
      int input_segment_index;
      RNScalar begin_timestamp, end_timestamp;  
      RNScalar x1, y1, x2, y2;
      if (fscanf(fp, "%s%d%lf%lf%lf%lf%lf%lf", 
                 input_run_name, &input_segment_index, 
                 &begin_timestamp, &end_timestamp, 
                 &x1, &y1, &x2, &y2) != (unsigned int) 8) {
        fprintf(stderr, "Error parsing command %d in %s\n", command_count, filename);
        return NULL;
      }

      // Find input run 
      GSVRun *input_run = input_scene->Run(input_run_name);
      if (!input_run) {
        fprintf(stderr, "Unable to find run with name matching %s\n", input_run_name);
        return 0;
      }
    
      // Check input segment index
      if ((input_segment_index < 0) || (input_segment_index >= input_run->NSegments())) {
        fprintf(stderr, "Unable to find segment %d in run %s\n", input_segment_index, input_run_name);
        return 0;
      }

      // Get input segment
      GSVSegment *input_segment = input_run->Segment(input_segment_index);
      if (!input_segment) {
        fprintf(stderr, "Unable to find segment %d in run %s\n", input_segment_index, input_run_name);
        return 0;
      }

      // Check timestamp
      if (begin_timestamp == -1) begin_timestamp = -FLT_MAX;
      if (end_timestamp == -1) end_timestamp = FLT_MAX;

    // Create/check bounding box
      R3Box bbox(-FLT_MAX, -FLT_MAX, -FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX);
      if (x1 != -1) bbox[0][0] = x1;
      if (y1 != -1) bbox[0][1] = y1;
      if (x2 != -1) bbox[1][0] = x2;
      if (y2 != -1) bbox[1][1] = y2;

      // Get/create output run
      GSVRun *output_run = output_scene->Run(input_run->Name());
      if (!output_run) {
        // Create run
        output_run = new GSVRun(input_run->Name());
        output_scene->InsertRun(output_run);

        // Create cameras
        for (int ic = 0; ic < input_run->NCameras(); ic++) {
          GSVCamera *input_camera = input_run->Camera(ic);
          GSVCamera *output_camera = new GSVCamera();
          output_camera->SetDistortionType(input_camera->DistortionType());
          output_camera->SetXFocal(input_camera->XFocal());
          output_camera->SetYFocal(input_camera->YFocal());
          output_camera->SetXCenter(input_camera->XCenter());
          output_camera->SetYCenter(input_camera->YCenter());
          output_camera->SetK1(input_camera->K1());
          output_camera->SetK2(input_camera->K2());
          output_camera->SetK3(input_camera->K3());
          output_camera->SetP1(input_camera->P1());
          output_camera->SetP2(input_camera->P2());
          output_camera->SetMaxFov(input_camera->MaxFov());
          output_run->InsertCamera(output_camera);
        }

        // Create lasers
        for (int ic = 0; ic < input_run->NLasers(); ic++) {
          GSVLaser *output_laser = new GSVLaser();
          output_run->InsertLaser(output_laser);
        }
      }

      // Create output segment
      GSVSegment *output_segment = new GSVSegment();
      output_run->InsertSegment(output_segment);

      // Create tapestries
      for (int it = 0; it < input_segment->NTapestries(); it++) {
        GSVTapestry *output_tapestry = new GSVTapestry();
        GSVCamera *output_camera = output_run->Camera(it);
        output_camera->InsertTapestry(output_tapestry);
        output_segment->InsertTapestry(output_tapestry);
      }      

      // Get input image directory
      char input_image_directory[4096];
      const char *input_raw_directory = input_run->Scene()->RawDataDirectoryName();
      sprintf(input_image_directory, "%s/%s/segment_%02d", input_raw_directory, input_run->Name(), input_segment->RunIndex());

      // Create output image directory
      char cmd[4096], output_image_directory[4096];
      const char *output_raw_directory = output_run->Scene()->RawDataDirectoryName();
      sprintf(output_image_directory, "%s/%s/segment_%02d", output_raw_directory, output_run->Name(), output_segment->RunIndex());
      sprintf(cmd, "mkdir -p %s", output_image_directory);
      system(cmd);

      // Copy images
      for (int ip = 0; ip < input_segment->NPanoramas(); ip++) {
        GSVPanorama *input_panorama = input_segment->Panorama(ip);

        // Check timestamp
        RNScalar timestamp = input_panorama->Timestamp();
        if (timestamp < begin_timestamp) continue;
        if (timestamp > end_timestamp) continue;

        // Check viewpoint
        const R3Point& viewpoint = input_panorama->Viewpoint();
        if (!bbox.IsEmpty() && !R3Contains(bbox, viewpoint)) continue;

        // Copy panorama
        GSVPanorama *output_panorama = new GSVPanorama();
        output_panorama->SetViewpoint(viewpoint);
        output_panorama->SetTimestamp(timestamp);
        output_segment->InsertPanorama(output_panorama);

        // Copy images
        for (int ii = 0; ii < input_panorama->NImages(); ii++) {
          GSVImage *input_image = input_panorama->Image(ii);
          GSVImage *output_image = new GSVImage();

          // Insert image into panorama
          output_panorama->InsertImage(output_image);
 
          // Insert image into tapestry
          GSVTapestry *output_tapestry = output_segment->Tapestry(ii);
          output_tapestry->InsertImage(output_image);

          // Copy image variables
          RNScalar timestamp0 = input_image->Timestamp(0);
          RNScalar timestamp1 = input_image->Timestamp(input_image->Width()-1);
          const GSVPose& pose0 = input_image->Pose(0);
          const GSVPose& pose1 = input_image->Pose(input_image->Width()-1);
          output_image->SetTimestamp(timestamp0, timestamp1);
          output_image->SetPose(pose0, pose1);

          // Copy image file
          GSVRun *input_run = input_segment->Run();
          char input_image_directory[4096], input_image_name[4096];
          char output_image_directory[4096], output_image_name[4096];
          sprintf(input_image_directory, "%s/%s/segment_%02d", input_raw_directory, input_run->Name(), input_segment->RunIndex());
          sprintf(output_image_directory, "%s/%s/segment_%02d", output_raw_directory, output_run->Name(), output_segment->RunIndex());
          sprintf(input_image_name,  "%s/unstitched_%06d_%02d.jpg", input_image_directory, input_panorama->RunIndex(), ii);
          sprintf(output_image_name, "%s/unstitched_%06d_%02d.jpg", output_image_directory, output_panorama->RunIndex(), ii);
          sprintf(cmd, "cp %s %s", input_image_name, output_image_name);
          system(cmd);
        }
      }

      // Copy scans
      for (int ia = 0; ia < input_segment->NScans(); ia++) {
        GSVScan *input_scan = input_segment->Scan(ia);
    
        // Create scan
        GSVScan *output_scan = new GSVScan();

        // Insert scan into segment
        output_segment->InsertScan(output_scan);

        // Insert scan into laser
        GSVLaser *output_laser = output_run->Laser(ia);
        output_laser->InsertScan(output_scan);

        // Read points
        input_scan->ReadPoints();
    
        // Copy scanlines
        for (int ie = 0; ie < input_scan->NScanlines(); ie++) {
          GSVScanline *input_scanline = input_scan->Scanline(ie);
          RNScalar timestamp = input_scanline->Timestamp();
          if (timestamp < begin_timestamp) continue;
          if (timestamp > end_timestamp) continue;
          const GSVPose& pose = input_scanline->Pose();
          const R3Point& viewpoint = pose.Viewpoint();
          if (!bbox.IsEmpty() && !R3Contains(bbox, viewpoint)) continue;
          int npoints = input_scanline->NPoints();
          R3Point *points = new R3Point [ npoints ];
          for (int k = 0; k < npoints; k++) points[k] = input_scanline->PointPosition(k);
          GSVScanline *output_scanline = new GSVScanline();
          output_scanline->SetPose(pose);
          output_scanline->SetTimestamp(timestamp);
          output_scanline->SetPointPositions(points, npoints);
          output_scan->InsertScanline(output_scanline);
          delete [] points;
        }
    
        // Release points
        input_scan->ReleasePoints();
      }
    }
  }

  // Close the file 
  fclose(fp);

  // Print statistics
  if (print_verbose) {
    printf("Copied scene ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Commands = %d\n", command_count);
    printf("  # Runs = %d\n", output_scene->NRuns());
    printf("  # Cameras = %d\n", output_scene->NCameras());
    printf("  # Lasers = %d\n", output_scene->NLasers());
    printf("  # Segments = %d\n", output_scene->NSegments());
    printf("  # Scans = %d\n", output_scene->NScans());
    printf("  # Tapestries = %d\n", output_scene->NTapestries());
    printf("  # Panoramas = %d\n", output_scene->NPanoramas());
    printf("  # Images = %d\n", output_scene->NImages());
    printf("  # Scanlines = %d\n", output_scene->NScanlines());
    fflush(stdout);
  }

  // Return output scene
  return output_scene;
}



static int 
SelectRun(GSVScene *scene, const char *selected_run_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Make array of runs (so don't change runs while traversing list)
  RNArray<GSVRun *> runs;
  for (int i = 0; i < scene->NRuns(); i++) {
    GSVRun *run = scene->Run(i);
    runs.Insert(run);
  }

  // Find run matching name
  GSVRun *selected_run = NULL;
  for (int i = 0; i < runs.NEntries(); i++) {
    GSVRun *run = runs.Kth(i);
    if (!strcmp(run->Name(), selected_run_name)) {
      selected_run = run;
      break;
    }
  }

  // Check whether any run was found matching name
  if (!selected_run) {
    fprintf(stderr, "Unable to find run with name matching %s\n", selected_run_name);
    return 0;
  }
    
  // Remove runs not matching name
  for (int i = 0; i < runs.NEntries(); i++) {
    GSVRun *run = runs.Kth(i);
    if (run != selected_run) {
      scene->RemoveRun(run);
      delete run;
    }
  }

  // Print statistics
  if (print_verbose) {
    printf("Selected run ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int 
FillCache(GSVScene *scene)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Access scene elements represented by files in gsv_data cache
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int it = 0; it < segment->NTapestries(); it++) {
        GSVTapestry *tapestry = segment->Tapestry(it);
        for (int ii = 0; ii < tapestry->NImages(); ii++) {
          GSVImage *image = tapestry->Image(ii);
          R2Image *undistorted_image = image->UndistortedImage();
          delete undistorted_image;
        }
      }
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        GSVMesh *mesh = scan->Mesh();
        delete mesh;
      }
    }
  }

  // Print statistics
  if (print_verbose) {
    printf("Filled cache ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int 
RemovePoints(GSVScene *scene)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int count = 0;

  // Remove points
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        for (int ie = 0; ie < scan->NScanlines(); ie++) {
          GSVScanline *scanline = scan->Scanline(ie);
          count += scanline->NPoints();
          scanline->SetPointPositions(NULL, 0);
        }
      }
    }
  }

  // Print statistics
  if (print_verbose) {
    printf("Removed points ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Points = %d\n", count);
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
      else if (!strcmp(*argv, "-fill_cache")) fill_cache = 1; 
      else if (!strcmp(*argv, "-remove_points")) remove_points = 1; 
      else if (!strcmp(*argv, "-select_run")) { argc--; argv++; selected_run_name = *argv; } 
      else if (!strcmp(*argv, "-commands")) { argc--; argv++; commands_name = *argv; } 
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
    else {
      if (!input_scene_name) input_scene_name = *argv;
      else if (!output_scene_name) output_scene_name = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check scene names
  if (!input_scene_name || !output_scene_name) {
    fprintf(stderr, "Usage: gsv2gsv inputscenefile outputscenefile [options]\n");
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

  // Select run
  if (selected_run_name) {
    if (!SelectRun(scene, selected_run_name)) exit(-1);
  }

  // Create new scene with commands
  if (commands_name) {
    GSVScene *output_scene = CopyScene(scene, commands_name);
    if (output_scene) { delete scene; scene = output_scene; }
    else exit(-1);
  }

  // Remove points
  if (remove_points) RemovePoints(scene);

  // Fill cache
  if (fill_cache) FillCache(scene);

  // Write scene
  if (!WriteScene(scene, output_scene_name)) exit(-1);

  // Return success 
  return 0;
}

















