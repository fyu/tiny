// Source file for GSV scene class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV.h"



////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

GSVScene::
GSVScene(void)
  : filename(NULL),
    raw_data_directory_name(strdup("raw_data")),
    cache_data_directory_name(strdup("gsv_data")),
    bbox(FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX),
    data(NULL)
{
}



GSVScene::
~GSVScene(void)
{
  // Delete runs
  while (NRuns() > 0) delete Run(0);

  // Delete filename
  if (filename) free(filename);

  // Delete raw data directory name
  if (raw_data_directory_name) free(raw_data_directory_name);

  // Delete cache data_directory name
  if (cache_data_directory_name) free(cache_data_directory_name);
}


////////////////////////////////////////////////////////////////////////
// Access functions
////////////////////////////////////////////////////////////////////////

GSVRun *GSVScene::
Run(const char *name) const
{
  // Return run with matching name
  for (int ir = 0; ir < NRuns(); ir++) {
    GSVRun *run = Run(ir);
    if (!strcmp(run->Name(), name)) return run;
  }

  // None found
  return NULL;
}



int GSVScene::
NSegments(void) const
{
  // Return number of segments in scene
  int count = 0;
  for (int ir = 0; ir < NRuns(); ir++) {
    GSVRun *run = Run(ir);
    count += run->NSegments();
  }
  return count;
}



int GSVScene::
NLasers(void) const
{
  // Return number of lasers in scene
  int count = 0;
  for (int ir = 0; ir < NRuns(); ir++) {
    GSVRun *run = Run(ir);
    count += run->NLasers();
  }
  return count;
}



int GSVScene::
NCameras(void) const
{
  // Return number of cameras in scene
  int count = 0;
  for (int ir = 0; ir < NRuns(); ir++) {
    GSVRun *run = Run(ir);
    count += run->NCameras();
  }
  return count;
}



int GSVScene::
NTapestries(void) const
{
  // Return number of tapestries in scene
  int count = 0;
  for (int i = 0; i < NRuns(); i++) {
    GSVRun *run = Run(i);
    count += run->NTapestries();
  }
  return count;
}



int GSVScene::
NPanoramas(void) const
{
  // Return number of panoramas in scene
  int count = 0;
  for (int i = 0; i < NRuns(); i++) {
    GSVRun *run = Run(i);
    count += run->NPanoramas();
  }
  return count;
}



int GSVScene::
NImages(void) const
{
  // Return number of images in scene
  int count = 0;
  for (int i = 0; i < NRuns(); i++) {
    GSVRun *run = Run(i);
    count += run->NImages();
  }
  return count;
}



int GSVScene::
NScans(void) const
{
  // Return number of scans in scene
  int count = 0;
  for (int i = 0; i < NRuns(); i++) {
    GSVRun *run = Run(i);
    count += run->NScans();
  }
  return count;
}



int GSVScene::
NScanlines(void) const
{
  // Return number of scanlines in scene
  int count = 0;
  for (int i = 0; i < NRuns(); i++) {
    GSVRun *run = Run(i);
    count += run->NScanlines();
  }
  return count;
}



////////////////////////////////////////////////////////////////////////
// Manipulation functions
////////////////////////////////////////////////////////////////////////

void GSVScene::
InsertRun(GSVRun *run)
{
  // Just checking
  assert(run->scene_index == -1);
  assert(run->scene == NULL);

  // Insert run
  run->scene = this;
  run->scene_index = runs.NEntries();
  runs.Insert(run);
}



void GSVScene::
RemoveRun(GSVRun *run)
{
  // Just checking
  assert(run->scene_index >= 0);
  assert(run->scene_index < runs.NEntries());
  assert(run->scene == this);

  // Remove run
  RNArrayEntry *entry = runs.KthEntry(run->scene_index);
  GSVRun *tail = runs.Tail();
  tail->scene_index = run->scene_index;
  runs.EntryContents(entry) = tail;
  runs.RemoveTail();
  run->scene_index = -1;
  run->scene = NULL;
}



////////////////////////////////////////////////////////////////////////
// Display functions
////////////////////////////////////////////////////////////////////////

void GSVScene::
Draw(RNFlags flags) const
{
  // Draw runs
  for (int i = 0; i < NRuns(); i++) {
    GSVRun *run = Run(i);
    run->Draw(flags);
  }
}



void GSVScene::
Print(FILE *fp, const char *prefix, const char *suffix) const
{
  // Check fp
  if (!fp) fp = stdout;

  // Print scene header
  if (prefix) fprintf(fp, "%s", prefix);
  fprintf(fp, "%s", "GSVScene");
  if (suffix) fprintf(fp, "%s", suffix);
  fprintf(fp, "\n");

  // Add indentation to prefix
  char indented_prefix[1024];
  sprintf(indented_prefix, "%s  ", prefix);

   // Print runs
  for (int i = 0; i < NRuns(); i++) {
    GSVRun *run = Run(i);
    run->Print(fp, indented_prefix, suffix);
  }
}



////////////////////////////////////////////////////////////////////////
// Data directory functions
////////////////////////////////////////////////////////////////////////

void GSVScene::
SetRawDataDirectoryName(const char *directory_name)
{
  // Set name of directory to store rawd data
  if (this->raw_data_directory_name) free(this->raw_data_directory_name);
  this->raw_data_directory_name = strdup(directory_name);
}


void GSVScene::
SetCacheDataDirectoryName(const char *directory_name)
{
  // Set name of directory to store cached data
  if (this->cache_data_directory_name) free(this->cache_data_directory_name);
  this->cache_data_directory_name = strdup(directory_name);
}



int GSVScene::
CreateCacheDataDirectory(const char *subdirectory, GSVRun *run) const
{
  // Create the cache directory
  if (subdirectory && (strlen(subdirectory) > 0)) {
    char cmd[4096];
    if (run->Name()) sprintf(cmd, "mkdir -p %s/%s/%s", cache_data_directory_name, subdirectory, run->Name());
    else sprintf(cmd, "mkdir -p %s/%s/Run%d", cache_data_directory_name, subdirectory, run->SceneIndex());
    system(cmd);
  }
  else if (run) {
    char cmd[4096];
    if (run->Name()) sprintf(cmd, "mkdir -p %s/run_data/%s", cache_data_directory_name, run->Name());
    else sprintf(cmd, "mkdir -p %s/run_data/Run%d", cache_data_directory_name, run->SceneIndex());
    system(cmd);
  }
  else {
    char cmd[4096];
    sprintf(cmd, "mkdir -p %s", cache_data_directory_name);
    system(cmd);
  }

  // Return success
  return 1;
}


 
////////////////////////////////////////////////////////////////////////
// Update functions
////////////////////////////////////////////////////////////////////////

void GSVScene::
UpdateBBox(void) 
{
  // Update bounding box
  bbox = R3null_box;
  for (int i = 0; i < NRuns(); i++) {
    bbox.Union(Run(i)->BBox());
  }
}



void GSVScene::
InvalidateBBox(void) 
{
  // Invalidate bounding box
  if (bbox.XMin() != FLT_MAX) {
    bbox = R3Box(FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX);
  }
}



#if 0

void GSVScene::
Update(void) 
{
  // Update lists
  for (int ir = 0; ir < NRuns(); ir++) {
    GSVRun *run = Run(i);
    for (int ic = 0; ic < run->NCameras(); ic++) {
      GSVCamera *camera = run->Camera(ic);
      camera->scene_index = cameras.NEntries();
      cameras.Insert(camera);
    }
    for (int il = 0; il < run->NLasers(); il++) {
      GSVLaser *laser = run->Laser(il);
      laser->scene_index = lasers.NEntries();
      lasers.Insert(laser);
    }
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      segments.Insert(segment);
      segment->scene_index = segments.NEntries();
      for (it = 0; it < segment->NTapestries; it++) {
        GSVTapestry *tapestry = segment->Tapestry(it);
        tapestry->scene_index = tapestries.NEntries();
        tapestries.Insert(tapestry);
      }
      for (int ip = 0; ip < segment->NPanoramas(); ip++) {
        GSVPanorama *panorama = segment->Panorama(ip);
        panorama->scene_index = panoramas.NEntries();
        panoramas.Insert(panorama);
      }
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        scan->scene_index = scans.NEntries();
        scans.Insert(scan);
        for (int ie = 0; ie < scan->Scanlines(); ie++) {
          GSVScanline *scanline = scan->Scanline(ie);
          scanline->scene_index = scanlines.NEntries();
          scanline.Insert(scanline);
        }
      }
    }
  }
}

#endif



////////////////////////////////////////////////////////////////////////
// Read/Write functions
////////////////////////////////////////////////////////////////////////

int GSVScene::
ReadFile(const char *filename, RNBoolean read_points)
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    fprintf(stderr, "Filename %s has no extension (e.g., .gsv)\n", filename);
    return 0;
  }

  // Read file of appropriate type
  int status = 0;
  if (!strncmp(extension, ".txt", 4)) status = ReadASCIIFile(filename, read_points);
  else if (!strncmp(extension, ".lsr", 4)) status = ReadLaserFile(filename, read_points);
  else if (!strncmp(extension, ".gsf", 4)) status = ReadGSFFile(filename, read_points);
  else if (!strncmp(extension, ".gsv", 4)) status = ReadGSVFile(filename, read_points);
  else fprintf(stderr, "Unrecognized file extension %s in %s\n", extension, filename); 
  if (status == 0) return 0;

  // Remember filename
  this->filename = strdup(filename);

  // Set raw directory name
  char buffer[2048];
  strcpy(buffer, filename);
  char *bufferp = strrchr(buffer, '/');
  if (bufferp) *bufferp = '\0';
  else strcpy(buffer, ".");
  strcat(buffer, "/raw_data");
  if (RNFileExists(buffer)) {
    SetRawDataDirectoryName(buffer);
  }

  // Set cache directory name
  strcpy(buffer, filename);
  bufferp = strrchr(buffer, '/');
  if (bufferp) *bufferp = '\0';
  else strcpy(buffer, ".");
  strcat(buffer, "/gsv_data");
  if (RNFileExists(buffer)) {
    SetCacheDataDirectoryName(buffer);
  }

  // Return success
  return 1;
}



int GSVScene::
WriteFile(const char *filename)
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    fprintf(stderr, "Filename %s has no extension (e.g., .gsv)\n", filename);
    return 0;
  }

  // Write file of appropriate type
  if (!strncmp(extension, ".gsv", 4)) return WriteGSVFile(filename);

  // Did not recognize extension
  fprintf(stderr, "Unable to write file %s (unrecognized extension: %s)\n", filename, extension);
  return 0;
}



////////////////////////////////////////////////////////////////////////
// ASCII Read/Write functions
////////////////////////////////////////////////////////////////////////

int GSVScene::
ReadASCIIFile(const char *filename, RNBoolean read_points)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open GSV scene file: %s\n", filename);
    return 0;
  }

  // Read run names from file
  int line_count = 0;
  const int buffer_size = 16*1024;
  char buffer[buffer_size];
  while (fgets(buffer, buffer_size, fp)) {
    line_count++;
    if (buffer[0] == '\0') continue;
    if (buffer[0] == '#') continue;
    char run_name[4096];
    if (sscanf(buffer, "%s", run_name) != (unsigned int) 1) {
      fprintf(stderr, "Unable to read run name from line %d of %s\n", line_count, filename);
      return 0;
    }

    // Create run
    GSVRun *run = new GSVRun(run_name);
    InsertRun(run);

    // Read camera info
    if (!run->ReadCameraInfoFile()) return 0;

    // Read image info
    RNArray<GSVPanorama *> panoramas;
    if (!run->ReadImageSegmentFile(panoramas)) return 0;
    if (!run->ReadImagePoseFile(panoramas)) return 0;

    // Read laser stuff
    const int max_laser_index = 2;
    for (int laser_index = 0; laser_index <= max_laser_index; laser_index++) {
      if (laser_index == 1) continue;
      RNArray<GSVScanline *> scanlines;
      RNArray<GSVScan *> scanline_scans;
      if (!run->ReadLaserObjFile(laser_index, scanlines, scanline_scans, read_points)) return 0;
      if (!run->ReadLaserPoseFile(laser_index, scanlines, scanline_scans)) return 0;
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Laser Read/Write functions
////////////////////////////////////////////////////////////////////////

int GSVScene::
ReadLaserFile(const char *filename, RNBoolean read_points)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open GSV scene file: %s\n", filename);
    return 0;
  }

  // Read run names from file
  int line_count = 0;
  const int buffer_size = 16*1024;
  char buffer[buffer_size];
  while (fgets(buffer, buffer_size, fp)) {
    line_count++;
    if (buffer[0] == '\0') continue;
    if (buffer[0] == '#') continue;
    char run_name[4096];
    if (sscanf(buffer, "%s", run_name) != (unsigned int) 1) {
      fprintf(stderr, "Unable to read run name from line %d of %s\n", line_count, filename);
      return 0;
    }

    // Create run
    GSVRun *run = new GSVRun(run_name);
    InsertRun(run);

    // Read laser stuff
    const int max_laser_index = 2;
    for (int laser_index = 0; laser_index <= max_laser_index; laser_index++) {
      RNArray<GSVScanline *> scanlines;
      RNArray<GSVScan *> scanline_scans;
      if (!run->ReadLaserObjFile(laser_index, scanlines, scanline_scans, read_points)) return 0;
      if (!run->ReadLaserPoseFile(laser_index, scanlines, scanline_scans)) return 0;
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// GSF Read/Write functions
////////////////////////////////////////////////////////////////////////

int GSVScene::
ReadGSFFile(const char *filename, RNBoolean read_points)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open GSV scene file: %s\n", filename);
    return 0;
  }

  // Read run names from file
  int line_count = 0;
  const int buffer_size = 16*1024;
  char buffer[buffer_size];
  while (fgets(buffer, buffer_size, fp)) {
    line_count++;
    if (buffer[0] == '\0') continue;
    if (buffer[0] == '#') continue;

    // Parse command keyword
    char cmd[4096];
    if (sscanf(buffer, "%s", cmd) != (unsigned int) 1) {
      fprintf(stderr, "Unable to read command from line %d of %s\n", line_count, filename);
      return 0;
    }

    // Parse rest of command
    if (!strcmp(cmd, "raw_data")) {
      // Parse directory name
      char directory_name[1024];
      if (sscanf(buffer, "%s%s", cmd, directory_name) != (unsigned int) 2) {
        fprintf(stderr, "Unable to read scan_directory from line %d of %s\n", line_count, filename);
        return 0;
      }

      // Set directory name
      SetRawDataDirectoryName(directory_name);
    }
    else if (!strcmp(cmd, "cache_data")) {
      // Parse directory name
      char directory_name[1024];
      if (sscanf(buffer, "%s%s", cmd, directory_name) != (unsigned int) 2) {
        fprintf(stderr, "Unable to read scan_directory from line %d of %s\n", line_count, filename);
        return 0;
      }

      // Set directory name
      SetCacheDataDirectoryName(directory_name);
    }
    else if (!strcmp(cmd, "run")) {
      // Parse run name
      char run_name[4096];
      if (sscanf(buffer, "%s%s", cmd, run_name) != (unsigned int) 2) {
        fprintf(stderr, "Unable to read run name from line %d of %s\n", line_count, filename);
        return 0;
      }

      // Create run
      GSVRun *run = new GSVRun(run_name);
      InsertRun(run);
      
      // Read camera info
      if (!run->ReadCameraInfoFile()) return 0;
      
      // Read image info
      RNArray<GSVPanorama *> panoramas;
      if (!run->ReadImageSegmentFile(panoramas)) return 0;
      if (!run->ReadImagePoseFile(panoramas)) return 0;

      // Read laser stuff
      const int max_laser_index = 2;
      for (int laser_index = 0; laser_index <= max_laser_index; laser_index++) {
        RNArray<GSVScanline *> scanlines;
        RNArray<GSVScan *> scanline_scans;
        if (!run->ReadLaserObjFile(laser_index, scanlines, scanline_scans, read_points)) return 0;
        if (!run->ReadLaserPoseFile(laser_index, scanlines, scanline_scans)) return 0;
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// GSV Read/Write functions
////////////////////////////////////////////////////////////////////////

int GSVScene::
ReadGSVFile(const char *filename, RNBoolean read_points)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open GSV scene file: %s\n", filename);
    return 0;
  }

  // Convenient variables
  int dummy = 0;

  // Read magic keyword
  char magic[16] = { '\0' };
  if (fread(magic, sizeof(char), 16, fp) != (unsigned int) 16) {
    fprintf(stderr, "Unable to read first bytes from GSV scene file\n");
    fclose(fp);
    return 0;
  }

  // Check magic keyword
  if (strcmp(magic, "GSV")) {
    fprintf(stderr, "Invalid file format for GSV scene file\n");
    fclose(fp);
    return 0;
  }

  // Read header
  int endian_test;
  int major_version;
  int minor_version;
  unsigned long long scene_file_offset = 0;
  fread(&dummy, sizeof(int), 1, fp);
  fread(&dummy, sizeof(int), 1, fp);
  fread(&scene_file_offset, sizeof(unsigned long long), 1, fp);
  fread(&endian_test, sizeof(int), 1, fp);
  fread(&major_version, sizeof(int), 1, fp);
  fread(&minor_version, sizeof(int), 1, fp);
  if ((major_version != 0) || (minor_version != 1)) {
    for (int j = 0; j < 53; j++) {
      fread(&dummy, sizeof(int), 1, fp);
    }
  }

  // Seek to beginning of scene stuff
  assert(((major_version != 0) || (minor_version != 1)) || (scene_file_offset == 0));
  assert(((major_version == 0) && (minor_version == 1)) || (scene_file_offset > 0));
  unsigned long long point_file_offset = RNFileTell(fp);
  if (scene_file_offset > 0) {
    if (!RNFileSeek(fp, scene_file_offset, RN_FILE_SEEK_SET)) {
      fprintf(stderr, "Unable to seek to scene file offset\n");
      fclose(fp);
      return 0;
    }
  }

  // Read scene stuff
  int scene_nruns;
  R3Box scene_bbox;
  fread(&scene_nruns, sizeof(int), 1, fp);
  fread(&(scene_bbox[0][0]), sizeof(RNCoord), 6, fp);
  for (int j = 0; j < 121; j++) {
    fread(&dummy, sizeof(int), 1, fp);
  }
  
  // Read runs
  for (int ir = 0; ir < scene_nruns; ir++) {
    // Create run
    GSVRun *run = new GSVRun();

    // Read run stuff
    int run_ncameras;
    int run_nlasers;
    int run_nsegments;
    R3Box run_bbox;
    char run_name[256];
    fread(run_name, sizeof(char), 256, fp);
    fread(&run_ncameras, sizeof(int), 1, fp);
    fread(&run_nlasers, sizeof(int), 1, fp);
    fread(&run_nsegments, sizeof(int), 1, fp);
    fread(&(run_bbox[0][0]), sizeof(RNCoord), 6, fp);
    for (int j = 0; j < 23; j++) {
      fread(&dummy, sizeof(int), 1, fp);
    }

    // Read cameras
    for (int ic = 0; ic < run_ncameras; ic++) {
      GSVCamera *camera = new GSVCamera();

      // Read camera stuff
      int camera_ntapestries;
      int camera_distortion_type;
      RNLength camera_focal_x;
      RNLength camera_focal_y;
      RNCoord camera_center_x;
      RNCoord camera_center_y;
      RNScalar camera_k1;
      RNScalar camera_k2;
      RNScalar camera_k3;
      RNScalar camera_p1;
      RNScalar camera_p2;
      RNAngle camera_fov_max;
      R3Box camera_bbox;
      fread(&camera_ntapestries, sizeof(int), 1, fp);
      fread(&camera_distortion_type, sizeof(int), 1, fp);
      fread(&camera_focal_x, sizeof(RNLength), 1, fp);
      fread(&camera_focal_y, sizeof(RNLength), 1, fp);
      fread(&camera_center_x, sizeof(RNCoord), 1, fp);
      fread(&camera_center_y, sizeof(RNCoord), 1, fp);
      fread(&camera_k1, sizeof(RNScalar), 1, fp);
      fread(&camera_k2, sizeof(RNScalar), 1, fp);
      fread(&camera_k3, sizeof(RNScalar), 1, fp);
      fread(&camera_p1, sizeof(RNScalar), 1, fp);
      fread(&camera_p2, sizeof(RNScalar), 1, fp);
      fread(&camera_fov_max, sizeof(RNAngle), 1, fp);
      fread(&(camera_bbox[0][0]), sizeof(RNCoord), 6, fp);
      for (int j = 0; j < 15; j++) {
        fread(&dummy, sizeof(int), 1, fp);
      }

      // Set camera stuff
      camera->SetDistortionType(camera_distortion_type);
      camera->SetXFocal(camera_focal_x);
      camera->SetYFocal(camera_focal_y);
      camera->SetXCenter(camera_center_x);
      camera->SetYCenter(camera_center_y);
      camera->SetK1(camera_k1);
      camera->SetK2(camera_k2);
      camera->SetK3(camera_k3);
      camera->SetP1(camera_p1);
      camera->SetP2(camera_p2);
      camera->SetMaxFov(camera_fov_max);
      camera->SetBBox(camera_bbox);
      run->InsertCamera(camera);
    }

    // Read lasers
    for (int il = 0; il < run_nlasers; il++) {
      GSVLaser *laser = new GSVLaser();

      // Read laser stuff
      int laser_nscans;
      R3Box laser_bbox;
      fread(&laser_nscans, sizeof(int), 1, fp);
      fread(&(laser_bbox[0][0]), sizeof(RNCoord), 6, fp);
      for (int j = 0; j < 25; j++) {
        fread(&dummy, sizeof(int), 1, fp);
      }

      // Set laser stuff
      laser->SetBBox(laser_bbox);
      run->InsertLaser(laser);
    }

    // Read segments
    for (int is = 0; is < run_nsegments; is++) {
      GSVSegment *segment = new GSVSegment();

      // Read segment stuff
      int segment_nscans;
      int segment_ntapestries;
      int segment_npanoramas;
      R3Box segment_bbox;
      fread(&segment_nscans, sizeof(int), 1, fp);
      fread(&segment_ntapestries, sizeof(int), 1, fp);
      fread(&segment_npanoramas, sizeof(int), 1, fp);
      fread(&(segment_bbox[0][0]), sizeof(RNCoord), 6, fp);
      for (int j = 0; j < 23; j++) {
        fread(&dummy, sizeof(int), 1, fp);
      }

      // Read scans
      for (int ia = 0; ia < segment_nscans; ia++) {
        GSVScan *scan = new GSVScan();

        // Read scan stuff
        int scan_laser_index;
        int scan_nscanlines;
        R3Box scan_bbox;
        fread(&scan_laser_index, sizeof(int), 1, fp);
        fread(&scan_nscanlines, sizeof(int), 1, fp);
        fread(&(scan_bbox[0][0]), sizeof(RNCoord), 6, fp);
        for (int j = 0; j < 24; j++) {
          fread(&dummy, sizeof(int), 1, fp);
        }

        // Read scanlines
        for (int ie = 0; ie < scan_nscanlines; ie++) {
          GSVScanline *scanline = new GSVScanline();

          // Read scanline stuff
          int scanline_npoints;
          R3Point scanline_viewpoint;
          R3Quaternion scanline_orientation;
          RNScalar scanline_timestamp;
          R3Box scanline_bbox;
          unsigned long long scanline_file_offset;
          fread(&(scanline_npoints), sizeof(int), 1, fp);
          fread(&(scanline_viewpoint[0]), sizeof(RNCoord), 3, fp);
          fread(&(scanline_orientation[0]), sizeof(RNCoord), 4, fp);
          fread(&scanline_timestamp, sizeof(RNScalar), 1, fp);
          fread(&(scanline_bbox[0][0]), sizeof(RNCoord), 6, fp);
          fread(&(scanline_file_offset), sizeof(unsigned long long), 1, fp);
          for (int j = 0; j < 15; j++) {
            fread(&dummy, sizeof(int), 1, fp);
          }

          // Read points
          if ((major_version == 0) && (minor_version == 1)) {
            scanline->ReadPoints(fp, FALSE);
          }

          // Set scanline stuff
          GSVPose scanline_pose(scanline_viewpoint, scanline_orientation);
          scanline->SetPose(scanline_pose);
          scanline->SetTimestamp(scanline_timestamp);
          scanline->SetBBox(scanline_bbox);
          scanline->SetNPoints(scanline_npoints);
          scanline->SetFileOffset(scanline_file_offset);
          scan->InsertScanline(scanline);
        }

        // Set scan stuff
        scan->SetBBox(scan_bbox);
        segment->InsertScan(scan);

        // Insert scan into laser
        if (scan_laser_index >= 0) {
          GSVLaser *scan_laser = run->Laser(scan_laser_index);
          scan_laser->InsertScan(scan);
        }
      }

      // Read tapestries
      for (int it = 0; it < segment_ntapestries; it++) {
        GSVTapestry *tapestry = new GSVTapestry();

        // Read tapestry stuff
        int tapestry_camera_index;
        int tapestry_nimages;
        R3Box tapestry_bbox;
        fread(&tapestry_camera_index, sizeof(int), 1, fp);
        fread(&tapestry_nimages, sizeof(int), 1, fp);
        fread(&(tapestry_bbox[0][0]), sizeof(RNCoord), 6, fp);
        for (int j = 0; j < 24; j++) {
          fread(&dummy, sizeof(int), 1, fp);
        }

        // Set tapestry stuff
        tapestry->SetBBox(tapestry_bbox);
        segment->InsertTapestry(tapestry);
        if (tapestry_camera_index >= 0) {
          GSVCamera *camera = run->Camera(tapestry_camera_index);
          camera->InsertTapestry(tapestry);
        }
      }

      // Read panoramas
      for (int ip = 0; ip < segment_npanoramas; ip++) {
        GSVPanorama *panorama = new GSVPanorama();

        // Read panorama stuff
        int panorama_nimages = panorama->NImages();
        R3Point panorama_viewpoint = panorama->Viewpoint();
        RNScalar panorama_timestamp = panorama->Timestamp();
        fread(&panorama_nimages, sizeof(int), 1, fp);
        fread(&(panorama_viewpoint[0]), sizeof(RNCoord), 3, fp);
        fread(&panorama_timestamp, sizeof(RNScalar), 1, fp);
        for (int j = 0; j < 22; j++) {
          fread(&dummy, sizeof(int), 1, fp);
        }

        // Read images
        for (int ii = 0; ii < panorama_nimages; ii++) {
          GSVImage *image = new GSVImage();

          // Read image stuff
          int image_tapestry_index, image_ncolumns;
          R3Point image_viewpoint0, image_viewpoint1;
          R3Quaternion image_orientation0, image_orientation1;
          RNScalar image_timestamp0, image_timestamp1;
          fread(&image_tapestry_index, sizeof(int), 1, fp);
          fread(&image_ncolumns, sizeof(int), 1, fp);
          fread(&(image_viewpoint0[0]), sizeof(RNCoord), 3, fp);
          fread(&(image_viewpoint1[0]), sizeof(RNCoord), 3, fp);
          fread(&(image_orientation0[0]), sizeof(RNCoord), 4, fp);
          fread(&(image_orientation1[0]), sizeof(RNCoord), 4, fp);
          fread(&image_timestamp0, sizeof(RNScalar), 1, fp);
          fread(&image_timestamp1, sizeof(RNScalar), 1, fp);
          for (int j = 0; j < 14; j++) {
            fread(&dummy, sizeof(int), 1, fp);
          }

          // Set image stuff
          GSVPose image_pose0(image_viewpoint0, image_orientation0);
          GSVPose image_pose1(image_viewpoint1, image_orientation1);
          image->SetPose(image_pose0, image_pose1);
          image->SetTimestamp(image_timestamp0, image_timestamp1);
          panorama->InsertImage(image);

          // Insert image into tapestry
          if (image_tapestry_index >= 0) {
            GSVTapestry *tapestry = segment->Tapestry(image_tapestry_index);
            tapestry->InsertImage(image);
          }
        }

        // Set panorama stuff
        panorama->SetViewpoint(panorama_viewpoint);
        panorama->SetTimestamp(panorama_timestamp);
        segment->InsertPanorama(panorama);
      }

      // Set segment stuff
      segment->SetBBox(segment_bbox);
      run->InsertSegment(segment);
    }

    // Set run stuff
    run->SetName(run_name);
    run->SetBBox(run_bbox);
    InsertRun(run);
  }

  // Read points
  if (read_points) {
    if ((major_version != 0) || (minor_version != 1)) {
      if (!RNFileSeek(fp, point_file_offset, RN_FILE_SEEK_SET)) {
        fprintf(stderr, "Unable to seek to point file offset\n");
        fclose(fp);
        return 0;
      }
      for (int ir = 0; ir < NRuns(); ir++) {
        GSVRun *run = Run(ir);
        for (int is = 0; is < run->NSegments(); is++) {
          GSVSegment *segment = run->Segment(is);
          for (int ia = 0; ia < segment->NScans(); ia++) {
            GSVScan *scan = segment->Scan(ia);
            for (int ie = 0; ie < scan->NScanlines(); ie++) {
              GSVScanline *scanline = scan->Scanline(ie);
              scanline->ReadPoints(fp, FALSE);
            }
          }
        }
      }
    }
  }

  // Set scene stuff
  SetBBox(scene_bbox);

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int GSVScene::
WriteGSVFile(const char *filename)
{
  // Convenient variables
  int dummy = 0;

  // Make sure everything is uptodate
  R3Box foobar = BBox();
  if (foobar[0][0] == FLT_MAX) printf("Ouch\n");

  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open GSV scene file: %s\n", filename);
    return 0;
  }

  // Write magic keyword
  char magic[16] = { '\0' };
  strncpy(magic, "GSV", 16);
  if (fwrite(magic, sizeof(char), 16, fp) != (unsigned int) 16) {
    fprintf(stderr, "Unable to write GSV scene file: %s\n", filename);
    return 0;
  }

  // Write header
  int endian_test = 1;
  int major_version = 0;
  int minor_version = 2;
  unsigned long long scene_file_offset = 0;
  fwrite(&dummy, sizeof(int), 1, fp);
  fwrite(&dummy, sizeof(int), 1, fp);
  fwrite(&scene_file_offset, sizeof(unsigned long long), 1, fp);
  fwrite(&endian_test, sizeof(int), 1, fp);
  fwrite(&major_version, sizeof(int), 1, fp);
  fwrite(&minor_version, sizeof(int), 1, fp);
  for (int j = 0; j < 53; j++) {
    fwrite(&dummy, sizeof(int), 1, fp);
  }

  // Write points
  for (int ir = 0; ir < NRuns(); ir++) {
    GSVRun *run = Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        scan->ReadPoints();
        for (int ie = 0; ie < scan->NScanlines(); ie++) {
          GSVScanline *scanline = scan->Scanline(ie);
          scanline->WritePoints(fp, FALSE);
        }
        scan->ReleasePoints();
      }
    }
  }

  // Get scene file offset
  scene_file_offset = RNFileTell(fp);

  // Write scene 
  int scene_nruns = NRuns();
  R3Box scene_bbox = BBox();
  fwrite(&scene_nruns, sizeof(int), 1, fp);
  fwrite(&(scene_bbox[0][0]), sizeof(RNCoord), 6, fp);
  for (int j = 0; j < 121; j++) {
    fwrite(&dummy, sizeof(int), 1, fp);
  }

  // Write runs
  for (int ir = 0; ir < NRuns(); ir++) {
    GSVRun *run = Run(ir);
    int run_ncameras = run->NCameras();
    int run_nlasers = run->NLasers();
    int run_nsegments = run->NSegments();
    R3Box run_bbox = run->BBox();
    char run_name[256] = { '\0' };
    strncpy(run_name, run->Name(), 255);
    fwrite(run_name, sizeof(char), 256, fp);
    fwrite(&run_ncameras, sizeof(int), 1, fp);
    fwrite(&run_nlasers, sizeof(int), 1, fp);
    fwrite(&run_nsegments, sizeof(int), 1, fp);
    fwrite(&(run_bbox[0][0]), sizeof(RNCoord), 6, fp);
    for (int j = 0; j < 23; j++) {
      fwrite(&dummy, sizeof(int), 1, fp);
    }

    // Write cameras
    for (int ic = 0; ic < run->NCameras(); ic++) {
      GSVCamera *camera = run->Camera(ic);
      int camera_ntapestries = camera->NTapestries();
      int camera_distortion_type = camera->DistortionType();
      RNLength camera_focal_x = camera->XFocal();
      RNLength camera_focal_y = camera->YFocal();
      RNCoord camera_center_x = camera->XCenter();
      RNCoord camera_center_y = camera->YCenter();
      RNScalar camera_k1 = camera->K1();
      RNScalar camera_k2 = camera->K2();
      RNScalar camera_k3 = camera->K3();
      RNScalar camera_p1 = camera->P1();
      RNScalar camera_p2 = camera->P2();
      RNAngle camera_fov_max = camera->MaxFov();
      R3Box camera_bbox = camera->BBox();
      fwrite(&camera_ntapestries, sizeof(int), 1, fp);
      fwrite(&camera_distortion_type, sizeof(int), 1, fp);
      fwrite(&camera_focal_x, sizeof(RNLength), 1, fp);
      fwrite(&camera_focal_y, sizeof(RNLength), 1, fp);
      fwrite(&camera_center_x, sizeof(RNCoord), 1, fp);
      fwrite(&camera_center_y, sizeof(RNCoord), 1, fp);
      fwrite(&camera_k1, sizeof(RNScalar), 1, fp);
      fwrite(&camera_k2, sizeof(RNScalar), 1, fp);
      fwrite(&camera_k3, sizeof(RNScalar), 1, fp);
      fwrite(&camera_p1, sizeof(RNScalar), 1, fp);
      fwrite(&camera_p2, sizeof(RNScalar), 1, fp);
      fwrite(&camera_fov_max, sizeof(RNAngle), 1, fp);
      fwrite(&(camera_bbox[0][0]), sizeof(RNCoord), 6, fp);
      for (int j = 0; j < 15; j++) {
        fwrite(&dummy, sizeof(int), 1, fp);
      }
    }

    // Write lasers
    for (int il = 0; il < run->NLasers(); il++) {
      GSVLaser *laser = run->Laser(il);
      int laser_nscans = laser->NScans();
      R3Box laser_bbox = laser->BBox();
      fwrite(&laser_nscans, sizeof(int), 1, fp);
      fwrite(&(laser_bbox[0][0]), sizeof(RNCoord), 6, fp);
      for (int j = 0; j < 25; j++) {
        fwrite(&dummy, sizeof(int), 1, fp);
      }
    }

    // Write segments
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      int segment_nscans = segment->NScans();
      int segment_ntapestries = segment->NTapestries();
      int segment_npanoramas = segment->NPanoramas();
      R3Box segment_bbox = segment->BBox();
      fwrite(&segment_nscans, sizeof(int), 1, fp);
      fwrite(&segment_ntapestries, sizeof(int), 1, fp);
      fwrite(&segment_npanoramas, sizeof(int), 1, fp);
      fwrite(&(segment_bbox[0][0]), sizeof(RNCoord), 6, fp);
      for (int j = 0; j < 23; j++) {
        fwrite(&dummy, sizeof(int), 1, fp);
      }

      // Write scans
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        int scan_laser_index = (scan->Laser()) ? scan->Laser()->RunIndex() : -1;
        int scan_nscanlines = scan->NScanlines();
        R3Box scan_bbox = scan->BBox();
        fwrite(&scan_laser_index, sizeof(int), 1, fp);
        fwrite(&scan_nscanlines, sizeof(int), 1, fp);
        fwrite(&(scan_bbox[0][0]), sizeof(RNCoord), 6, fp);
        for (int j = 0; j < 24; j++) {
          fwrite(&dummy, sizeof(int), 1, fp);
        }

        // Write scanlines
        for (int ie = 0; ie < scan->NScanlines(); ie++) {
          GSVScanline *scanline = scan->Scanline(ie);
          int scanline_npoints = scanline->NPoints();
          R3Point scanline_viewpoint = scanline->Pose().Viewpoint();
          R3Quaternion scanline_orientation = scanline->Pose().Orientation();
          RNScalar scanline_timestamp = scanline->Timestamp();
          R3Box scanline_bbox = scanline->BBox();
          unsigned long long scanline_file_offset = scanline->FileOffset();
          fwrite(&(scanline_npoints), sizeof(int), 1, fp);
          fwrite(&(scanline_viewpoint[0]), sizeof(RNCoord), 3, fp);
          fwrite(&(scanline_orientation[0]), sizeof(RNCoord), 4, fp);
          fwrite(&scanline_timestamp, sizeof(RNScalar), 1, fp);
          fwrite(&(scanline_bbox[0][0]), sizeof(RNCoord), 6, fp);
          fwrite(&(scanline_file_offset), sizeof(unsigned long long), 1, fp);
          for (int j = 0; j < 15; j++) {
            fwrite(&dummy, sizeof(int), 1, fp);
          }
        }
      }

      // Write tapestries
      for (int it = 0; it < segment->NTapestries(); it++) {
        GSVTapestry *tapestry = segment->Tapestry(it);
        int tapestry_camera_index = (tapestry->Camera()) ? tapestry->Camera()->RunIndex() : -1;
        int tapestry_nimages = tapestry->NImages();
        R3Box tapestry_bbox = segment->BBox();
        fwrite(&tapestry_camera_index, sizeof(int), 1, fp);
        fwrite(&tapestry_nimages, sizeof(int), 1, fp);
        fwrite(&(tapestry_bbox[0][0]), sizeof(RNCoord), 6, fp);
        for (int j = 0; j < 24; j++) {
          fwrite(&dummy, sizeof(int), 1, fp);
        }
      }

      // Write panoramas
      for (int ip = 0; ip < segment->NPanoramas(); ip++) {
        GSVPanorama *panorama = segment->Panorama(ip);
        int panorama_nimages = panorama->NImages();
        R3Point panorama_viewpoint = panorama->Viewpoint();
        RNScalar panorama_timestamp = panorama->Timestamp();
        fwrite(&panorama_nimages, sizeof(int), 1, fp);
        fwrite(&(panorama_viewpoint[0]), sizeof(RNCoord), 3, fp);
        fwrite(&panorama_timestamp, sizeof(RNScalar), 1, fp);
        for (int j = 0; j < 22; j++) {
          fwrite(&dummy, sizeof(int), 1, fp);
        }

        // Write images
        for (int ii = 0; ii < panorama->NImages(); ii++) {
          GSVImage *image = panorama->Image(ii);
          int image_tapestry_index = (image->Tapestry()) ? image->Tapestry()->SegmentIndex() : -1;
          int image_ncolumns = image->Width();
          R3Point image_viewpoint0 = image->Pose(0).Viewpoint();
          R3Point image_viewpoint1 = image->Pose(image->Width()-1).Viewpoint();
          R3Quaternion image_orientation0 = image->Pose(0).Orientation();
          R3Quaternion image_orientation1 = image->Pose(image->Width()-1).Orientation();
          RNScalar image_timestamp0 = image->Timestamp(0);
          RNScalar image_timestamp1 = image->Timestamp(image->Width()-1);
          fwrite(&image_tapestry_index, sizeof(int), 1, fp);
          fwrite(&image_ncolumns, sizeof(int), 1, fp);
          fwrite(&(image_viewpoint0[0]), sizeof(RNCoord), 3, fp);
          fwrite(&(image_viewpoint1[0]), sizeof(RNCoord), 3, fp);
          fwrite(&(image_orientation0[0]), sizeof(RNCoord), 4, fp);
          fwrite(&(image_orientation1[0]), sizeof(RNCoord), 4, fp);
          fwrite(&image_timestamp0, sizeof(RNScalar), 1, fp);
          fwrite(&image_timestamp1, sizeof(RNScalar), 1, fp);
          for (int j = 0; j < 14; j++) {
            fwrite(&dummy, sizeof(int), 1, fp);
          }
        }
      }
    }
  }

  // Re-write header with updated scene file offset
  RNFileSeek(fp, 0, RN_FILE_SEEK_SET);
  fwrite(magic, sizeof(char), 16, fp);
  fwrite(&dummy, sizeof(int), 1, fp);
  fwrite(&dummy, sizeof(int), 1, fp);
  fwrite(&scene_file_offset, sizeof(unsigned long long), 1, fp);
  fwrite(&endian_test, sizeof(int), 1, fp);
  fwrite(&major_version, sizeof(int), 1, fp);
  fwrite(&minor_version, sizeof(int), 1, fp);

  // Close file
  fclose(fp);

  // Return success
  return 1;
}






