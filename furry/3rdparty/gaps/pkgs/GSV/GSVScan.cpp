// Source file for GSV laser class


////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV.h"



////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

GSVScan::
GSVScan(void)
  : segment(NULL),
    segment_index(-1),
    laser(NULL),
    laser_index(-1),
    read_count(0),
    scanlines(),
    travel_distances(NULL),
    bbox(FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX),
    data(NULL)
{
}



GSVScan::
~GSVScan(void)
{
  // Remove scan from segment
  if (segment) segment->RemoveScan(this);

  // Remove scan from laser
  if (laser) laser->RemoveScan(this);

  // Delete travel distances
  if (travel_distances) delete [] travel_distances;

  // Delete scanlines
  while (NScanlines()) delete Scanline(0);
}



////////////////////////////////////////////////////////////////////////
// Access functions
////////////////////////////////////////////////////////////////////////

GSVScene *GSVScan::
Scene(void) const
{
  // Return scene
  if (!segment) return NULL;
  return segment->Scene();
}



GSVRun *GSVScan::
Run(void) const
{
  // Return run
  if (!segment) return NULL;
  return segment->Run();
}



int GSVScan::
SceneIndex(void) const
{
  // Get convenient variables
  GSVRun *run = Run();
  if (!run) return -1;
  GSVScene *scene = run->Scene();
  if (!scene) return -1;

  // Compute scene index
  int scene_index = RunIndex();
  for (int ir = 0; ir < run->SceneIndex(); ir++) {
    GSVRun *r = scene->Run(ir);
    scene_index += r->NScans();
  }

  // Return scene index
  return scene_index;
}



int GSVScan::
RunIndex(void) const
{
  // Get convenient variables
  GSVRun *run = Run();
  if (!run) return -1;

  // Compute run index
  int run_index = SegmentIndex();
  for (int i = 0; i < segment->RunIndex(); i++) {
    GSVSegment *s = run->Segment(i);
    run_index += s->NScans();
  }

  // Return run index
  return run_index;
}



////////////////////////////////////////////////////////////////////////
// Property functions
////////////////////////////////////////////////////////////////////////

GSVScanline *GSVScan::
FindScanlineBeforeTimestamp(RNScalar timestamp, int imin, int imax) const
{
  // Binary search
  int i = (imin + imax) / 2;
  if (i == imin) return Scanline(imin);
  assert(i < imax);
  GSVScanline *scanline = Scanline(i);
  RNScalar t = scanline->Timestamp();
  if (t > timestamp) return FindScanlineBeforeTimestamp(timestamp, imin, i);
  else if (t < timestamp) return FindScanlineBeforeTimestamp(timestamp, i, imax);
  else return scanline;
}



GSVScanline *GSVScan::
FindScanlineBeforeTimestamp(RNScalar timestamp) const
{
  // Binary search
  if (NScanlines() == 0) return NULL;
  if (timestamp <= Scanline(0)->Timestamp()) return Scanline(0);
  if (timestamp >= Scanline(NScanlines()-1)->Timestamp()) return Scanline(NScanlines()-1);
  return FindScanlineBeforeTimestamp(timestamp, 0, NScanlines()-1);
}



GSVPose GSVScan::
Pose(RNScalar timestamp) const
{
  // Find scanlines before and after timestamp
  GSVScanline *scanline1 = FindScanlineBeforeTimestamp(timestamp);
  if (!scanline1) return GSVPose(R3zero_point, R3zero_quaternion);
  int index1 = scanline1->ScanIndex();
  if (index1 < 0) return Scanline(0)->Pose();
  if (index1 >= NScanlines()-1) return Scanline(NScanlines()-1)->Pose();
  int index2 = index1 + 1;
  GSVScanline *scanline2 = Scanline(index2);

  // Compute interpolation parameter (0 <= t <= 1)
  RNScalar timestamp1 = scanline1->Timestamp();
  RNScalar timestamp2 = scanline2->Timestamp();
  RNScalar delta_timestamp = timestamp2 - timestamp1;
  if (delta_timestamp == 0) return scanline1->Pose();
  RNScalar t = (timestamp - timestamp1) / delta_timestamp;

  // Interpolate poses (would be better with quaternions)
  const GSVPose& pose1 = scanline1->Pose();
  const GSVPose& pose2 = scanline2->Pose();
  const R3Point& viewpoint1 = pose1.Viewpoint();
  const R3Point& viewpoint2 = pose2.Viewpoint();
  const R3Quaternion& orientation1 = pose1.Orientation();
  const R3Quaternion& orientation2 = pose2.Orientation();
  R3Point viewpoint = (1-t)*viewpoint1 + t*viewpoint2;
  R3Quaternion orientation = R3QuaternionSlerp(orientation1, orientation2, t);

  // Return interpolated pose
  return GSVPose(viewpoint, orientation);
}



GSVMesh *GSVScan::
Mesh(void) const
{
  // Get convenient variables (parameters)
  GSVScan *scan = (GSVScan *) this;
  if (NScanlines() == 0) return NULL;
  int scan_index = SegmentIndex();
  GSVSegment *segment = Segment();
  if (!segment) return NULL;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return NULL;
  char run_name[4096];
  if (run->Name()) sprintf(run_name, "%s", run->Name());
  else sprintf(run_name, "Run%d", run->SceneIndex());
  GSVScene *scene = run->Scene();
  if (!scene) return NULL;
  const char *cache_name = scene->CacheDataDirectoryName();
  if (!cache_name) return NULL;

  // Create mesh
  GSVMesh *mesh = new GSVMesh();
  if (!mesh) {
    fprintf(stderr, "Unable to allocate mesh\n");
    return NULL;
  }

  // Construct mesh filename
  char mesh_name[4096];
  sprintf(mesh_name, "%s/laser_meshes/%s/%02d_%02d_Mesh.ply", 
    cache_name, run_name, segment_index, scan_index);

  // Check if mesh file exists
  if (RNFileExists(mesh_name)) {
    // Read file 
    if (!mesh->ReadPlyFile(scan, mesh_name)) {
      fprintf(stderr, "Unable to read mesh file %s\n", mesh_name);
      delete mesh;
      return NULL;
    }
  }
  else {
    // Compute mesh 
    if (!mesh->LoadScan(scan)) {
      fprintf(stderr, "Unable to create mesh\n");
      delete mesh;
      return NULL;
    }

    // Write mesh to file
    if (scene->CreateCacheDataDirectory("laser_meshes", run)) {
      fprintf(stderr, "    Creating %s\n", mesh_name);
      if (!mesh->WritePlyFile(scan, mesh_name, TRUE)) {
        fprintf(stderr, "Unable to write mesh file %s\n", mesh_name);
      }
    }
  }

  // Return mesh
  return mesh;
}



////////////////////////////////////////////////////////////////////////
//Input/Output functions
////////////////////////////////////////////////////////////////////////

int GSVScan::
ReadPoints(void) 
{
  // Check read counter
  if (++read_count > 1) return 1;

  // Get the scene filename
  if (NScanlines() == 0) return 1;
  GSVSegment *segment = Segment();
  if (!segment) return 0;
  GSVRun *run = segment->Run();
  if (!run) return 0;
  GSVScene *scene = run->Scene();
  if (!scene) return 0;
  const char *filename = scene->Filename();
  if (!filename) return 0;

  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open GSV scene file: %s\n", filename);
    return 0;
  }

  // Read scanlines
  for (int i = 0; i < NScanlines(); i++) {
    GSVScanline *scanline = Scanline(i);
    if (!scanline->ReadPoints(fp)) return 0;
  }

  // Close file 
  fclose(fp);

  // Return success
  return 1;
}



int GSVScan::
ReleasePoints(void) 
{
  // Check read counter
  if (--read_count > 0) return 1;

  // Release points
  for (int i = 0; i < NScanlines(); i++) {
    GSVScanline *scanline = Scanline(i);
    if (!scanline->ReleasePoints()) return 0;
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Manipulation functions
////////////////////////////////////////////////////////////////////////

void GSVScan::
InsertScanline(GSVScanline *scanline)
{
  // Just checking
  assert(scanline->scan_index == -1);
  assert(scanline->scan == NULL);

  // Insert scanline
  scanline->scan = this;
  scanline->scan_index = scanlines.NEntries();
  scanlines.Insert(scanline);

  // Invalidate bounding box
  InvalidateBBox();
}



void GSVScan::
RemoveScanline(GSVScanline *scanline)
{
  // Just checking
  assert(scanline->scan_index >= 0);
  assert(scanline->scan_index < scanlines.NEntries());
  assert(scanline->scan == this);

  // Remove scanline
  RNArrayEntry *entry = scanlines.KthEntry(scanline->scan_index);
  GSVScanline *tail = scanlines.Tail();
  tail->scan_index = scanline->scan_index;
  scanlines.EntryContents(entry) = tail;
  scanlines.RemoveTail();
  scanline->scan_index = -1;
  scanline->scan = NULL;

  // Invalidate bounding box
  InvalidateBBox();
}



////////////////////////////////////////////////////////////////////////
// Display functions
////////////////////////////////////////////////////////////////////////

void GSVScan::
Draw(RNFlags flags) const
{
  if (flags & GSV_DRAW_POINTS_WITH_LASER_INDEX_COLOR) {
    // Draw points  colored by laser index
    const int ncolors = 3;
    GLfloat color[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
    glColor3fv(color[segment_index % ncolors]);
    for (int i = 0; i < NScanlines(); i++) {
      GSVScanline *scanline = Scanline(i);
      scanline->Draw(flags);
    }
  }
  else if (flags & GSV_DRAW_MESHES_WITHOUT_COLOR) {
    glEnable(GL_CULL_FACE);
    static GSVMesh *last_mesh = NULL;
    static const GSVScan *last_scan = NULL;
    if (this != last_scan) {
      if (last_mesh) delete last_mesh;
      last_mesh = Mesh();
      last_scan = this;
    }
    if (last_mesh) last_mesh->Draw();
    glDisable(GL_CULL_FACE);
  }
  else {
    // Draw scanlines
    for (int i = 0; i < NScanlines(); i++) {
      GSVScanline *scanline = Scanline(i);
      scanline->Draw(flags);
    }
  }
}



void GSVScan::
Print(FILE *fp, const char *prefix, const char *suffix) const
{
  // Check fp
  if (!fp) fp = stdout;

  // Print scan header
  if (prefix) fprintf(fp, "%s", prefix);
  fprintf(fp, "Scan %d:", segment_index);
  if (suffix) fprintf(fp, "%s", suffix);
  fprintf(fp, "\n");

  // Add indentation to prefix
  char indented_prefix[1024];
  sprintf(indented_prefix, "%s  ", prefix);

   // Print scanlines
  for (int i = 0; i < NScanlines(); i++) {
    GSVScanline *scanline = Scanline(i);
    scanline->Print(fp, indented_prefix, suffix);
  }
}



////////////////////////////////////////////////////////////////////////
// Update functions
////////////////////////////////////////////////////////////////////////

void GSVScan::
UpdateBBox(void) 
{
  // Read points (assumes all scanlines are in memory or not) 
  RNBoolean read_points = FALSE;
  if (NScanlines() > 0) {
    GSVScanline *scanline = Scanline(0);
    if (scanline->DoesBBoxNeedUpdate()) {
      ReadPoints();
      read_points = TRUE;
    }
  }

  // Update bounding box
  bbox = R3null_box;
  for (int i = 0; i < NScanlines(); i++) {
    bbox.Union(Scanline(i)->BBox());
  }

  // Release points
  if (read_points) {
    ReleasePoints();
  }
}



void GSVScan::
InvalidateBBox(void) 
{
  // Invalidate bounding box
  if (bbox.XMin() != FLT_MAX) {
    bbox = R3Box(FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX);
    if (segment) segment->InvalidateBBox();
    if (laser) laser->InvalidateBBox();
  }
}



void GSVScan::
UpdateTravelDistances(void) 
{
  // Check stuff
  if (travel_distances) return;
  if (NScanlines() == 0) return;

  // Allocate travel distances
  travel_distances = new RNLength [ NScanlines() ];
  if (!travel_distances) {
    fprintf(stderr, "Unable to allocate travel distances\n");
    return;
  }

  // Compute travel distances
  RNLength travel_distance = 0;
  R3Point prev_viewpoint = Scanline(0)->Pose().Viewpoint();
  for (int i = 0; i < NScanlines(); i++) {
    GSVScanline *scanline = Scanline(i);
    const R3Point& viewpoint = scanline->Pose().Viewpoint();
    travel_distance += R3Distance(viewpoint, prev_viewpoint);
    travel_distances[i] = travel_distance;
    prev_viewpoint = viewpoint;
  }
}



void GSVScan::
InvalidateTravelDistances(void) 
{
  // Invalidate travel distances
  if (!travel_distances) return;
  delete [] travel_distances;
  travel_distances = NULL;
}



