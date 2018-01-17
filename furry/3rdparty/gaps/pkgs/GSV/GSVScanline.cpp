// Source file for GSV scanline class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV.h"



////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

GSVScanline::
GSVScanline(void)
  : scan(NULL),
    scan_index(-1),
    file_offset(0),
    read_count(0),
    points(NULL),
    npoints(0),
    pose(0,0,0,0,0,0,0),
    timestamp(-1), 
    bbox(FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX),
    data(NULL)
{
}



GSVScanline::
GSVScanline(const GSVPose& pose, RNScalar timestamp)
  : scan(NULL),
    scan_index(-1),
    file_offset(0),
    read_count(0),
    points(NULL),
    npoints(0),
    pose(pose),
    timestamp(timestamp),
    bbox(FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX),
    data(NULL)
{
}



GSVScanline::
~GSVScanline(void)
{
  // Delete points
  if (points) delete [] points;

  // Remove scanline from laser scan
  if (scan) scan->RemoveScanline(this);
}



////////////////////////////////////////////////////////////////////////
// Access functions
////////////////////////////////////////////////////////////////////////

GSVScene *GSVScanline::
Scene(void) const
{
  // Return scene
  if (!scan) return NULL;
  return scan->Scene();
}



GSVRun *GSVScanline::
Run(void) const
{
  // Return run
  if (!scan) return NULL;
  return scan->Run();
}



GSVSegment *GSVScanline::
Segment(void) const
{
  // Return segment
  if (!scan) return NULL;
  return scan->Segment();
}



int GSVScanline::
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
    scene_index += r->NScanlines();
  }

  // Return scene index
  return scene_index;
}



int GSVScanline::
RunIndex(void) const
{
  // Get convenient variables
  GSVSegment *segment = Segment();
  if (!segment) return -1;
  GSVRun *run = segment->Run();
  if (!run) return -1;

  // Compute run index
  int run_index = SegmentIndex();
  for (int i = 0; i < segment->RunIndex(); i++) {
    GSVSegment *s = run->Segment(i);
    run_index += s->NScanlines();
  }

  // Return runindex
  return run_index;
}



int GSVScanline::
SegmentIndex(void) const
{
  // Get convenient variables
  GSVSegment *segment = Segment();
  if (!segment) return -1;

  // Compute segment index
  int segment_index = ScanIndex();
  for (int i = 0; i < scan->SegmentIndex(); i++) {
    GSVScan *s = segment->Scan(i);
    segment_index += s->NScanlines();
  }

  // Return segment index
  return segment_index;
}



////////////////////////////////////////////////////////////////////////
// Property functions
////////////////////////////////////////////////////////////////////////

RNLength GSVScanline::
TravelDistance(void) const
{
  // Return distance viewpoint of laser of scan travelled up to scanline
  if (!scan) return 0.0;
  return scan->TravelDistance(scan_index);
}



RNCoord GSVScanline::
EstimatedGroundZ(void) const
{
  // Just checking
  if (NPoints() == 0) return RN_UNKNOWN;
  if (!points) return RN_UNKNOWN;

  // Parameters
  const int scan_ground_index = 25;
  const int radius = 1;

  // Check number of points
  if (NPoints() < scan_ground_index+radius) return RN_UNKNOWN;

  // Get three points around scan_ground_index
  R3Point p0 = PointPosition(NPoints() - scan_ground_index - radius);
  R3Point p1 = PointPosition(NPoints() - scan_ground_index);
  R3Point p2 = PointPosition(NPoints() - scan_ground_index + radius);

  // Find median Z coordinate
  RNCoord z;
  if (p0.Z() < p1.Z()) {
    if (p1.Z() < p2.Z()) {
      z = p1.Z();
    }
    else {
      if (p0.Z() < p2.Z()) z = p2.Z();
      else z = p0.Z();
    }
  }
  else {
    if (p0.Z() < p2.Z()) {
      z = p0.Z();
    }
    else {
      if (p1.Z() < p2.Z()) z = p2.Z();
      else z = p1.Z();
    }
  }

  // Return estimate of ground Z coordinate
  return z;
}



RNAngle GSVScanline::
PointAngle(int point_index) const
{
  // Return angle of laser when point was sampled (0 - PI)
  const R3Point position = PointPosition(point_index);
  const R3Point& viewpoint = pose.Viewpoint();
  const R3Vector& towards = pose.Towards();
  const R3Vector& up = pose.Up();
  R3Vector right = towards % up;
  R3Vector v = position - viewpoint;
  RNLength distance = v.Length();
  if (RNIsZero(distance)) return 0;
  v /= distance;
  RNScalar dot = -(v.Dot(up));
  RNAngle angle = acos(dot);
  if (v.Dot(right) > 0) angle = RN_PI_OVER_TWO + angle;
  else angle = RN_PI_OVER_TWO - angle;
  return angle;
}



////////////////////////////////////////////////////////////////////////
// Manipulation functions
////////////////////////////////////////////////////////////////////////

void GSVScanline::
SetPointPositions(const R3Point *positions, int npositions) 
{
  // Delete previous points
  if (points) delete [] points;

  // Reset everything
  npoints = npositions;
  points = NULL;
  bbox = R3null_box;

  // Copy points
  if (positions && (npositions > 0)) {
    // Allocate points 
    points = new R3Point [ npoints ];
    if (!points) {
      fprintf(stderr, "Unable to allocate points for scanline\n");
      abort();
    }
    
    // Copy points
    for (int i = 0; i < npoints; i++) {
      points[i] = positions[i];
    }

    // Update bounding box
    for (int i = 0; i < npositions; i++) {
      bbox.Union(positions[i]);
    }
  }

  // Set read count
  read_count = 1;
}



void GSVScanline::
SetPointPosition(int point_index, const R3Point& position)
{
  // Just checking
  assert(points);

  // Set point position
  points[point_index] = position;

  // Update bounding box
  bbox.Union(position);
}



void GSVScanline::
SetPose(const GSVPose& pose)
{
  // Update points
  if ((NPoints() > 0) && (points) && (!(this->pose.Orientation().IsZero()))) {
    // Determine transformation of points
    R4Matrix matrix = R4identity_matrix;
    matrix.Translate(pose.Viewpoint().Vector());
    matrix.Multiply(pose.Orientation().Matrix());
    matrix.Multiply(this->pose.Orientation().Inverse().Matrix());
    matrix.Translate(-this->pose.Viewpoint().Vector());

    // Transform points
    for (int i = 0; i < npoints; i++) {
      points[i] = matrix * points[i];
    }
  }

  // Set pose
  this->pose = pose;

  // Invalidate bounding box
  InvalidateBBox();
}



////////////////////////////////////////////////////////////////////////
// Display functions
////////////////////////////////////////////////////////////////////////

void GSVScanline::
Draw(RNFlags flags) const
{
  // Just checking
  if (NPoints() == 0) return;
  assert(points);

  if (flags & GSV_DRAW_POINTS_WITH_VIEWPOINT_DISTANCE_COLOR) {
    // Get viewpoint
    R3Point v = Pose().Viewpoint();

    // Draw points
    glBegin(GL_POINTS);
    for (int i = 0; i < NPoints(); i++) {
      const R3Point& p = PointPosition(i);
      double dx = v[0] - p[0];
      double dy = v[1] - p[1];
      double dz = v[2] - p[2];
      double d2 = dx*dx + dy*dy + dz*dz;
      // double value = 0.0025 * d2;
      double value = 0.001 * d2;
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
      RNLoadRgb(r, g, b);
      R3LoadPoint(PointPosition(i));
    }
    glEnd();
  }
  else if (flags & GSV_DRAW_POINTS_WITH_HEIGHT_COLOR) {
    // Get height
    RNCoord ground_z = EstimatedGroundZ();

    // Draw points
    glBegin(GL_POINTS);
    for (int i = 0; i < NPoints(); i++) {
      const R3Point& p = PointPosition(i);
      double dz = p[2] - ground_z;
      double value = 0.25 * dz;
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
      RNLoadRgb(r, g, b);
      R3LoadPoint(PointPosition(i));
    }
    glEnd();
  }
  else if (flags & GSV_DRAW_SCANLINES_WITH_POINT_INDEX_COLOR) {
    const int ncolors = 3;
    GLfloat color[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
    RNLength max_edge_length_squared = 1;
    glBegin(GL_LINES);
    R3Point prev_point = PointPosition(0);
    for (int ik = 1; ik < NPoints(); ik++) {
      const R3Point& current_point = PointPosition(ik);
      RNScalar dx = prev_point[0] - current_point[0];
      RNScalar dy = prev_point[1] - current_point[1];
      RNScalar dd = dx*dx + dy*dy;
      if (dd < max_edge_length_squared) {
        glColor3fv(color[(int) (0.1 * (NPoints()-ik)) % ncolors]);
        R3LoadPoint(prev_point);
        R3LoadPoint(current_point);
      }
      prev_point = current_point;
    }
    glEnd();
  }
  else {
    // Draw points
    glBegin(GL_POINTS);
    for (int i = 0; i < NPoints(); i++) 
      R3LoadPoint(PointPosition(i));
    glEnd();
  }
}



void GSVScanline::
Print(FILE *fp, const char *prefix, const char *suffix) const
{
  // Check fp
  if (!fp) fp = stdout;

  // Print scanline header
  if (prefix) fprintf(fp, "%s", prefix);
  fprintf(fp, "Scanline %d:", scan_index);
  if (suffix) fprintf(fp, "%s", suffix);
  fprintf(fp, "\n");

  // Add indentation to prefix
  char indented_prefix[1024];
  sprintf(indented_prefix, "%s  ", prefix);

  // Print points
  // ???
}



////////////////////////////////////////////////////////////////////////
// Update functions
////////////////////////////////////////////////////////////////////////

RNBoolean GSVScanline::
DoesBBoxNeedUpdate(void) const
{
  // Return whether bounding box needs update
  if (bbox.XMin() == FLT_MAX) return TRUE;
  else return FALSE;
}



void GSVScanline::
UpdateBBox(void) 
{
  // Just checking
  assert(read_count > 0);
  assert(points);

  // Update bounding box
  bbox = R3null_box;
  for (int i = 0; i < NPoints(); i++) {
    bbox.Union(PointPosition(i));
  }
}



void GSVScanline::
InvalidateBBox(void) 
{
  // Invalidate bounding box
  if (bbox.XMin() != FLT_MAX) {
    bbox = R3Box(FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX);
    if (scan) scan->InvalidateBBox();
  }
}



////////////////////////////////////////////////////////////////////////
//Input/Output functions
////////////////////////////////////////////////////////////////////////

int GSVScanline::
ReadPoints(FILE *fp, RNBoolean seek) 
{
  // Check read counter
  if (++read_count > 1) return 1;
  assert(!points);

  // Seek within file
  if (seek && (file_offset > 0)) {
    if (!RNFileSeek(fp, file_offset, RN_FILE_SEEK_SET)) {
      fprintf(stderr, "Unable to seek to scanline points\n");
      return 0;
    }
  }
  else {
    // Update file offset
    file_offset = RNFileTell(fp);
  }

  // Read number of points 
  if (fread(&npoints, sizeof(int), 1, fp) != (unsigned int) 1) {
    fprintf(stderr, "Unable to read number of scanline points\n");
    return 0;
  }

  // Check number of points
  if (npoints > 0) {
    // Allocate points 
    points = new R3Point [ npoints ];
    if (!points) {
      fprintf(stderr, "Unable to allocate scanline points\n");
      return 0;
    }
    
    // Read points
    if (fread(points, sizeof(R3Point), npoints, fp) != (unsigned int) npoints) {
      fprintf(stderr, "Unable to read scanline points\n");
      return 0;
    }
    
    // Update bounding box
    bbox = R3null_box;
    for (int i = 0; i < npoints; i++) {
      bbox.Union(points[i]);
    }
  }
  else {
    // No points
    points = NULL;
    bbox = R3null_box;
  }

  // Return success
  return 1;
}



int GSVScanline::
WritePoints(FILE *fp, RNBoolean seek) 
{
  // Seek within file
  if (seek && (file_offset > 0)) {
    // Seek to file offset
    if (!RNFileSeek(fp, file_offset, RN_FILE_SEEK_SET)) {
      fprintf(stderr, "Unable to seek to scanline points\n");
      return 0;
    }
  }
  else {
    // Update file offset
    file_offset = RNFileTell(fp);
  }

  // Write number of points 
  if (fwrite(&npoints, sizeof(int), 1, fp) != (unsigned int) 1) {
    fprintf(stderr, "Unable to write number of scanline points\n");
    return 0;
  }

  // Write points
  if (npoints > 0) {
    if (fwrite(points, sizeof(R3Point), npoints, fp) != (unsigned int) npoints) {
      fprintf(stderr, "Unable to write scanline points\n");
      return 0;
    }
  }

  // Return success
  return 1;
}



int GSVScanline::
ReleasePoints(void) 
{
  // Check read counter
  if (--read_count > 0) return 1;
  assert(points);

  // Release points
  delete [] points;
  points = NULL;

  // Return success
  return 1;
}



