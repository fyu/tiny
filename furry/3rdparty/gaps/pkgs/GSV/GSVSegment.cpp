// Source file for GSV segment class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV.h"



////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

GSVSegment::
GSVSegment(void)
  : scene_index(-1),
    run(NULL),
    run_index(-1),
    scans(), 
    tapestries(), 
    panoramas(), 
    bbox(FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX),
    data(NULL)
{
}



GSVSegment::
~GSVSegment(void)
{
  // Remove segment from scene
  assert(scene_index == -1);

  // Remove segment from run
  if (run) run->RemoveSegment(this);

  // Delete laser tapestries
  while (NScans() > 0) delete Scan(0);

  // Delete camera tapestries
  while (NTapestries() > 0) delete Tapestry(0);

  // Delete camera panoramas
  while (NPanoramas() > 0) delete Panorama(0);
}



////////////////////////////////////////////////////////////////////////
// Access functions
////////////////////////////////////////////////////////////////////////

GSVScene *GSVSegment::
Scene(void) const
{
  // Return scene
  if (!run) return NULL;
  return run->Scene();
}



int GSVSegment::
SceneIndex(void) const
{
  // Get convenient variables
  GSVScene *scene = Scene();
  if (!scene) return -1;

  // Compute scene index
  int scene_index = RunIndex();
  for (int i = 0; i < run->SceneIndex(); i++) {
    GSVRun *r = scene->Run(i);
    scene_index += r->NSegments();
  }

  // Return scene index
  return scene_index;
}



int GSVSegment::
NImages(void) const
{
  // Return number of images in segment
  int count = 0;
  for (int i = 0; i < NTapestries(); i++) {
    GSVTapestry *tapestry = Tapestry(i);
    count += tapestry->NImages();
  }
  return count;
}



GSVImage *GSVSegment::
Image(int image_index) const
{
  // Return image
  for (int i = 0; i < NTapestries(); i++) {
    GSVTapestry *tapestry = Tapestry(i);
    if (image_index >= tapestry->NImages()) image_index -= tapestry->NImages();
    else return tapestry->Image(image_index);
  }

  // Index too high
  return NULL;
}




int GSVSegment::
NScanlines(void) const
{
  // Return number of scanlines in segment
  int count = 0;
  for (int i = 0; i < NScans(); i++) {
    GSVScan *scan = Scan(i);
    count += scan->NScanlines();
  }
  return count;
}



GSVScanline *GSVSegment::
Scanline(int scanline_index) const
{
  // Return scanline
  for (int i = 0; i < NScans(); i++) {
    GSVScan *scan = Scan(i);
    if (scanline_index >= scan->NScanlines()) scanline_index -= scan->NScanlines();
    else return scan->Scanline(scanline_index);
  }

  // Index too high
  return NULL;
}




GSVPanorama *GSVSegment::
FindPanoramaBeforeTimestamp(RNScalar timestamp, int imin, int imax) const
{
  // Binary search
  int i = (imin + imax) / 2;
  if (i == imin) return Panorama(imin);
  assert(i < imax);
  GSVPanorama *panorama = Panorama(i);
  RNScalar t = panorama->Timestamp();
  if (t > timestamp) return FindPanoramaBeforeTimestamp(timestamp, imin, i);
  else if (t < timestamp) return FindPanoramaBeforeTimestamp(timestamp, i, imax);
  else return panorama;
}



GSVPanorama *GSVSegment::
FindPanoramaBeforeTimestamp(RNScalar timestamp) const
{
  // Binary search
  if (NPanoramas() == 0) return NULL;
  if (timestamp <= Panorama(0)->Timestamp()) return Panorama(0);
  if (timestamp >= Panorama(NPanoramas()-1)->Timestamp()) return Panorama(NPanoramas()-1);
  return FindPanoramaBeforeTimestamp(timestamp, 0, NPanoramas()-1);
}



////////////////////////////////////////////////////////////////////////
// Manipulation functions
////////////////////////////////////////////////////////////////////////

void GSVSegment::
InsertScan(GSVScan *scan)
{
  // Just checking
  assert(scan->segment_index == -1);
  assert(scan->segment == NULL);

  // Insert laser scan
  scan->segment = this;
  scan->segment_index = scans.NEntries();
  scans.Insert(scan);

  // Invalidate bounding box
  InvalidateBBox();
}



void GSVSegment::
RemoveScan(GSVScan *scan)
{
  // Just checking
  assert(scan->segment_index >= 0);
  assert(scan->segment_index < scans.NEntries());
  assert(scan->segment == this);

  // Remove laser scan
  RNArrayEntry *entry = scans.KthEntry(scan->segment_index);
  GSVScan *tail = scans.Tail();
  tail->segment_index = scan->segment_index;
  scans.EntryContents(entry) = tail;
  scans.RemoveTail();
  scan->segment_index = -1;
  scan->segment = NULL;

  // Invalidate bounding box
  InvalidateBBox();
}



void GSVSegment::
InsertTapestry(GSVTapestry *tapestry)
{
  // Just checking
  assert(tapestry->segment_index == -1);
  assert(tapestry->segment == NULL);

  // Insert camera tapestry
  tapestry->segment = this;
  tapestry->segment_index = tapestries.NEntries();
  tapestries.Insert(tapestry);

  // Invalidate bounding box
  InvalidateBBox();
}



void GSVSegment::
RemoveTapestry(GSVTapestry *tapestry)
{
  // Just checking
  assert(tapestry->segment_index >= 0);
  assert(tapestry->segment_index < tapestries.NEntries());
  assert(tapestry->segment == this);

  // Remove camera tapestry
  RNArrayEntry *entry = tapestries.KthEntry(tapestry->segment_index);
  GSVTapestry *tail = tapestries.Tail();
  tail->segment_index = tapestry->segment_index;
  tapestries.EntryContents(entry) = tail;
  tapestries.RemoveTail();
  tapestry->segment_index = -1;
  tapestry->segment = NULL;

  // Invalidate bounding box
  InvalidateBBox();
}



void GSVSegment::
InsertPanorama(GSVPanorama *panorama)
{
  // Just checking
  assert(panorama->segment_index == -1);
  assert(panorama->segment == NULL);

  // Insert camera panorama
  panorama->segment = this;
  panorama->segment_index = panoramas.NEntries();
  panoramas.Insert(panorama);

  // Invalidate bounding box
  InvalidateBBox();
}



void GSVSegment::
RemovePanorama(GSVPanorama *panorama)
{
  // Just checking
  assert(panorama->segment_index >= 0);
  assert(panorama->segment_index < panoramas.NEntries());
  assert(panorama->segment == this);

  // Remove camera panorama
  RNArrayEntry *entry = panoramas.KthEntry(panorama->segment_index);
  GSVPanorama *tail = panoramas.Tail();
  tail->segment_index = panorama->segment_index;
  panoramas.EntryContents(entry) = tail;
  panoramas.RemoveTail();
  panorama->segment_index = -1;
  panorama->segment = NULL;

  // Invalidate bounding box
  InvalidateBBox();
}



////////////////////////////////////////////////////////////////////////
// Display functions
////////////////////////////////////////////////////////////////////////

void GSVSegment::
Draw(RNFlags flags) const
{
  // Draw scans
  for (int i = 0; i < NScans(); i++) {
    GSVScan *scan = Scan(i);
    scan->Draw(flags);
  }

  // Draw tapestries
  for (int i = 0; i < NTapestries(); i++) {
    GSVTapestry *tapestry = Tapestry(i);
    tapestry->Draw(flags);
  }
}



void GSVSegment::
Print(FILE *fp, const char *prefix, const char *suffix) const
{
  // Check fp
  if (!fp) fp = stdout;

  // Print segment 
  if (prefix) fprintf(fp, "%s", prefix);
  fprintf(fp, "Segment %d:", run_index);
  if (suffix) fprintf(fp, "%s", suffix);
  fprintf(fp, "\n");

  // Add indentation to prefix
  char indented_prefix[1024];
  sprintf(indented_prefix, "%s  ", prefix);

  // Print laser tapestries
  for (int i = 0; i < NScans(); i++) {
    GSVScan *scan = Scan(i);
    scan->Print(fp, indented_prefix, suffix);
  }

  // Print camera tapestries
  for (int i = 0; i < NTapestries(); i++) {
    GSVTapestry *tapestry = Tapestry(i);
    tapestry->Print(fp, indented_prefix, suffix);
  }

  // Print camera panoramas
  for (int i = 0; i < NPanoramas(); i++) {
    GSVPanorama *panorama = Panorama(i);
    panorama->Print(fp, indented_prefix, suffix);
  }
}



////////////////////////////////////////////////////////////////////////
// Update functions
////////////////////////////////////////////////////////////////////////

void GSVSegment::
UpdateBBox(void) 
{
  // Update bounding box
  bbox = R3null_box;
  for (int i = 0; i < NScans(); i++) bbox.Union(Scan(i)->BBox());
  for (int i = 0; i < NTapestries(); i++) bbox.Union(Tapestry(i)->BBox());
  for (int i = 0; i < NPanoramas(); i++) bbox.Union(Panorama(i)->Viewpoint());
}



void GSVSegment::
InvalidateBBox(void) 
{
  // Invalidate bounding box
  if (bbox.XMin() != FLT_MAX) {
    bbox = R3Box(FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX);
    if (run) run->InvalidateBBox();
  }
}




