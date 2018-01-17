// Source file for GSV laser class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV.h"



////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

GSVLaser::
GSVLaser(void)
  : run(NULL),
    run_index(-1),
    scans(),
    bbox(FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX),
    data(NULL)
{
}



GSVLaser::
~GSVLaser(void)
{
  // Remove laser from run
  if (run) run->RemoveLaser(this);

  // Remove scans from this laser
  while (NScans() > 0) RemoveScan(Scan(0));
}



////////////////////////////////////////////////////////////////////////
// Access functions
////////////////////////////////////////////////////////////////////////

GSVScene *GSVLaser::
Scene(void) const
{
  // Return scene
  if (!run) return NULL;
  return run->Scene();
}



int GSVLaser::
SceneIndex(void) const
{
  // Get convenient variables
  GSVScene *scene = Scene();
  if (!scene) return -1;

  // Compute scene index
  int scene_index = RunIndex();
  for (int i = 0; i < run->SceneIndex(); i++) {
    GSVRun *r = scene->Run(i);
    scene_index += r->NLasers();
  }

  // Return scene index
  return scene_index;
}



////////////////////////////////////////////////////////////////////////
// Manipulation functions
////////////////////////////////////////////////////////////////////////

void GSVLaser::
InsertScan(GSVScan *scan)
{
  // Just checking
  assert(scan->laser_index == -1);
  assert(scan->laser == NULL);

  // Insert scan
  scan->laser = this;
  scan->laser_index = scans.NEntries();
  scans.Insert(scan);

  // Invalidate bounding box
  InvalidateBBox();
}



void GSVLaser::
RemoveScan(GSVScan *scan)
{
  // Just checking
  assert(scan->laser_index >= 0);
  assert(scan->laser_index < scans.NEntries());
  assert(scan->laser == this);

  // Remove scan
  RNArrayEntry *entry = scans.KthEntry(scan->laser_index);
  GSVScan *tail = scans.Tail();
  tail->laser_index = scan->laser_index;
  scans.EntryContents(entry) = tail;
  scans.RemoveTail();
  scan->laser_index = -1;
  scan->laser = NULL;

  // Invalidate bounding box
  InvalidateBBox();
}



////////////////////////////////////////////////////////////////////////
// Display functions
////////////////////////////////////////////////////////////////////////

void GSVLaser::
Draw(RNFlags flags) const
{
  // Draw tapestries
  for (int i = 0; i < NScans(); i++) {
    GSVScan *scan = Scan(i);
    scan->Draw(flags);
  }
}



void GSVLaser::
Print(FILE *fp, const char *prefix, const char *suffix) const
{
  // Check fp
  if (!fp) fp = stdout;

  // Print laser header
  if (prefix) fprintf(fp, "%s", prefix);
  fprintf(fp, "Laser %d:", run_index);
  if (suffix) fprintf(fp, "%s", suffix);
  fprintf(fp, "\n");

  // Add indentation to prefix
  char indented_prefix[1024];
  sprintf(indented_prefix, "%s  ", prefix);

   // Print tapestries
  for (int i = 0; i < NScans(); i++) {
    GSVScan *scan = Scan(i);
    scan->Print(fp, indented_prefix, suffix);
  }
}



////////////////////////////////////////////////////////////////////////
// Update functions
////////////////////////////////////////////////////////////////////////

void GSVLaser::
UpdateBBox(void) 
{
  // Update bounding box
  bbox = R3null_box;
  for (int i = 0; i < NScans(); i++) {
    bbox.Union(Scan(i)->BBox());
  }
}



void GSVLaser::
InvalidateBBox(void) 
{
  // Invalidate bounding box
  if (bbox.XMin() != FLT_MAX) {
    bbox = R3Box(FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX);
    if (run) run->InvalidateBBox();
  }
}



