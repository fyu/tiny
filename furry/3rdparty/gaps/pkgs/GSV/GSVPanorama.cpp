// Source file for GSV panorama class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV.h"



////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

GSVPanorama::
GSVPanorama(void)
  : segment(NULL),
    segment_index(-1),
    images(),
    viewpoint(FLT_MAX, FLT_MAX, FLT_MAX),
    timestamp(-1), 
    data(NULL)
{
}



GSVPanorama::
~GSVPanorama(void)
{
  // Remove panorama from segment
  if (segment) segment->RemovePanorama(this);

  // Delete images 
  while (NImages() > 0) delete Image(0);
}



////////////////////////////////////////////////////////////////////////
// Access functions
////////////////////////////////////////////////////////////////////////

GSVScene *GSVPanorama::
Scene(void) const
{
  // Return scene
  if (!segment) return NULL;
  return segment->Scene();
}



GSVRun *GSVPanorama::
Run(void) const
{
  // Return run
  if (!segment) return NULL;
  return segment->Run();
}



int GSVPanorama::
SceneIndex(void) const
{
  // Get convenient variables
  GSVRun *run = Run();
  if (!run) return -1;
  GSVScene *scene = run->Scene();
  if (!scene) return -1;

  // Compute scene index
  int scene_index = RunIndex();
  for (int i = 0; i < run->SceneIndex(); i++) {
    GSVRun *r = scene->Run(i);
    scene_index += r->NPanoramas();
  }

  // Return scene index
  return scene_index;
}



int GSVPanorama::
RunIndex(void) const
{
  // Get convenient variables
  GSVRun *run = Run();
  if (!run) return -1;

  // Compute run index
  int run_index = SegmentIndex();
  for (int i = 0; i < segment->RunIndex(); i++) {
    GSVSegment *s = run->Segment(i);
    run_index += s->NPanoramas();
  }

  // Return run index
  return run_index;
}



////////////////////////////////////////////////////////////////////////
// Manipulation functions
////////////////////////////////////////////////////////////////////////

void GSVPanorama::
InsertImage(GSVImage *image)
{
  // Just checking
  assert(image->panorama_index == -1);
  assert(image->panorama == NULL);

  // Insert image
  image->panorama = this;
  image->panorama_index = images.NEntries();
  images.Insert(image);

  // Invalidate stuff
  viewpoint.Reset(FLT_MAX, FLT_MAX, FLT_MAX);
  timestamp = -1;
}



void GSVPanorama::
RemoveImage(GSVImage *image)
{
  // Just checking
  assert(image->panorama_index >= 0);
  assert(image->panorama_index < images.NEntries());
  assert(image->panorama == this);

  // Remove image
  RNArrayEntry *entry = images.KthEntry(image->panorama_index);
  GSVImage *tail = images.Tail();
  tail->panorama_index = image->panorama_index;
  images.EntryContents(entry) = tail;
  images.RemoveTail();
  image->panorama_index = -1;
  image->panorama = NULL;

  // Invalidate stuff
  viewpoint.Reset(FLT_MAX, FLT_MAX, FLT_MAX);
  timestamp = -1;
}



////////////////////////////////////////////////////////////////////////
// Display functions
////////////////////////////////////////////////////////////////////////

void GSVPanorama::
Draw(RNFlags flags) const
{
  // Draw images
  for (int i = 0; i < NImages(); i++) {
    GSVImage *image = Image(i);
    image->Draw(flags);
  }
}



void GSVPanorama::
Print(FILE *fp, const char *prefix, const char *suffix) const
{
  // Check fp
  if (!fp) fp = stdout;

  // Print panorama header
  if (prefix) fprintf(fp, "%s", prefix);
  fprintf(fp, "Panorama %d:", segment_index);
  if (suffix) fprintf(fp, "%s", suffix);
  fprintf(fp, "\n");

  // Add indentation to prefix
  char indented_prefix[1024];
  sprintf(indented_prefix, "%s  ", prefix);

   // Print images
  for (int i = 0; i < NImages(); i++) {
    GSVImage *image = Image(i);
    image->Print(fp, indented_prefix, suffix);
  }
}



////////////////////////////////////////////////////////////////////////
// Update functions
////////////////////////////////////////////////////////////////////////

void GSVPanorama::
UpdateViewpoint(void) 
{
  // Update viewpoint
  viewpoint = R3zero_point;
  if (NImages() > 0) {
    for (int i = 0; i < NImages(); i++) 
      viewpoint += Image(i)->Pose().Viewpoint();
    viewpoint /= NImages();
  }
}


void GSVPanorama::
UpdateTimestamp(void) 
{
  // Update timestamp
  timestamp = 0;
  if (NImages() > 0) {
    for (int i = 0; i < NImages(); i++) 
      timestamp += Image(i)->Timestamp();
    timestamp /= NImages();
  }
}



