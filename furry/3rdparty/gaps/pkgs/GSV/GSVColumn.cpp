// Source file for GSV column class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV.h"



////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

GSVColumn::
GSVColumn(void)
  : image(NULL),
    image_index(-1),
    pose(0,0,0,0,0,0,0),
    timestamp(-1),
    data(NULL)
{
}



GSVColumn::
GSVColumn(const GSVPose& pose, RNScalar timestamp)
  : image(NULL),
    image_index(-1),
    pose(pose),
    timestamp(timestamp),
    data(NULL)
{
}



GSVColumn::
~GSVColumn(void)
{
  // Remove column from camera image
  if (image) image->RemoveColumn(this);
}



////////////////////////////////////////////////////////////////////////
// Access functions
////////////////////////////////////////////////////////////////////////

GSVScene *GSVColumn::
Scene(void) const
{
  // Return scene
  if (!image) return NULL;
  return image->Scene();
}



GSVRun *GSVColumn::
Run(void) const
{
  // Return run
  if (!image) return NULL;
  return image->Run();
}



GSVSegment *GSVColumn::
Segment(void) const
{
  // Return segment
  if (!image) return NULL;
  return image->Segment();
}



GSVCamera *GSVColumn::
Camera(void) const
{
  // Return segment
  if (!image) return NULL;
  return image->Camera();
}



////////////////////////////////////////////////////////////////////////
// Display functions
////////////////////////////////////////////////////////////////////////

void GSVColumn::
Draw(RNFlags flags) const
{
  // Draw pixels
  // ???
}



void GSVColumn::
Print(FILE *fp, const char *prefix, const char *suffix) const
{
  // Check fp
  if (!fp) fp = stdout;

  // Print column header
  if (prefix) fprintf(fp, "%s", prefix);
  fprintf(fp, "Column %d:", image_index);
  if (suffix) fprintf(fp, "%s", suffix);
  fprintf(fp, "\n");

  // Add indentation to prefix
  char indented_prefix[1024];
  sprintf(indented_prefix, "%s  ", prefix);

  // Print pixels
  // ???
}




