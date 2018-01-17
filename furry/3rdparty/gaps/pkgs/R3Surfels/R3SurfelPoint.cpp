/* Source file for the R3 surfel point class */



/* Include files */

#include "R3Surfels/R3Surfels.h"



/* Public functions */

R3SurfelPoint::
R3SurfelPoint(void)
  : block(NULL),
    surfel(NULL)
{
}



R3SurfelPoint::
R3SurfelPoint(const R3SurfelPoint& point)
  : block(point.block),
    surfel(point.surfel)
{
  // Update block reference count
  if (block && block->Database()) block->Database()->ReadBlock(block);  
}



R3SurfelPoint::
R3SurfelPoint(R3SurfelBlock *block, const R3Surfel *surfel)
  : block(block),
    surfel(surfel)
{
  // Update block reference count
  if (block && block->Database()) block->Database()->ReadBlock(block);  
}



R3SurfelPoint::
~R3SurfelPoint(void)
{
  // Update block reference count
  if (block && block->Database()) block->Database()->ReleaseBlock(block);  
}



void R3SurfelPoint::
Copy(const R3SurfelPoint *point)
{
  // Update block reference counts
  if (this->block != point->block) {
    if (this->block && this->block->Database()) this->block->Database()->ReleaseBlock(this->block);
    if (point->block && point->block->Database()) point->block->Database()->ReadBlock(point->block);
  }

  // Copy block and surfel
  this->block = point->block;
  this->surfel = point->surfel;
}



void R3SurfelPoint::
Reset(R3SurfelBlock *block, const R3Surfel *surfel)
{
  // Update block reference counts
  if (this->block != block) {
    if (this->block && this->block->Database()) this->block->Database()->ReleaseBlock(this->block);
    if (block && block->Database()) block->Database()->ReadBlock(block);
  }

  // Copy block and surfel
  this->block = block;
  this->surfel = surfel;
}



R3SurfelPoint& R3SurfelPoint::
operator=(const R3SurfelPoint& point)
{
  // Copy the point
  Copy(&point);
  return *this;
}



void R3SurfelPoint::
SetMark(RNBoolean mark)
{
  // Set mark 
  R3Surfel *s = (R3Surfel *) this->surfel;
  s->SetMark(mark);
}



void R3SurfelPoint::
Draw(RNFlags flags) const
{
  // Draw surfel
  glBegin(GL_POINTS);
  if (flags[R3_SURFEL_COLOR_DRAW_FLAG]) glColor3ubv(Color());
  glVertex3dv(Position().Coords());
  glEnd();
}










