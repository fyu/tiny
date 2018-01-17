/* Source file for the R3 surfel class */



/* Include files */

#include "R3Surfels/R3Surfels.h"



// Flag definitions

#define R3_SURFEL_AERIAL_FLAG  0x02
#define R3_SURFEL_MARKED_FLAG  0x80



/* Public functions */

R3Surfel::
R3Surfel(void)
  : flags(0)
{
  // Set everything
  SetCoords(0, 0, 0);
  SetColor(0, 0, 0);
}



R3Surfel::
R3Surfel(float x, float y, float z, unsigned char r, unsigned char g, unsigned char b, RNBoolean aerial)
  : flags(0)
{
  // Set everything
  SetCoords(x, y, z);
  SetColor(r, g, b);
  SetAerial(aerial);
}



RNBoolean R3Surfel::
IsAerial(void) const
{
  // Return whether point was captured with aerial scanner
  return flags & R3_SURFEL_AERIAL_FLAG;
}



void R3Surfel::
SetAerial(RNBoolean aerial)
{
  // Set whether point was captured with aerial scanner
  if (aerial) flags |= R3_SURFEL_AERIAL_FLAG;
  else flags &= ~R3_SURFEL_AERIAL_FLAG;
}



RNBoolean R3Surfel::
IsMarked(void) const
{
  // Return whether surfel is marked (useful for set and traversal operations)
  return flags & R3_SURFEL_MARKED_FLAG;
}



void R3Surfel::
SetMark(RNBoolean mark)
{
  // Set whether surfel is marked (useful for set and traversal operations)
  if (mark) flags |= R3_SURFEL_MARKED_FLAG;
  else flags &= ~R3_SURFEL_MARKED_FLAG;
}



