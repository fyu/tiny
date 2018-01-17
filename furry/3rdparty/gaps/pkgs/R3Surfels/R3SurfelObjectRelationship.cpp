/* Source file for the R3 surfel object relationship class */



////////////////////////////////////////////////////////////////////////
// INCLUDE FILES
////////////////////////////////////////////////////////////////////////

#include "R3Surfels/R3Surfels.h"



////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS/DESTRUCTORS
////////////////////////////////////////////////////////////////////////

R3SurfelObjectRelationship::
R3SurfelObjectRelationship(int type, const RNArray<R3SurfelObject *>& objects, RNScalar *operands)
  : scene(NULL),
    scene_index(-1),
    objects(objects),
    operands(NULL),
    type(type)
{
  // Copy operands
  int noperands = NOperands();
  if (noperands > 0) {
    this->operands = new RNScalar [ noperands ];
    for (int i = 0; i < noperands; i++) {
      this->operands[i] = operands[i];
    }
  }
}



R3SurfelObjectRelationship::
R3SurfelObjectRelationship(int type, R3SurfelObject *object0, R3SurfelObject *object1, RNScalar *operands)
  : scene(NULL),
    scene_index(-1),
    objects(),
    operands(NULL),
    type(type)
{
  // Insert objects
  objects.Insert(object0);
  objects.Insert(object1);

  // Copy operands
  int noperands = NOperands();
  if (noperands > 0) {
    this->operands = new RNScalar [ noperands ];
    for (int i = 0; i < noperands; i++) {
      this->operands[i] = operands[i];
    }
  }
}



R3SurfelObjectRelationship::
R3SurfelObjectRelationship(const R3SurfelObjectRelationship& relationship)
  : scene(NULL),
    scene_index(-1),
    objects(relationship.objects),
    operands(NULL),
    type(relationship.type)
{
  // Copy operands
  int noperands = NOperands();
  if (noperands > 0) {
    this->operands = new RNScalar [ noperands ];
    for (int i = 0; i < noperands; i++) {
      this->operands[i] = relationship.operands[i];
    }
  }
}



R3SurfelObjectRelationship::
~R3SurfelObjectRelationship(void)
{
  // Remove from scene
  if (scene) scene->RemoveObjectRelationship(this);

  // Delete operands
  if (operands) delete [] operands;
}






////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCDTIONS
////////////////////////////////////////////////////////////////////////

int R3SurfelObjectRelationship::
NOperands(int type)
{
  // Return number of operands
  switch (type) {
  case R3_SURFEL_OBJECT_SIMILARITY_RELATIONSHIP: 
    // operands[0] = similarity (1 = similar, 0 = dissimilar)
    return 1;
  }

  // Default case
  return 0;
}



