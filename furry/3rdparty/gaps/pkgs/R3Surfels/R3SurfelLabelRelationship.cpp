/* Source file for the R3 surfel label relationship class */



////////////////////////////////////////////////////////////////////////
// INCLUDE FILES
////////////////////////////////////////////////////////////////////////

#include "R3Surfels/R3Surfels.h"



////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS/DESTRUCTORS
////////////////////////////////////////////////////////////////////////

R3SurfelLabelRelationship::
R3SurfelLabelRelationship(int type, const RNArray<R3SurfelLabel *>& labels, RNScalar *operands)
  : scene(NULL),
    scene_index(-1),
    labels(labels),
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



R3SurfelLabelRelationship::
R3SurfelLabelRelationship(int type, R3SurfelLabel *label0, R3SurfelLabel *label1, RNScalar *operands)
  : scene(NULL),
    scene_index(-1),
    labels(),
    operands(NULL),
    type(type)
{
  // Insert labels
  labels.Insert(label0);
  labels.Insert(label1);

  // Copy operands
  int noperands = NOperands();
  if (noperands > 0) {
    this->operands = new RNScalar [ noperands ];
    for (int i = 0; i < noperands; i++) {
      this->operands[i] = operands[i];
    }
  }
}



R3SurfelLabelRelationship::
R3SurfelLabelRelationship(const R3SurfelLabelRelationship& relationship)
  : scene(NULL),
    scene_index(-1),
    labels(relationship.labels),
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



R3SurfelLabelRelationship::
~R3SurfelLabelRelationship(void)
{
  // Remove from scene
  if (scene) scene->RemoveLabelRelationship(this);

  // Delete operands
  if (operands) delete [] operands;
}






////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCDTIONS
////////////////////////////////////////////////////////////////////////

int R3SurfelLabelRelationship::
NOperands(int type) 
{
  // Return number of operands
  return R3SurfelObjectRelationship::NOperands(type);
}



