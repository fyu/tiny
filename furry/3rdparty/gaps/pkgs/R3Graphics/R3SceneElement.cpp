/* Source file for the R3 shape class */



/* Include files */

#include "R3Graphics.h"



/* Member functions */

R3SceneElement::
R3SceneElement(R3Material *material)
  : node(NULL),
    material(material),
    shapes(),
    opengl_id(0),
    bbox(FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX)
{
}



R3SceneElement::
~R3SceneElement(void)
{
  // Delete display list
  if (opengl_id > 0) glDeleteLists(opengl_id, 1); 

  // Remove from node
  if (node) node->RemoveElement(this);
}



void R3SceneElement::
SetMaterial(R3Material *material) 
{
  // Set material
  this->material = material;
}



void R3SceneElement::
InsertShape(R3Shape *shape) 
{
  // Insert shape
  shapes.Insert(shape);

  // Invalidate bounding box
  InvalidateBBox();
}



void R3SceneElement::
RemoveShape(R3Shape *shape) 
{
  // Remove shape
  shapes.Remove(shape);

  // Invalidate bounding box
  InvalidateBBox();
}



RNBoolean R3SceneElement::
Intersects(const R3Ray& ray, R3Shape **hit_shape,
  R3Point *hit_point, R3Vector *hit_normal, RNScalar *hit_t) const
{
  // Variables
  RNScalar closest_t = FLT_MAX;
  R3Point point;
  R3Vector normal;
  RNScalar t;

  // Intersect with shapes
  for (int i = 0; i < NShapes(); i++) {
    R3Shape *shape = Shape(i);
    if (shape->Intersects(ray, &point, &normal, &t)) {
      if (t < closest_t) {
        if (hit_shape) *hit_shape = shape;
        if (hit_point) *hit_point = point;
        if (hit_normal) *hit_normal = normal;
        if (hit_t) *hit_t = t;
        closest_t = t;
      }
    }
  }

  // Return whether hit any shape
  return (closest_t == FLT_MAX) ? FALSE : TRUE;
}



void R3SceneElement::
Draw(const R3DrawFlags draw_flags) const
{
  // Check shapes
  if (NShapes() == 0) return;

  // Draw material
  if (material) material->Draw();

  // Draw shapes
  if (0 && (draw_flags == R3_DEFAULT_DRAW_FLAGS)) {
    // Create display list
    if (opengl_id == 0) {
      // Begin display list
      R3SceneElement *element = (R3SceneElement *) this;
      element->opengl_id = glGenLists(1);
      glNewList(opengl_id, GL_COMPILE);
      
      // Draw shapes
      for (int i = 0; i < NShapes(); i++) {
        R3Shape *shape = Shape(i);
        shape->Draw(draw_flags);
      }

      // End display list
      glEndList();
    }

    // Call display list
    glCallList(opengl_id);
  }
  else {
    // Draw shapes with unusual parameters
    for (int i = 0; i < NShapes(); i++) {
      R3Shape *shape = Shape(i);
      shape->Draw(draw_flags);
    }
  }
}



void R3SceneElement::
UpdateBBox(void)
{
  // Update bounding box
  bbox = R3null_box;
  for (int i = 0; i < NShapes(); i++) {
    R3Shape *shape = Shape(i);
    bbox.Union(shape->BBox());
  }
}



void R3SceneElement::
InvalidateBBox(void)
{
  // Invalidate bounding box
  bbox[0][0] = FLT_MAX;

  // Invalidate node bounding box
  if (node) node->InvalidateBBox();
}

