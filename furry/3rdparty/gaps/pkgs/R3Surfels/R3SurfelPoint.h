/* Include file for the R3 surfel point class */



/* Class definition */

class R3SurfelPoint {
public:
  // Constructor functions
  R3SurfelPoint(void);
  R3SurfelPoint(const R3SurfelPoint& surfel);
  R3SurfelPoint(R3SurfelBlock *block, const R3Surfel *surfel);
  ~R3SurfelPoint(void);

  // Position property functions
  RNCoord X(void) const;
  RNCoord Y(void) const;
  RNCoord Z(void) const;
  RNCoord Coord(int dimension) const;
  R3Point Position(void) const;

  // Color property functions
  unsigned char R(void) const;
  unsigned char G(void) const;
  unsigned char B(void) const;
  const unsigned char *Color(void) const;
  RNRgb Rgb(void) const;

  // Other property functions
  RNBoolean IsAerial(void) const;
  RNBoolean IsMarked(void) const;

  // Access functions
  R3SurfelBlock *Block(void) const;
  int BlockIndex(void) const;
  const R3Surfel *Surfel(void) const;

  // Copy and reset functions
  void Copy(const R3SurfelPoint *point);
  void Reset(R3SurfelBlock *block, const R3Surfel *surfel);
  R3SurfelPoint& operator=(const R3SurfelPoint& point);

  // Manipulation functions
  void SetMark(RNBoolean mark);

  // Draw functions
  void Draw(RNFlags flags = R3_SURFEL_DEFAULT_DRAW_FLAGS) const;

private:
  R3SurfelBlock *block;
  const R3Surfel *surfel;
};



/* Inline functions */

inline RNCoord R3SurfelPoint::
X(void) const
{
  // Return X coordinate of surfel point in global coordinate system
  assert(block && surfel);
 return surfel->X() + block->Origin().X();
}



inline RNCoord R3SurfelPoint::
Y(void) const
{
  // Return Y coordinate of surfel point in global coordinate system
  assert(block && surfel);
  return surfel->Y() + block->Origin().Y();
}



inline RNCoord R3SurfelPoint::
Z(void) const
{
  // Return Z coordinate of surfel point in global coordinate system
  assert(block && surfel);
  return surfel->Z() + block->Origin().Z();
}



inline RNCoord R3SurfelPoint::
Coord(int dimension) const
{
  // Return coordinate of surfel point in global coordinate system
  assert(block && surfel);
  return surfel->Coord(dimension) + block->Origin().Coord(dimension);
}



inline R3Point R3SurfelPoint::
Position(void) const
{
  // Return position of surfel point in global coordinate system
  return R3Point(X(), Y(), Z());
}



inline unsigned char R3SurfelPoint::
R(void) const
{
  // Return red component of color
  return surfel->R();
}



inline unsigned char R3SurfelPoint::
G(void) const
{
  // Return green component of color
  return surfel->G();
}



inline unsigned char R3SurfelPoint::
B(void) const
{
  // Return blue component of color
  return surfel->B();
}



inline const unsigned char *R3SurfelPoint::
Color(void) const
{
  // Return pointer to color
  return surfel->Color();
}



inline RNRgb R3SurfelPoint::
Rgb(void) const
{
  // Return RGB
  return surfel->Rgb();
}



inline RNBoolean R3SurfelPoint::
IsAerial(void) const
{
  // Return whether point was captured with aerial scanner
  return surfel->IsAerial();
}



inline RNBoolean R3SurfelPoint::
IsMarked(void) const
{
  // Return whether point is marked
  return surfel->IsMarked();
}



inline R3SurfelBlock *R3SurfelPoint::
Block(void) const
{
  // Return block
  return block;
}



inline int R3SurfelPoint::
BlockIndex(void) const
{
  // Return index of surfel within block
  return surfel - block->Surfels();
}



inline const R3Surfel *R3SurfelPoint::
Surfel(void) const
{
  // Return surfel 
  return surfel;
}



