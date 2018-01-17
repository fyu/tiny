/* Include file for the R3 surfel class */



/* Class definition */

class R3Surfel {
public:
  // Constructor functions
  R3Surfel(void);
  R3Surfel(float x, float y, float z, unsigned char r = 0, unsigned char g = 0, unsigned char b = 0, RNBoolean aerial = FALSE);

  // Position property functions
  // NOTE THAT THESE COORDINATES ARE RELATIVE TO THE BLOCK ORIGIN
  // TO GET THE WORLD COORDINATES, YOU MUST ADD THE BLOCK ORIGIN
  float X(void) const;
  float Y(void) const;
  float Z(void) const;
  float Coord(int dimension) const;
  const float *Coords(void) const;
 
  // Color property functions
  unsigned char R(void) const;
  unsigned char G(void) const;
  unsigned char B(void) const;
  const unsigned char *Color(void) const;
  RNRgb Rgb(void) const;

  // Other property functions
  RNBoolean IsMarked(void) const;
  RNBoolean IsAerial(void) const;
  RNBoolean IsTerrestrial(void) const;

  // Manipulation functions
  void SetCoords(float x, float y, float z);
  void SetCoords(const float *xyz);
  void SetColor(unsigned char r, unsigned char g, unsigned char b);
  void SetColor(const unsigned char *rgb);
  void SetColor(const RNRgb& rgb);
  void SetAerial(RNBoolean aerial = TRUE);
  void SetMark(RNBoolean mark = TRUE);

  // Draw functions
  void Draw(RNFlags flags = R3_SURFEL_DEFAULT_DRAW_FLAGS) const;


////////////////////////////////////////////////////////////////////////
// INTERNAL STUFF BELOW HERE
////////////////////////////////////////////////////////////////////////

private:
  // Internal data
  float position[3];
  unsigned char color[3];
  unsigned char flags;
};



////////////////////////////////////////////////////////////////////////
// INLINE FUNCTION DEFINITIONS
////////////////////////////////////////////////////////////////////////

inline float R3Surfel::
X(void) const
{
  // Return X coordinate
  return position[0];
}



inline float R3Surfel::
Y(void) const
{
  // Return Y coordinate
  return position[1];
}



inline float R3Surfel::
Z(void) const
{
  // Return Z coordinate
  return position[2];
}



inline float R3Surfel::
Coord(int dimension) const
{
  // Return coordinate
  return position[dimension];
}



inline const float *R3Surfel::
Coords(void) const
{
  // Return pointer to position
  return position;
}



inline unsigned char R3Surfel::
R(void) const
{
  // Return red component of color
  return color[0];
}



inline unsigned char R3Surfel::
G(void) const
{
  // Return green component of color
  return color[1];
}



inline unsigned char R3Surfel::
B(void) const
{
  // Return blue component of color
  return color[2];
}



inline const unsigned char *R3Surfel::
Color(void) const
{
  // Return pointer to color
  return color;
}



inline RNRgb R3Surfel::
Rgb(void) const
{
  // Return RGB
  return RNRgb(color[0]/255.0, color[1]/255.0, color[2]/255.0);
}



inline RNBoolean R3Surfel::
IsTerrestrial(void) const
{
  // Return whether point was captured with terrestrial scanner
  return (!IsAerial());
}



inline void R3Surfel::
SetCoords(float x, float y, float z)
{
  // Set position
  position[0] = x;
  position[1] = y;
  position[2] = z;
}



inline void R3Surfel::
SetCoords(const float *xyz)
{
  // Set position
  position[0] = xyz[0];
  position[1] = xyz[1];
  position[2] = xyz[2];
}



inline void R3Surfel::
SetColor(unsigned char r, unsigned char g, unsigned char b)
{
  // Set color
  color[0] = r;
  color[1] = g;
  color[2] = b;
}



inline void R3Surfel::
SetColor(const unsigned char *rgb)
{
  // Set color
  color[0] = rgb[0];
  color[1] = rgb[1];
  color[2] = rgb[2];
}



inline void R3Surfel::
SetColor(const RNRgb& rgb)
{
  // Set color
  color[0] = (unsigned char) (255.0 * rgb.R());
  color[1] = (unsigned char) (255.0 * rgb.G());
  color[2] = (unsigned char) (255.0 * rgb.B());
}



inline void R3Surfel::
Draw(RNFlags flags) const
{
  // Draw point at surfel
  glBegin(GL_POINTS);
  if (flags[R3_SURFEL_COLOR_DRAW_FLAG]) glColor3ubv(color);
  glVertex3fv(position);
  glEnd();
}



