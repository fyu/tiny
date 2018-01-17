/* Include file for the R3 surfel block class */



////////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
////////////////////////////////////////////////////////////////////////

class R3SurfelBlock {
public:
  //////////////////////////////////////////
  //// CONSTRUCTOR/DESTRUCTOR FUNCTIONS ////
  //////////////////////////////////////////

  // Constructor functions
  R3SurfelBlock(void);
  R3SurfelBlock(const R3SurfelBlock& block);
  R3SurfelBlock(const R3SurfelPointSet *set);
  R3SurfelBlock(const R3SurfelPointSet *set, const R3Point& origin);
  R3SurfelBlock(const R3Surfel *surfels, int nsurfels, const R3Point& origin = R3zero_point);
  R3SurfelBlock(const RNArray<const R3Surfel *>& surfels, const R3Point& origin = R3zero_point);
  R3SurfelBlock(const R3Point *points, int npoints);
  R3SurfelBlock(const RNArray<R3Point *>& points);

  // Destructor function
  ~R3SurfelBlock(void);


  //////////////////////////
  //// ACCESS FUNCTIONS ////
  //////////////////////////

  // Node access functions
  R3SurfelNode *Node(void) const;

  // Database access functions
  R3SurfelDatabase *Database(void) const;
  int DatabaseIndex(void) const;

  // Surfel access functions
  int NSurfels(void) const;
  const R3Surfel *Surfels(void) const;
  const R3Surfel *Surfel(int k) const;
  const R3Surfel *operator[](int k) const;


  //////////////////////////////////
  //// BLOCK PROPERTY FUNCTIONS ////
  //////////////////////////////////

  // Geometric property functions
  const R3Point& Origin(void) const;
  const R3Box& BBox(void) const;
  R3Point Centroid(void) const;
  RNLength Resolution(void) const;

  // Scanner source property functions
  RNBoolean HasAerial(void) const;
  RNBoolean HasTerrestrial(void) const;


  ///////////////////////////////////
  //// SURFEL PROPERTY FUNCTIONS ////
  ///////////////////////////////////

  // Surfel property functions
  R3Point SurfelPosition(int surfel_index);
  RNRgb SurfelColor(int surfel_index);
  RNBoolean IsSurfelMarked(int surfel_index);
  RNBoolean IsSurfelAerial(int surfel_index);
  RNBoolean IsSurfelTerrestrial(int surfel_index);


  //////////////////////////////////////////
  //// BLOCK MANIPULATION FUNCTIONS ////
  //////////////////////////////////////////

  // Property manipulation functions
  void SetOrigin(const R3Point& position);

  // Surfel mark manipulation functions
  void SetMarks(RNBoolean mark = TRUE);


  ///////////////////////////////////////
  //// SURFEL MANIPULATION FUNCTIONS ////
  ///////////////////////////////////////

  // Surfel manipulation functions
  void SetSurfelPosition(int surfel_index, const R3Point& position);
  void SetSurfelColor(int surfel_index, const RNRgb& color);
  void SetSurfelAerial(int surfel_index, RNBoolean aerial = TRUE);
  void SetSurfelMark(int surfel_index, RNBoolean mark = TRUE);


  ///////////////////////////
  //// DISPLAY FUNCTIONS ////
  ///////////////////////////

  // Draw function
  void Draw(RNFlags flags = R3_SURFEL_DEFAULT_DRAW_FLAGS) const;

  // Print function
  void Print(FILE *fp = NULL, const char *prefix = NULL, const char *suffix = NULL) const;

  //////////////////////////////////////////
  //// I/O functions
  //////////////////////////////////////////

  // File I/O
  int ReadFile(const char *filename);
  int ReadXYZFile(const char *filename);
  int ReadBinaryFile(const char *filename);
  int ReadUPCFile(const char *filename);
  int ReadOBJFile(const char *filename);
  int ReadXYZ(FILE *fp);
  int ReadBinary(FILE *fp);
  int ReadUPC(FILE *fp);
  int ReadOBJ(FILE *fp);
  int WriteFile(const char *filename) const;
  int WriteXYZFile(const char *filename) const;
  int WriteBinaryFile(const char *filename) const;
  int WriteXYZ(FILE *fp) const;
  int WriteBinary(FILE *fp) const;

  ////////////////////////////////////////////////////////////////////////
  // INTERNAL STUFF BELOW HERE
  ////////////////////////////////////////////////////////////////////////

  // Manipulation functions
  RNBoolean IsDirty(void) const;
  void SetDirty(RNBoolean dirty = TRUE);

  // Temporary
  int ReadCount(void) const { return file_read_count; }

public:
  // Update functions
  void UpdateProperties(void);

protected:
  // Database update functions
  void UpdateAfterInsert(R3SurfelDatabase *database);
  void UpdateBeforeRemove(R3SurfelDatabase *database);

  // Tree node update functions
  void UpdateAfterInsert(R3SurfelNode *node);
  void UpdateBeforeRemove(R3SurfelNode *node);

private:
  // Manipulation functions
  void SetDatabase(R3SurfelDatabase *database);

  // Update functions
  void UpdateBBox(void);
  void UpdateResolution(void);
  void UpdateScannerSources(void);

private:
  // Surfel data
  R3Surfel *surfels;
  int nsurfels;

  // Property data
  R3Point origin;
  R3Box bbox;
  RNLength resolution;
  RNFlags flags;

  // Database data
  friend class R3SurfelDatabase;
  friend class R3SurfelPointSet;
  class R3SurfelDatabase *database;
  int database_index;
  unsigned long long file_surfels_offset;
  unsigned int file_surfels_count;
  unsigned int file_read_count;

  // Node data
  friend class R3SurfelNode;
  R3SurfelNode *node;

  // Display data
  GLuint opengl_id;
};



////////////////////////////////////////////////////////////////////////
// BLOCK FLAGS
////////////////////////////////////////////////////////////////////////

#define R3_SURFEL_BLOCK_PROPERTY_FLAGS                 0x00FF
#define R3_SURFEL_BLOCK_BBOX_UPTODATE_FLAG             0x0001
#define R3_SURFEL_BLOCK_RESOLUTION_UPTODATE_FLAG       0x0002
#define R3_SURFEL_BLOCK_SCANNER_SOURCES_UPTODATE_FLAG  0x0004

#define R3_SURFEL_BLOCK_HAS_AERIAL_FLAG                0x0010
#define R3_SURFEL_BLOCK_HAS_TERRESTRIAL_FLAG           0x0020

#define R3_SURFEL_BLOCK_DATABASE_FLAGS                 0xFF00
#define R3_SURFEL_BLOCK_DIRTY_FLAG                     0x0100
#define R3_SURFEL_BLOCK_DELETE_PENDING_FLAG            0x0200



////////////////////////////////////////////////////////////////////////
// INLINE FUNCTION DEFINITIONS
////////////////////////////////////////////////////////////////////////

inline const R3Point& R3SurfelBlock::
Origin(void) const
{
  // Return position of origin
  // Surfels coordinates are relative to this position
  return origin;
}



inline R3Point R3SurfelBlock::
Centroid(void) const
{
  // Return centroid of bounding box
  return BBox().Centroid();
}



inline R3SurfelNode *R3SurfelBlock::
Node(void) const
{
  // Return tree node this block is part of
  return node;
}



inline int R3SurfelBlock::
NSurfels(void) const
{
  // Return number of surfels
  return nsurfels;
}



inline const R3Surfel *R3SurfelBlock::
Surfels(void) const
{
  // Return pointer to block of surfels
  return surfels;
}



inline const R3Surfel *R3SurfelBlock::
Surfel(int k) const
{
  // Return kth surfel
  assert((k >= 0) && (k < nsurfels));
  return &surfels[k];
}



inline const R3Surfel *R3SurfelBlock::
operator[](int k) const
{
  // Return kth surfel
  return Surfel(k);
}



inline R3SurfelDatabase *R3SurfelBlock::
Database(void) const
{
  // Return database
  return database;
}



inline int R3SurfelBlock::
DatabaseIndex(void) const
{
  // Return index of block in database
  return database_index;
}



inline R3Point R3SurfelBlock::
SurfelPosition(int surfel_index)
{
  // Return position of kth surfel
  R3Surfel& surfel = surfels[surfel_index];
  return R3Point(origin.X() + surfel.X(), origin.Y() + surfel.Y(), origin.Z() + surfel.Z());
}



inline RNRgb R3SurfelBlock::
SurfelColor(int surfel_index)
{
  // Return color of kth surfel
  return surfels[surfel_index].Rgb();
}



inline RNBoolean R3SurfelBlock::
IsSurfelMarked(int surfel_index)
{
  // Return whether kth surfel is marked
  return surfels[surfel_index].IsMarked();
}



inline RNBoolean R3SurfelBlock::
IsSurfelAerial(int surfel_index)
{
  // Return whether kth surfel is aerial
  return surfels[surfel_index].IsAerial();
}



inline RNBoolean R3SurfelBlock::
IsSurfelTerrestrial(int surfel_index)
{
  // Return whether kth surfel is aerial
  return surfels[surfel_index].IsTerrestrial();
}



inline void R3SurfelBlock::
UpdateAfterInsert(R3SurfelDatabase *database)
{
}



inline void R3SurfelBlock::
UpdateBeforeRemove(R3SurfelDatabase *database)
{
}



