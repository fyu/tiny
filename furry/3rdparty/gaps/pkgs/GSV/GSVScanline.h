// Include file for the GSV scanline class 



////////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
////////////////////////////////////////////////////////////////////////

class GSVScanline {
public:
  //////////////////////////////////////////
  //// CONSTRUCTOR/DESTRUCTOR FUNCTIONS ////
  //////////////////////////////////////////

  // Constructor functions
  GSVScanline(void);
  GSVScanline(const GSVPose& pose, RNScalar timestamp);
  virtual ~GSVScanline(void);


  //////////////////////////
  //// ACCESS FUNCTIONS ////
  //////////////////////////

  // Scene access functions
  GSVScene *Scene(void) const;
  int SceneIndex(void) const;

  // Run access functions
  GSVRun *Run(void) const;
  int RunIndex(void) const;

  // Segment access functions
  GSVSegment *Segment(void) const;
  int SegmentIndex(void) const;

  // Scan access functions
  GSVScan *Scan(void) const;
  int ScanIndex(void) const;

  // Point access functions
  int NPoints(void) const;
  const R3Point& PointPosition(int point_index) const;
  RNAngle PointAngle(int point_index) const;


  ////////////////////////////
  //// PROPERTY FUNCTIONS ////
  ////////////////////////////

  // Extrinsic sensor property functions
  const GSVPose& Pose(void) const;

  // Timing property functions
  RNScalar Timestamp(void) const;

  // Spatial property functions
  const R3Box& BBox(void) const;

  // Estimated property functions
  RNCoord EstimatedGroundZ(void) const;

  // Distance laser travelled along scan
  RNLength TravelDistance(void) const;
  
  // Get user data 
  void *Data(void) const;


  ////////////////////////////////
  //// MANIPULATION FUNCTIONS ////
  ////////////////////////////////

  // Point manipulation functions
  void SetPointPositions(const R3Point *positions, int npositions);
  void SetPointPosition(int point_index, const R3Point& position);

  // Extrinsic sensor property functions
  void SetPose(const GSVPose& pose);

  // Timing property functions
  void SetTimestamp(RNScalar timestamp);

  // Spatial property manipulation
  void SetBBox(const R3Box& box);

  // Set user data 
  void SetData(void *data);
  

  ///////////////////////////
  //// DISPLAY FUNCTIONS ////
  ///////////////////////////

  // Draw function
  virtual void Draw(RNFlags flags = GSV_DEFAULT_DRAW_FLAGS) const;

  // Print function
  virtual void Print(FILE *fp = NULL, const char *prefix = NULL, const char *suffix = NULL) const;


////////////////////////////////////////////////////////////////////////
// INTERNAL STUFF BELOW HERE
////////////////////////////////////////////////////////////////////////

public:
  // Update functions 
  void UpdateBBox(void);
  void InvalidateBBox(void);
  RNBoolean DoesBBoxNeedUpdate(void) const;
  
  // I/O functions
  RNBoolean ArePointsResident(void) const;
  int ReadPoints(FILE *fp, RNBoolean seek = TRUE);
  int WritePoints(FILE *fp, RNBoolean seek = TRUE);
  int ReleasePoints(void);

  // File functions
  unsigned long long FileOffset(void) const;
  void SetFileOffset(unsigned long long file_offset);
  void SetNPoints(int npoints);

protected:
  friend class GSVScan;
  GSVScan *scan;
  int scan_index;
  unsigned long long file_offset;
  unsigned int read_count;
  R3Point *points;
  int npoints;
  GSVPose pose;
  RNScalar timestamp;
  R3Box bbox;
  void *data;
};



////////////////////////////////////////////////////////////////////////
// INLINE FUNCTION DEFINITIONS
////////////////////////////////////////////////////////////////////////

inline GSVScan *GSVScanline::
Scan(void) const
{
  // Return laser scan
  return scan;
}



inline int GSVScanline::
ScanIndex(void) const
{
  // Return index of this scanline in scan
  return scan_index;
}



inline int GSVScanline::
NPoints(void) const
{
  // Return number of points
  return npoints;
}



inline const R3Point& GSVScanline::
PointPosition(int point_index) const
{
  // Return point position
  return points[point_index];
}



inline const GSVPose& GSVScanline::
Pose(void) const
{
  // Return sensor pose
  return pose;
}



inline RNScalar GSVScanline::
Timestamp(void) const
{
  // Return timestamp 
  return timestamp;
}



inline const R3Box& GSVScanline::
BBox(void) const
{
  // Return bounding box
  if (DoesBBoxNeedUpdate()) 
    ((GSVScanline *) this)->UpdateBBox();
  return bbox;
}



inline void *GSVScanline::
Data(void) const
{
  // Return user data
  return data;
}



inline void GSVScanline::
SetTimestamp(RNScalar timestamp)
{
  // Set timestamp
  this->timestamp = timestamp;
}



inline void GSVScanline::
SetBBox(const R3Box& box)
{
  // Set bbox
  this->bbox = box;
}



inline void GSVScanline::
SetData(void *data)
{
  // Set user data
  this->data = data;
}



inline RNBoolean GSVScanline::
ArePointsResident(void) const
{
  // Return whether points are resident
  return (points) ? TRUE: FALSE;
}



inline unsigned long long GSVScanline::
FileOffset(void) const
{
  // Return offset of scanline points in GSV file
  return file_offset;
}



inline void GSVScanline::
SetFileOffset(unsigned long long file_offset)
{
  // Set offset of scanline points in GSV file
  this->file_offset = file_offset;
}



inline void GSVScanline::
SetNPoints(int npoints)
{
  // Set number of points
  this->npoints = npoints;
}



