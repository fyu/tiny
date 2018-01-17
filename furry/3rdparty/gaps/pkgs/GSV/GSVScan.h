// Include file for the GSV laser class 



////////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
////////////////////////////////////////////////////////////////////////

class GSVScan {
public:
  //////////////////////////////////////////
  //// CONSTRUCTOR/DESTRUCTOR FUNCTIONS ////
  //////////////////////////////////////////

  // Constructor/destructor functions
  GSVScan(void);
  virtual ~GSVScan(void);


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

  // Laser access functions
  GSVLaser *Laser(void) const;
  int LaserIndex(void) const;

  // Scanline access functions
  int NScanlines(void) const;
  GSVScanline *Scanline(int scanline_index) const;
  GSVScanline *FindScanlineBeforeTimestamp(RNScalar timestamp) const;


  ////////////////////////////
  //// PROPERTY FUNCTIONS ////
  ///////////////////////////

  // Spatial property functions
  const R3Box& BBox(void) const;

  // Surface property functions
  GSVMesh *Mesh(void) const;

  // Pose property functions
  GSVPose Pose(RNScalar timestamp) const;

  // Get user data 
  void *Data(void) const;


  ////////////////////////////////
  //// MANIPULATION FUNCTIONS ////
  ////////////////////////////////

  // Scanline manipulation
  void InsertScanline(GSVScanline *scanline);
  void RemoveScanline(GSVScanline *scanline);

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
  // I/O functions
  int ReadPoints(void);
  int ReadPoints(FILE *fp);
  int ReleasePoints(void);

  // Travel distances
  RNLength TravelDistance(int scanline_index) const;

  // Update functions
  void UpdateBBox(void);
  void InvalidateBBox(void);
  void UpdateTravelDistances(void);
  void InvalidateTravelDistances(void);

  // Memory management functions
  void EmptyCache(void);

protected:
  // Internal support functions
  GSVScanline *FindScanlineBeforeTimestamp(RNScalar timestamp, int imin, int imax) const;

protected:
  friend class GSVSegment;
  friend class GSVLaser;
  GSVSegment *segment;
  int segment_index;
  GSVLaser *laser;
  int laser_index;
  unsigned int read_count;
  RNArray<GSVScanline *> scanlines;
  RNLength *travel_distances;
  R3Box bbox;
  void *data;
};



////////////////////////////////////////////////////////////////////////
// Boundary types 
////////////////////////////////////////////////////////////////////////

#define GSV_LEFT_SILHOUETTE_BOUNDARY   (0b0000000000001000)
#define GSV_RIGHT_SILHOUETTE_BOUNDARY  (0b0000000000000100)
#define GSV_DOWN_SILHOUETTE_BOUNDARY   (0b0000000000000010)
#define GSV_UP_SILHOUETTE_BOUNDARY     (0b0000000000000001)
#define GSV_LEFT_SHADOW_BOUNDARY       (0b0000000010000000)
#define GSV_RIGHT_SHADOW_BOUNDARY      (0b0000000001000000)
#define GSV_DOWN_SHADOW_BOUNDARY       (0b0000000000100000)
#define GSV_UP_SHADOW_BOUNDARY         (0b0000000000010000)
#define GSV_LEFT_UNKNOWN_BOUNDARY      (0b0000100000000000)
#define GSV_RIGHT_UNKNOWN_BOUNDARY     (0b0000010000000000)
#define GSV_DOWN_UNKNOWN_BOUNDARY      (0b0000001000000000)
#define GSV_UP_UNKNOWN_BOUNDARY        (0b0000000100000000)
#define GSV_HORIZONTAL_RIDGE_BOUNDARY  (0b1000000000000000)
#define GSV_HORIZONTAL_VALLEY_BOUNDARY (0b0100000000000000)
#define GSV_VERTICAL_RIDGE_BOUNDARY    (0b0010000000000000)
#define GSV_VERTICAL_VALLEY_BOUNDARY   (0b0001000000000000)
#define GSV_SILHOUETTE_BOUNDARIES      (0b0000000000001111)
#define GSV_SILHOUETTE_BOUNDARIES      (0b0000000000001111)
#define GSV_SHADOW_BOUNDARIES          (0b0000000011110000)
#define GSV_GSV_UNKNOWN_BOUNDARIES     (0b0000111100000000)
#define RIDGE_BOUNDARIES               (0b1010000000000000)
#define GSV_VALLEY_BOUNDARIES          (0b0101000000000000)
#define GSV_ALL_BOUNDARIES             (0b1111111111111111)



////////////////////////////////////////////////////////////////////////
// INLINE FUNCTION DEFINITIONS
////////////////////////////////////////////////////////////////////////

inline GSVSegment *GSVScan::
Segment(void) const
{
  // Return segment
  return segment;
}



inline int GSVScan::
SegmentIndex(void) const
{
  // Return index of this scan in segment
  return segment_index;
}



inline GSVLaser *GSVScan::
Laser(void) const
{
  // Return laser
  return laser;
}



inline int GSVScan::
LaserIndex(void) const
{
  // Return index of this scan in laser
  return laser_index;
}



inline int GSVScan::
NScanlines(void) const
{
  // Return number of scanlines
  return scanlines.NEntries();
}



inline GSVScanline *GSVScan::
Scanline(int scanline_index) const
{
  // Return kth scanline
  return scanlines.Kth(scanline_index);
}



inline const R3Box& GSVScan::
BBox(void) const
{
  // Return bounding box
  if (bbox.XMin() == FLT_MAX) 
    ((GSVScan *) this)->UpdateBBox();
  return bbox;
}



inline RNLength GSVScan::
TravelDistance(int scanline_index) const
{
  // Return distance car laser viewpoint traveled before arriving at scanline
  assert((scanline_index >= 0) && (scanline_index < NScanlines()));
  if (!travel_distances) ((GSVScan *) this)->UpdateTravelDistances();
  return travel_distances[scanline_index];
}



inline void *GSVScan::
Data(void) const
{
  // Return user data
  return data;
}



inline void GSVScan::
SetBBox(const R3Box& box)
{
  // Set bbox
  this->bbox = box;
}



inline void GSVScan::
SetData(void *data)
{
  // Set user data
  this->data = data;
}



