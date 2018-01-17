// Include file for the GSV panorama class 



////////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
////////////////////////////////////////////////////////////////////////

class GSVPanorama {
public:
  //////////////////////////////////////////
  //// CONSTRUCTOR/DESTRUCTOR FUNCTIONS ////
  //////////////////////////////////////////

  // Constructor/destructor functions
  GSVPanorama(void);
  virtual ~GSVPanorama(void);


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

  // Image access functions
  int NImages(void) const;
  GSVImage *Image(int image_index) const;


  ///////////////////////////
  //// PROPERTY FUNCTIONS ////
  ///////////////////////////

  // Viewpoint property functions
  const R3Point& Viewpoint(void) const;

  // Timestamp property functins
  RNScalar Timestamp(void) const;

  // Get user data 
  void *Data(void) const;


  ////////////////////////////////
  //// MANIPULATION FUNCTIONS ////
  ////////////////////////////////

  // Camera image manipulation
  void InsertImage(GSVImage *image);
  void RemoveImage(GSVImage *image);

  // Extrinsic sensor property functions
  void SetViewpoint(const R3Point& viewpoint);

  // Timing property functions
  void SetTimestamp(RNScalar timestamp);

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
  void UpdateViewpoint(void);
  void UpdateTimestamp(void);

protected:
  friend class GSVSegment;
  GSVSegment *segment;
  int segment_index;
  RNArray<GSVImage *> images;
  R3Point viewpoint;
  RNScalar timestamp;
  void *data;
};



////////////////////////////////////////////////////////////////////////
// INLINE FUNCTION DEFINITIONS
////////////////////////////////////////////////////////////////////////

inline GSVSegment *GSVPanorama::
Segment(void) const
{
  // Return segment
  return segment;
}



inline int GSVPanorama::
SegmentIndex(void) const
{
  // Return index of this panaorama in segment
  return segment_index;
}



inline int GSVPanorama::
NImages(void) const
{
  // Return number of images
  return images.NEntries();
}



inline GSVImage *GSVPanorama::
Image(int k) const
{
  // Return kth image
  return images.Kth(k);
}



inline const R3Point& GSVPanorama::
Viewpoint(void) const
{
  // Return (average) viewpoint
  if (viewpoint.X() == FLT_MAX) 
    ((GSVPanorama *) this)->UpdateViewpoint();
  return viewpoint;
}



inline RNScalar GSVPanorama::
Timestamp(void) const
{
  // Return (average) timestamp 
  if (timestamp == -1) 
    ((GSVPanorama *) this)->UpdateTimestamp();
  return timestamp;
}



inline void *GSVPanorama::
Data(void) const
{
  // Return user data
  return data;
}



inline void GSVPanorama::
SetViewpoint(const R3Point& viewpoint)
{
  // Set viewpoint
  this->viewpoint = viewpoint;
}



inline void GSVPanorama::
SetTimestamp(RNScalar timestamp)
{
  // Set timestamp
  this->timestamp = timestamp;
}



inline void GSVPanorama::
SetData(void *data)
{
  // Set user data
  this->data = data;
}


