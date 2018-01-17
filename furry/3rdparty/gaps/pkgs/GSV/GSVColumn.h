// Include file for the GSV column class 



////////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
////////////////////////////////////////////////////////////////////////

class GSVColumn {
public:
  //////////////////////////////////////////
  //// CONSTRUCTOR/DESTRUCTOR FUNCTIONS ////
  //////////////////////////////////////////

  // Constructor functions
  GSVColumn(void);
  GSVColumn(const GSVPose& pose, RNScalar timestamp);
  virtual ~GSVColumn(void);


  //////////////////////////
  //// ACCESS FUNCTIONS ////
  //////////////////////////

  // Scene access functions
  GSVScene *Scene(void) const;

  // Run access functions
  GSVRun *Run(void) const;

  // Segment access functions
  GSVSegment *Segment(void) const;

  // Camera access functions
  GSVCamera *Camera(void) const;

  // Camera image access functions
  GSVImage *Image(void) const;
  int ImageIndex(void) const;


  //////////////////////////
  //// PROPERTY FUNCTIONS ////
  //////////////////////////

  // Extrinsic sensor properties
  const GSVPose& Pose(void) const;

  // Timing properties
  RNScalar Timestamp(void) const;


  /////////////////////////////
  //// USER DATA FUNCTIONS ////
  /////////////////////////////

  // Get user data 
  void *Data(void) const;

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
  // Extrinsic sensor property functions
  void SetPose(const GSVPose& pose);

  // Timing property functions
  void SetTimestamp(RNScalar timestamp);

  // Update functions
  void Update(void);

protected:
  friend class GSVImage;
  GSVImage *image;
  int image_index;
  GSVPose pose;
  RNScalar timestamp;
  void *data;
};



////////////////////////////////////////////////////////////////////////
// INLINE FUNCTION DEFINITIONS
////////////////////////////////////////////////////////////////////////

inline GSVImage *GSVColumn::
Image(void) const
{
  // Return camera image
  return image;
}



inline int GSVColumn::
ImageIndex(void) const
{
  // Return index of this column in camera image
  return image_index;
}



inline const GSVPose& GSVColumn::
Pose(void) const
{
  // Return sensor pose
  return pose;
}



inline RNScalar GSVColumn::
Timestamp(void) const
{
  // Return timestamp 
  return timestamp;
}



inline void *GSVColumn::
Data(void) const
{
  // Return user data
  return data;
}



inline void GSVColumn::
Update(void)
{
}


inline void GSVColumn::
SetPose(const GSVPose& pose)
{
  // Set pose
  this->pose = pose;
}



inline void GSVColumn::
SetTimestamp(RNScalar timestamp)
{
  // Set timestamp
  this->timestamp = timestamp;
}



inline void GSVColumn::
SetData(void *data)
{
  // Set user data
  this->data = data;
}



