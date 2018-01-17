// Include file for the GSV segment class 



////////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
////////////////////////////////////////////////////////////////////////

class GSVSegment {
public:
  //////////////////////////////////////////
  //// CONSTRUCTOR/DESTRUCTOR FUNCTIONS ////
  //////////////////////////////////////////

  // Constructor/destructor functions
  GSVSegment(void);
  virtual ~GSVSegment(void);


  //////////////////////////
  //// ACCESS FUNCTIONS ////
  //////////////////////////

  // Scene access functions
  GSVScene *Scene(void) const;
  int SceneIndex(void) const;

  // Run access functions
  GSVRun *Run(void) const;
  int RunIndex(void) const;

  // Scan access functions
  int NScans(void) const;
  GSVScan *Scan(int scan_index) const;

  // Tapestry access functions
  int NTapestries(void) const;
  GSVTapestry *Tapestry(int tapestry_index) const;

  // Panorama access functions
  int NPanoramas(void) const;
  GSVPanorama *Panorama(int panorama_index) const;
  GSVPanorama *FindPanoramaBeforeTimestamp(RNScalar timestamp) const;

  // Image access functions (indirect - i.e., slow) 
  int NImages(void) const;
  GSVImage *Image(int image_index) const;

  // Scanline access functions (indirect - i.e., slow) 
  int NScanlines(void) const;
  GSVScanline *Scanline(int scanline_index) const;


  //////////////////////////
  //// PROPERTY FUNCTIONS ////
  //////////////////////////

  // Spatial property functions
  const R3Box& BBox(void) const;

  // Get user data 
  void *Data(void) const;


  ////////////////////////////////
  //// MANIPULATION FUNCTIONS ////
  ////////////////////////////////

  // Laser tapestry manipulation functions
  void InsertScan(GSVScan *scan);
  void RemoveScan(GSVScan *scan);

  // Camera tapestry manipulation functions
  void InsertTapestry(GSVTapestry *tapestry);
  void RemoveTapestry(GSVTapestry *tapestry);

  // Camera panorama manipulation functions
  void InsertPanorama(GSVPanorama *panorama);
  void RemovePanorama(GSVPanorama *panorama);

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

protected:
  // Internal support functions
  GSVPanorama *FindPanoramaBeforeTimestamp(RNScalar timestamp, int imin, int imax) const;

protected:
  friend class GSVScene;
  friend class GSVRun;
  int scene_index;
  GSVRun *run;
  int run_index;
  RNArray<GSVScan *> scans;
  RNArray<GSVTapestry *> tapestries;
  RNArray<GSVPanorama *> panoramas;
  R3Box bbox;
  void *data;
};



////////////////////////////////////////////////////////////////////////
// INLINE FUNCTION DEFINITIONS
////////////////////////////////////////////////////////////////////////

inline GSVRun *GSVSegment::
Run(void) const
{
  // Return run
  return run;
}



inline int GSVSegment::
RunIndex(void) const
{
  // Return index of this segment in run
  return run_index;
}



inline int GSVSegment::
NScans(void) const
{
  // Return number of laser tapestries
  return scans.NEntries();
}



inline GSVScan *GSVSegment::
Scan(int scan_index) const
{
  // Return kth laser tapestry
  return scans.Kth(scan_index);
}



inline int GSVSegment::
NTapestries(void) const
{
  // Return number of camera tapestries
  return tapestries.NEntries();
}



inline GSVTapestry *GSVSegment::
Tapestry(int tapestry_index) const
{
  // Return kth camera tapestry
  return tapestries.Kth(tapestry_index);
}



inline int GSVSegment::
NPanoramas(void) const
{
  // Return number of camera panoramas
  return panoramas.NEntries();
}



inline GSVPanorama *GSVSegment::
Panorama(int panorama_index) const
{
  // Return kth camera panorama
  return panoramas.Kth(panorama_index);
}



inline const R3Box& GSVSegment::
BBox(void) const
{
  // Return bounding box
  if (bbox.XMin() == FLT_MAX) 
    ((GSVSegment *) this)->UpdateBBox();
  return bbox;
}



inline void *GSVSegment::
Data(void) const
{
  // Return user data
  return data;
}



inline void GSVSegment::
SetBBox(const R3Box& box)
{
  // Set bbox
  this->bbox = box;
}



inline void GSVSegment::
SetData(void *data)
{
  // Set user data
  this->data = data;
}


