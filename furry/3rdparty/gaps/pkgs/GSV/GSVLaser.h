// Include file for the GSV laser class 



////////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
////////////////////////////////////////////////////////////////////////

class GSVLaser {
public:
  //////////////////////////////////////////
  //// CONSTRUCTOR/DESTRUCTOR FUNCTIONS ////
  //////////////////////////////////////////

  // Constructor/destructor functions
  GSVLaser(void);
  virtual ~GSVLaser(void);


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


  ///////////////////////////
  //// PROPERTY FUNCTIONS ////
  ///////////////////////////

  // Spatial property functions
  const R3Box& BBox(void) const;

  // Get user data 
  void *Data(void) const;


  ////////////////////////////////
  //// MANIPULATION FUNCTIONS ////
  ////////////////////////////////

  // Scan manipulation
  void InsertScan(GSVScan *scan);
  void RemoveScan(GSVScan *scan);

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
  friend class GSVRun;
  GSVRun *run;
  int run_index;
  RNArray<GSVScan *> scans;
  R3Box bbox;
  void *data;
};



////////////////////////////////////////////////////////////////////////
// INLINE FUNCTION DEFINITIONS
////////////////////////////////////////////////////////////////////////

inline GSVRun *GSVLaser::
Run(void) const
{
  // Return run
  return run;
}



inline int GSVLaser::
RunIndex(void) const
{
  // Return index of this laser in run
  return run_index;
}



inline int GSVLaser::
NScans(void) const
{
  // Return number of scans
  return scans.NEntries();
}



inline GSVScan *GSVLaser::
Scan(int scan_index) const
{
  // Return kth scan
  return scans.Kth(scan_index);
}



inline const R3Box& GSVLaser::
BBox(void) const
{
  // Return bounding box
  if (bbox.XMin() == FLT_MAX) 
    ((GSVLaser *) this)->UpdateBBox();
  return bbox;
}



inline void *GSVLaser::
Data(void) const
{
  // Return user data
  return data;
}



inline void GSVLaser::
SetBBox(const R3Box& box)
{
  // Set bbox
  this->bbox = box;
}



inline void GSVLaser::
SetData(void *data)
{
  // Set user data
  this->data = data;
}


