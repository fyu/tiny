// Include file for pose optimization 



////////////////////////////////////////////////////////////////////////
// Dependencies
////////////////////////////////////////////////////////////////////////

#include "RNMath/RNMath.h"



////////////////////////////////////////////////////////////////////////
// Class declarations
////////////////////////////////////////////////////////////////////////

struct GSVDescriptor {
public:
  GSVDescriptor(int descriptor_type = 0, RNScalar *values = NULL, int nvalues = 0);
  GSVDescriptor(const GSVDescriptor& descriptor);
  ~GSVDescriptor(void);
  RNScalar SquaredDistance(const GSVDescriptor& descriptor) const;
public:
  int descriptor_type;
  RNScalar *values; 
  int nvalues;
};

struct GSVFeature {
public:
  GSVFeature(void);
  GSVFeature(int feature_type, GSVScanline *scanline, int scan_point_index, 
    const R3Point& scan_position, const R3Vector& scan_direction, const R3Vector& scan_normal, 
    RNLength scan_scale, RNScalar score, int group = 0);
  GSVFeature(int feature_type, GSVImage *image, 
    const R2Point& image_position, const R2Vector& image_direction, 
    RNLength image_scale, RNScalar score, int group = 0);
  GSVFeature(int feature_type, GSVImage *image, GSVScanline *scanline, int scan_point_index, 
    const R3Point& scan_position, const R3Vector& scan_direction, const R3Vector& scan_normal, RNLength scan_scale,
    const R2Point& image_position, const R2Vector& image_direction, RNLength image_scale, 
    RNScalar image_t, RNScalar score, int group = 0);
  void Draw(void) const;
public:
  int feature_type;
  GSVImage *image;
  GSVScanline *scanline;
  int scan_point_index;
  R3Point scan_position;
  R3Vector scan_direction;
  R3Vector scan_normal;
  RNScalar scan_scale;
  R2Point image_position;
  R2Vector image_direction;
  RNScalar image_scale;
  RNScalar image_t;
  RNScalar score;
  GSVDescriptor descriptor;
  int group;
  int index;
};

struct GSVFeaturePair {
public:
  GSVFeaturePair(void);
  GSVFeaturePair(const GSVFeaturePair& pair);
  GSVFeaturePair(const struct GSVFeatureCorrespondence& correspondence, int dummy);
  GSVFeaturePair(GSVFeature *feature0, GSVFeature *feature1, RNScalar score = 1);
  GSVFeature *Feature(int k) const;
  void Draw(void) const;
public:
  GSVFeature *features[2];
  RNScalar score;
  int index;
};

struct GSVFeatureCorrespondence {
public:
  GSVFeatureCorrespondence(void);
  GSVFeatureCorrespondence(const GSVFeaturePair& pair, int dummy);
  GSVFeatureCorrespondence(const GSVFeatureCorrespondence& correspondence);
  GSVFeatureCorrespondence(GSVFeature *feature0, GSVFeature *feature1, RNScalar score = 1);
  GSVFeature *Feature(int k) const;
  void Draw(void) const;
public:
  GSVFeature *features[2];
  RNScalar score;
  int index;
};

struct GSVFeatureCluster {
public:
  GSVFeatureCluster(void);
  GSVFeatureCluster(int feature_type, const R3Point& scan_position, 
    const R3Vector& scan_direction, const R3Vector& scan_normal, 
    RNLength scan_scale, RNScalar score);
  GSVFeatureCluster(const GSVFeatureCluster& cluster);
  int NFeatures(void) const;
  GSVFeature *Feature(int k) const;
  void InsertFeature(GSVFeature *feature);
  void RemoveFeature(GSVFeature *feature);
  int NClusters(void) const;
  GSVFeatureCluster *Cluster(int k) const;
  void InsertCluster(GSVFeatureCluster *cluster);
  void RemoveCluster(GSVFeatureCluster *cluster);
  void Update(struct GSVPoseOptimization *optimization = NULL);
  void Draw(void) const;
public:
  int feature_type;
  R3Point scan_position;
  R3Vector scan_direction;
  R3Vector scan_normal;
  RNScalar scan_scale;
  RNScalar score;
  R3Vector translation;
  R3Vector rotation;
  RNArray<GSVFeature *> features;
  RNArray<GSVFeatureCluster *> clusters;
  RNScalar inertia[6];
  int index;
};

class GSVPathVertex {
public:
  GSVPathVertex(void);
public:
  struct GSVPath *path;
  GSVPose pose;
  R3Vector translation;
  R3Vector rotation;
  RNScalar parameter;
  RNScalar timestamp;
  RNScalar inertia[6];
  int index;
};

struct GSVPath {
public:
  GSVPath(GSVSegment *segment = NULL, RNLength max_vertex_spacing = 2);
  ~GSVPath(void);
  GSVPose TransformedPose(RNScalar u) const;
  GSVPose Pose(RNScalar u) const;
  R3Vector Translation(RNScalar u) const;
  R3Vector Rotation(RNScalar u) const;
  RNScalar Parameter(RNScalar timestamp) const;
  RNScalar Timestamp(RNScalar u) const;
  RNScalar Inertia(RNScalar u, int variable_index) const;
  RNInterval ParameterRange(void) const;
  RNInterval TimestampRange(void) const;
  void Draw(void) const;
public:
  friend struct GSVPoseOptimization;
  int NVertices(void) const;
  GSVPathVertex *Vertex(int vertex_index) const;
  const GSVPose& VertexPose(int vertex_index) const;
  const R3Point& VertexPosition(int vertex_index) const;
  const R3Vector& VertexTranslation(int vertex_index) const;
  const R3Vector& VertexRotation(int vertex_index) const;
  RNScalar VertexParameter(int vertex_index) const;
  RNScalar VertexTimestamp(int vertex_index) const;
  RNScalar VertexInertia(int vertex_index, int variable_index) const;
  GSVPose VertexTransformedPose(int vertex_index) const;
  void SetVertexTranslation(int vertex_index, const R3Vector& rotation);
  void SetVertexRotation(int vertex_index, const R3Vector& rotation);
  int FindVertexIndexBeforeTimestamp(RNScalar timestamp) const;
  int FindVertexIndexBeforeTimestamp(RNScalar timestamp, int imin, int imax) const;
public:
  GSVSegment *segment;
  R3CatmullRomSpline *spline;
  RNArray<class GSVPathVertex *> vertices;
  int index;
};

struct GSVPoseOptimization {
public:
  // Constructor/destructor
  GSVPoseOptimization(GSVScene *scene, RNLength path_vertex_spacing = 1);
  ~GSVPoseOptimization(void);

  // Access functions
  int NFeatures(void) const;
  GSVFeature *Feature(int k) const;
  int NPairs(void) const;
  GSVFeaturePair *Pair(int k) const;
  int NCorrespondences(void) const;
  GSVFeatureCorrespondence *Correspondence(int k) const;
  int NClusters(void) const;
  GSVFeatureCluster *Cluster(int k) const;

  // Property functions
  RNScalar Score(void);
  RNScalar WorldRMSD(void);
  RNScalar ImageRMSD(void);

  // Manipulation functions
  void InsertFeature(GSVFeature *feature);
  void RemoveFeature(GSVFeature *feature);
  int TruncateFeatures(int max_features = 0);
  void InsertPair(GSVFeaturePair *pair);
  void RemovePair(GSVFeaturePair *pair);
  int TruncatePairs(int max_pairs = 0);
  void InsertCorrespondence(GSVFeatureCorrespondence *correspondence);
  void RemoveCorrespondence(GSVFeatureCorrespondence *correspondence);
  int TruncateCorrespondences(int max_correspondences = 0);
  void InsertCluster(GSVFeatureCluster *cluster);
  void RemoveCluster(GSVFeatureCluster *cluster);
  int TruncateClusters(int max_clusters = 0);

  // Feature creation functions
  int CreateImageCornerFeatures(void);
  int CreateImageSiftFeatures(void);
  int CreateImageLineFeatures(void);
  int CreateScanPoleFeatures(void);
  int CreateScanCurbFeatures(void);
  int CreateScanRidgeAndValleyFeatures(void);
  int CreateScanEdgeFeatures(void);
  int CreateScanPlaneFeatures(void);

  // Pair creation functions (from features)
  int CreateAllPairs(RNBoolean check_duplicates = TRUE);
  int CreateFeasiblePairs(RNBoolean check_duplicates = TRUE);
  int CreateClosePairs(RNBoolean check_duplicates = TRUE);
  int CreateCoplanarPairs(RNBoolean check_duplicates = TRUE);
  int CreatePixelPixelPlanePairs(RNBoolean check_duplicates = TRUE);

  // Correspondence creation functions (from pairs)
  int CreateAllCorrespondences(RNBoolean check_duplicates = TRUE);
  int CreateFeasibleCorrespondences(RNBoolean check_duplicates = TRUE);
  int CreateCloseCorrespondences(RNBoolean check_duplicates = TRUE);
  int CreateCoplanarCorrespondences(RNBoolean check_duplicates = TRUE);
  int CreateICPCorrespondences(RNBoolean check_duplicates = TRUE);

  // Correspondence creation functions (from other correspondences)
  int CreatePixelPixelPlaneCorrespondences(RNBoolean check_duplicates = TRUE);
  int CreateTransitiveClosureCorrespondences(RNBoolean check_duplicates = TRUE);

  // Cluster creation functions (from correspondences)
  int CreateAllClusters(RNBoolean check_duplicates = TRUE);
  int CreateCloseClusters(RNBoolean check_duplicates = TRUE);
  int CreateICPClusters(RNBoolean check_duplicates = TRUE);
  int CreateHierarchicalClusters(RNBoolean check_duplicates = TRUE);
  int CreateMergedClusters(RNBoolean check_duplicates = TRUE);
  int CreateCoplanarClusters(RNBoolean check_duplicates = TRUE);
  int CreateCorrespondenceClusters(RNBoolean check_duplicates = TRUE);

  // Input/output functions
  int ReadFile(const char *filename);
  int WriteFile(const char *filename, RNBoolean apply_pose_transformations = FALSE) const;
  int ReadAsciiFile(const char *filename);
  int WriteAsciiFile(const char *filename, RNBoolean apply_pose_transformations = FALSE) const;
  int ReadAscii(FILE *fp = NULL);
  int WriteAscii(FILE *fp = NULL, RNBoolean apply_pose_transformations = FALSE) const;
  int ReadBinaryFile(const char *filename);
  int WriteBinaryFile(const char *filename, RNBoolean apply_pose_transformations = FALSE) const;
  int ReadBinary(FILE *fp = NULL);
  int WriteBinary(FILE *fp = NULL, RNBoolean apply_pose_transformations = FALSE) const;
  int ReadSiftFile(GSVImage *image, const char *filename);
  int ReadSiftFile(const char *filename);
  int ReadLineFile(GSVImage *image, const char *filename);
  int ReadLineFile(const char *filename);
  int ReadTransformationsFile(const char *filename);
  int WriteTransformationsFile(const char *filename) const;
  int ReadInertiaFile(const char *filename);
  int WriteInertiaFile(const char *filename) const;
  int ReadOptionsFile(const char *filename);
  int WriteOptionsFile(const char *filename) const;
  int ReadOptions(FILE *fp = NULL);
  int WriteOptions(FILE *fp = NULL) const;
  int ParseOption(const char *keyword, const char *value);

  // Scene update function
  int ApplyOptimizedTransformationsToScene(void);

  // Optimization functions
  int Solve(
    RNBoolean update_laser_translations = TRUE, RNBoolean update_laser_rotations = TRUE, 
    RNBoolean update_camera_translations = TRUE, RNBoolean update_camera_rotations = TRUE, 
    RNBoolean update_path_translations = TRUE, RNBoolean update_path_rotations = TRUE, 
    RNBoolean update_cluster_positions = TRUE, RNBoolean update_cluster_rotations = TRUE,
    RNBoolean update_pixel_depths = TRUE, 
    RNBoolean include_correspondence_equations = TRUE, 
    RNBoolean include_cluster_equations = TRUE,
    RNScalar rigidity = 0.5, int solver = RN_CERES_SOLVER);

public:
  // Pixel position updates
  int UpdatePixelPositionsFromCorrespondenceSolutions(void);
  int UpdatePixelPositionsFromRayIntersections(void);
  int UpdatePixelPositionsFromMeshIntersections(void);
  int ClearPixelPositions(void);

  // Cluster updates
  int UpdateClusterPositionsFromFeaturePositions(void);
  int UpdateClusterInertiasFromFeatureInertias(void);

  // Transformation utility functions
  R3Affine OptimizedTransformation(const GSVFeatureCluster *cluster) const;
  R3Affine OptimizedTransformation(const GSVFeature *feature) const;
  R3Affine OptimizedTransformation(const GSVScanline *scanline) const;
  R3Affine OptimizedTransformation(const GSVImage *image, const R2Point& point) const;

  // Feature utility functions
  R3Point WorldPosition(const GSVFeatureCluster *cluster, RNBoolean apply_pose_transformation = TRUE) const;
  R3Vector WorldDirection(const GSVFeatureCluster *cluster, RNBoolean apply_pose_transformation = TRUE) const;
  R3Vector WorldNormal(const GSVFeatureCluster *cluster, RNBoolean apply_pose_transformation = TRUE) const;
  R3Point WorldPosition(const GSVFeature *feature, RNBoolean apply_pose_transformation = TRUE) const;
  R3Vector WorldDirection(const GSVFeature *feature, RNBoolean apply_pose_transformation = TRUE) const;
  R3Vector WorldNormal(const GSVFeature *feature, RNBoolean apply_pose_transformation = TRUE) const;
  R2Point ImagePosition(const GSVFeature *feature, const GSVFeature *image_feature, RNBoolean apply_pose_transformation = TRUE) const;
  R2Vector ImageDirection(const GSVFeature *feature, const GSVFeature *image_feature, RNBoolean apply_pose_transformation = TRUE) const;

  // Distance utility functions
  RNScalar SquaredWorldDistance(const GSVFeature *feature1, const GSVFeature *feature2, RNBoolean apply_pose_transformation = TRUE) const;
  RNScalar SquaredImageDistance(const GSVFeature *feature1, const GSVFeature *feature2, RNBoolean apply_pose_transformation = TRUE) const;
  RNScalar SquaredWorldDistance(const GSVFeatureCluster *cluster, const GSVFeature *feature, RNBoolean apply_pose_transformation = TRUE) const;
  RNScalar SquaredImageDistance(const GSVFeatureCluster *cluster, const GSVFeature *feature, RNBoolean apply_pose_transformation = TRUE) const;
  RNScalar SquaredWorldDistance(const GSVFeatureCluster *cluster1, const GSVFeatureCluster *cluster2, RNBoolean apply_pose_transformation = TRUE) const;
  RNScalar SquaredDescriptorDistance(const GSVFeature *feature1, const GSVFeature *feature2) const;

  // Correspondence scoring
  RNScalar PairScore(const GSVFeature *feature1, const GSVFeature *feature2) const;
  RNScalar CorrespondenceScore(const GSVFeature *feature1, const GSVFeature *feature2) const;
  RNScalar ClusterScore(const GSVFeatureCluster *cluster1, const GSVFeatureCluster *cluster2) const;
  RNScalar Score(const GSVFeature *feature0, const GSVFeature *feature1, 
    RNScalar max_world_distance, RNScalar max_world_distance_ratio, 
    RNScalar max_world_direction_angle, RNScalar max_world_normal_angle,
    RNScalar max_image_distance, RNScalar max_image_distance_ratio, RNScalar max_image_direction_angle, 
    RNScalar max_spin_image_descriptor_distance, RNScalar max_shape_context_descriptor_distance,
    RNScalar max_sift_descriptor_distance, RNScalar max_line_descriptor_distance,
    RNScalar max_descriptor_distance_ratio, RNScalar min_path_parameter_difference) const;

public:
  GSVScene *scene;
  RNArray<GSVFeature *> features;
  RNArray<GSVFeaturePair *> pairs;
  RNArray<GSVFeatureCorrespondence *> correspondences;
  RNArray<GSVFeatureCluster *> clusters;
  RNArray<GSVLaser *> lasers;
  RNArray<GSVCamera *> cameras;
  RNArray<GSVSegment *> segments;
  RNArray<GSVScanline *> scanlines;
  RNArray<GSVImage *> images;
  RNArray<GSVPathVertex *> vertices;
  RNScalar score;

  // Path parameters
  RNScalar max_path_vertex_spacing;

  // Feature parameters
  RNScalar sift_image_scale;

  // Pair parameters
  RNLength max_pair_world_distance;
  RNScalar max_pair_world_distance_ratio;
  RNScalar max_pair_world_direction_angle;
  RNScalar max_pair_world_normal_angle;
  RNLength max_pair_image_distance;
  RNScalar max_pair_image_distance_ratio;
  RNScalar max_pair_image_direction_angle;
  RNScalar max_pair_spin_image_descriptor_distance;
  RNScalar max_pair_shape_context_descriptor_distance;
  RNScalar max_pair_sift_descriptor_distance;
  RNScalar max_pair_line_descriptor_distance;
  RNScalar max_pair_descriptor_distance_ratio;
  RNScalar min_pair_path_parameter_difference;

  // Correspondence parameters
  RNScalar max_correspondence_world_distance;
  RNScalar max_correspondence_world_distance_ratio;
  RNScalar max_correspondence_world_direction_angle;
  RNScalar max_correspondence_world_normal_angle;
  RNScalar max_correspondence_image_distance;
  RNScalar max_correspondence_image_distance_ratio;
  RNScalar max_correspondence_image_direction_angle;
  RNScalar max_correspondence_spin_image_descriptor_distance;
  RNScalar max_correspondence_shape_context_descriptor_distance;
  RNScalar max_correspondence_sift_descriptor_distance;
  RNScalar max_correspondence_line_descriptor_distance;
  RNScalar max_correspondence_descriptor_distance_ratio;
  RNScalar min_correspondence_path_parameter_difference;

  // Correspondence search parameters
  RNScalar icp_max_world_distance_start;
  RNScalar icp_max_world_distance_end;
  RNScalar icp_max_image_distance_start;
  RNScalar icp_max_image_distance_end;

  // Optimization variable indicies
  static const int LASER_NV = 6;
  static const int CAMERA_NV = 6;
  static const int VERTEX_NV = 6;
  static const int TX = 0;
  static const int TY = 1;
  static const int TZ = 2;
  static const int RX = 3;
  static const int RY = 4;
  static const int RZ = 5;
  static const int PIXEL_FEATURE_NV = 1;
  static const int T = 0;
  static const int CLUSTER_NV = 6;
  // static const int TX = 0;
  // static const int TY = 1;
  // static const int TZ = 2;
  // static const int RX = 3;
  // static const int RY = 4;
  // static const int RZ = 5;

  // Optimization variables
  RNArray<GSVFeature *> pixel_features;
  int *pixel_feature_indices;
  int laser_v[LASER_NV];
  int laser_nv;
  int camera_v[CAMERA_NV];
  int camera_nv;
  int vertex_v[VERTEX_NV];
  int vertex_nv;
  int pixel_feature_v[PIXEL_FEATURE_NV];
  int pixel_feature_nv;
  int cluster_v[CLUSTER_NV];
  int cluster_nv;
  int expression_type;

  // Optimization parameters
  RNScalar laser_inertia_weights[LASER_NV];
  RNScalar camera_inertia_weights[CAMERA_NV];
  RNScalar vertex_inertia_weights[VERTEX_NV];
  RNScalar pixel_feature_inertia_weights[PIXEL_FEATURE_NV];
  RNScalar cluster_inertia_weights[CLUSTER_NV];
  int path_rigidity_radius;
  RNScalar path_rigidity_sigma;
  RNScalar path_rigidity_weight;
  RNScalar camera_camera_rigidity_weight;
  RNScalar laser_laser_rigidity_weight;
  RNScalar camera_laser_rigidity_weight;
  RNScalar scan_point_scan_point_correspondence_weight;
  RNScalar scan_point_scan_line_correspondence_weight;
  RNScalar scan_point_scan_plane_correspondence_weight;
  RNScalar scan_point_image_point_correspondence_weight;
  RNScalar scan_point_image_line_correspondence_weight;
  RNScalar scan_line_scan_line_correspondence_weight;
  RNScalar scan_line_scan_plane_correspondence_weight;
  RNScalar scan_line_image_point_correspondence_weight;
  RNScalar scan_line_image_line_correspondence_weight;
  RNScalar scan_plane_image_point_correspondence_weight;
  RNScalar image_point_image_point_correspondence_weight;
  RNScalar image_point_image_line_correspondence_weight;
  RNScalar image_line_image_line_correspondence_weight;
  RNScalar scan_point_image_point_correspondence_reprojection_weight;
  RNScalar scan_point_image_line_correspondence_reprojection_weight;
  RNScalar scan_line_image_point_correspondence_reprojection_weight;
  RNScalar scan_line_image_line_correspondence_reprojection_weight;
  RNScalar scan_plane_image_point_correspondence_reprojection_weight;
  RNScalar image_point_image_point_correspondence_reprojection_weight;
  RNScalar image_point_image_line_correspondence_reprojection_weight;
  RNScalar image_line_image_line_correspondence_reprojection_weight;
  RNScalar scan_direction_scan_direction_correspondence_weight;
  RNScalar image_residual_threshold;
  RNScalar scan_residual_threshold;
  RNScalar dot_product_residual_threshold;

public:
  // OPTIMIZATION STUFF

  // Optimization functions
  int ClearOptimizationVariables(void);
  int InitializeOptimizationVariables(double *x);
  int ExtractOptimizationVariables(double *x);
  int ExecuteOptimization(RNSystemOfEquations *system, double *x, int solver);

  // Inertial equations
  void AddInertiaEquations(RNSystemOfEquations *system);

  // Rigidity equations
  void AddRigidityEquations(RNSystemOfEquations *system, RNScalar rigidity); 
  void AddPathRigidityEquations(RNSystemOfEquations *system, 
    int spline_index0, int spline_index1, const R3Point& v0, const R3Point& v1, RNScalar w);
  void AddLaserRigidityEquations(RNSystemOfEquations *system, 
    GSVScanline *scanline0, GSVScanline *scanline1, RNScalar w);
  void AddCameraRigidityEquations(RNSystemOfEquations *system, 
    GSVImage *image0, GSVImage *image1, RNScalar w);

  // Correspondence equations
  void AddClusterEquations(RNSystemOfEquations *system);
  void AddCorrespondenceEquations(RNSystemOfEquations *system);
  void AddCorrespondenceEquations(RNSystemOfEquations *system, 
    GSVFeature *feature0, GSVFeature *feature1, RNScalar weight);
  void AddClusterEquations(RNSystemOfEquations *system, 
    GSVFeatureCluster *cluster0, GSVFeature *feature1, RNScalar weight);
  void AddClusterEquations(RNSystemOfEquations *system, 
    GSVFeatureCluster *cluster0, GSVFeatureCluster *cluster1, RNScalar weight);

  // Mid-level system of equation utility functions 
  void AddPointPointDistanceEquations(RNSystemOfEquations *system, 
    RNAlgebraic *px0, RNAlgebraic *py0, RNAlgebraic *pz0, 
    RNAlgebraic *px1, RNAlgebraic *py1, RNAlgebraic *pz1,
    RNScalar weight);
  void AddPointLineDistanceEquations(RNSystemOfEquations *system, 
    RNAlgebraic *px0, RNAlgebraic *py0, RNAlgebraic *pz0, 
    RNAlgebraic *px1, RNAlgebraic *py1, RNAlgebraic *pz1,
    RNAlgebraic *dx1, RNAlgebraic *dy1, RNAlgebraic *dz1,
    RNScalar weight);
  void AddPointPlaneDistanceEquations(RNSystemOfEquations *system, 
    RNAlgebraic *px0, RNAlgebraic *py0, RNAlgebraic *pz0, 
    RNAlgebraic *px1, RNAlgebraic *py1, RNAlgebraic *pz1,
    RNAlgebraic *nx1, RNAlgebraic *ny1, RNAlgebraic *nz1,
    RNScalar weight);
  void AddPointPointDistanceEquations(RNSystemOfEquations *system, 
    RNAlgebraic *px0, RNAlgebraic *py0, 
    const R2Point& p1,RNScalar weight);
  void AddPointLineDistanceEquations(RNSystemOfEquations *system, 
    RNAlgebraic *px0, RNAlgebraic *py0, 
    const R2Point& p1, const R2Vector& d1, RNScalar weight);
  void AddVectorVectorDotProductEquations(RNSystemOfEquations *system, 
    RNAlgebraic *dx0, RNAlgebraic *dy0, RNAlgebraic *dz0, 
    RNAlgebraic *dx1, RNAlgebraic *dy1, RNAlgebraic *dz1,
    RNScalar w);

  // Low-level system of equation utility functions
  int ComputeTransformedPointCoordinates(GSVFeature *feature,
    RNAlgebraic *& px, RNAlgebraic *& py, RNAlgebraic *& pz);
  int ComputeTransformedPointCoordinates(GSVFeatureCluster *cluster, const R3Point& point,
    RNAlgebraic *& px, RNAlgebraic *& py, RNAlgebraic *& pz);
  int ComputeTransformedPointCoordinates(GSVScanline *scanline, const R3Point& point,
    RNAlgebraic *& px, RNAlgebraic *& py, RNAlgebraic *& pz);
  int ComputeTransformedPointCoordinates(GSVImage *image, const R3Point& point,
    RNAlgebraic *& px, RNAlgebraic *& py, RNAlgebraic *& pz);
  int ComputeTransformedPointCoordinates(GSVImage *image, const R2Point& point, int pixel_feature_index,
    RNAlgebraic *& px, RNAlgebraic *& py, RNAlgebraic *& pz);

  int ComputeTransformedVectorCoordinates(GSVFeatureCluster *cluster, const R3Vector& vector,
    RNAlgebraic *& px, RNAlgebraic *& py, RNAlgebraic *& pz);
  int ComputeTransformedVectorCoordinates(GSVScanline *scanline, const R3Vector& vector,
    RNAlgebraic *& px, RNAlgebraic *& py, RNAlgebraic *& pz);
  int ComputeTransformedVectorCoordinates(GSVImage *image, const R3Vector& vector,
    RNAlgebraic *& px, RNAlgebraic *& py, RNAlgebraic *& pz);
  int ComputeTransformedPixelDepth(GSVImage *image, int pixel_feature_index,
    RNAlgebraic *& t);

  // Image-projection utility functions 
  int ComputeTransformedImageCoordinates(GSVImage *image, int column_index,
    RNAlgebraic *px, RNAlgebraic *py, RNAlgebraic *pz,
    RNAlgebraic *&ix, RNAlgebraic *&iy);
  int ComputeTransformedImageCoordinates(GSVFeature *image_feature, GSVFeature *feature,
    RNAlgebraic *& ix, RNAlgebraic *& iy);
  int ComputeTransformedImageCoordinates(GSVFeature *image_feature, GSVFeatureCluster *cluster,
    RNAlgebraic *& ix, RNAlgebraic *& iy);

  // Variable equation utilities
  void AddLaserVariableToEquation(RNPolynomial *equation, int laser_index, int variable_index, 
    double a = 1.0, double e = 1.0, RNBoolean already_unique = FALSE);
  void AddCameraVariableToEquation(RNPolynomial *equation, int camera_index, int variable_index, 
    double a = 1.0, double e = 1.0, RNBoolean already_unique = FALSE);
  void AddPathVertexVariableToEquation(RNPolynomial *equation, int vertex_index, int variable_index, 
    double a = 1.0, double e = 1.0, RNBoolean already_unique = FALSE);
  void AddPathVertexVariableToEquation(RNPolynomial *equation, GSVPath *path, RNScalar path_parameter, int variable_index, 
    double a = 1.0, double e = 1.0, RNBoolean already_unique = FALSE);
  void AddScanlineVariableToEquation(RNPolynomial *equation, int scanline_index, int variable_index, 
    double a = 1.0, double e = 1.0, RNBoolean already_unique = FALSE);
  void AddImageVariableToEquation(RNPolynomial *equation, int image_index, int variable_index, 
    double a = 1.0, double e = 1.0, RNBoolean already_unique = FALSE);
  void AddPixelFeatureVariableToEquation(RNPolynomial *equation, int pixel_feature_index, int variable_index, 
    double a = 1.0, double e = 1.0, RNBoolean already_unique = FALSE);
  void AddClusterVariableToEquation(RNPolynomial *equation, int cluster_index, int variable_index, 
    double a = 1.0, double e = 1.0, RNBoolean already_unique = FALSE);

  // Variable index utilities
  int LaserVariableIndex(int laser_index, int variable_index);
  int CameraVariableIndex(int camera_index, int variable_index);
  int VertexVariableIndex(int vertex_index, int variable_index);
  int PixelFeatureVariableIndex(int pixel_feature_index, int variable_index);
  int ClusterVariableIndex(int cluster_index, int variable_index);

  // Variable value utilities
  RNScalar LaserVariableValue(int laser_index, int variable_index, RNScalar *x = NULL);
  RNScalar CameraVariableValue(int camera_index, int variable_index, RNScalar *x = NULL);
  RNScalar VertexVariableValue(int vertex_index, int variable_index, RNScalar *x = NULL);
  RNScalar PixelFeatureVariableValue(int pixel_feature_index, int variable_index, RNScalar *x = NULL);
  RNScalar ClusterVariableValue(int cluster_index, int variable_index, RNScalar *x = NULL);

  // Variable assignment utilities
  int SetLaserVariableValue(int laser_index, int variable_index, RNScalar *x);
  int SetCameraVariableValue(int camera_index, int variable_index, RNScalar *x);
  int SetVertexVariableValue(int vertex_index, int variable_index, RNScalar *x);
  int SetPixelFeatureVariableValue(int pixel_feature_index, int variable_index, RNScalar *x);
  int SetClusterVariableValue(int cluster_index, int variable_index, RNScalar *x);

  // Error reporting
  void PrintErrors(RNSystemOfEquations *equations, double *x, double rigidity, int solver, 
    int print_details = 0, int print_values = 0, int print_residuals = 0, 
    int print_partial_derivatives = 0, int print_equations = 0);
};



////////////////////////////////////////////////////////////////////////
// Other type defs???
////////////////////////////////////////////////////////////////////////

struct LaserData {
  GSVLaser *laser;
  R3Vector translation;
  R3Vector rotation;
  RNScalar inertia[6];
  int index;
};

struct CameraData {
  GSVCamera *camera;
  R3Vector translation;
  R3Vector rotation;
  RNScalar inertia[6];
  int index;
};

struct SegmentData {
  GSVSegment *segment;
  GSVPath *path;
  int index;
};

struct ScanData {
  GSVScan *scan;
  GSVMesh *mesh;
  R2Grid scanline_grid;
  R2Grid position_x_grid;
  R2Grid position_y_grid;
  R2Grid position_z_grid;
  R2Grid timestamp_grid;
  R2Grid segment_id_grid;
  RNArray<R3PlanarGrid *> planar_grids;
  RNScalar timestamp;
  int index;
};

struct ImageData {
  GSVImage *image;
  RNScalar path_parameter;
  int index;
};

struct ScanlineData {
  GSVScanline *scanline;
  RNScalar path_parameter;
  int index;
};



////////////////////////////////////////////////////////////////////////
// Constants
////////////////////////////////////////////////////////////////////////

enum {
  GSV_NULL_DESCRIPTOR_TYPE,
  GSV_SPIN_IMAGE_DESCRIPTOR_TYPE,
  GSV_SIFT_DESCRIPTOR_TYPE,
  GSV_SHAPE_CONTEXT_DESCRIPTOR_TYPE,
  GSV_LINE_DESCRIPTOR_TYPE,
  GSV_NUM_DESCRIPTOR_TYPES
};

enum {
  GSV_NULL_FEATURE_TYPE,
  GSV_SCAN_POINT_FEATURE_TYPE,
  GSV_SCAN_PLANE_FEATURE_TYPE,
  GSV_SCAN_LINE_FEATURE_TYPE,
  GSV_IMAGE_POINT_FEATURE_TYPE,
  GSV_IMAGE_LINE_FEATURE_TYPE,
  GSV_NUM_FEATURE_TYPES
};



////////////////////////////////////////////////////////////////////////
// Inline functions
////////////////////////////////////////////////////////////////////////

inline GSVFeature *GSVFeaturePair::
Feature(int k) const
{
  // Return kth feature
  assert((k == 0) || (k == 1));
  return features[k];
}



inline GSVFeature *GSVFeatureCorrespondence::
Feature(int k) const
{
  // Return kth feature
  assert((k == 0) || (k == 1));
  return features[k];
}



inline int GSVFeatureCluster::
NFeatures(void) const
{
  // Return number of features
  return features.NEntries();
}



inline GSVFeature *GSVFeatureCluster::
Feature(int k) const
{
  // Return kth feature
  return features[k];
}



inline int GSVFeatureCluster::
NClusters(void) const
{
  // Return number of clusters
  return clusters.NEntries();
}



inline GSVFeatureCluster *GSVFeatureCluster::
Cluster(int k) const
{
  // Return kth cluster
  return clusters[k];
}



inline int GSVPoseOptimization::
NFeatures(void) const
{
  // Return number of features
  return features.NEntries();
}



inline GSVFeature *GSVPoseOptimization::
Feature(int k) const
{
  // Return kth feature
  return features.Kth(k);
}



inline int GSVPoseOptimization::
NPairs(void) const
{
  // Return number of correspodnences
  return pairs.NEntries();
}



inline GSVFeaturePair *GSVPoseOptimization::
Pair(int k) const
{
  // Return kth correspodnence
  return pairs.Kth(k);
}



inline int GSVPoseOptimization::
NCorrespondences(void) const
{
  // Return number of correspondences
  return correspondences.NEntries();
}



inline GSVFeatureCorrespondence *GSVPoseOptimization::
Correspondence(int k) const
{
  // Return kth correspondence
  return correspondences.Kth(k);
}



inline int GSVPoseOptimization::
NClusters(void) const
{
  // Return number of clusters
  return clusters.NEntries();
}



inline GSVFeatureCluster *GSVPoseOptimization::
Cluster(int k) const
{
  // Return kth cluster
  return clusters.Kth(k);
}



