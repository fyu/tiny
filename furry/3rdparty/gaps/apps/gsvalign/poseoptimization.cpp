// Source file for GSV pose optimization



////////////////////////////////////////////////////////////////////////
// Select solvers (guides compilation in RNMath/RNSystemOfEquations.h)
////////////////////////////////////////////////////////////////////////

// #define RN_USE_SPLM
// #define RN_USE_CERES
// #define RN_USE_MINPACK



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV/GSV.h"
#include "RNMath/RNMath.h"
#include "R3CatmullRomSpline.h"
#include "poseoptimization.h"



/////////////////////////////////////////////////////////////////////////
// GSVDescriptor utility functions
////////////////////////////////////////////////////////////////////////

GSVDescriptor::
GSVDescriptor(int descriptor_type, RNScalar *values, int nvalues)
  : descriptor_type(descriptor_type),
    values(NULL),
    nvalues(nvalues)
{
  // Copy values
  if (nvalues > 0) {
    this->values = new RNScalar [ nvalues ];
    for (int i = 0; i < nvalues; i++) {
      this->values[i] = values[i];
    }
  }
}



GSVDescriptor::
GSVDescriptor(const GSVDescriptor& descriptor)
  : descriptor_type(descriptor.descriptor_type),
    values(NULL),
    nvalues(descriptor.nvalues)
{
  // Copy values
  if (nvalues > 0) {
    values = new RNScalar [ nvalues ];
    for (int i = 0; i < nvalues; i++) {
      values[i] = descriptor.values[i];
    }
  }
}



GSVDescriptor::
~GSVDescriptor(void)
{
  // Delete values
  if (values) delete [] values;
}



RNScalar GSVDescriptor::
SquaredDistance(const GSVDescriptor& descriptor) const
{
  // Check descriptor type
  if (descriptor_type != descriptor.descriptor_type) return FLT_MAX;

  // Compute squared distance
  RNScalar sum = 0;
  for (int i = 0; i < nvalues; i++) {
    RNScalar delta = values[i] - descriptor.values[i];
    sum += delta * delta;
  }

  // Return sum of squared differences
  return sum;
}



static int
ComputeSpinImageDescriptor(GSVDescriptor& descriptor,
  const R2Grid& dh_position_x_grid, const R2Grid& dh_position_y_grid, int ix, int iy, 
  double height = 4, double radius = 2, int nstacks = 4, int nshells = 2)
{
  // Initialize descriptor
  if (descriptor.values) delete [] descriptor.values;
  descriptor.descriptor_type = GSV_NULL_DESCRIPTOR_TYPE;
  descriptor.values = NULL;
  descriptor.nvalues = 0;

  // Get spin image center position
  RNScalar cx = dh_position_x_grid.GridValue(ix, iy);
  if (cx == R2_GRID_UNKNOWN_VALUE) return 0;
  RNScalar cy = dh_position_y_grid.GridValue(ix, iy);
  if (cy == R2_GRID_UNKNOWN_VALUE) return 0;

  // Compute convenient variables
  int ixradius = (int) (radius * dh_position_x_grid.WorldToGridScaleFactor() + 0.5);
  int iyradius = (int) (height * dh_position_x_grid.WorldToGridScaleFactor() + 0.5);
  if (iyradius >= dh_position_y_grid.YResolution()) iyradius = dh_position_y_grid.YResolution();
  double radius_squared = radius * radius;

  // Allocate temporary values
  int nvalues = nshells * nstacks;
  RNScalar *values = new RNScalar [ nvalues ];
  for (int i = 0; i < nvalues; i++) values[i] = 0;

  // Count number of points in every spin image bin
  int count = 0;
  for (int i = -ixradius; i <= ixradius; i++) {
    if (ix + i < 0) continue;
    if (ix + i >= dh_position_x_grid.XResolution()) continue;
    for (int j = 0; j < iyradius; j++) {
      RNScalar x = dh_position_x_grid.GridValue(ix + i, j);
      if (x == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar y = dh_position_y_grid.GridValue(ix + i, j);
      if (y == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar dx = x - cx;
      RNScalar dy = y - cy;
      RNScalar rr = dx*dx + dy*dy;
      if (rr >= radius_squared) continue;
      RNScalar r = sqrt(rr);
      if (r >= radius) continue;
      int stack = nstacks * j / iyradius;
      int shell = (int) (nshells * (r / radius));
      int bin = stack*nshells + shell;
      assert((bin >= 0) && (bin < nvalues));
      values[bin] += 1;
      count++;
    }
  }

  // Normalize number of points in spin image bins
  if (count > 0) {
    for (int i = 0; i < nvalues; i++) {
      values[i] /= count;
    }
  }

  // Update descriptor
  descriptor.descriptor_type = GSV_SPIN_IMAGE_DESCRIPTOR_TYPE;
  descriptor.nvalues = nvalues;
  descriptor.values = values;

  // Return success
  return 1;
}



static int
ComputeShapeContextDescriptor(GSVDescriptor& descriptor, const R2Image *image, 
  RNScalar cx, RNScalar cy, int radius = 16, int nsectors = 4, int nshells = 2)
{
  // Initialize descriptor
  if (descriptor.values) delete [] descriptor.values; 
  descriptor.descriptor_type = GSV_NULL_DESCRIPTOR_TYPE;
  descriptor.values = NULL;
  descriptor.nvalues = 0;

  // Allocate temporary values
  int nvalues = nsectors * nshells;
  if (nvalues == 0) return 0;
  RNScalar *values = new RNScalar [ nvalues ];
  for (int i = 0; i < nvalues; i++) values[i] = 0;

  // Sum luminance in every shape context bin
  double sum = 0;
  int radius_squared = radius * radius;
  for (int dx = -radius; dx <= radius; dx++) {
    int ix = (int) (cx + dx + 0.5);
    if ((ix < 0) || (ix >= image->Width())) continue;
    for (int dy = -radius; dy <= radius; dy++) {
      int iy = (int) (cy + dy + 0.5);
      if ((iy < 0) || (iy >= image->Width())) continue;
      int rr = dx*dx + dy*dy;
      if (rr >= radius_squared) continue;
      RNAngle angle = atan2(dy, dx) + RN_PI;
      int sector = (int) (nsectors * angle / RN_TWO_PI);
      RNScalar r = sqrt(rr);
      if (r >= radius) continue;
      int shell = (int) (nshells * (r / radius));
      if (sector < 0) sector = 0;
      if (sector >= nsectors) sector = nsectors-1;
      if (shell < 0) shell = 0;
      if (shell >= nshells) shell = nshells-1;
      RNRgb pixel = image->PixelRGB(ix, iy);
      RNScalar weight = pixel.Luminance();
      int bin = sector*nshells + shell;
      assert((bin >= 0) && (bin < nvalues));
      values[bin] += weight;
      sum += weight;
    }
  }

  // Normalize shape context bins
  if (sum > 0) {
    for (int i = 0; i < nvalues; i++) {
      values[i] /= sum;
    }
  }

  // Update descriptor
  descriptor.descriptor_type = GSV_SHAPE_CONTEXT_DESCRIPTOR_TYPE;
  descriptor.nvalues = nvalues;
  descriptor.values = values;

  // Return success
  return 1;
}



/////////////////////////////////////////////////////////////////////////
// GSVFeature member functions
////////////////////////////////////////////////////////////////////////

GSVFeature::
GSVFeature(void)
  : feature_type(GSV_NULL_FEATURE_TYPE),
    image(NULL),
    scanline(NULL),
    scan_point_index(-1),
    scan_position(RN_UNKNOWN, RN_UNKNOWN, RN_UNKNOWN),
    scan_direction(RN_UNKNOWN, RN_UNKNOWN, RN_UNKNOWN),
    scan_normal(RN_UNKNOWN, RN_UNKNOWN, RN_UNKNOWN),
    scan_scale(RN_UNKNOWN),
    image_position(RN_UNKNOWN, RN_UNKNOWN),
    image_direction(RN_UNKNOWN, RN_UNKNOWN),
    image_scale(RN_UNKNOWN),
    image_t(RN_UNKNOWN),
    score(RN_UNKNOWN),
    descriptor(),
    group(0),
    index(-1)
{
}



GSVFeature::
GSVFeature(int feature_type, 
  GSVScanline *scanline, int scan_point_index, 
  const R3Point& scan_position, const R3Vector& scan_direction, const R3Vector& scan_normal, 
  RNScalar scan_scale, RNScalar score, int group)
  : feature_type(feature_type),
    image(NULL),
    scanline(scanline),
    scan_point_index(scan_point_index),
    scan_position(scan_position),
    scan_direction(scan_direction),
    scan_normal(scan_normal),
    scan_scale(scan_scale),
    image_position(RN_UNKNOWN, RN_UNKNOWN),
    image_direction(RN_UNKNOWN, RN_UNKNOWN),
    image_scale(RN_UNKNOWN),
    image_t(RN_UNKNOWN),
    score(score),
    descriptor(),
    group(group),
    index(-1)
{
}



GSVFeature::
GSVFeature(int feature_type, GSVImage *image, 
  const R2Point& image_position, const R2Vector& image_direction, 
  RNScalar image_scale, RNScalar score, int group)
  : feature_type(feature_type),
    image(image),
    scanline(NULL),
    scan_point_index(-1),
    scan_position(RN_UNKNOWN, RN_UNKNOWN, RN_UNKNOWN),
    scan_direction(RN_UNKNOWN, RN_UNKNOWN, RN_UNKNOWN),
    scan_normal(RN_UNKNOWN, RN_UNKNOWN, RN_UNKNOWN),
    scan_scale(RN_UNKNOWN),
    image_position(image_position),
    image_direction(image_direction),
    image_scale(image_scale),
    image_t(RN_UNKNOWN),
    score(score),
    descriptor(),
    group(group),
    index(-1)
{
  // Useful variables
  const R2Point unknown_point2(RN_UNKNOWN, RN_UNKNOWN);
  const R2Vector unknown_vector2(RN_UNKNOWN, RN_UNKNOWN);
  const double default_image_t = 10;

  // Compute scan stuff from image stuff
  if (image) {
    // Compute scanline
    if (!scanline) {
      RNScalar timestamp = image->Timestamp();
      GSVSegment *segment = image->Segment();
      if (!segment) return;
      if (segment->NScans() == 0) return;
      GSVScan *scan = segment->Scan(0);
      scanline = scan->FindScanlineBeforeTimestamp(timestamp);
    }

    // Compute scan position
    if (image_position != unknown_point2) {
      R3Ray ray = image->RayThroughUndistortedPosition(image_position);
      this->scan_position = ray.Point(default_image_t);
      this->image_t = 0;
    }

    // Compute scan direction
    if (image_direction != unknown_vector2) {
      R3Vector image_plane_direction = R3zero_vector;
      image_plane_direction += image_direction.X() * image->Pose().Right();
      image_plane_direction += image_direction.Y() * image->Pose().Up();
      this->scan_direction = image_plane_direction;
    }

    // Compute scan scale
    if (image_scale >= 0) {
      R2Point image_position2 = image_position + image_scale * image_direction;
      R3Ray ray2 = image->RayThroughUndistortedPosition(image_position2);
      R3Point scan_position2 = ray2.Point(default_image_t);
      this->scan_scale = R3Distance(this->scan_position, scan_position2);
    }
  }
}



GSVFeature::
GSVFeature(int feature_type, GSVImage *image, GSVScanline *scanline, int scan_point_index, 
  const R3Point& scan_position, const R3Vector& scan_direction, const R3Vector& scan_normal, RNScalar scan_scale,
  const R2Point& image_position, const R2Vector& image_direction, RNScalar image_scale, RNScalar image_t, 
  RNScalar score, int group)
  : feature_type(feature_type),
    image(image),
    scanline(scanline),
    scan_point_index(scan_point_index),
    scan_position(scan_position),
    scan_direction(scan_direction),
    scan_normal(scan_normal),
    scan_scale(scan_scale),
    image_position(image_position),
    image_direction(image_direction),
    image_scale(image_scale),
    image_t(image_t),
    score(score),
    descriptor(),
    group(group),
    index(-1)
{
  // Useful variables
  const R3Point unknown_point3(RN_UNKNOWN, RN_UNKNOWN, RN_UNKNOWN);
  const R3Vector unknown_vector3(RN_UNKNOWN, RN_UNKNOWN, RN_UNKNOWN);
  const R2Point unknown_point2(RN_UNKNOWN, RN_UNKNOWN);
  const R2Vector unknown_vector2(RN_UNKNOWN, RN_UNKNOWN);
  const double default_image_t = 10;

  // Compute scan stuff from image stuff
  if (image) {
    // Compute scanline
    if (!scanline) {
      RNScalar timestamp = image->Timestamp();
      GSVSegment *segment = image->Segment();
      if (!segment) return;
      if (segment->NScans() == 0) return;
      GSVScan *scan = segment->Scan(0);
      this->scanline = scan->FindScanlineBeforeTimestamp(timestamp);
    }

    // Compute scan position
    if ((scan_position == unknown_point3) && (image_position != unknown_point2)) {
      R3Ray ray = image->RayThroughUndistortedPosition(image_position);
      this->scan_position = ray.Point(default_image_t);
      this->image_t = 0;
    }

    // Compute scan direction
    if ((scan_direction == unknown_vector3) && (image_direction != unknown_vector2)) {
      R3Vector image_plane_direction = R3zero_vector;
      image_plane_direction += image_direction.X() * image->Pose().Right();
      image_plane_direction += image_direction.Y() * image->Pose().Up();
      this->scan_direction = image_plane_direction;
    }

    // Compute scan scale
    if ((scan_scale == RN_UNKNOWN) && (image_scale >= 0) && 
        (image_position != unknown_point2) && (image_direction != unknown_vector2))  {
      R2Point image_position2 = image_position + image_scale * image_direction;
      R3Ray ray2 = image->RayThroughUndistortedPosition(image_position2);
      R3Point scan_position2 = ray2.Point(default_image_t);
      this->scan_scale = R3Distance(this->scan_position, scan_position2);
    }
  }
}



void GSVFeature::
Draw(void) const
{
  // Parameters
  double radius = 0.1;

  // Check feature type
  if (feature_type == GSV_SCAN_POINT_FEATURE_TYPE) {
    R3Sphere(scan_position, radius).Draw();
  }
  else if (feature_type == GSV_SCAN_LINE_FEATURE_TYPE) {
    // Draw point
    glPointSize(5);
    glBegin(GL_POINTS);
    R3LoadPoint(scan_position);
    glEnd();
    glPointSize(1);

    // Draw line
    glLineWidth(3);
    glBegin(GL_LINES);
    R3LoadPoint(scan_position - 0.5 * scan_scale * scan_direction);
    R3LoadPoint(scan_position + 0.5 * scan_scale * scan_direction);
    glEnd();
    glLineWidth(1);
  }
  else if (feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) {
    const R3Vector& normal = scan_normal;
    int dim = normal.MinDimension();
    R3Vector axis1 = R3xyz_triad[dim] % normal;
    axis1.Normalize();
    R3Vector axis2 = axis1 % normal;
    axis2.Normalize();
    glBegin(GL_POLYGON);
    R3LoadPoint(scan_position - 1.5 * radius * axis1 - 1.5 * radius * axis2);
    R3LoadPoint(scan_position + 1.5 * radius * axis1 - 1.5 * radius * axis2);
    R3LoadPoint(scan_position + 1.5 * radius * axis1 + 1.5 * radius * axis2);
    R3LoadPoint(scan_position - 1.5 * radius * axis1 + 1.5 * radius * axis2);
    glEnd();
    glBegin(GL_POLYGON);
    R3LoadPoint(scan_position - 1.5 * radius * axis1 - 1.5 * radius * axis2);
    R3LoadPoint(scan_position - 1.5 * radius * axis1 + 1.5 * radius * axis2);
    R3LoadPoint(scan_position + 1.5 * radius * axis1 + 1.5 * radius * axis2);
    R3LoadPoint(scan_position + 1.5 * radius * axis1 - 1.5 * radius * axis2);
    glEnd();
    glLineWidth(5);
    glBegin(GL_LINES);
    R3LoadPoint(scan_position);
    R3LoadPoint(scan_position + normal);
    glEnd();
    glLineWidth(1);
  }
  else if (feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) {
    // Draw point
    glPointSize(5);
    glBegin(GL_POINTS);
    R3LoadPoint(scan_position);
    glEnd();
    glPointSize(1);
  }
  else if (feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) {
    // Draw point
    glPointSize(5);
    glBegin(GL_POINTS);
    R3LoadPoint(scan_position);
    glEnd();
    glPointSize(1);
 
    // Draw line
    glLineWidth(3);
    glBegin(GL_LINES);
    R3LoadPoint(scan_position - 0.5*scan_scale*scan_direction);
    R3LoadPoint(scan_position + 0.5*scan_scale*scan_direction);
    glEnd();
    glLineWidth(1);
  }
}



////////////////////////////////////////////////////////////////////////
// GSVFeaturePair member functions
////////////////////////////////////////////////////////////////////////

GSVFeaturePair::
GSVFeaturePair(void)
  : score(RN_UNKNOWN),
    index(-1)
{
  // Assign features
  features[0] = NULL;
  features[1] = NULL;
}



GSVFeaturePair::
GSVFeaturePair(const GSVFeaturePair& pair)
  : score(pair.score),
    index(-1)
{
  // Assign features
  features[0] = pair.features[0];
  features[1] = pair.features[1];
}



GSVFeaturePair::
GSVFeaturePair(const GSVFeatureCorrespondence& correspondence, int dummy)
  : score(correspondence.score),
    index(-1)
{
  // Assign features
  features[0] = correspondence.Feature(0);
  features[1] = correspondence.Feature(1);
}



GSVFeaturePair::
GSVFeaturePair(GSVFeature *feature0, GSVFeature *feature1, RNScalar score)
  : score(score),
    index(-1)
{
  // Assign features
  features[0] = feature0;
  features[1] = feature1;
}



void GSVFeaturePair::
Draw(void) const
{
  // Check features
  if (!features[0] || !features[1]) return;

  // Draw features
  features[0]->Draw();
  features[1]->Draw();
}



#if 0
static int
GSVComparePairs(const void *data1, const void *data2)
{
  GSVFeaturePair *pair1 = *((GSVFeaturePair **) data1);
  GSVFeaturePair *pair2 = *((GSVFeaturePair **) data2);
  if (pair1->score > pair2->score) return -1;
  else if (pair1->score < pair2->score) return -1;
  else return 0;
}
#endif



////////////////////////////////////////////////////////////////////////
// GSVFeatureCorrespondence member functions
////////////////////////////////////////////////////////////////////////

GSVFeatureCorrespondence::
GSVFeatureCorrespondence(void)
  : score(RN_UNKNOWN),
    index(-1)
{
  // Assign features
  features[0] = NULL;
  features[1] = NULL;
}



GSVFeatureCorrespondence::
GSVFeatureCorrespondence(const GSVFeaturePair& pair, int dummy)
  : score(pair.score),
    index(-1)
{
  // Assign features
  features[0] = pair.Feature(0);
  features[1] = pair.Feature(1);
}



GSVFeatureCorrespondence::
GSVFeatureCorrespondence(const GSVFeatureCorrespondence& correspondence)
  : score(correspondence.score),
    index(-1)
{
  // Assign features
  features[0] = correspondence.Feature(0);
  features[1] = correspondence.Feature(1);
}



GSVFeatureCorrespondence::
GSVFeatureCorrespondence(GSVFeature *feature0, GSVFeature *feature1, RNScalar score)
  : score(score),
    index(-1)
{
  // Assign features
  features[0] = feature0;
  features[1] = feature1;
}



void GSVFeatureCorrespondence::
Draw(void) const
{
  // Check features
  if (!features[0] || !features[1]) return;

  // Draw features
  features[0]->Draw();
  features[1]->Draw();
}



#if 0
static int
GSVCompareCorrespondences(const void *data1, const void *data2)
{
  GSVFeatureCorrespondence *correspondence1 = *((GSVFeatureCorrespondence **) data1);
  GSVFeatureCorrespondence *correspondence2 = *((GSVFeatureCorrespondence **) data2);
  if (correspondence1->score > correspondence2->score) return -1;
  else if (correspondence1->score < correspondence2->score) return -1;
  else return 0;
}
#endif



////////////////////////////////////////////////////////////////////////
// GSVFeatureCluster member functions
////////////////////////////////////////////////////////////////////////

GSVFeatureCluster::
GSVFeatureCluster(void)
  : feature_type(GSV_NULL_FEATURE_TYPE),
    scan_position(0, 0, 0),
    scan_direction(0, 0, 0),
    scan_normal(0, 0, 0),
    scan_scale(0),
    score(1),
    translation(0, 0, 0),
    rotation(0, 0, 0),
    features(),
    clusters(),
    index(-1)
{
  // Initialize inertias
  for (int i = 0; i < 6; i++) inertia[i] = 0;
}



GSVFeatureCluster::
GSVFeatureCluster(int feature_type, 
  const R3Point& scan_position, const R3Vector& scan_direction, 
  const R3Vector& scan_normal, RNScalar scan_scale, RNScalar score)
  : feature_type(feature_type),
    scan_position(scan_position),
    scan_direction(scan_direction),
    scan_normal(scan_normal),
    scan_scale(scan_scale),
    score(score),
    translation(0,0,0),
    rotation(0,0,0),
    features(),
    clusters(),
    index(-1)
{
  // Initialize inertias
  for (int i = 0; i < 6; i++) inertia[i] = 0;
}



GSVFeatureCluster::
GSVFeatureCluster(const GSVFeatureCluster& cluster)
  : feature_type(cluster.feature_type),
    scan_position(cluster.scan_position),
    scan_direction(cluster.scan_direction),
    scan_normal(cluster.scan_normal),
    scan_scale(cluster.scan_scale),
    score(cluster.score),
    translation(cluster.translation),
    rotation(cluster.rotation),
    features(cluster.features),
    clusters(cluster.clusters),
    index(-1)
{
  // Initialize inertias
  for (int i = 0; i < 6; i++) inertia[i] = cluster.inertia[i];
}



void GSVFeatureCluster::
InsertFeature(GSVFeature *feature)
{
  // Insert feature
  features.Insert(feature);
}



void GSVFeatureCluster::
RemoveFeature(GSVFeature *feature)
{
  // Remove feature
  features.Remove(feature);
}



void GSVFeatureCluster::
InsertCluster(GSVFeatureCluster *cluster)
{
  // Insert cluster
  clusters.Insert(cluster);
}



void GSVFeatureCluster::
RemoveCluster(GSVFeatureCluster *cluster)
{
  // Remove cluster
  clusters.Remove(cluster);
}



void GSVFeatureCluster::
Update(GSVPoseOptimization *optimization)
{
  // Initialize everything
  scan_position = R3zero_point;
  scan_direction = R3zero_vector;
  scan_normal = R3zero_vector;
  scan_scale = 0;
  score = 0;

  // Initialize weights
  RNScalar position_weight = 0;
  RNScalar direction_weight = 0;
  RNScalar normal_weight = 0;  

  // Aggregate feature info
  for (int j = 0; j < NFeatures(); j++) {
    GSVFeature *feature = Feature(j);
    RNScalar weight = (feature->score >= 0) ? feature->score : 1;

    if (feature->scan_position != R3unknown_point) {
      R3Point position = feature->scan_position;
      if (optimization) position.Transform(optimization->OptimizedTransformation(feature));
      scan_position += weight * position;
      position_weight += weight;
    }
    if (feature->scan_direction != R3unknown_vector) {
      R3Vector direction = feature->scan_direction;
      if (optimization) direction.Transform(optimization->OptimizedTransformation(feature));
      if ((direction_weight > 0) && (scan_direction.Dot(direction) < 0)) direction.Flip();
      scan_direction += weight * direction;
      direction_weight += weight;
    }
    if (feature->scan_normal != R3unknown_vector) {
      R3Vector normal = feature->scan_normal;
      if (optimization) normal.Transform(optimization->OptimizedTransformation(feature));
      if ((normal_weight > 0) && (scan_normal.Dot(normal) < 0)) normal.Flip();
      scan_normal += weight * normal;
      normal_weight += weight;
    }
    if (feature->scan_scale > this->scan_scale) {
      this->scan_scale = feature->scan_scale;
    }
    if (feature->score > 0) {
      this->score += feature->score;
    }
  }

  // Aggregate subcluster info
  for (int j = 0; j < NClusters(); j++) {
    GSVFeatureCluster *cluster = Cluster(j);
    RNScalar weight = (cluster->score >= 0) ? cluster->score : 1;

    if (cluster->scan_position != R3unknown_point) {
      R3Point position = cluster->scan_position;
      if (optimization) position.Transform(optimization->OptimizedTransformation(cluster));
      scan_position += weight * position;
      position_weight += weight;
    }
    if (cluster->scan_direction != R3unknown_vector) {
      R3Vector direction = cluster->scan_direction;
      if (optimization) direction.Transform(optimization->OptimizedTransformation(cluster));
      if ((direction_weight > 0) && (scan_direction.Dot(direction) < 0)) direction.Flip();
      scan_direction += weight * direction;
      direction_weight += weight;
    }
    if (cluster->scan_normal != R3unknown_vector) {
      R3Vector normal = cluster->scan_normal;
      if (optimization) normal.Transform(optimization->OptimizedTransformation(cluster));
      if ((normal_weight > 0) && (scan_normal.Dot(normal) < 0)) normal.Flip();
      scan_normal += weight * normal;
      normal_weight += weight;
    }
    if (cluster->scan_scale > this->scan_scale) {
      this->scan_scale = cluster->scan_scale;
    }
    if (cluster->score > 0) {
      this->score += cluster->score;
    }
  }
    
  // Divide by weights and apply inverse transform
  if (position_weight > 0) {
    scan_position /= position_weight;
    // if (optimization) scan_position.InverseTransform(optimization->OptimizedTransformation(this));
  }
  if (direction_weight > 0) {
    scan_direction.Normalize();
    // if (optimization) scan_direction.InverseTransform(optimization->OptimizedTransformation(this));
    scan_direction.Normalize();
  }
  if (normal_weight > 0) {
    scan_normal.Normalize();
    // if (optimization) scan_normal.InverseTransform(optimization->OptimizedTransformation(this));
    scan_normal.Normalize();
  }

  // Zero the transformation
  translation = R3zero_vector;
  rotation = R3zero_vector;

  // Update feature type
  if (feature_type == GSV_NULL_FEATURE_TYPE) {
    feature_type = GSV_SCAN_POINT_FEATURE_TYPE;
    for (int j = 0; j < NFeatures(); j++) {
      GSVFeature *feature = Feature(j);
      if (feature->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) {
        feature_type = GSV_SCAN_LINE_FEATURE_TYPE;
      }
      else if (feature->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) {
        feature_type = GSV_SCAN_LINE_FEATURE_TYPE;
      }
      else if (feature->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) {
        feature_type = GSV_SCAN_PLANE_FEATURE_TYPE;
        break;
      }
    }
    if (feature_type != GSV_SCAN_PLANE_FEATURE_TYPE) {
      for (int j = 0; j < NClusters(); j++) {
        GSVFeatureCluster *cluster = Cluster(j);
        if (cluster->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) {
          feature_type = GSV_SCAN_LINE_FEATURE_TYPE;
        }
        else if (cluster->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) {
          feature_type = GSV_SCAN_PLANE_FEATURE_TYPE;
          break;
        }
      }
    }
  }
}



void GSVFeatureCluster::
Draw(void) const
{
  // Parameters
  double radius = 0.1;

  // Check feature type
  if (feature_type == GSV_SCAN_POINT_FEATURE_TYPE) {
    R3Sphere(scan_position, radius).Draw();
  }
  else if (feature_type == GSV_SCAN_LINE_FEATURE_TYPE) {
    // Draw point
    glPointSize(5);
    glBegin(GL_POINTS);
    R3LoadPoint(scan_position);
    glEnd();
    glPointSize(1);

    // Draw line
    glLineWidth(3);
    glBegin(GL_LINES);
    R3LoadPoint(scan_position - 0.5 * scan_scale * scan_direction);
    R3LoadPoint(scan_position + 0.5 * scan_scale * scan_direction);
    glEnd();
    glLineWidth(1);
  }
  else if (feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) {
    const R3Vector& normal = scan_normal;
    int dim = normal.MinDimension();
    R3Vector axis1 = R3xyz_triad[dim] % normal;
    axis1.Normalize();
    R3Vector axis2 = axis1 % normal;
    axis2.Normalize();
    glBegin(GL_POLYGON);
    R3LoadPoint(scan_position - 1.5 * radius * axis1 - 1.5 * radius * axis2);
    R3LoadPoint(scan_position + 1.5 * radius * axis1 - 1.5 * radius * axis2);
    R3LoadPoint(scan_position + 1.5 * radius * axis1 + 1.5 * radius * axis2);
    R3LoadPoint(scan_position - 1.5 * radius * axis1 + 1.5 * radius * axis2);
    glEnd();
    glBegin(GL_POLYGON);
    R3LoadPoint(scan_position - 1.5 * radius * axis1 - 1.5 * radius * axis2);
    R3LoadPoint(scan_position - 1.5 * radius * axis1 + 1.5 * radius * axis2);
    R3LoadPoint(scan_position + 1.5 * radius * axis1 + 1.5 * radius * axis2);
    R3LoadPoint(scan_position + 1.5 * radius * axis1 - 1.5 * radius * axis2);
    glEnd();
    glLineWidth(5);
    glBegin(GL_LINES);
    R3LoadPoint(scan_position);
    R3LoadPoint(scan_position + normal);
    glEnd();
    glLineWidth(1);
  }
}



#if 0
static int
GSVCompareClusters(const void *data1, const void *data2)
{
  GSVFeatureCluster *cluster1 = *((GSVFeatureCluster **) data1);
  GSVFeatureCluster *cluster2 = *((GSVFeatureCluster **) data2);
  if (cluster1->score > cluster2->score) return -1;
  else if (cluster1->score < cluster2->score) return -1;
  else return 0;
}
#endif



////////////////////////////////////////////////////////////////////////
// GSVPath member functions
////////////////////////////////////////////////////////////////////////

GSVPath::
GSVPath(GSVSegment *segment, RNLength max_vertex_spacing)
  : segment(segment),
    spline(NULL),
    vertices(),
    index(-1)
{
  // Get scans
  if (segment->NScans() < 3) return;
  GSVScan *scan0 = segment->Scan(0);
  if (scan0->NScanlines() == 0) return;
  GSVScan *scan2 = segment->Scan(2);
  if (scan2->NScanlines() == 0) return;

  // Create array of poses
  GSVPose *poses = new GSVPose [ scan0->NScanlines() ];
  for (int ie = 0; ie < scan0->NScanlines(); ie++) {
    GSVScanline *scanline0 = scan0->Scanline(ie);
    RNScalar timestamp = scanline0->Timestamp();
    GSVPose pose0 = scanline0->Pose();
    GSVPose pose2 = scan2->Pose(timestamp);
    R3Point viewpoint = 0.5 * (pose0.Viewpoint() + pose2.Viewpoint());
    R3Quaternion orientation0 = R3Quaternion(pose0.Up(), -RN_PI_OVER_TWO) * pose0.Orientation();
    R3Quaternion orientation2 = R3Quaternion(pose2.Up(),  RN_PI_OVER_TWO) * pose2.Orientation();
    R3Quaternion orientation = R3QuaternionSlerp(orientation0, orientation2, 0.5);
    poses[ie].SetViewpoint(viewpoint);
    poses[ie].SetOrientation(orientation);
  }    

  // Count total travel distance
  RNScalar total_travel_distance = 0;
  R3Point previous_viewpoint = poses[0].Viewpoint();
  for (int ie = 0; ie < scan0->NScanlines(); ie++) {
    const R3Point& viewpoint = poses[ie].Viewpoint();
    total_travel_distance += R3Distance(viewpoint, previous_viewpoint);
    previous_viewpoint = viewpoint;
  }      
  
  // Compute spline control point parameters
  int num_vertices = 2 + (int) (total_travel_distance / max_vertex_spacing);
  if (num_vertices > scan0->NScanlines()) num_vertices = scan0->NScanlines();
  RNScalar vertex_spacing = total_travel_distance / (num_vertices - 1);

  // Create spline vertices
  RNScalar travel_distance = 0;
  previous_viewpoint = poses[0].Viewpoint();
  int previous_vertex_scanline_index = 0;
  RNScalar previous_vertex_travel_distance = 0;
  RNScalar next_vertex_travel_distance = vertex_spacing;
  RNArray<R3Point *> vertex_positions;
  for (int ie = 0; ie < scan0->NScanlines(); ie++) {
    GSVScanline *scanline0 = scan0->Scanline(ie);
    const R3Point& viewpoint = poses[ie].Viewpoint();
    travel_distance += R3Distance(viewpoint, previous_viewpoint);
    previous_viewpoint = viewpoint;
    if ((ie == 0) || (ie == scan0->NScanlines()-1) || 
        ((travel_distance >= next_vertex_travel_distance) && 
         (total_travel_distance - travel_distance > 1.5 * vertex_spacing))) {
      // Add spline vertex
      GSVPathVertex *vertex = new GSVPathVertex();
      vertex->path = this;
      vertex->pose = poses[ie];
      vertex->translation = R3zero_vector;
      vertex->rotation = R3zero_vector;
      vertex->parameter = vertex_positions.NEntries();
      vertex->timestamp = scanline0->Timestamp();
      vertices.Insert(vertex);
      vertex_positions.Insert((R3Point *) &(poses[ie].Viewpoint()));
      next_vertex_travel_distance = travel_distance + vertex_spacing;
      previous_vertex_travel_distance = travel_distance;
      previous_vertex_scanline_index = ie;
    }
  }

  // Create spline
  spline = new R3CatmullRomSpline(vertex_positions);

  // Delete array of poses
  delete [] poses;
}



GSVPath::
~GSVPath(void)
{
  // Delete spline
  if (spline) delete spline;
}



GSVPose GSVPath::
TransformedPose(RNScalar u) const
{
  // Check if fewer than two vertices
  if (NVertices() == 0) return GSVnull_pose;
  if (NVertices() == 1) return VertexPose(0);

  // Get translation and rotation
  int iu0 = (int) u;
  if (iu0 < 0) iu0 = 0;
  if (iu0 >= NVertices()-1) iu0 = NVertices()-2;
  int iu1 = iu0 + 1;
  RNScalar t = u - iu0;
  R3Vector translation0 = VertexTranslation(iu0);
  R3Vector translation1 = VertexTranslation(iu1);
  R3Vector rotation0 = VertexRotation(iu0);
  R3Vector rotation1 = VertexRotation(iu1);
  R3Vector translation = (1-t) * translation0 + t * translation1;
  R3Vector rotation = (1-t) * rotation0 + t * rotation1;

  // Return transformed pose at parameter value
  GSVPose pose = Pose(u);
  pose.Translate(translation);
  pose.Rotate(rotation);
  return pose;
}



GSVPose GSVPath::
Pose(RNScalar u) const
{
  // Check if fewer than two vertices
  if (NVertices() == 0) return GSVnull_pose;
  if (NVertices() == 1) return VertexPose(0);

  // Get translation and rotation
  int iu0 = (int) u;
  if (iu0 < 0) iu0 = 0;
  if (iu0 >= NVertices()-1) iu0 = NVertices()-2;
  int iu1 = iu0 + 1;
  RNScalar t = u - iu0;
  GSVPose pose0 = VertexPose(iu0);
  GSVPose pose1 = VertexPose(iu1);
  return GSVInterpolatedPose(pose0, pose1, t);
}



R3Vector GSVPath::
Translation(RNScalar u) const
{
  // Check if fewer than two vertices
  if (NVertices() == 0) return R3zero_vector;
  if (NVertices() == 1) return VertexTranslation(0);

  // Get translation and rotation
  int iu0 = (int) u;
  if (iu0 < 0) iu0 = 0;
  if (iu0 >= NVertices()-1) iu0 = NVertices()-2;
  int iu1 = iu0 + 1;
  RNScalar t = u - iu0;
  R3Vector translation0 = VertexTranslation(iu0);
  R3Vector translation1 = VertexTranslation(iu1);
  return (1-t) * translation0 + t * translation1;
}



R3Vector GSVPath::
Rotation(RNScalar u) const
{
  // Check if fewer than two vertices
  if (NVertices() == 0) return R3zero_vector;
  if (NVertices() == 1) return VertexRotation(0);

  // Get translation and rotation
  int iu0 = (int) u;
  if (iu0 < 0) iu0 = 0;
  if (iu0 >= NVertices()-1) iu0 = NVertices()-2;
  int iu1 = iu0 + 1;
  RNScalar t = u - iu0;
  R3Vector rotation0 = VertexRotation(iu0);
  R3Vector rotation1 = VertexRotation(iu1);
  return (1-t) * rotation0 + t * rotation1;
}



RNScalar GSVPath::
Timestamp(RNScalar u) const
{
  // Check if fewer than two vertices
  if (NVertices() == 0) return 0.0;
  if (NVertices() == 1) return VertexTimestamp(0);

  // Get translation and rotation
  int iu0 = (int) u;
  if (iu0 < 0) iu0 = 0;
  if (iu0 >= NVertices()-1) iu0 = NVertices()-2;
  int iu1 = iu0 + 1;
  RNScalar t = u - iu0;
  RNScalar timestamp0 = VertexTimestamp(iu0);
  RNScalar timestamp1 = VertexTimestamp(iu1);
  return (1-t) * timestamp0 + t * timestamp1;
}



RNScalar GSVPath::
Inertia(RNScalar u, int variable_index) const
{
  // Check if fewer than two vertices
  if (NVertices() == 0) return 0.0;
  if (NVertices() == 1) return VertexInertia(0, variable_index);

  // Get translation and rotation
  int iu0 = (int) u;
  if (iu0 < 0) iu0 = 0;
  if (iu0 >= NVertices()-1) iu0 = NVertices()-2;
  int iu1 = iu0 + 1;
  RNScalar t = u - iu0;
  RNScalar inertia0 = VertexInertia(iu0, variable_index);
  RNScalar inertia1 = VertexInertia(iu1, variable_index);
  return (1-t) * inertia0 + t * inertia1;
}



RNScalar GSVPath::
Parameter(RNScalar timestamp) const
{
  // Check if fewer than two vertices
  if (NVertices() == 0) return RNScalar(0.0);
  if (NVertices() == 1) return VertexParameter(0);

  // Return parameter u at timestamp
  int i0 = FindVertexIndexBeforeTimestamp(timestamp);
  assert((i0 >= 0) && (i0 < NVertices()-1));
  RNScalar timestamp0 = VertexTimestamp(i0);
  RNScalar timestamp1 = VertexTimestamp(i0 + 1);
  RNScalar denom = timestamp1 - timestamp0;
  if (RNIsZero(denom)) return 0.5 * (timestamp0 + timestamp1);
  RNScalar t = (timestamp - timestamp0) / denom;
  RNScalar u0 = VertexParameter(i0);
  RNScalar u1 = VertexParameter(i0 + 1);
  return RNScalar((1-t)*u0 + t*u1);
}



RNInterval GSVPath::
ParameterRange(void) const
{
  // Return range of parameter u
  GSVPathVertex *tail = vertices.Tail();
  return RNInterval(0.0, tail->parameter);
}



RNInterval GSVPath::
TimestampRange(void) const
{
  // Return range of timestamps
  GSVPathVertex *tail = vertices.Tail();
  return RNInterval(0.0, tail->timestamp);
}



int GSVPath::
NVertices(void) const
{
  // Return number of vertices
  return vertices.NEntries();
}



GSVPathVertex *GSVPath::
Vertex(int k) const
{
  // Return kth vertex
  return vertices.Kth(k);
}



const GSVPose& GSVPath::
VertexPose(int vertex_index) const
{
  // Return pose at kth vertex
  return vertices[vertex_index]->pose;
}



const R3Point& GSVPath::
VertexPosition(int vertex_index) const
{
  // Return viewpoint position at kth vertex
  return vertices[vertex_index]->pose.Viewpoint();
}



const R3Vector& GSVPath::
VertexTranslation(int vertex_index) const
{
  // Return translation at kth vertex
  return vertices[vertex_index]->translation;
}



const R3Vector& GSVPath::
VertexRotation(int vertex_index) const
{
  // Return rotation at kth vertex
  return vertices[vertex_index]->rotation;
}



RNScalar GSVPath::
VertexParameter(int vertex_index) const
{
  // Return parameter u at kth vertex
  return vertices[vertex_index]->parameter;
}



RNScalar GSVPath::
VertexTimestamp(int vertex_index) const
{
  // Return timestamp at kth vertex
  return vertices[vertex_index]->timestamp;
}



RNScalar GSVPath::
VertexInertia(int vertex_index, int variable_index) const
{
  // Return timestamp at kth vertex
  return vertices[vertex_index]->inertia[variable_index];
}



GSVPose GSVPath::
VertexTransformedPose(int vertex_index) const
{
  // Return transformed pose 
  R3Vector translation = VertexTranslation(vertex_index);
  R3Vector rotation = VertexRotation(vertex_index);
  GSVPose pose = VertexPose(vertex_index);
  pose.Translate(translation);
  pose.Rotate(rotation);
  return pose;
}



void GSVPath::
Draw(void) const
{
  // Draw spline
  if (!spline) return;
  spline->Draw();
}


int GSVPath::
FindVertexIndexBeforeTimestamp(RNScalar timestamp) const
{
  // Binary search
  if (NVertices() == 0) return -1;
  if (timestamp <= VertexTimestamp(0)) return 0;
  if (timestamp >= VertexTimestamp(NVertices()-2)) return NVertices()-2;
  return FindVertexIndexBeforeTimestamp(timestamp, 0, NVertices()-2);
}



int GSVPath::
FindVertexIndexBeforeTimestamp(RNScalar timestamp, int imin, int imax) const
{
  // Binary search
  int i = (imin + imax) / 2;
  if (i == imin) return imin;
  assert((i >= 0) && (i < imax));
  RNScalar t = VertexTimestamp(i);
  if (t > timestamp) return FindVertexIndexBeforeTimestamp(timestamp, imin, i);
  else if (t < timestamp) return FindVertexIndexBeforeTimestamp(timestamp, i, imax);
  else return i;
}



GSVPathVertex:: 
GSVPathVertex(void)
  : path(NULL),
    pose(0,0,0,0,0,0,0),
    translation(0,0,0),
    rotation(0,0,0),
    parameter(-1),
    timestamp(-1),
    index(-1)
{
  // Initialize inertias
  for (int i = 0; i < 6; i++) inertia[i] = 0;
}



////////////////////////////////////////////////////////////////////////
// GSVPoseOptimization member functions
////////////////////////////////////////////////////////////////////////

GSVPoseOptimization::
GSVPoseOptimization(GSVScene *scene, RNLength path_vertex_spacing)
  : scene(scene),
    features(),
    pairs(),
    correspondences(),
    lasers(),
    cameras(),
    segments(),
    scanlines(),
    images(),
    vertices(),
    score(-1),
    pixel_features(),
    pixel_feature_indices(NULL),
    laser_nv(0),
    camera_nv(0),
    vertex_nv(0),
    pixel_feature_nv(0),
    cluster_nv(0),
    expression_type(1),
    image_residual_threshold(0),
    scan_residual_threshold(0),
    dot_product_residual_threshold(0)
{
  // Path parameters
  max_path_vertex_spacing = path_vertex_spacing;

  // Set feature extraction parameters
  sift_image_scale = 0.125;

  // Set pair parameters
  max_pair_world_distance = 20;
  max_pair_world_distance_ratio = 2;
  max_pair_world_direction_angle = RN_PI / 8.0;
  max_pair_world_normal_angle = RN_PI / 8.0;
  max_pair_image_distance = 0;
  max_pair_image_distance_ratio = 0;
  max_pair_image_direction_angle = RN_PI / 8.0;
  max_pair_spin_image_descriptor_distance = 0.25;
  max_pair_shape_context_descriptor_distance = 0;
  max_pair_sift_descriptor_distance = 128;
  max_pair_line_descriptor_distance = 2.5;
  max_pair_descriptor_distance_ratio = 0;
  min_pair_path_parameter_difference = 50;

  // Set correspondence parameters
  max_correspondence_world_distance = 1;
  max_correspondence_world_distance_ratio = 1;
  max_correspondence_world_direction_angle = RN_PI / 8.0;
  max_correspondence_world_normal_angle = RN_PI / 8.0;
  max_correspondence_image_distance = 0;
  max_correspondence_image_distance_ratio = 0;
  max_correspondence_image_direction_angle = RN_PI / 8.0;
  max_correspondence_spin_image_descriptor_distance = 0.25;
  max_correspondence_shape_context_descriptor_distance = 0;
  max_correspondence_sift_descriptor_distance = 128;
  max_correspondence_line_descriptor_distance = 10;
  max_correspondence_descriptor_distance_ratio = 0;
  min_correspondence_path_parameter_difference = 50;

  // ICP parameters
  icp_max_world_distance_start = 20;
  icp_max_world_distance_end = 1.0;
  icp_max_image_distance_start = 0;
  icp_max_image_distance_end = 0;

  // Set optimization parameters
  RNScalar laser_translation_inertia = 1E3;
  RNScalar laser_rotation_inertia = 1E1;
  RNScalar camera_translation_inertia = 1E6;
  RNScalar camera_rotation_inertia = 1E3;
  RNScalar vertex_translation_inertia = 1E-3;
  RNScalar vertex_rotation_inertia = 1E2 * vertex_translation_inertia; // 100 meters ~ 1 radian ~ 60 degrees
  RNScalar cluster_translation_inertia = vertex_translation_inertia;
  RNScalar cluster_rotation_inertia = vertex_rotation_inertia;
  laser_inertia_weights[TX] = laser_translation_inertia;
  laser_inertia_weights[TY] = laser_translation_inertia;
  laser_inertia_weights[TZ] = laser_translation_inertia;
  laser_inertia_weights[RX] = laser_rotation_inertia;
  laser_inertia_weights[RY] = laser_rotation_inertia;
  laser_inertia_weights[RZ] = laser_rotation_inertia;
  camera_inertia_weights[TX] = camera_translation_inertia;
  camera_inertia_weights[TY] = camera_translation_inertia;
  camera_inertia_weights[TZ] = camera_translation_inertia;
  camera_inertia_weights[RX] = camera_rotation_inertia;
  camera_inertia_weights[RY] = camera_rotation_inertia;
  camera_inertia_weights[RZ] = camera_rotation_inertia;
  vertex_inertia_weights[TX] = vertex_translation_inertia;
  vertex_inertia_weights[TY] = vertex_translation_inertia;
  vertex_inertia_weights[TZ] = vertex_translation_inertia;
  vertex_inertia_weights[RX] = 1E4 * vertex_rotation_inertia;  // don't tilt
  vertex_inertia_weights[RY] = 1E4 * vertex_rotation_inertia;  // don't tilt
  vertex_inertia_weights[RZ] = vertex_rotation_inertia;
  pixel_feature_inertia_weights[T] = 1E-9;
  cluster_inertia_weights[TX] = cluster_translation_inertia;
  cluster_inertia_weights[TY] = cluster_translation_inertia;
  cluster_inertia_weights[TZ] = cluster_translation_inertia;
  cluster_inertia_weights[RX] = cluster_rotation_inertia;
  cluster_inertia_weights[RY] = cluster_rotation_inertia;
  cluster_inertia_weights[RZ] = cluster_rotation_inertia;
  path_rigidity_weight = 1;
  camera_camera_rigidity_weight = 1E6;
  laser_laser_rigidity_weight = 1E3;
  camera_laser_rigidity_weight = 1E3;
  path_rigidity_radius = (int) (4 * path_vertex_spacing + 0.5);
  path_rigidity_sigma = 0.5 * path_rigidity_radius;
  scan_point_scan_point_correspondence_weight = (2*path_rigidity_radius + 1) * path_rigidity_weight;
  scan_point_scan_line_correspondence_weight = scan_point_scan_point_correspondence_weight;
  scan_point_scan_plane_correspondence_weight = scan_point_scan_point_correspondence_weight;
  scan_line_scan_line_correspondence_weight = scan_point_scan_point_correspondence_weight;
  scan_line_scan_plane_correspondence_weight = scan_point_scan_point_correspondence_weight;

  scan_point_image_point_correspondence_weight = 0;
  scan_point_image_line_correspondence_weight = scan_point_image_point_correspondence_weight;
  scan_line_image_point_correspondence_weight = scan_point_image_point_correspondence_weight;
  scan_line_image_line_correspondence_weight = scan_point_image_point_correspondence_weight;
  scan_plane_image_point_correspondence_weight = scan_point_image_point_correspondence_weight;
  image_point_image_point_correspondence_weight = scan_point_image_point_correspondence_weight;
  image_point_image_line_correspondence_weight = image_point_image_point_correspondence_weight;
  image_line_image_line_correspondence_weight = image_point_image_point_correspondence_weight;

  RNScalar reprojection_scale = 0.01; // 1 meter = 100 pixels
  scan_point_image_point_correspondence_reprojection_weight = reprojection_scale * scan_point_scan_point_correspondence_weight;
  scan_point_image_line_correspondence_reprojection_weight = scan_point_image_point_correspondence_reprojection_weight;
  scan_line_image_point_correspondence_reprojection_weight = scan_point_image_point_correspondence_reprojection_weight;
  scan_line_image_line_correspondence_reprojection_weight = scan_point_image_point_correspondence_reprojection_weight;
  scan_plane_image_point_correspondence_reprojection_weight = scan_point_image_point_correspondence_reprojection_weight;
  image_point_image_point_correspondence_reprojection_weight = scan_point_image_point_correspondence_reprojection_weight;
  image_point_image_line_correspondence_reprojection_weight = scan_point_image_point_correspondence_reprojection_weight;
  image_line_image_line_correspondence_reprojection_weight = scan_point_image_point_correspondence_reprojection_weight;

  // RNScalar direction_scale = 50; // 1 meter ~ 10 degrees 
  RNScalar direction_scale = 0; // XXX TEMPORARY XXX
  scan_direction_scan_direction_correspondence_weight = direction_scale * scan_point_scan_point_correspondence_weight;

  // Initialize optimization variable selection
  for (int i = 0; i < LASER_NV; i++) laser_v[i] = -1;
  for (int i = 0; i < CAMERA_NV; i++) camera_v[i] = -1;
  for (int i = 0; i < VERTEX_NV; i++) vertex_v[i] = -1;
  for (int i = 0; i < PIXEL_FEATURE_NV; i++) pixel_feature_v[i] = -1;
  for (int i = 0; i < CLUSTER_NV; i++) cluster_v[i] = -1;

  // Check scene
  if (scene) {
    // Create data
    for (int ir = 0; ir < scene->NRuns(); ir++) {
      GSVRun *run = scene->Run(ir);

      // Create data for every laser
      for (int is = 0; is < run->NLasers(); is++) {
        GSVLaser *laser = run->Laser(is);
        assert(!laser->Data());
        LaserData *laser_data = new LaserData();
        laser_data->laser = laser;
        laser_data->translation = R3zero_vector;
        laser_data->rotation = R3zero_vector;
        for (int k = 0; k < LASER_NV; k++) laser_data->inertia[k] = 0;
        laser_data->index = lasers.NEntries();
        laser->SetData(laser_data);
        lasers.Insert(laser);
      }

      // Create data for every camera
      for (int is = 0; is < run->NCameras(); is++) {
        GSVCamera *camera = run->Camera(is);
        assert(!camera->Data());
        CameraData *camera_data = new CameraData();
        camera_data->camera = camera;
        camera_data->translation = R3zero_vector;
        camera_data->rotation = R3zero_vector;
        for (int k = 0; k < CAMERA_NV; k++) camera_data->inertia[k] = 0;
        camera_data->index = cameras.NEntries();
        camera->SetData(camera_data);
        cameras.Insert(camera);
      }

      // Create data for every segment
      for (int is = 0; is < run->NSegments(); is++) {
        GSVSegment *segment = run->Segment(is);
        
        // Create path
        GSVPath *path = new GSVPath(segment, max_path_vertex_spacing);
        if (!path || (path->NVertices() == 0)) {
          fprintf(stderr, "Unable to create path for segment %d\n", is);
          continue;
        }

        // Insert vertices from path
        for (int i = 0; i < path->NVertices(); i++) {
          GSVPathVertex *vertex = path->Vertex(i);
          vertex->index = vertices.NEntries();
          vertices.Insert(vertex);
        }

        // Create segment data
        assert(!segment->Data());
        SegmentData *segment_data = new SegmentData();
        segment_data->segment = segment;
        segment_data->path = path;
        segment_data->index = segments.NEntries();
        segment->SetData((void *) segment_data);
        segments.Insert(segment);

        // Create scanline data 
        for (int ia = 0; ia < segment->NScans(); ia++) {
          GSVScan *scan = segment->Scan(ia);
          if (scan->NScanlines() == 0) continue;
          for (int ie = 0; ie < scan->NScanlines(); ie++) {
            GSVScanline *scanline = scan->Scanline(ie);
            R3Point viewpoint = scanline->Pose().Viewpoint();
            RNScalar timestamp = scanline->Timestamp();
            RNScalar path_parameter = path->Parameter(timestamp);
            R3Point path_viewpoint = path->Pose(path_parameter).Viewpoint();
            assert(!scanline->Data());
            ScanlineData *scanline_data = new ScanlineData();
            scanline_data->scanline = scanline;
            scanline_data->path_parameter = path_parameter;
            scanline_data->index = scanlines.NEntries();
            scanline->SetData((void *) scanline_data);
            scanlines.Insert(scanline);
          }
        }

        // Create image data
        for (int ip = 0; ip < segment->NPanoramas(); ip++) {
          GSVPanorama *panorama = segment->Panorama(ip);
          for (int ii = 0; ii < panorama->NImages(); ii++) {
            GSVImage *image = panorama->Image(ii);
            R3Point viewpoint = image->Pose().Viewpoint();
            RNScalar timestamp = image->Timestamp();
            RNScalar path_parameter = path->Parameter(timestamp);
            R3Point path_viewpoint = path->Pose(path_parameter).Viewpoint();
            assert(!image->Data());
            ImageData *image_data = new ImageData();
            image_data->image = image;
            image_data->path_parameter = path_parameter;
            image_data->index = images.NEntries();
            image->SetData((void *) image_data);
            images.Insert(image);
          }
        }
      }
    }
  }
}



GSVPoseOptimization::
~GSVPoseOptimization(void)
{
  // XXX NEED TO FILL THIS IN XXX
}



////////////////////////////////////////////////////////////////////////
// Manipulation functions
////////////////////////////////////////////////////////////////////////

void GSVPoseOptimization::
InsertFeature(GSVFeature *feature)
{
  // Insert into list of features
  assert(feature->index == -1);
  feature->index = features.NEntries();
  features.Insert(feature);
}



void GSVPoseOptimization::
RemoveFeature(GSVFeature *feature)
{
  // Remove correspondences using feature
  RNArray<GSVFeatureCorrespondence *> corr = correspondences;
  for (int i = 0; i < corr.NEntries(); i++) {
    GSVFeatureCorrespondence *correspondence = corr.Kth(i);
    if ((correspondence->Feature(0) == feature) || (correspondence->Feature(1) == feature)) {
      RemoveCorrespondence(correspondence);
      delete correspondence;
    }
  }

  // Remove from list of features
  assert(feature->index >= 0);
  RNArrayEntry *entry = features.KthEntry(feature->index);
  GSVFeature *tail = features.Tail();
  tail->index = feature->index;
  features.EntryContents(entry) = tail;
  features.RemoveTail();
  feature->index = -1;
}



void GSVPoseOptimization::
InsertPair(GSVFeaturePair *pair)
{
  // Insert pair
  assert(pair->index == -1);
  pair->index = pairs.NEntries();
  pairs.Insert(pair);
}


void GSVPoseOptimization::
RemovePair(GSVFeaturePair *pair)
{
  // Remove pair
  RNArrayEntry *entry = pairs.KthEntry(pair->index);
  GSVFeaturePair *tail = pairs.Tail();
  tail->index = pair->index;
  pairs.EntryContents(entry) = tail;
  pairs.RemoveTail();
  pair->index = -1;
}



void GSVPoseOptimization::
InsertCorrespondence(GSVFeatureCorrespondence *correspondence)
{
  // Insert correspondence
  assert(correspondence->index == -1);
  correspondence->index = correspondences.NEntries();
  correspondences.Insert(correspondence);
}


void GSVPoseOptimization::
RemoveCorrespondence(GSVFeatureCorrespondence *correspondence)
{
  // Remove correspondence
  RNArrayEntry *entry = correspondences.KthEntry(correspondence->index);
  GSVFeatureCorrespondence *tail = correspondences.Tail();
  tail->index = correspondence->index;
  correspondences.EntryContents(entry) = tail;
  correspondences.RemoveTail();
  correspondence->index = -1;
}



void GSVPoseOptimization::
InsertCluster(GSVFeatureCluster *cluster)
{
  // Insert cluster
  assert(cluster->index == -1);
  cluster->index = clusters.NEntries();
  clusters.Insert(cluster);
}


void GSVPoseOptimization::
RemoveCluster(GSVFeatureCluster *cluster)
{
  // Remove cluster
  RNArrayEntry *entry = clusters.KthEntry(cluster->index);
  GSVFeatureCluster *tail = clusters.Tail();
  tail->index = cluster->index;
  clusters.EntryContents(entry) = tail;
  clusters.RemoveTail();
  cluster->index = -1;
}



int GSVPoseOptimization::
TruncateFeatures(int max_features) 
{
  // Delete all features beyond the first max_features
  while (features.NEntries() > max_features) {
    GSVFeature *feature = features.Tail();
    features.RemoveTail();
    delete feature;
  }

  // Return success
  return 1;
}



int GSVPoseOptimization::
TruncatePairs(int max_pairs) 
{
  // Delete all pairs beyond the first max_pairs
  while (pairs.NEntries() > max_pairs) {
    GSVFeaturePair *pair = pairs.Tail();
    pairs.RemoveTail();
    delete pair;
  }

  // Return success
  return 1;
}



int GSVPoseOptimization::
TruncateCorrespondences(int max_correspondences) 
{
  // Delete all correspondences beyond the first max_correspondences
  while (correspondences.NEntries() > max_correspondences) {
    GSVFeatureCorrespondence *correspondence = correspondences.Tail();
    correspondences.RemoveTail();
    delete correspondence;
  }

  // Return success
  return 1;
}



int GSVPoseOptimization::
TruncateClusters(int max_clusters) 
{
  // Delete all clusters beyond the first max_clusters
  while (clusters.NEntries() > max_clusters) {
    GSVFeatureCluster *cluster = clusters.Tail();
    clusters.RemoveTail();
    delete cluster;
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Harris corner feature creation functions
////////////////////////////////////////////////////////////////////////

#include "klt/klt.h"

static int 
CreateCornerFeatures(GSVPoseOptimization *optimization, GSVImage *image)
{
  // Parameters
 int nfeatures = 1024;
 RNScalar border_fraction = 0.05;
 RNScalar min_spacing_fraction = 0.01;
  
  // Get image
  R2Image *img = image->UndistortedImage();
  if (!img) return 0;

  // Create image in KLT format
  int width = img->Width();
  int height = img->Height();
  unsigned char *pixels = new unsigned char [ width * height ];
  unsigned char *pixelp = pixels;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      RNRgb rgb = img->PixelRGB(i, j);
      int value = (int) (rgb.Luminance() * 256.0);
      if (value > 255) value = 255;
      else if (value < 0) value = 0;
      (*pixelp++) = value;
    }
  }

  // Create KLT context
  KLTSetVerbosity(0);
  KLT_TrackingContext tc = KLTCreateTrackingContext();
  tc->writeInternalImages = FALSE;
  tc->affineConsistencyCheck = -1;  /* set this to 2 to turn on affine consistency check */
  tc->borderx = (int) (border_fraction * width);
  tc->bordery = (int) (border_fraction * height);
  tc->mindist = (int) (min_spacing_fraction * width);

  // Create KLT feature list
  KLT_FeatureList fl = KLTCreateFeatureList(nfeatures);
  KLTSelectGoodFeatures(tc, pixels, width, height, fl);

  // Find maximum score
  RNScalar max_score = 0;
  for (int i = 0; i < fl->nFeatures; i++) {
    RNScalar score = fl->feature[i]->val;
    if (score > max_score) max_score = score;
  }

  // Create features
  if (max_score > 0) {
    for (int i = 0; i < fl->nFeatures; i++) {
      if (fl->feature[i]->val <= 0) continue;
      RNScalar x = fl->feature[i]->x;
      RNScalar y = fl->feature[i]->y;
      RNScalar score = fl->feature[i]->val / max_score;
      GSVFeature *feature = new GSVFeature(GSV_IMAGE_POINT_FEATURE_TYPE, image, R2Point(x,y), R2zero_vector, 1, score);
      if (!feature) { fprintf(stderr, "Unable to create klt corner feature\n"); return 0; }
      ComputeShapeContextDescriptor(feature->descriptor, img, x, y);
      optimization->InsertFeature(feature);
    }
  }

  // Delete KLT stuff
  KLTFreeFeatureList(fl);
  KLTFreeTrackingContext(tc);
  delete pixels;

  // Delete image
  delete img;

  // Return success
  return 1;
}



int GSVPoseOptimization::
CreateImageCornerFeatures(void)
{
  // Make temporary directory
  char tmp_directory[4096], mkdir_cmd[4096];
  sprintf(tmp_directory, "%s/image_corner_features", scene->CacheDataDirectoryName());
  sprintf(mkdir_cmd, "mkdir -p %s", tmp_directory); 
  system(mkdir_cmd);

  // Make run directories
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    char run_directory[4096], mkdir_cmd[4096];
    sprintf(run_directory, "%s/%s", tmp_directory, run->Name());
    sprintf(mkdir_cmd, "mkdir -p %s", run_directory); 
    system(mkdir_cmd);
  }

  // Create corner features for every image in scene
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    char run_directory[4096];
    sprintf(run_directory, "%s/%s", tmp_directory, run->Name());
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ip = 0; ip < segment->NPanoramas(); ip++) {
        GSVPanorama *panorama = segment->Panorama(ip);
        for (int ii = 0; ii < panorama->NImages(); ii++) {
          GSVImage *image = panorama->Image(ii);
          int saved_nfeatures = NFeatures();
          if (!CreateCornerFeatures(this, image)) return 0;
          printf("%d %d %d %d : %d\n", ir, is, ip, ii, NFeatures() - saved_nfeatures);
        }
      }
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Sift feature creation functions
////////////////////////////////////////////////////////////////////////

static int 
CreateSiftFeatures(GSVPoseOptimization *optimization, GSVImage *image, const char *directory_name)
{
  // Get filenames
  char pgm_name[4096], sift_name[4096];
  int image_index = image->PanoramaIndex();
  GSVPanorama *panorama = image->Panorama();
  int panorama_index = panorama->SegmentIndex();
  GSVSegment *segment = panorama->Segment();
  int segment_index = segment->RunIndex();
  sprintf(pgm_name, "%s/%02d_%06d_%02d.pgm", directory_name, segment_index, panorama_index, image_index);
  sprintf(sift_name, "%s/%02d_%06d_%02d.sft", directory_name, segment_index, panorama_index, image_index);

  // Check if sift file already exists 
  if (!RNFileExists(sift_name)) {
    // Get image
    R2Image *img = image->UndistortedImage();
    if (!img) return 0;
    
    // Scale image (this could be better)
    RNScalar scale = optimization->sift_image_scale;
    int width = scale * img->Width();
    int height = scale * img->Height();
    R2Image scaled_img(width, height, 3);
    for (int i = 0; i < width; i++) {
      for (int j = 0; j < height; j++) {
        int x = (int) (i / scale + 0.5);
        int y = (int) (j / scale + 0.5);
        if ((x < 0) || (x >= img->Width())) continue;
        if ((y < 0) || (y >= img->Height())) continue;
        RNRgb pixel = img->PixelRGB(x, y);
        scaled_img.SetPixelRGB(i, j, pixel);
      }
    }
    

    // Write pgm image
    if (!scaled_img.Write(pgm_name)) {
      delete img;
      return 0;
    }

    // Run program to extract sift features from image
    char sift_command[4096];
    sprintf(sift_command, "sift < %s > %s", pgm_name, sift_name);
    system(sift_command);
    sprintf(sift_command, "rm -f %s", pgm_name);
    system(sift_command);
    
    // Delete image
    delete img;
  }

  // Read sift file
  return optimization->ReadSiftFile(image, sift_name);

  // Return success
  return 1;
}



int GSVPoseOptimization::
CreateImageSiftFeatures(void)
{
  // Make temporary directory
  char tmp_directory[4096], mkdir_cmd[4096];
  sprintf(tmp_directory, "%s/image_sift_features", scene->CacheDataDirectoryName());
  sprintf(mkdir_cmd, "mkdir -p %s", tmp_directory); 
  system(mkdir_cmd);

  // Make run directories
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    char run_directory[4096], mkdir_cmd[4096];
    sprintf(run_directory, "%s/%s", tmp_directory, run->Name());
    sprintf(mkdir_cmd, "mkdir -p %s", run_directory); 
    system(mkdir_cmd);
  }

  // Create sift features for every image in scene
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    char run_directory[4096];
    sprintf(run_directory, "%s/%s", tmp_directory, run->Name());
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ip = 0; ip < segment->NPanoramas(); ip++) {
        GSVPanorama *panorama = segment->Panorama(ip);
        for (int ii = 0; ii < panorama->NImages(); ii++) {
          GSVImage *image = panorama->Image(ii);
          if (!CreateSiftFeatures(this, image, run_directory)) return 0;
        }
      }
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Line feature creation functions
////////////////////////////////////////////////////////////////////////

static int 
CreateLineFeatures(GSVPoseOptimization *optimization, GSVImage *image, const char *directory_name)
{
  // Get filenames
  char image_name[4096], line_name[4096];
  int image_index = image->PanoramaIndex();
  GSVPanorama *panorama = image->Panorama();
  int panorama_index = panorama->SegmentIndex();
  GSVSegment *segment = panorama->Segment();
  int segment_index = segment->RunIndex();
  sprintf(image_name, "%s/%02d_%06d_%02d.pgm", directory_name, segment_index, panorama_index, image_index);
  sprintf(line_name, "%s/%02d_%06d_%02d.lin", directory_name, segment_index, panorama_index, image_index);

  // Check if line file already exists 
  if (!RNFileExists(line_name)) {
    // Get image
    R2Image *img = image->UndistortedImage();
    if (!img) return 0;

    // Write image to temporary directory
    if (!img->Write(image_name)) {
      delete img;
      return 0;
    }

    // Run program to extract line features from image
    char line_command[4096];
    sprintf(line_command, "pfm2lin %s %s", image_name, line_name);
    system(line_command);
    sprintf(line_command, "rm -f %s", image_name);
    system(line_command);
    
    // Delete image
    delete img;
  }

  // Read line file
  return optimization->ReadLineFile(image, line_name);

  // Return success
  return 1;
}



int GSVPoseOptimization::
CreateImageLineFeatures(void)
{
  // Make temporary directory
  char tmp_directory[4096], mkdir_cmd[4096];
  sprintf(tmp_directory, "%s/image_line_features", scene->CacheDataDirectoryName());
  sprintf(mkdir_cmd, "mkdir -p %s", tmp_directory); 
  system(mkdir_cmd);

  // Make run directories
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    char run_directory[4096], mkdir_cmd[4096];
    sprintf(run_directory, "%s/%s", tmp_directory, run->Name());
    sprintf(mkdir_cmd, "mkdir -p %s", run_directory); 
    system(mkdir_cmd);
  }

  // Create line features for every image in scene
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    char run_directory[4096];
    sprintf(run_directory, "%s/%s", tmp_directory, run->Name());
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ip = 0; ip < segment->NPanoramas(); ip++) {
        GSVPanorama *panorama = segment->Panorama(ip);
        for (int ii = 0; ii < panorama->NImages(); ii++) {
          GSVImage *image = panorama->Image(ii);
          if (!CreateLineFeatures(this, image, run_directory)) return 0;
        }
      }
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Pole feature creation functions
////////////////////////////////////////////////////////////////////////

static int 
CompareScalarPointers(const void *data1, const void *data2)
{
  const RNScalar *ptr1 = *((const RNScalar **) data1);
  const RNScalar *ptr2 = *((const RNScalar **) data2);
  RNScalar value1 = *ptr1;
  RNScalar value2 = *ptr2;
  if (value1 > value2) return -1;
  else if (value2 > value1) return 1;
  else return 0;
}



static void 
FilterVotes(R2Grid& vote_grid, const R2Grid& depth_grid, 
  int min_iy, int max_iy, RNScalar y_sigma, RNScalar depth_sigma)
{
  // Make copy of vote grid
  R2Grid vote_copy(vote_grid);
  vote_grid.Clear(R2_GRID_UNKNOWN_VALUE);
  
  // Get convenient variables
  const RNScalar sqrt_two_pi = sqrt(RN_TWO_PI);
  double y_fac = 1.0 / (sqrt_two_pi * y_sigma);
  double y_denom = -2.0 * y_sigma * y_sigma;
  double d_fac = 1.0 / (sqrt_two_pi * depth_sigma);
  double d_denom = -2.0 * depth_sigma * depth_sigma;
  int y_radius = (int) (3.0*y_sigma) + 1;
  RNScalar depth_radius = 3 * depth_sigma;

  // Set every vote to be bilateral filter of vertical region
  for (int ix = 0; ix < vote_grid.XResolution(); ix++) {
    for (int iy = min_iy; iy <= max_iy; iy++) {
      RNScalar vote = vote_copy.GridValue(ix, iy);
      if (vote == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar depth = depth_grid.GridValue(ix, iy);
      if (depth == R2_GRID_UNKNOWN_VALUE) continue;

      // Compute blurred vote value
      RNScalar sum = 0;
      RNScalar weight = 0;
      for (int y = iy - y_radius; y <= iy + y_radius; y++) {
        if ((y < 0) || (y >= vote_grid.YResolution())) continue;
        RNScalar vote_sample = vote_copy.GridValue(ix, y);
        if (vote_sample == R2_GRID_UNKNOWN_VALUE) continue;
        RNScalar depth_sample = depth_grid.GridValue(ix, y);
        if (depth_sample == R2_GRID_UNKNOWN_VALUE) continue;
        int dy = y - iy;
        RNScalar dd = depth - depth_sample;
        if (dd > depth_radius) continue;
        RNScalar w = y_fac * exp(dy*dy/y_denom) * d_fac * exp(dd*dd/d_denom);
        sum += w * vote_sample;
        weight += w;
      }
      
      // Set grid value
      // if (weight > 0) vote_grid.SetGridValue(ix, iy, sum / weight);
      vote_grid.SetGridValue(ix, iy, sum);
    }
  }
}



// THIS IS A BETTER WAY TO FIND POLES
// pfm2pfm 00_00_DA_Density.grd mask.grd -threshold 5 0 1
// pfm2pfm 00_00_DA_Lambda21.grd foo.grd -v 
//   -negate -add 1 
//   -multiply_grid 00_00_DA_PrincipleAxis1Z.grd 
//   -threshold 0.8 0 keep 
//   -mask_grid mask.grd 



int GSVPoseOptimization::
CreateScanPoleFeatures(void)
{
  // Parameters
  const RNScalar min_height = 0.25;
  const RNScalar max_height = 4.0;
  const RNScalar max_depth = 10;
  const RNScalar min_depth_difference = 0.5;
  const RNScalar expected_pole_radius = 0.25;
  const RNScalar expected_pole_height = 2.0;
  const RNScalar min_feature_spacing = 1;
  const RNScalar min_coverage = 0.5;
  const RNScalar max_travel_distance_for_duplicates = 50;

  // Check scene
  if (!scene) return 0;

  // Get image directory name
  char laser_image_directory[4096];
  sprintf(laser_image_directory, "%s/laser_images", scene->CacheDataDirectoryName());

  // Create pole features
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (scan->NScanlines() == 0) continue;
        if (ia == 1) continue;

        // Read points (needed for EstimatedGroundZ)
        if (!scan->ReadPoints()) {
          fprintf(stderr, "Unable to read points for scan\n");
          return 0;
        }

        // Read DH grids
        char image_name[4096];
        R2Grid dh_viewpoint_depth_grid;
        R2Grid dh_scanline_index_grid;
        R2Grid dh_point_index_grid;
        R2Grid dh_position_x_grid;
        R2Grid dh_position_y_grid;
        sprintf(image_name, "%s/%s/%02d_%02d_DH_Scanline.grd", laser_image_directory, run->Name(), is, ia);
        if (!dh_scanline_index_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DH_PointIndex.grd", laser_image_directory, run->Name(), is, ia);
        if (!dh_point_index_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DH_ViewpointDepth.grd", laser_image_directory, run->Name(), is, ia);
        if (!dh_viewpoint_depth_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DH_PositionX.grd", laser_image_directory, run->Name(), is, ia);
        if (!dh_position_x_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DH_PositionY.grd", laser_image_directory, run->Name(), is, ia);
        if (!dh_position_y_grid.Read(image_name)) return 0;

        // Create DH grid containing vertical valley indicator function
        int pole_radius = (int) (dh_viewpoint_depth_grid.WorldToGridScaleFactor() * expected_pole_radius) + 1;
        R2Grid dh_vertical_valley_grid(dh_viewpoint_depth_grid);
        dh_vertical_valley_grid.Clear(-1);
        for (int i = 0; i < dh_viewpoint_depth_grid.XResolution(); i++) {
          for (int j = 0; j < dh_viewpoint_depth_grid.YResolution(); j++) {
            int left_ix = i - pole_radius;
            int right_ix = i + pole_radius;
            if (left_ix < 0) left_ix = 0;
            if (right_ix > dh_viewpoint_depth_grid.XResolution()-1) 
              right_ix = dh_viewpoint_depth_grid.XResolution()-1;
            RNScalar value = dh_viewpoint_depth_grid.GridValue(i, j);
            if (value == R2_GRID_UNKNOWN_VALUE) continue;
            if (value > max_depth) continue;
            RNScalar left_value = dh_viewpoint_depth_grid.GridValue(left_ix, j);
            if (left_value == R2_GRID_UNKNOWN_VALUE) continue;
            RNScalar right_value = dh_viewpoint_depth_grid.GridValue(right_ix, j);
            if (right_value == R2_GRID_UNKNOWN_VALUE) continue;
            if (left_value - value < min_depth_difference) continue;
            if (right_value - value < min_depth_difference) continue;
            dh_vertical_valley_grid.SetGridValue(i, j, 1);
          }
        }

        // Create DH vote grid
        R2Grid dh_vote_grid(dh_vertical_valley_grid);
        int y1 = (int) (dh_vote_grid.WorldToGridScaleFactor() * min_height);
        int y2 = (int) (dh_vote_grid.WorldToGridScaleFactor() * max_height);
        if (y1 >= dh_vote_grid.YResolution()) y1 = dh_vote_grid.YResolution()-1;
        if (y2 >= dh_vote_grid.YResolution()) y2 = dh_vote_grid.YResolution()-1;
        double y_sigma = 0.5 * dh_vote_grid.WorldToGridScaleFactor() * expected_pole_height;
        double depth_sigma = expected_pole_radius;
        FilterVotes(dh_vote_grid, dh_viewpoint_depth_grid, y1, y2, y_sigma, depth_sigma);

        // Find local maxima of dh_vote_grid
        RNArray<const RNScalar *> local_maxima;
        R2Grid dh_maxima_grid(dh_vote_grid);
        dh_maxima_grid.MaskNonMaxima(2 * expected_pole_radius * dh_maxima_grid.WorldToGridScaleFactor());
        dh_maxima_grid.Threshold(min_coverage, R2_GRID_UNKNOWN_VALUE, R2_GRID_KEEP_VALUE);
        for (int i = 0; i < dh_maxima_grid.XResolution(); i++) {
          // Find maximum value in column
          int maximum_j = -1;
          RNScalar maximum_value = -FLT_MAX;
          for (int j = y1; j <= y2; j++) {
            RNScalar value = dh_maxima_grid.GridValue(i, j);
            if (value == R2_GRID_UNKNOWN_VALUE) continue;
            if (value == 0) continue;
            if (value < maximum_value) continue;
            maximum_value = value;
            maximum_j = j;
          }

          // Insert local maximum
          if (maximum_j >= 0) {
            const RNScalar *local_maximum = &dh_maxima_grid(i, maximum_j);
            local_maxima.Insert(local_maximum);
          }
        }

        // Sort local maxima
        local_maxima.Sort(CompareScalarPointers);

        // Create features at local maxima with minimum spacing
        RNArray<GSVFeature *> created_features;
        RNLength min_feature_spacing_squared = min_feature_spacing * min_feature_spacing;
        for (int k = 0; k < local_maxima.NEntries(); k++) {
          // Get local maximum info
          const RNScalar *local_maximum = local_maxima.Kth(k);
          int grid_index = local_maximum - dh_maxima_grid.GridValues();
          assert((grid_index >= 0) && (grid_index < dh_maxima_grid.NEntries()));

          // Determine the x,y position
          RNScalar x = dh_position_x_grid.GridValue(grid_index);
          RNScalar y = dh_position_y_grid.GridValue(grid_index);

          // Determine travel distance
          int scanline_index = (int) (dh_scanline_index_grid.GridValue(grid_index) + 0.5);
          GSVScanline *scanline = scan->Scanline(scanline_index);
          RNLength travel_distance = scanline->TravelDistance();

          // Check if there is a previously created (stronger) feature within minimum spacing
          RNBoolean well_spaced = TRUE;
          R2Point local_maximum_position(x, y);
          for (int j = 0; j < created_features.NEntries(); j++) {
            GSVFeature *created_feature = created_features.Kth(j);
            GSVScanline *created_feature_scanline = created_feature->scanline;
            RNLength created_feature_travel_distance = created_feature_scanline->TravelDistance();
            RNScalar travel_distance_delta = fabs(created_feature_travel_distance - travel_distance);
            if (travel_distance_delta < max_travel_distance_for_duplicates) {
              R2Point created_feature_position(created_feature->scan_position.X(), created_feature->scan_position.Y());
              RNScalar distance_squared = R2SquaredDistance(local_maximum_position, created_feature_position);
             if (distance_squared < min_feature_spacing_squared) { well_spaced = FALSE; break; }
            }
          }

          // Create feature if not too close to previous (stronger) one
          if (well_spaced) {
            int grid_i = -1, grid_j = -1;
            dh_maxima_grid.IndexToIndices(grid_index, grid_i, grid_j);
            int scanline_index = (int) (dh_scanline_index_grid.GridValue(grid_index) + 0.5);
            int point_index = (int) (dh_point_index_grid.GridValue(grid_index) + 0.5);
            GSVScanline *scanline = scan->Scanline(scanline_index);
            RNScalar z = scanline->EstimatedGroundZ();
            if (z == RN_UNKNOWN) continue;
            R3Point position(x, y, z); 
            RNScalar score = dh_maxima_grid.GridValue(grid_index);
            GSVFeature *feature = new GSVFeature(GSV_SCAN_POINT_FEATURE_TYPE, scanline, point_index, position, R3zero_vector, R3zero_vector, 1, score);
            ComputeSpinImageDescriptor(feature->descriptor, dh_position_x_grid, dh_position_y_grid, grid_i, grid_j);
            created_features.Insert(feature);
            InsertFeature(feature);
          }
        }

        // Release points
        if (!scan->ReleasePoints()) {
          fprintf(stderr, "Unable to release points for scan\n");
          return 0;
        }

#if 0
        // Write grids for debugging
        char buffer[1024];
        sprintf(buffer, "tmp/%d_%d_%d_dh_vertical_valley.grd", ir, is, ia);
        dh_vertical_valley_grid.WriteFile(buffer);
        sprintf(buffer, "tmp/%d_%d_%d_dh_vote.grd", ir, is, ia);
        dh_vote_grid.WriteFile(buffer);
        sprintf(buffer, "tmp/%d_%d_%d_dh_maxima.grd", ir, is, ia);
        dh_maxima_grid.WriteFile(buffer);
#endif
      }
    }
  }

  // Return success
  return 1;
}



#if 0

int GSVPoseOptimization::
CreateScanPoleFeatures(void)
{
  // THIS IS AN ALTERNATIVE STRATEGY -- NOT FULLY IMPLEMENTED

  // Check scene
  if (!scene) return 0;

  // Get image directory name
  char laser_image_directory[4096];
  sprintf(laser_image_directory, "%s/laser_images", scene->CacheDataDirectoryName());

  // Create pole features
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (scan->NScanlines() == 0) continue;
        if (ia == 1) continue;

        // Read grids
        char image_name[4096];
        R2Grid density_grid;
        R2Grid lambda21_grid;
        R2Grid principle_axis_3_z_grid;
        R2Grid normal_z_grid;
        R2Grid viewpoint_depth_grid;
        sprintf(image_name, "%s/%s/%02d_%02d_DH_Density.grd", laser_image_directory, run->Name(), is, ia);
        if (!density_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DH_Lambda21.grd", laser_image_directory, run->Name(), is, ia);
        if (!lambda21_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DH_PrincipleAxis3Z.grd", laser_image_directory, run->Name(), is, ia);
        if (!principle_axis_3_z_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DH_NormalZ.grd", laser_image_directory, run->Name(), is, ia);
        if (!normal_z_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DH_ViewpointDepth.grd", laser_image_directory, run->Name(), is, ia);
        if (!viewpoint_depth_grid.Read(image_name)) return 0;

        // Create estimate of vertical poles
        R2Grid density_mask(density_grid);
        density_mask.Threshold(5, 0, 1);
        R2Grid principle_axis_3_z_mask(principle_axis_3_z_grid);
        principle_axis_3_z_mask.Abs();
        principle_axis_3_z_mask.Threshold(0.1, 1, 0);
        R2Grid viewpoint_depth_mask(viewpoint_depth_grid);
        viewpoint_depth_mask.Blur(RN_X, 2);
        viewpoint_depth_mask.Subtract(viewpoint_depth_grid);
        viewpoint_depth_mask.Threshold(0.5, 0, 1);
        R2Grid pole_grid(lambda21_grid);
        pole_grid.Negate();
        pole_grid.Add(1.0);
        pole_grid.Threshold(0.8, 0, R2_GRID_KEEP_VALUE);
        pole_grid.Mask(density_mask);
        pole_grid.Mask(principle_axis_3_z_mask);
        pole_grid.Mask(viewpoint_depth_mask);
        // pole_grid.Subtract(0.8);
        // pole_grid.Multiply(5.0);
        // pole_grid.MaskNonMaxima(10);
        // pole_grid.Substitute(R2_GRID_UNKNOWN_VALUE, 0);

        pole_grid.Multiply(10.0);
        pole_grid.Add(lambda21_grid);

#if 0
        // Write grids for debugging
        char buffer[1024];
        sprintf(buffer, "tmp/%d_%d_%d_da_pole.grd", ir, is, ia);
        pole_grid.WriteFile(buffer);
        sprintf(buffer, "tmp/%d_%d_%d_da_viewpoint_depth_mask.grd", ir, is, ia);
        viewpoint_depth_mask.WriteFile(buffer);
#endif

      }
    }
  }

  // Return success
  return 1;
}

#endif



////////////////////////////////////////////////////////////////////////
// Curb feature creation functions
////////////////////////////////////////////////////////////////////////

static int
ReadCurbFeatures(GSVPoseOptimization *optimization, GSVScan *scan, const char *filename)
{
  // Parameters
  const RNAngle max_curvature = 1;

  // Open curb path file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open curb path file %s\n", filename);
    return 0;
  }

  // Read XYZ positions of curb path vertices
  RNScalar score[3];
  R3Point position[3];
  int scanline_index[3], point_index[3], count;
  while (fscanf(fp, "%lf%lf%lf%lf%d%d%d\n", &position[2][0], &position[2][1], &position[2][2], 
    &score[2], &scanline_index[2], &point_index[2], &count) == (unsigned int) 7) {
    // Create feature
    if (count >= 2) {
      R3Vector va = position[0] - position[1];
      R3Vector vb = position[2] - position[1];
      RNLength lena = va.Length();
      RNLength lenb = vb.Length();
      if (RNIsNotZero(lena) && RNIsNotZero(lenb)) {
        RNScalar angle = R3InteriorAngle(va, vb);
        RNScalar curvature = (RN_PI - angle) / (lena + lenb);
        if (fabs(curvature) <= max_curvature) {
          R3Vector direction = position[2] - position[0];
          RNLength length = direction.Length();
          if (RNIsNotZero(length)) {
            RNScalar scale = 0.5 * length;
            direction /= length;
            R3Vector normal = direction % R3posz_vector;
            normal.Normalize();
            if (scan->SegmentIndex() == 0) normal.Flip();
            GSVScanline *scanline = scan->Scanline(scanline_index[1]);
            GSVFeature *feature = new GSVFeature(GSV_SCAN_PLANE_FEATURE_TYPE, scanline, point_index[1],  
              position[1], direction, normal, scale, score[1]);
            optimization->InsertFeature(feature);
          }
        }
      }
    }

    // Remember stuff
    position[0] = position[1];
    scanline_index[0] = scanline_index[1];
    point_index[0] = point_index[1];
    score[0] = score[1];
    position[1] = position[2];
    scanline_index[1] = scanline_index[2];
    point_index[1] = point_index[2];
    score[1] = score[2];
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int GSVPoseOptimization::
CreateScanCurbFeatures(void)
{
  // Make temporary directory
  char tmp_directory[4096], mkdir_cmd[4096];
  sprintf(tmp_directory, "%s/scan_curb_features", scene->CacheDataDirectoryName());
  sprintf(mkdir_cmd, "mkdir -p %s", tmp_directory); 
  system(mkdir_cmd);

  // Run program to extract curbs from scans
  char gsv2map_command[4096];
  sprintf(gsv2map_command, "gsv2map %s %s", scene->Filename(), tmp_directory);
  system(gsv2map_command);
    
  // Create curb features for every scan in scene
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    char run_directory[4096];
    sprintf(run_directory, "%s/%s", tmp_directory, run->Name());
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (scan->NScanlines() == 0) continue;
        if (ia == 1) continue;
        char curb_name[1024];
        sprintf(curb_name, "%s/%02d_%02d_XYZ_CurbPath.txt", run_directory, is, ia);
        if (!ReadCurbFeatures(this, scan, curb_name)) return 0;
      }
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Ridge and valley feature creation functions
////////////////////////////////////////////////////////////////////////

struct SegmentBoundary {
  int segments[2];
  RNArray<const RNScalar *> points[2];
  RNScalar score;
};


int GSVPoseOptimization::
CreateScanRidgeAndValleyFeatures(void)
{
  // Parameters
  int min_segment_points = 16;
  int min_boundary_points = 8;
  RNAngle min_dihedral_angle = RN_PI / 4.0;
  RNAngle max_dihedral_angle = 3.0 * RN_PI / 4.0;
  RNLength min_length = 1;
  RNLength max_max_distance = 1;
  RNLength max_mean_distance = 1;
  RNScalar min_score = RN_EPSILON;

  // Check scene
  if (!scene) return 0;

  // Get image directory name
  char laser_image_directory[4096];
  sprintf(laser_image_directory, "%s/laser_images", scene->CacheDataDirectoryName());

  // Create ridge and valley features
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (scan->NScanlines() == 0) continue;
        if (ia == 1) continue;

        // Read position and info grids
        char image_name[4096];
        R2Grid position_x_grid;
        R2Grid position_y_grid;
        R2Grid position_z_grid;
        R2Grid scanline_grid;
        R2Grid point_index_grid;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionX.grd", laser_image_directory, run->Name(), is, ia);
        if (!position_x_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionY.grd", laser_image_directory, run->Name(), is, ia);
        if (!position_y_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionZ.grd", laser_image_directory, run->Name(), is, ia);
        if (!position_z_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_Scanline.grd", laser_image_directory, run->Name(), is, ia);
        if (!scanline_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PointIndex.grd", laser_image_directory, run->Name(), is, ia);
        if (!point_index_grid.Read(image_name)) return 0;

        // Read plane segment grids 
        R2Grid plane_segment_a_grid;
        R2Grid plane_segment_b_grid;
        R2Grid plane_segment_c_grid;
        R2Grid plane_segment_d_grid;
        R2Grid plane_segment_id_grid;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_SmallPlaneA.grd", laser_image_directory, run->Name(), is, ia);
        if (!plane_segment_a_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_SmallPlaneB.grd", laser_image_directory, run->Name(), is, ia);
        if (!plane_segment_b_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_SmallPlaneC.grd", laser_image_directory, run->Name(), is, ia);
        if (!plane_segment_c_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_SmallPlaneD.grd", laser_image_directory, run->Name(), is, ia);
        if (!plane_segment_d_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_SmallPlaneId.grd", laser_image_directory, run->Name(), is, ia);
        if (!plane_segment_id_grid.Read(image_name)) return 0;
        
        // Parse segments from plane segment grids
        RNArray<RNArray<const RNScalar *> *> segments;
        const RNScalar *grid_values = plane_segment_id_grid.GridValues();
        for (int i = 0; i < plane_segment_id_grid.NEntries(); i++) {
          int id_value = plane_segment_id_grid.GridValue(i);
          if (id_value == R2_GRID_UNKNOWN_VALUE) continue;
          int id = (int) (id_value + 0.5);
          while (id >= segments.NEntries()) segments.Insert(new RNArray<const RNScalar *>());
          RNArray<const RNScalar *> *segment = segments.Kth(id);
          segment->Insert(&grid_values[i]);
        }

        printf("HEREA %d\n", segments.NEntries());

        // Create segment boundaries
        RNArray<SegmentBoundary *> boundaries;
        for (int ix0 = 0; ix0 < plane_segment_id_grid.XResolution(); ix0++) {
          for (int iy0 = 0; iy0 < plane_segment_id_grid.YResolution(); iy0++) {
            // Get segment id0
            int index0;
            plane_segment_id_grid.IndicesToIndex(ix0, iy0, index0);
            int id0_value = plane_segment_id_grid.GridValue(index0);
            if (id0_value == R2_GRID_UNKNOWN_VALUE) continue;
            int id0 = (int) (id0_value + 0.5);
            if ((id0 < 0) || (id0 >= segments.NEntries())) continue;
            if ((min_segment_points > 0) && (segments[id0]->NEntries() < min_segment_points)) continue;

            // Check neighbors 
            for (int s = 0; s < 4; s++) {
              int ix1 = ix0;
              int iy1 = iy0;
              if (s == 0) ix1++;
              else if (s == 1) ix1--;
              else if (s == 2) iy1++;
              else if (s == 3) iy1--;

              // Check if neighbor is outside grid
              if ((ix1 < 0) || (ix1 >= plane_segment_id_grid.XResolution())) continue;
              if ((iy1 < 0) || (iy1 >= plane_segment_id_grid.YResolution())) continue;

              // Get neighbor segment id
              int index1;
              plane_segment_id_grid.IndicesToIndex(ix1, iy1, index1);
              int id1_value = plane_segment_id_grid.GridValue(index1);
              if (id1_value == R2_GRID_UNKNOWN_VALUE) continue;
              int id1 = (int) (id1_value + 0.5);
              if ((id1 < 0) || (id1 >= segments.NEntries())) continue;
              if ((min_segment_points > 0) && (segments[id1]->NEntries() < min_segment_points)) continue;

              // Check if not boundary
              if (id0 == id1) continue;

              // Find boundary
              SegmentBoundary *boundary = NULL;
              for (int k = 0; k < boundaries.NEntries(); k++) {
                SegmentBoundary *b = boundaries.Kth(k);
                if ((b->segments[0] == id0) || (b->segments[1] == id1)) { boundary = b; break; }
                if ((b->segments[1] == id0) || (b->segments[0] == id1)) { boundary = b; break; }
              }

              // Create boundary
              if (!boundary) {
                boundary = new SegmentBoundary();
                boundary->segments[0] = id0;
                boundary->segments[1] = id1;
                boundary->score = sqrt(segments[id0]->NEntries() * segments[id1]->NEntries());
                boundaries.Insert(boundary);
              }

              // Add points to boundary
              boundary->points[0].Insert(grid_values + index0);
              boundary->points[1].Insert(grid_values + index1);
            }
          }
        }

        printf("HEREB %d\n", boundaries.NEntries());

        // Process boundaries
        int feature_count = 0;
        for (int i = 0; i < boundaries.NEntries(); i++) {
          SegmentBoundary *boundary = boundaries.Kth(i);

          // Check boundary score
          if ((min_score > 0) && (boundary->score < min_score)) continue;

          // Check boundary size
          if ((min_boundary_points > 0) && (boundary->points[0].NEntries() < min_boundary_points)) continue;
          if ((min_boundary_points > 0) && (boundary->points[1].NEntries() < min_boundary_points)) continue;
          
          // Get planes
          R3Plane planes[2];
          for (int j = 0; j < 2; j++) {
            const RNScalar *grid_value = boundary->points[j].Kth(0);
            int grid_index = grid_value - grid_values;
            RNScalar a = plane_segment_a_grid.GridValue(grid_index);
            if (a == R2_GRID_UNKNOWN_VALUE) continue;
            RNScalar b = plane_segment_b_grid.GridValue(grid_index);
            if (b == R2_GRID_UNKNOWN_VALUE) continue;
            RNScalar c = plane_segment_c_grid.GridValue(grid_index);
            if (c == R2_GRID_UNKNOWN_VALUE) continue;
            RNScalar d = plane_segment_d_grid.GridValue(grid_index);
            if (d == R2_GRID_UNKNOWN_VALUE) continue;
            R3Vector normal(a, b, c);
            RNLength normal_length = normal.Length();
            if (RNIsZero(normal_length)) continue;
            normal /= normal_length;
            planes[j] = R3Plane(normal, d);
          }

          // Check dihedral angle
          RNScalar dot = fabs(planes[0].Normal().Dot(planes[1].Normal()));
          RNAngle dihedral_angle = (dot < 1.0) ? acos(dot) : 0.0;
          if ((min_dihedral_angle > 0) && (dihedral_angle < min_dihedral_angle)) continue;
          if ((max_dihedral_angle > 0) && (dihedral_angle > max_dihedral_angle)) continue;

          // Make array of boundary points
          int npoints = 0;
          int max_points = boundary->points[0].NEntries() + boundary->points[1].NEntries();
          R3Point *points = new R3Point [ max_points ];
          for (int j = 0; j < 2; j++) {
            for (int k = 0; k < boundary->points[j].NEntries(); k++) {
              // Get grid index
              const RNScalar *grid_value = boundary->points[j].Kth(k);
              int grid_index = grid_value - grid_values;

              // Get point position
              RNScalar x = position_x_grid.GridValue(grid_index);
              if (x == R2_GRID_UNKNOWN_VALUE) continue;
              RNScalar y = position_y_grid.GridValue(grid_index);
              if (y == R2_GRID_UNKNOWN_VALUE) continue;
              RNScalar z = position_z_grid.GridValue(grid_index);
              if (z == R2_GRID_UNKNOWN_VALUE) continue;
              points[npoints++].Reset(x, y, z);
            }
          }

          // Get line (ray) throughh boundary points
          RNScalar variances[3] = { 0 };
          R3Point centroid = R3Centroid(npoints, points);
          R3Triad triad = R3PrincipleAxes(centroid, npoints, points, NULL, variances);
          R3Ray ray(centroid, triad.Axis(0));

          // Project all points onto ray
          int point_count = 0;
          RNLength max_distance = 0;
          RNLength total_distance = 0;
          RNScalar endpoint_t[2] = { FLT_MAX, -FLT_MAX };
          int endpoint_grid_index[2] = { -1, -1 };
          for (int j = 0; j < 2; j++) {
            for (int k = 0; k < boundary->points[j].NEntries(); k++) {
              // Get grid index
              const RNScalar *grid_value = boundary->points[j].Kth(k);
              int grid_index = grid_value - grid_values;

              // Get point position
              RNScalar x = position_x_grid.GridValue(grid_index);
              if (x == R2_GRID_UNKNOWN_VALUE) continue;
              RNScalar y = position_y_grid.GridValue(grid_index);
              if (y == R2_GRID_UNKNOWN_VALUE) continue;
              RNScalar z = position_z_grid.GridValue(grid_index);
              if (z == R2_GRID_UNKNOWN_VALUE) continue;
              R3Point position(x, y, z);

              // Project point onto ray
              RNScalar t = ray.T(position);
              R3Point projection = ray.Point(t);

              // Check distance
              RNLength distance = R3Distance(position, projection);
              if (distance > max_distance) max_distance = distance;
              total_distance += distance;

              // Update endpoints
              if (t < endpoint_t[0]) {
                endpoint_grid_index[0] = grid_index;
                endpoint_t[0] = t;
              }
              if (t > endpoint_t[1]) {
                endpoint_grid_index[1] = grid_index;
                endpoint_t[1] = t;
              }

              // Update count
              point_count++;
            }
          }

          // printf("  HEREB1 %d : %g %g : %g\n", point_count, 
          //   endpoint_t[0], endpoint_t[1], max_distance);

          // Check boundary
          if ((min_boundary_points > 0) && (point_count < min_boundary_points)) continue;
          if ((max_max_distance > 0) && (max_distance > max_max_distance)) continue;
          RNScalar mean_distance = total_distance / point_count;
          if ((max_mean_distance > 0) && (mean_distance > max_mean_distance)) continue;

          // Get stuff for feature
          R3Vector direction = ray.Vector();
          R3Vector normal = R3zero_vector;
          RNScalar score = boundary->score;
          RNScalar length = endpoint_t[1] - endpoint_t[0];
          if ((min_length > 0) && (length < min_length)) continue;

          // printf("  HEREB2 %d %d %g %g\n", 
          //   endpoint_grid_index[0], endpoint_grid_index[1], length, score);

          // Create first line feature
          if (endpoint_grid_index[0] >= 0) {
            // RNScalar x = position_x_grid.GridValue(endpoint_grid_index[0]);
            // if (x == R2_GRID_UNKNOWN_VALUE) continue;
            // RNScalar y = position_y_grid.GridValue(endpoint_grid_index[0]);
            // if (y == R2_GRID_UNKNOWN_VALUE) continue;
            // RNScalar z = position_z_grid.GridValue(endpoint_grid_index[0]);
            // if (z == R2_GRID_UNKNOWN_VALUE) continue;
            // R3Point position(x, y, z);
            R3Point position = ray.Point(endpoint_t[0]);
            RNScalar scanline_value = scanline_grid.GridValue(endpoint_grid_index[0]);
            if (scanline_value == R2_GRID_UNKNOWN_VALUE) continue;
            RNScalar point_index_value = point_index_grid.GridValue(endpoint_grid_index[0]);
            if (point_index_value == R2_GRID_UNKNOWN_VALUE) continue;
            int scanline_index = (int) (scanline_value + 0.5);
            if ((scanline_index < 0) || (scanline_index >= scan->NScanlines())) continue;
            GSVScanline *scanline = scan->Scanline(scanline_index);
            int point_index = (int) (point_index_value + 0.5);
            GSVFeature *feature = new GSVFeature(GSV_SCAN_LINE_FEATURE_TYPE,
              scanline, point_index, position, direction, normal, length, score);
            InsertFeature(feature);
            feature_count++;
          }

          // Create second line feature
          if (endpoint_grid_index[1] >= 0) {
            // RNScalar x = position_x_grid.GridValue(endpoint_grid_index[1]);
            // if (x == R2_GRID_UNKNOWN_VALUE) continue;
            // RNScalar y = position_y_grid.GridValue(endpoint_grid_index[1]);
            // if (y == R2_GRID_UNKNOWN_VALUE) continue;
            // RNScalar z = position_z_grid.GridValue(endpoint_grid_index[1]);
            // if (z == R2_GRID_UNKNOWN_VALUE) continue;
            // R3Point position(x, y, z);
            R3Point position = ray.Point(endpoint_t[1]);
            RNScalar scanline_value = scanline_grid.GridValue(endpoint_grid_index[1]);
            if (scanline_value == R2_GRID_UNKNOWN_VALUE) continue;
            RNScalar point_index_value = point_index_grid.GridValue(endpoint_grid_index[1]);
            if (point_index_value == R2_GRID_UNKNOWN_VALUE) continue;
            int scanline_index = (int) (scanline_value + 0.5);
            if ((scanline_index < 0) || (scanline_index >= scan->NScanlines())) continue;
            GSVScanline *scanline = scan->Scanline(scanline_index);
            int point_index = (int) (point_index_value + 0.5);
            GSVFeature *feature = new GSVFeature(GSV_SCAN_LINE_FEATURE_TYPE,
              scanline, point_index, position, -direction, normal, length, score);
            InsertFeature(feature);
            feature_count++;
          }
        }

        // Delete segments
        for (int i = 0; i < segments.NEntries(); i++) {
          RNArray<const RNScalar *> *segment = segments.Kth(i);
          delete segment;
        }

        // Delete boundaries
        for (int i = 0; i < boundaries.NEntries(); i++) {
          SegmentBoundary *boundary = boundaries.Kth(i);
          delete boundary;
        }

        printf("  HEREC %d\n", feature_count);
      }
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Edge feature creation functions
////////////////////////////////////////////////////////////////////////

struct FeatureCluster {
  FeatureCluster(void) : points(), span(0,0,0,0,0,0), id(-1)
    { variances[0] = variances[1] = variances[2] = 0.0; };
  RNArray<const RNScalar *> points;
  RNScalar variances[3];
  R3Span span;
  int id;
};



static int 
CreateEdgeFeatures(GSVPoseOptimization *optimization, GSVScan *scan)
{
  // Parameters
  int min_points = 4;
  int boundary_type_mask = 0;
  RNScalar min_lambda1 = 1E-2;
  RNScalar max_lambda2 = 1E-2;
  RNScalar max_lambda21 = 0.25;
  RNLength min_length = 0.5;
  RNScalar min_density = 1;
  RNScalar axis_alignment = cos(RN_PI/16.0);
  const RNScalar min_feature_spacing = 5;

  // Get useful variables
  int scan_index = scan->SegmentIndex();
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  int segment_index = segment->RunIndex();
  GSVRun *run = segment->Run();
  if (!run) return 0;
  const char *run_name = run->Name();
  GSVScene *scene = run->Scene();
  if (!scene) return 0;
  char laser_image_directory[4096];
  sprintf(laser_image_directory, "%s/laser_images", scene->CacheDataDirectoryName());

  // Read grids
  char image_name[4096];
  R2Grid position_x_grid;
  R2Grid position_y_grid;
  R2Grid position_z_grid;
  R2Grid scanline_index_grid;
  R2Grid point_index_grid;
  R2Grid ground_z_grid;
  R2Grid normal_z_grid;
  R2Grid boundary_type_grid;
  R2Grid boundary_id_grid;
  R2Grid boundary_size_grid;
  sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionX.grd", laser_image_directory, run_name, segment_index, scan_index);
  if (!position_x_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionY.grd", laser_image_directory, run_name, segment_index, scan_index);
  if (!position_y_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionZ.grd", laser_image_directory, run_name, segment_index, scan_index);
  if (!position_z_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_DA_Scanline.grd", laser_image_directory, run_name, segment_index, scan_index);
  if (!scanline_index_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_DA_PointIndex.grd", laser_image_directory, run_name, segment_index, scan_index);
  if (!point_index_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_DA_GroundZ.grd", laser_image_directory, run_name, segment_index, scan_index);
  if (!ground_z_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_DA_BoundaryType.grd", laser_image_directory, run_name, segment_index, scan_index);
  if (!boundary_type_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_DA_BoundaryLineId.grd", laser_image_directory, run_name, segment_index, scan_index);
  if (!boundary_id_grid.Read(image_name)) return 0;
  sprintf(image_name, "%s/%s/%02d_%02d_DA_BoundaryLineSize.grd", laser_image_directory, run_name, segment_index, scan_index);
  if (!boundary_size_grid.Read(image_name)) return 0;

  // Allocate temporary memory
  R3Point *positions = new R3Point [ boundary_id_grid.NEntries() ];

  // Parse line clusters
  RNArray<FeatureCluster *> clusters;
  const RNScalar *grid_values = boundary_id_grid.GridValues();
  for (int i = 0; i < boundary_id_grid.NEntries(); i++) {
    RNScalar boundary_size_value = boundary_size_grid.GridValue(i);
    if (boundary_size_value == R2_GRID_UNKNOWN_VALUE) continue;
    if ((min_points > 0) && (boundary_size_value < min_points - 0.5)) continue;
    RNScalar boundary_id_value = boundary_id_grid.GridValue(i);
    if (boundary_id_value == R2_GRID_UNKNOWN_VALUE) continue;
    RNScalar boundary_type_value = boundary_type_grid.GridValue(i);
    if (boundary_type_value == R2_GRID_UNKNOWN_VALUE) continue;
    int boundary_type = (int) (boundary_type_value + 0.5);
    if ((boundary_type_mask != 0) && (!(boundary_type & boundary_type_mask))) continue;
    RNScalar x = position_x_grid.GridValue(i);
    if (x == R2_GRID_UNKNOWN_VALUE) continue;
    RNScalar y = position_y_grid.GridValue(i);
    if (y == R2_GRID_UNKNOWN_VALUE) continue;
    RNScalar z = position_z_grid.GridValue(i);
    if (z == R2_GRID_UNKNOWN_VALUE) continue;
    RNScalar ground_z = ground_z_grid.GridValue(i);
    if (ground_z == R2_GRID_UNKNOWN_VALUE) continue;
    RNScalar scanline_index_value = scanline_index_grid.GridValue(i);
    if (scanline_index_value == R2_GRID_UNKNOWN_VALUE) continue;
    RNScalar point_index_value = point_index_grid.GridValue(i);
    if (point_index_value == R2_GRID_UNKNOWN_VALUE) continue;
    int boundary_id = (int) (boundary_id_value + 0.5);
    while (boundary_id >= clusters.NEntries()) 
      clusters.Insert(new FeatureCluster());
    FeatureCluster *cluster = clusters.Kth(boundary_id);
    cluster->points.Insert(&grid_values[i]);
    cluster->id = boundary_id;
  }

  // Create line features
  for (int i = 0; i < clusters.NEntries(); i++) {
    FeatureCluster *cluster = clusters.Kth(i);

    // Check number of points
    if ((min_points > 0) && (cluster->points.NEntries() < min_points)) continue;

    // Create array of point positions
    int npositions = 0;
    for (int j = 0; j < cluster->points.NEntries(); j++) {
      const RNScalar *point = cluster->points.Kth(j);
      int grid_index = point - grid_values;
      RNScalar x = position_x_grid.GridValue(grid_index);
      assert(x != R2_GRID_UNKNOWN_VALUE);
      RNScalar y = position_y_grid.GridValue(grid_index);
      assert(y != R2_GRID_UNKNOWN_VALUE);
      RNScalar z = position_z_grid.GridValue(grid_index);
      assert(z != R2_GRID_UNKNOWN_VALUE);
      positions[npositions++].Reset(x, y, z);
    }

    // Check number of point positions
    if ((min_points > 0) && (npositions < min_points)) continue;
    if (npositions < 3) continue;

    // Compute centroid and principal axes 
    RNScalar variances[3] = { 0, 0, 0 };
    R3Point centroid = R3Centroid(npositions, positions);
    R3Triad triad = R3PrincipleAxes(centroid, npositions, positions, NULL, variances);
    R3Vector principal_axis = triad.Axis(0);
    if (principal_axis.Z() < 0) principal_axis.Flip();

    // Check variances
    if ((min_lambda1 > 0) && (variances[0] < min_lambda1)) continue;
    if ((max_lambda2 > 0) && (variances[1] > max_lambda2)) continue;
    if ((max_lambda21 > 0) && (variances[0] > 0) && (variances[1] / variances[0] > max_lambda21)) continue; 

    // Check if axis aligned
    if (axis_alignment > 0) {
      // Should be sqrt(1.0 - axis_alignment * axis_alignment)
      if ((principal_axis.Z() < axis_alignment) && (principal_axis.Z() > 1.0 - axis_alignment)) continue;
    }

    // Determine range of positions on ray
    R3Ray ray(centroid, principal_axis);
    RNInterval range( FLT_MAX, -FLT_MAX );
    for (int j = 0; j < npositions; j++) {
      RNScalar t = ray.T(positions[j]);
      range.Union(t);
    }

    // Determine span
    R3Point endpoint0 = ray.Point(range.Min());
    R3Point endpoint1 = ray.Point(range.Max());
    R3Span span(endpoint0, endpoint1);

    // Check length
    if ((min_length > 0) && (span.Length() < min_length)) continue;
    if (span.Length() == 0) continue;

    // Check density
    RNScalar density = npositions / span.Length();
    if ((min_density > 0) && (density < min_density)) continue;

    // Determine scanline and point index
    int point_index = -1;
    int scanline_index = -1;
    RNScalar min_d = FLT_MAX;
    RNScalar ground_z = 0;
    int boundary_types[6] = { 0, 0, 0, 0, 0, 0 };
    for (int j = 0; j < cluster->points.NEntries(); j++) {
      const RNScalar *point = cluster->points.Kth(j);
      int grid_index = point - grid_values;
      RNScalar x = position_x_grid.GridValue(grid_index);
      assert(x != R2_GRID_UNKNOWN_VALUE);
      RNScalar y = position_y_grid.GridValue(i);
      assert(y != R2_GRID_UNKNOWN_VALUE);
      RNScalar z = position_z_grid.GridValue(i);
      assert(z != R2_GRID_UNKNOWN_VALUE);
      RNScalar boundary_type_value = boundary_type_grid.GridValue(grid_index);
      assert(boundary_type_value != R2_GRID_UNKNOWN_VALUE);
      int boundary_type = (int) (boundary_type_value + 0.5);
      if (boundary_type & GSV_LEFT_SILHOUETTE_BOUNDARY) boundary_types[0]++;
      else if (boundary_type & GSV_RIGHT_SILHOUETTE_BOUNDARY) boundary_types[1]++;
      else if (boundary_type & GSV_DOWN_SILHOUETTE_BOUNDARY) boundary_types[2]++;
      else if (boundary_type & GSV_UP_SILHOUETTE_BOUNDARY) boundary_types[3]++;
      else if (boundary_type & GSV_HORIZONTAL_RIDGE_BOUNDARY) boundary_types[4]++;
      else if (boundary_type & GSV_VERTICAL_RIDGE_BOUNDARY) boundary_types[5]++;
      R3Point p(x, y, z);
      RNScalar d = R3SquaredDistance(p, centroid);
      if (d < min_d) {
        RNScalar scanline_index_value = scanline_index_grid.GridValue(grid_index);
        assert(scanline_index_value != R2_GRID_UNKNOWN_VALUE);
        RNScalar point_index_value = point_index_grid.GridValue(grid_index);
        assert(point_index_value != R2_GRID_UNKNOWN_VALUE);
        ground_z = ground_z_grid.GridValue(grid_index);
        assert(ground_z != R2_GRID_UNKNOWN_VALUE);
        scanline_index = (int) (scanline_index_value + 0.5);
        point_index = (int) (point_index_value + 0.5);
        min_d = d;
      }
    }
    
    // Check scanline and point index
    if (scanline_index == -1) continue;
    if (point_index == -1) continue;

    // Compute boundary type
    int boundary_type = 0;
    int boundary_count = 0;
    for (int i = 0; i < 6; i++) {
      if (boundary_types[i] > boundary_count) { 
        boundary_count = boundary_types[i];
        boundary_type = i;
      }
    }

    // Check boundary type
    if (boundary_count < 0.5 * cluster->points.NEntries()) continue;

    // Temporary
    R3Vector normal = span.Vector() % R3posx_vector;
    normal.Normalize();

    // Create unique group number
    int group = (int) (RNRandomScalar() * INT_MAX);

    // Create line features
    int nfeatures = (int) (span.Length() / min_feature_spacing) + 1;
    RNScalar step = span.Length() / nfeatures;
    for (int i = 0; i < nfeatures; i++) {
      // Create line feature
      R3Point position = span.Point((i + 0.5) * span.Length() / nfeatures);
      GSVFeature *feature = new GSVFeature(GSV_SCAN_LINE_FEATURE_TYPE, 
        scan->Scanline(scanline_index), point_index, 
        position, span.Vector(), normal, step, 1, group);

      // Create line shape descriptor
      GSVDescriptor& descriptor = feature->descriptor;
      descriptor.descriptor_type = GSV_LINE_DESCRIPTOR_TYPE;
      descriptor.nvalues = 2;
      descriptor.values = new RNScalar [ descriptor.nvalues ];
      descriptor.values[0] = 100 * boundary_type;
      descriptor.values[1] = position.Z() - ground_z;

      // Insert feature
      optimization->InsertFeature(feature);
    }
  }

  // Delete temporary memory
  delete [] positions;

  // Return success
  return 1;
}



int GSVPoseOptimization::
CreateScanEdgeFeatures(void)
{
  // Check scene
  if (!scene) return 0;

  // Create edge features
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (scan->NScanlines() == 0) continue;
        if (ia == 1) continue;
        if (!CreateEdgeFeatures(this, scan)) return 0;
      }
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Plane feature creation functions
////////////////////////////////////////////////////////////////////////

static RNScalar *
CreatePlaneSegmentAreas(int nsegments, const R2Grid& id_grid, 
  const R2Grid& position_x_grid, const R2Grid& position_y_grid, const R2Grid& position_z_grid)
{
  // Allocate areas
  RNScalar *areas = new RNScalar [ nsegments ];
  if (!areas) {
    fprintf(stderr, "Unable to allocate plane segment areas\n");
    return NULL;
  }

  // Initialize areas
  for (int i = 0; i < nsegments; i++) areas[i] = 0;

  // Compute areas
  for (int iy = 0; iy < id_grid.YResolution()-1; iy++) {
    for (int ix = 0; ix < id_grid.XResolution()-1; ix++) {
      // Get segment id
      RNScalar id0_value = id_grid.GridValue(ix, iy);
      if (id0_value == R2_GRID_UNKNOWN_VALUE) continue;
      int id0 = (int) (id0_value + 0.5);
      assert((id0 >= 0) && (id0 < nsegments));

      // Determine position
      RNScalar x0 = position_x_grid.GridValue(ix, iy);
      RNScalar y0 = position_y_grid.GridValue(ix, iy);
      RNScalar z0 = position_z_grid.GridValue(ix, iy);
      if (x0 == R2_GRID_UNKNOWN_VALUE) continue;
      if (y0 == R2_GRID_UNKNOWN_VALUE) continue;
      if (z0 == R2_GRID_UNKNOWN_VALUE) continue;
      R3Point p0(x0, y0, z0);
      
      // Determine dx
      RNScalar dx = 0;
      if (1) {
        // Check segment id
        RNScalar id1_value = id_grid.GridValue(ix+1, iy);
        if (id1_value == R2_GRID_UNKNOWN_VALUE) continue;
        int id1 = (int) (id1_value + 0.5);
        assert((id1 >= 0) && (id1 < nsegments));
        if (id1 != id0) continue;

        // Get position
        RNScalar x1 = position_x_grid.GridValue(ix+1, iy);
        RNScalar y1 = position_y_grid.GridValue(ix+1, iy);
        RNScalar z1 = position_z_grid.GridValue(ix+1, iy);
        if (x1 == R2_GRID_UNKNOWN_VALUE) continue;
        if (y1 == R2_GRID_UNKNOWN_VALUE) continue;
        if (z1 == R2_GRID_UNKNOWN_VALUE) continue;
        R3Point p1(x1, y1, z1);

        // Compute distance
        dx = R3Distance(p0, p1);
      }

      // Determine dy
      RNScalar dy = 0;
      if (1) {
        // Check segment id
        RNScalar id1_value = id_grid.GridValue(ix, iy+1);
        if (id1_value == R2_GRID_UNKNOWN_VALUE) continue;
        int id1 = (int) (id1_value + 0.5);
        assert((id1 >= 0) && (id1 < nsegments));
        if (id1 != id0) continue;

        // Get position
        RNScalar x1 = position_x_grid.GridValue(ix, iy+1);
        RNScalar y1 = position_y_grid.GridValue(ix, iy+1);
        RNScalar z1 = position_z_grid.GridValue(ix, iy+1);
        if (x1 == R2_GRID_UNKNOWN_VALUE) continue;
        if (y1 == R2_GRID_UNKNOWN_VALUE) continue;
        if (z1 == R2_GRID_UNKNOWN_VALUE) continue;
        R3Point p1(x1, y1, z1);

        // Compute distance
        dy = R3Distance(p0, p1);
      }

      // Check diagonal point
      if (1) {
        // Check segment id
        RNScalar id1_value = id_grid.GridValue(ix+1, iy+1);
        if (id1_value == R2_GRID_UNKNOWN_VALUE) continue;
        int id1 = (int) (id1_value + 0.5);
        assert((id1 >= 0) && (id1 < nsegments));
        if (id1 != id0) continue;
      }

      // Sum area of point sample
      areas[id0] += dx * dy;
    }
  }

  // Return areas
  return areas;
}



static RNScalar *
CreatePlaneSegmentConnectivities(int nsegments, const R2Grid& id_grid)
{
  // Allocate connectivies
  RNScalar *connectivities = new RNScalar [ nsegments ];
  if (!connectivities) {
    fprintf(stderr, "Unable to allocate plane segment connectivities\n");
    return NULL;
  }

  // Allocate temporary data
  int *counts = new int [ nsegments ];

  // Initialize counts
  for (int i = 0; i < nsegments; i++) {
    connectivities[i] = 0.0;
    counts[i] = 0;
  }

  // Compute counts of neighbors in same segment
  for (int iy = 1; iy < id_grid.YResolution()-1; iy++) {
    for (int ix = 1; ix < id_grid.XResolution()-1; ix++) {
      // Get segment id
      RNScalar id0_value = id_grid.GridValue(ix, iy);
      if (id0_value == R2_GRID_UNKNOWN_VALUE) continue;
      int id0 = (int) (id0_value + 0.5);
      assert((id0 >= 0) && (id0 < nsegments));
      for (int s = -1; s <= 1; s++) {
        for (int t = -1; t <= 1; t++) {
          if ((s == 0) && (t == 0)) continue;
          RNScalar id1_value = id_grid.GridValue(ix+s, iy+t);
          if (id1_value == R2_GRID_UNKNOWN_VALUE) continue;
          int id1 = (int) (id1_value + 0.5);
          assert((id1 >= 0) && (id1 < nsegments));
          if (id0 == id1) connectivities[id0] += 1.0;
        }
      }
      counts[id0] += 8;
    }
  }

  // Convert counts to fractions
  for (int i = 0; i < nsegments; i++) {
    if (counts[i] == 0) continue;
    connectivities[i] /= counts[i];
  }

  // Delete temporary data
  delete [] counts;
  
  // Return connectivities
  return connectivities;
}



int GSVPoseOptimization::
CreateScanPlaneFeatures(void)
{
  // Parameters
  const RNScalar min_feature_spacing = 5;
  const RNScalar min_segment_area = 10;
  // const RNScalar min_segment_connectivity = 0.75;
  const RNScalar min_segment_connectivity = 0.5;
  const RNScalar target_segment_area = 200;
  const int min_segment_nentries = 100;
  RNLength min_feature_spacing_squared = min_feature_spacing * min_feature_spacing;

  // Check scene
  if (!scene) return 0;

  // Get image directory name
  char laser_image_directory[4096];
  sprintf(laser_image_directory, "%s/laser_images", scene->CacheDataDirectoryName());

  // Create plane features
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (scan->NScanlines() == 0) continue;
        if (ia == 1) continue;

        // Read DA grids
        char image_name[4096];
        R2Grid da_scanline_index_grid;
        R2Grid da_point_index_grid;
        R2Grid da_position_x_grid;
        R2Grid da_position_y_grid;
        R2Grid da_position_z_grid;
        R2Grid da_plane_segment_A_grid;
        R2Grid da_plane_segment_B_grid;
        R2Grid da_plane_segment_C_grid;
        R2Grid da_plane_segment_D_grid;
        R2Grid da_plane_segment_id_grid;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_Scanline.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_scanline_index_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PointIndex.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_point_index_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionX.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_position_x_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionY.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_position_y_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionZ.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_position_z_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_SmallPlaneA.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_plane_segment_A_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_SmallPlaneB.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_plane_segment_B_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_SmallPlaneC.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_plane_segment_C_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_SmallPlaneD.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_plane_segment_D_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_SmallPlaneId.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_plane_segment_id_grid.Read(image_name)) return 0;

        // Check the grid resolution
        if (da_plane_segment_id_grid.YResolution() < 25) continue;

        // Mask scan points on GSV car
        for (int iy = 0; iy < 25; iy++) {
          for (int ix = 0; ix < da_plane_segment_id_grid.XResolution(); ix++) {
            da_scanline_index_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_point_index_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_position_x_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_position_y_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_position_z_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_plane_segment_A_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_plane_segment_B_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_plane_segment_C_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_plane_segment_D_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_plane_segment_id_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
          }
        }

        // Allocate temporary segment data 
        int nsegments = (int) (da_plane_segment_id_grid.Maximum() + 0.5) + 1;
        RNScalar *segment_connectivities = CreatePlaneSegmentConnectivities(nsegments, da_plane_segment_id_grid); 
        RNScalar *segment_areas = CreatePlaneSegmentAreas(nsegments, da_plane_segment_id_grid, 
          da_position_x_grid, da_position_y_grid, da_position_z_grid);
        RNArray<const RNScalar *> *segment_entries = new RNArray<const RNScalar *> [ nsegments ];
        for (int i = 0; i < da_plane_segment_id_grid.NEntries(); i++) {
          RNScalar segment_id_value = da_plane_segment_id_grid.GridValue(i);
          if (segment_id_value == R2_GRID_UNKNOWN_VALUE) continue;
          int segment_id = (int) (segment_id_value + 0.5);
          assert((segment_id >= 0) && (segment_id < nsegments));
          const RNScalar *segment_entry = &da_plane_segment_id_grid(i);
          segment_entries[segment_id].Insert(segment_entry);
        }

        // Randomize the order of the segment entries
        for (int segment_id = 0; segment_id < nsegments; segment_id++) {
          for (int k1 = 0; k1 < segment_entries[segment_id].NEntries(); k1++) {
            int k2 = (int) (RNRandomScalar() * segment_entries[segment_id].NEntries());
            segment_entries[segment_id].Swap(k1, k2);
          }
        }

        // Create features at points sampled from within plane segment
        for (int segment_id = 0; segment_id < nsegments; segment_id++) {
          // Check segment area
          RNScalar segment_area = segment_areas[segment_id];
          if (segment_area < min_segment_area) continue;

          // Check segment connectivity
          RNScalar segment_connectivity = segment_connectivities[segment_id];
          if (segment_connectivity < min_segment_connectivity) continue;
          
          // Check segment entries
          int segment_nentries = segment_entries[segment_id].NEntries();
          if (segment_nentries < min_segment_nentries) continue;

          // Create unique group number
          int group = (int) (RNRandomScalar() * INT_MAX);

          // printf("%6d %6d %9.3g %9.3g\n", segment_id, segment_entries[segment_id].NEntries(), segment_area, segment_connectivity);

          // Create features on segment with minimum spacing
          RNArray<GSVFeature *> created_features;
          for (int i = 0; i < segment_entries[segment_id].NEntries(); i++) {
            int grid_ix, grid_iy;
            const RNScalar *segment_entry = segment_entries[segment_id].Kth(i);
            int grid_index = segment_entry - da_plane_segment_id_grid.GridValues();
            da_plane_segment_id_grid.IndexToIndices(grid_index, grid_ix, grid_iy);

            // Determine the position
            RNScalar x = da_position_x_grid.GridValue(grid_index);
            RNScalar y = da_position_y_grid.GridValue(grid_index);
            RNScalar z = da_position_z_grid.GridValue(grid_index);
            if (x == R2_GRID_UNKNOWN_VALUE) continue;
            if (y == R2_GRID_UNKNOWN_VALUE) continue;
            if (z == R2_GRID_UNKNOWN_VALUE) continue;
            R3Point position(x, y, z);

            // Check if there is a previously created feature within minimum spacing
            RNBoolean well_spaced = TRUE;
            for (int j = 0; j < created_features.NEntries(); j++) {
              GSVFeature *created_feature = created_features.Kth(j);
              RNScalar distance_squared = R3SquaredDistance(position, created_feature->scan_position);
              if (distance_squared < min_feature_spacing_squared) { well_spaced = FALSE; break; }
            }

             // Create feature if not too close to previously created one
            if (well_spaced) {
              // Get feature information
              RNScalar scanline_index_value = da_scanline_index_grid.GridValue(grid_index);
              RNScalar point_index_value = da_point_index_grid.GridValue(grid_index);
              RNScalar a = da_plane_segment_A_grid.GridValue(grid_index);
              RNScalar b = da_plane_segment_B_grid.GridValue(grid_index);
              RNScalar c = da_plane_segment_C_grid.GridValue(grid_index);
              if (scanline_index_value == R2_GRID_UNKNOWN_VALUE) continue;
              if (point_index_value == R2_GRID_UNKNOWN_VALUE) continue;
              if (a == R2_GRID_UNKNOWN_VALUE) continue;
              if (b == R2_GRID_UNKNOWN_VALUE) continue;
              if (c == R2_GRID_UNKNOWN_VALUE) continue;
              int scanline_index = (int) (scanline_index_value + 0.5);
              int point_index = (int) (point_index_value + 0.5);
              GSVScanline *scanline = scan->Scanline(scanline_index);
              R3Vector normal(a, b, c);
              normal.Normalize();
              RNScalar score = (segment_area < target_segment_area) ? segment_area / target_segment_area : 1.0;
              GSVFeature *feature = new GSVFeature(GSV_SCAN_PLANE_FEATURE_TYPE, scanline, point_index,  
                position, R3zero_vector, normal, min_feature_spacing, score, group);
              created_features.Insert(feature);
              InsertFeature(feature);
            }
          }
        }

        // Delete temporary segment data 
        delete [] segment_connectivities;
        delete [] segment_areas;
        delete [] segment_entries;
      }
    }
  }

  // Return success
  return 1;
}




////////////////////////////////////////////////////////////////////////
// Scan contour feature creation functions
////////////////////////////////////////////////////////////////////////

#if 0

int GSVPoseOptimization::
CreateScanContourFeatures(void)
{
  // Parameters
  const RNScalar min_feature_spacing = 5;
  const RNScalar min_segment_area = 10;
  const RNScalar min_segment_connectivity = 0.75;
  const RNScalar target_segment_area = 200;
  const int min_segment_nentries = 100;
  RNLength min_feature_spacing_squared = min_feature_spacing * min_feature_spacing;

  // Check scene
  if (!scene) return 0;

  // Get image directory name
  char laser_image_directory[4096];
  sprintf(laser_image_directory, "%s/laser_images", scene->CacheDataDirectoryName());

  // Create contour features
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (scan->NScanlines() == 0) continue;
        if (ia != 1) continue;

        // Read grids
        char image_name[4096];
        R2Grid scanline_index_grid;
        R2Grid point_index_grid;
        R2Grid viewpoint_depth_grid;
        R2Grid position_x_grid;
        R2Grid position_y_grid;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_Scanline.grd", laser_image_directory, run->Name(), is, ia);
        if (!scanline_index_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PointIndex.grd", laser_image_directory, run->Name(), is, ia);
        if (!point_index_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_ViewpointDepth.grd", laser_image_directory, run->Name(), is, ia);
        if (!viewpoint_depth_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionX.grd", laser_image_directory, run->Name(), is, ia);
        if (!position_x_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionY.grd", laser_image_directory, run->Name(), is, ia);
        if (!position_y_grid.Read(image_name)) return 0;

        // Mask points that are not local minima of viewpoint depth within same vertical scanline
        for (int ix = 0; ix < viewpoint_depth_grid.XResolution(); ix++) {
          viewpoint_depth_grid.SetGridValue(ix, 0, R2_GRID_UNKNOWN_VALUE);
          viewpoint_depth_grid.SetGridValue(ix, viewpoint_depth_grid.YResolution()-1, R2_GRID_UNKNOWN_VALUE);
          for (int iy = 1; iy < viewpoint_depth_grid.YResolution()-1; iy++) {
            RNScalar depth = viewpoint_depth_grid.GridValue(ix, iy);
            if (depth == R2_GRID_UNKNOWN_VALUE) continue;
            RNScalar depth0 = viewpoint_depth_grid.GridValue(ix, iy-1);
            if ((depth0 != R2_GRID_UNKNOWN_VALUE) && (depth0 < depth)) {
              viewpoint_depth_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
              continue;
            }
            RNScalar depth1 = viewpoint_depth_grid.GridValue(ix, iy+1);
            if ((depth1 != R2_GRID_UNKNOWN_VALUE) && (depth1 < depth)) {
              viewpoint_depth_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
              continue;
            }
          }
        }
            
        // Find mutually closest points
        for (int ix = 0; ix < scanline_index_grid.XResolution()-1; ix++) {
          for (int iy0 = 1; iy0 < scanline_index_grid.YResolution()-1; iy0++) {
            RNScalar depth0 = viewpoint_depth_grid.GridValue(ix, iy0);
            if (depth0 == R2_GRID_UNKNOWN_VALUE) continue;
            RNScalar x0 = position_x_grid.GridValue(ix, iy0);
            if (x0 == R2_GRID_UNKNOWN_VALUE) continue;
            RNScalar y0 = position_y_grid.GridValue(ix, iy0);
            if (y0 == R2_GRID_UNKNOWN_VALUE) continue;

            // Compute weight of point 0
            RNScalar depth0A = viewpoint_depth_grid.GridValue(ix, iy0-1);
            RNScalar delta0A = (depth0A == R2_GRID_UNKNOWN_VALUE) ? depth0A - depth0 : 1;
            if (delta0A <= 0) continue;
            RNScalar depth0B = viewpoint_depth_grid.GridValue(ix, iy0-1);
            RNScalar delta0B = (depth0B == R2_GRID_UNKNOWN_VALUE) ? depth0B - depth0 : 1;
            if (delta0B <= 0) continue;
            RNScalar weight0 = delta0A * delta0B;

            // Find closest point in next vertical scanline
            int closest_index = -1;
            RNScalar closest_x1 = -1;
            RNScalar closest_y1 = -1;
            RNScalar closest_dd = max_distance * max_distance;
            for (int iy1 = 1; iy1 < scanline_index_grid.YResolution()-1; iy1++) {
              RNScalar depth1 = viewpoint_depth_grid.GridValue(ix+1, iy1);
              if (depth1 == R2_GRID_UNKNOWN_VALUE) continue;
              RNScalar x1 = position_x_grid.GridValue(ix+1, iy1);
              if (x1 == R2_GRID_UNKNOWN_VALUE) continue;
              RNScalar y1 = position_y_grid.GridValue(ix+1, iy1);
              if (y1 == R2_GRID_UNKNOWN_VALUE) continue;
              RNScalar dx = x1 - x0;
              RNScalar dy = y1 - y0;
              RNScalar dd = dx*dx + dy*dy;
              if (dd < closest_dd) {
                closest_index = iy1;
                closest_dd = dd;
              }
            }

            // Compute weight of point 1
            RNScalar depth1A = viewpoint_depth_grid.GridValue(ix+1, iy1-1);
            RNScalar delta1A = (depth1A == R2_GRID_UNKNOWN_VALUE) ? depth1A - depth1 : 1;
            if (delta1A <= 0) continue;
            RNScalar depth1B = viewpoint_depth_grid.GridValue(ix+1, iy1+1);
            RNScalar delta1B = (depth1B == R2_GRID_UNKNOWN_VALUE) ? depth1B - depth1 : 1;
            if (delta1B <= 0) continue;
            RNScalar weight1 = delta1A * delta1B;

            // Check if mutually closest
            RNBoolean mutually_closest = TRUE;
            for (int iy2 = 1; iy2 < scanline_index_grid.YResolution()-1; iy2++) {
              if (iy2 == iy0) continue;
              RNScalar depth2 = viewpoint_depth_grid.GridValue(ix, iy2);
              if (depth2 == R2_GRID_UNKNOWN_VALUE) continue;
              RNScalar x2 = position_x_grid.GridValue(ix, iy2);
              if (x2 == R2_GRID_UNKNOWN_VALUE) continue;
              RNScalar y2 = position_y_grid.GridValue(ix, iy2);
              if (y2 == R2_GRID_UNKNOWN_VALUE) continue;
              RNScalar dx = x2 - closest_x1;
              RNScalar dy = y2 - closest_y1;
              RNScalar dd = dx*dx + dy*dy;
              if (dd < closest_dd) {
                mutually_closest = FALSE;
                break;
              }
            }

            // Create feature correspondence between mutually closest points
            if (mutually_closest) {
              // Create feature 0
              RNScalar feature_index0_value = feature_index_grid.GridValue(ix, iy0);
              if (feature_index0_value == R2_GRID_UNKNOWN_VALUE) {
                RNScalar scanline_index0_value = scanline_index_grid.GridValue(ix, iy0);
                if (scanline_index0_value == R2_GRID_UNKNOWN_VALUE) continue;
                RNScalar point_index0_value = point_index_grid.GridValue(ix, iy0);
                if (point_index0_value == R2_GRID_UNKNOWN_VALUE) continue;
                GSVScanline *scanline0 = scan->Scanline((int) (scanline_index0_value + 0.5));
                int point_index0 = (int) (point_index0_value + 0.5);
                R3Point position0(x0, y0, 0);
                feature_index0_value = NFeatures();
                feature_index_grid.SetGridValue(ix, iy0, feature_index0_value);
                GSVFeature *feature = new GSVFeature(GSV_SCAN_POINT_FEATURE_TYPE, 
                  scanline0, point_index0, position0, R3zero_vector, R3zero_vector, 1, weight0);
              }

              // Create feature 1
              RNScalar feature_index1_value = feature_index_grid.GridValue(ix, iy1);
              if (feature_index1_value == R2_GRID_UNKNOWN_VALUE) {
                RNScalar scanline_index1_value = scanline_index_grid.GridValue(ix, iy1);
                if (scanline_index1_value == R2_GRID_UNKNOWN_VALUE) continue;
                RNScalar point_index1_value = point_index_grid.GridValue(ix, iy1);
                if (point_index1_value == R2_GRID_UNKNOWN_VALUE) continue;
                GSVScanline *scanline1 = scan->Scanline((int) (scanline_index1_value + 0.5));
                int point_index1 = (int) (point_index1_value + 0.5);
                R3Point position1(x1, y1, 0);
                feature_index1_value = NFeatures();
                feature_index_grid.SetGridValue(ix, iy1, feature_index1_value);
                GSVFeature *feature = new GSVFeature(GSV_SCAN_POINT_FEATURE_TYPE, 
                  scanline1, point_index1, position1, R3zero_vector, R3zero_vector, 1, weight1);
              }

              // Create correspondence
              int feature_index0 = (int) (feature_index0_value + 0.5);
              int feature_index1 = (int) (feature_index1_value + 0.5);
              GSVFeature *feature0 = Feature(feature_index0);
              GSVFeature *feature1 = Feature(feature_index1);
              GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(feature0, feature1, 1);
              InsertCorrespondence(correspondence);
            }
          }
        }
      }
    }
  }

  // Return success
  return 1;
}

#endif



////////////////////////////////////////////////////////////////////////
// PoseOptimization correspondence creation functions
////////////////////////////////////////////////////////////////////////

RNScalar GSVPoseOptimization::
PairScore(const GSVFeature *feature0, const GSVFeature *feature1) const
{
  return Score(feature0, feature1, 
    max_pair_world_distance, max_pair_world_distance_ratio, 
    max_pair_world_direction_angle, max_pair_world_normal_angle,
    max_pair_image_distance, max_pair_image_distance_ratio, max_pair_image_direction_angle, 
    max_pair_spin_image_descriptor_distance, max_pair_shape_context_descriptor_distance, 
    max_pair_sift_descriptor_distance, max_pair_line_descriptor_distance, 
    max_pair_descriptor_distance_ratio, min_pair_path_parameter_difference);
}



RNScalar GSVPoseOptimization::
CorrespondenceScore(const GSVFeature *feature0, const GSVFeature *feature1) const
{
  return Score(feature0, feature1, 
    max_correspondence_world_distance, max_correspondence_world_distance_ratio, 
    max_correspondence_world_direction_angle, max_correspondence_world_normal_angle,
    max_correspondence_image_distance, max_correspondence_image_distance_ratio, max_correspondence_image_direction_angle, 
    max_correspondence_spin_image_descriptor_distance, max_correspondence_shape_context_descriptor_distance, 
    max_correspondence_sift_descriptor_distance, max_correspondence_line_descriptor_distance, 
    max_correspondence_descriptor_distance_ratio, min_correspondence_path_parameter_difference);
}



RNScalar GSVPoseOptimization::
Score(const GSVFeature *feature0, const GSVFeature *feature1,
  RNScalar max_world_distance,
  RNScalar max_world_distance_ratio,
  RNScalar max_world_direction_angle,
  RNScalar max_world_normal_angle,
  RNScalar max_image_distance,
  RNScalar max_image_distance_ratio,
  RNScalar max_image_direction_angle,
  RNScalar max_spin_image_descriptor_distance,
  RNScalar max_shape_context_descriptor_distance,
  RNScalar max_sift_descriptor_distance,
  RNScalar max_line_descriptor_distance,
  RNScalar max_descriptor_distance_ratio,
  RNScalar min_path_parameter_difference) const
{
  // Check features
  if (feature0 == feature1) return -1;

  // Check scanlines
  GSVScanline *scanline0 = feature0->scanline;
  GSVScanline *scanline1 = feature1->scanline;
  if (scanline0 && scanline1 && (scanline0 == scanline1)) return -1;

  // Check images 
  GSVImage *image0 = feature0->image;
  GSVImage *image1 = feature1->image;
  if (image0 && image1 && (image0 == image1)) return -1;

  // Check path parameters
  if (min_path_parameter_difference > 0) {
    if (scanline0 && scanline1) {
      if (!image0 && !image1) {
        GSVSegment *segment0 = scanline0->Scan()->Segment();
        GSVSegment *segment1 = scanline1->Scan()->Segment();
        if (segment0 && segment1 && (segment0 == segment1)) {
          ScanlineData *scanline_data0 = (ScanlineData *) scanline0->Data();
          ScanlineData *scanline_data1 = (ScanlineData *) scanline1->Data();
          if (scanline_data0 && scanline_data1) {
            RNScalar u0 = scanline_data0->path_parameter;
            RNScalar u1 = scanline_data1->path_parameter;
            RNScalar path_parameter_difference = fabs(u0 - u1);
            if (path_parameter_difference < min_path_parameter_difference) return -1;
          }
        }
      }
    }
  }

  // Check world positions
  RNScalar world_distance_score = 1;
  if (max_world_distance > 0) {
    RNScalar distance_squared = SquaredWorldDistance(feature0, feature1);
    RNScalar max_distance_squared = max_world_distance * max_world_distance;
    if (distance_squared > max_distance_squared) return -1;
    if ((feature0->scan_scale > 0) || (feature1->scan_scale > 0)) {
      R3Point position0 = WorldPosition(feature0);
      R3Point position1 = WorldPosition(feature1);
      RNScalar position_distance_squared = R3SquaredDistance(position0, position1);
      if (feature0->scan_scale > 0) {
        RNScalar scale0 = feature0->scan_scale + max_world_distance;
        RNScalar scale0_squared = scale0 * scale0;
        if (position_distance_squared > scale0_squared) return -1;
      }
      if (feature1->scan_scale > 0) {
        RNScalar scale1 = feature1->scan_scale + max_world_distance;
        RNScalar scale1_squared = scale1 * scale1;
        if (position_distance_squared > scale1_squared) return -1;
      }
    }
    world_distance_score = exp(-distance_squared / max_distance_squared);
  }

  // Check world normals
  RNScalar world_normal_angle_score = 1;
  if (max_world_normal_angle > 0) {
    if ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) {
      RNScalar max_world_normal_angle_squared = max_world_normal_angle * max_world_normal_angle;
      R3Vector world_normal0 = WorldNormal(feature0, TRUE);
      R3Vector world_normal1 = WorldNormal(feature1, TRUE);
      if (RNIsZero(world_normal0.Dot(world_normal0))) return -1;
      if (RNIsZero(world_normal1.Dot(world_normal1))) return -1;
      RNScalar dot = fabs(world_normal0.Dot(world_normal1));
      RNAngle world_normal_angle = (dot < 1) ? acos(dot) : 0;
      if (world_normal_angle > max_world_normal_angle) return -1;
      RNScalar world_normal_angle_squared = world_normal_angle * world_normal_angle;
      world_normal_angle_score = exp(-world_normal_angle_squared / max_world_normal_angle_squared);
    }
  }

  // Check world directions
  RNScalar world_direction_angle_score = 1;
  if (max_world_direction_angle > 0) {
    if (((feature0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE))) {
      RNScalar max_world_direction_angle_squared = max_world_direction_angle * max_world_direction_angle;
      R3Vector world_direction0 = WorldDirection(feature0, TRUE);
      R3Vector world_direction1 = WorldDirection(feature1, TRUE);
      if (RNIsZero(world_direction0.Dot(world_direction0))) return -1;
      if (RNIsZero(world_direction1.Dot(world_direction1))) return -1;
      RNScalar dot = fabs(world_direction0.Dot(world_direction1));
      RNAngle world_direction_angle = (dot < 1) ? acos(dot) : 0;
      if (world_direction_angle > max_world_direction_angle) return -1;
      RNScalar world_direction_angle_squared = world_direction_angle * world_direction_angle;
      world_direction_angle_score = exp(-world_direction_angle_squared / max_world_direction_angle_squared);
    }
  }

  // Check image positions
  RNScalar image_distance_score = 1;
  if (max_image_distance > 0) {
    if (feature0->image) {
      R2Point position0 = ImagePosition(feature0, feature0);
      R2Point position1 = ImagePosition(feature1, feature0);
      RNScalar distance_squared = R2SquaredDistance(position0, position1);
      RNScalar max_distance = max_image_distance;
      if ((feature0->image_scale != RN_UNKNOWN) && (max_distance < feature0->image_scale)) max_distance = feature0->image_scale;
      if ((feature1->image_scale != RN_UNKNOWN) && (max_distance < feature1->image_scale)) max_distance = feature1->image_scale;
      RNScalar max_distance_squared = max_distance * max_distance;
      if (distance_squared > max_distance_squared) return -1;
      image_distance_score = exp(-distance_squared / max_distance_squared);
    }
    if (feature1->image) {
      R2Point position0 = ImagePosition(feature0, feature1);
      R2Point position1 = ImagePosition(feature1, feature1);
      RNScalar distance_squared = R2SquaredDistance(position0, position1);
      RNScalar max_distance = max_image_distance;
      if ((feature0->image_scale != RN_UNKNOWN) && (max_distance < feature0->image_scale)) max_distance = feature0->image_scale;
      if ((feature1->image_scale != RN_UNKNOWN) && (max_distance < feature1->image_scale)) max_distance = feature1->image_scale;
      RNScalar max_distance_squared = max_distance * max_distance;
      if (distance_squared > max_distance_squared) return -1;
      image_distance_score = exp(-distance_squared / max_distance_squared);
    }
  }
  
  // Check image direction angles
  RNScalar image_direction_angle_score = 1;
  if (max_image_direction_angle > 0) {
    // Not implemented
  }

  // Check descriptors
  RNScalar descriptor_distance_score = 1;
  const GSVDescriptor& descriptor0 = feature0->descriptor;
  const GSVDescriptor& descriptor1 = feature1->descriptor;
  if (descriptor0.descriptor_type != descriptor1.descriptor_type) return -1;
  RNScalar max_d = 0;
  if (descriptor0.descriptor_type == GSV_SPIN_IMAGE_DESCRIPTOR_TYPE) 
    max_d = max_spin_image_descriptor_distance;
  else if (descriptor0.descriptor_type == GSV_SHAPE_CONTEXT_DESCRIPTOR_TYPE) 
    max_d = max_shape_context_descriptor_distance;
  else if (descriptor0.descriptor_type == GSV_SIFT_DESCRIPTOR_TYPE) 
    max_d = max_sift_descriptor_distance;
  else if (descriptor0.descriptor_type == GSV_LINE_DESCRIPTOR_TYPE) 
    max_d = max_line_descriptor_distance;
  RNScalar max_dd = max_d * max_d;
  if (max_dd > 0) {
    RNScalar dd = descriptor0.SquaredDistance(descriptor1);
    if (dd > max_dd) return -1;
    descriptor_distance_score = exp(-dd/max_dd);
  }

  // Compute score
  RNScalar score = feature0->score * feature1->score * 
    world_distance_score * world_normal_angle_score * world_direction_angle_score *
    image_distance_score * image_direction_angle_score * descriptor_distance_score;

  // Print message
  // printf("%9.6f : %5.3f %5.3f : %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", 
  //   score, feature0->score, feature1->score, 
  //   world_distance_score, world_normal_angle_score, world_direction_angle_score,
  //   image_distance_score, image_direction_angle_score, descriptor_distance_score);

  // Return score
  return score;
}



////////////////////////////////////////////////////////////////////////
// Optimized position and direction extraction
////////////////////////////////////////////////////////////////////////

R3Point GSVPoseOptimization::
WorldPosition(const GSVFeatureCluster *cluster, RNBoolean apply_pose_transformation) const
{
  // Return cluster world position (possibly with optimization transformation applied)
  R3Point result = cluster->scan_position;
  if (apply_pose_transformation) result.Transform(OptimizedTransformation(cluster));
  return result;
}



R3Vector GSVPoseOptimization::
WorldDirection(const GSVFeatureCluster *cluster, RNBoolean apply_pose_transformation) const
{
  // Return cluster world direction (possibly with optimization transformation applied)
  R3Vector result = cluster->scan_direction;
  if (apply_pose_transformation) result.Transform(OptimizedTransformation(cluster));
  return result;
}



R3Vector GSVPoseOptimization::
WorldNormal(const GSVFeatureCluster *cluster, RNBoolean apply_pose_transformation) const
{
  // Return cluster world normal (possibly with optimization transformation applied)
  R3Vector result = cluster->scan_normal;
  if (apply_pose_transformation) result.Transform(OptimizedTransformation(cluster));
  return result;
}



R3Point GSVPoseOptimization::
WorldPosition(const GSVFeature *feature, RNBoolean apply_pose_transformation) const
{
  // Return feature world position (possibly with optimization transformation applied)
  R3Point result = feature->scan_position;
  if (apply_pose_transformation) result.Transform(OptimizedTransformation(feature));
  return result;
}



R3Vector GSVPoseOptimization::
WorldDirection(const GSVFeature *feature, RNBoolean apply_pose_transformation) const
{
  // Return feature world direction (possibly with optimization transformation applied)
  R3Vector result = feature->scan_direction;
  if (apply_pose_transformation) result.Transform(OptimizedTransformation(feature));
  return result;
}



R3Vector GSVPoseOptimization::
WorldNormal(const GSVFeature *feature, RNBoolean apply_pose_transformation) const
{
  // Return feature world normal (possibly with optimization transformation applied)
  R3Vector result = feature->scan_normal;
  if (apply_pose_transformation) result.Transform(OptimizedTransformation(feature));
  return result;
}



R2Point GSVPoseOptimization::
ImagePosition(const GSVFeature *feature, const GSVFeature *image_feature, RNBoolean apply_pose_transformation) const
{
  // Return feature position in image (possibly with optimization transformation applied)
  if (feature->image == image_feature->image) return feature->image_position;
  R3Point world_position = feature->scan_position;
  if (apply_pose_transformation) {
    world_position.InverseTransform(OptimizedTransformation(image_feature));
    world_position.Transform(OptimizedTransformation(feature));
  }
  int column_index = (int) image_feature->image_position.X();
  assert((column_index >= 0) && (column_index < image_feature->image->Width()));
  return image_feature->image->UndistortedPosition(world_position, column_index);
}



R2Vector GSVPoseOptimization::
ImageDirection(const GSVFeature *feature, const GSVFeature *image_feature, RNBoolean apply_pose_transformation) const
{
  // Return feature direction (possibly with optimization transformation applied)
  if (feature->image == image_feature->image) return feature->image_direction;
  R3Point world_position = feature->scan_position;
  R3Vector world_direction = feature->scan_direction;
  if (apply_pose_transformation) {
    world_position.InverseTransform(OptimizedTransformation(image_feature));
    world_position.Transform(OptimizedTransformation(feature));
    world_direction.InverseTransform(OptimizedTransformation(image_feature));
    world_direction.Transform(OptimizedTransformation(feature));
  }
  int column_index = (int) image_feature->image_position.X();
  assert((column_index >= 0) && (column_index < image_feature->image->Width()));
  R2Point projectionA = image_feature->image->UndistortedPosition(world_position, column_index);
  R2Point projectionB = image_feature->image->UndistortedPosition(world_position + world_direction, column_index);
  return projectionB - projectionA;
}



////////////////////////////////////////////////////////////////////////
// Distance computation
////////////////////////////////////////////////////////////////////////

RNScalar GSVPoseOptimization::
SquaredWorldDistance(const GSVFeature *feature0, const GSVFeature *feature1, RNBoolean apply_pose_transformation) const
{
  // This function assumes that the scan_position of image features is uptodate

  // Get transformations
  R3Affine transformation0 = OptimizedTransformation(feature0);
  R3Affine transformation1 = OptimizedTransformation(feature1);

  // Compute squared distance between feature scan positions (in sq. meters)
  RNScalar sum = 0;
  int count = 0;

  if (((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE))) {
    R3Point point0 = feature0->scan_position;
    R3Point point1 = feature1->scan_position;
    if (apply_pose_transformation) point0.Transform(transformation0);
    if (apply_pose_transformation) point1.Transform(transformation1);
    sum += R3SquaredDistance(point0, point1);
    count++;
  }

  if (((feature0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE))) {
    R3Line line0(feature0->scan_position, feature0->scan_direction);
    R3Point point1 = feature1->scan_position;
    if (apply_pose_transformation) line0.Transform(transformation0);
    if (apply_pose_transformation) point1.Transform(transformation1);
    RNScalar distance = R3Distance(line0, point1);
    sum += distance * distance;
    count++;
  }    

  if (((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE))) {
    R3Plane plane0(feature0->scan_position, feature0->scan_normal);
    R3Point point1 = feature1->scan_position;
    if (apply_pose_transformation) plane0.Transform(transformation0);
    if (apply_pose_transformation) point1.Transform(transformation1);
    RNScalar distance = R3Distance(plane0, point1);
    sum += distance * distance;
    count++;
  }    

  if (((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE))) {
    R3Point point0 = feature0->scan_position;
    R3Point point1 = feature1->scan_position;
    R2Point point0A = feature0->image_position - feature0->image_scale * feature0->image_direction;
    R2Point point0B = feature0->image_position + feature0->image_scale * feature0->image_direction;
    R3Ray ray0A = feature0->image->RayThroughUndistortedPosition(point0A);
    R3Ray ray0B = feature0->image->RayThroughUndistortedPosition(point0B);
    R3Vector normal0 = ray0A.Vector() % ray0B.Vector(); normal0.Normalize();
    R3Plane plane0(point0, normal0);
    if (apply_pose_transformation) plane0.Transform(transformation0);
    if (apply_pose_transformation) point1.Transform(transformation1);
    RNLength distance = R3Distance(plane0, point1);
    sum += distance * distance;
    count++;
  }      

  if (((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE))) {
    R3Point point0 = feature0->scan_position;
    R3Line line1(feature1->scan_position, feature1->scan_direction);
    if (apply_pose_transformation) point0.Transform(transformation0);
    if (apply_pose_transformation) line1.Transform(transformation1);
    RNScalar distance = R3Distance(line1, point0);
    sum += distance * distance;
    count++;
  }      

  if (((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE))) {
    R3Point point0 = feature0->scan_position;
    R3Plane plane1(feature1->scan_position, feature1->scan_normal);
    if (apply_pose_transformation) point0.Transform(transformation0);
    if (apply_pose_transformation) plane1.Transform(transformation1);
    RNScalar distance = R3Distance(plane1, point0);
    sum += distance * distance;
    count++;
  }      

  if (((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE))) {
    R3Point point0 = feature0->scan_position;
    R3Point point1 = feature1->scan_position;
    R2Point point1A = feature1->image_position - feature1->image_scale * feature1->image_direction;
    R2Point point1B = feature1->image_position + feature1->image_scale * feature1->image_direction;
    R3Ray ray1A = feature1->image->RayThroughUndistortedPosition(point1A);
    R3Ray ray1B = feature1->image->RayThroughUndistortedPosition(point1B);
    R3Vector normal1 = ray1A.Vector() % ray1B.Vector(); normal1.Normalize();
    R3Plane plane1(point1, normal1);
    if (apply_pose_transformation) point0.Transform(transformation0);
    if (apply_pose_transformation) plane1.Transform(transformation1);
    RNLength distance = R3Distance(point0, plane1);
    sum += distance * distance;
    count++;
  }      

  // Just checking
  if (count == 0) RNAbort("Unrecognized feature types in SquaredWorldDistance");

  // Return average
  return sum / count;
}



RNScalar GSVPoseOptimization::
SquaredImageDistance(const GSVFeature *feature0, const GSVFeature *feature1, RNBoolean apply_pose_transformation) const
{
  // This function assumes that the scan_position of image features is uptodate

  // Check feature types
  if (!feature0->image && !feature1->image) return RN_INFINITY;
  if (feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) return RN_INFINITY;
  if (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) return RN_INFINITY;

  // Set maximum error (for projected features that lie outside viewport)
  RNScalar max_squared_distance = 1000 * 1000;

  // Get transformations
  R3Affine transformation_0_to_1 = R3identity_affine;
  R3Affine transformation_1_to_0 = R3identity_affine;
  if (feature0->image || feature1->image) {
    R3Affine transformation0 = OptimizedTransformation(feature0);
    R3Affine transformation1 = OptimizedTransformation(feature1);
    transformation_0_to_1.InverseTransform(transformation1);
    transformation_0_to_1.Transform(transformation0);
    transformation_1_to_0.InverseTransform(transformation0);
    transformation_1_to_0.Transform(transformation1);
  }

  // Compute squared distance between feature image positions (in sq. pixels)
  RNScalar sum = 0;
  int count = 0;

  if (((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE))) {
    R2Point projection0 = feature0->image_position;
    R3Point point1 = feature1->scan_position;
    if (apply_pose_transformation) point1.Transform(transformation_1_to_0);
    int column_index = (int) feature0->image_position.X();
    assert((column_index >= 0) && (column_index < feature0->image->Width()));
    R2Point projection1 = feature0->image->UndistortedPosition(point1, column_index);
    if ((projection1.X() < 0) || (projection1.Y() < 0)) return max_squared_distance;
    sum += R2SquaredDistance(projection0, projection1);
    count++;
  }

  if (((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE))) {
    R2Line projection0(feature0->image_position, feature0->image_direction);
    R3Point point1 = feature1->scan_position;
    if (apply_pose_transformation) point1.Transform(transformation_1_to_0);
    int column_index = (int) feature0->image_position.X();
    assert((column_index >= 0) && (column_index < feature0->image->Width()));
    R2Point projection1 = feature0->image->UndistortedPosition(point1, column_index);
    if ((projection1.X() < 0) || (projection1.Y() < 0)) return max_squared_distance;
    RNLength distance = R2Distance(projection0, projection1);
    sum += distance * distance;
    count++;
  }      

  if (((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)  && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)   && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)  && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)  && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE))) {
    R3Point point0 = feature0->scan_position;
    if (apply_pose_transformation) point0.Transform(transformation_0_to_1);
    int column_index = (int) feature1->image_position.X();
    assert((column_index >= 0) && (column_index < feature1->image->Width()));
    R2Point projection0 = feature1->image->UndistortedPosition(point0, column_index);
    if ((projection0.X() < 0) || (projection0.Y() < 0)) return max_squared_distance;
    R2Point projection1 = feature1->image_position;
    sum += R2SquaredDistance(projection0, projection1);
    count++;
  }

  if (((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)  && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)   && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)  && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) ||
      ((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)  && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE))) {
    R3Point point0 = feature0->scan_position;
    if (apply_pose_transformation) point0.Transform(transformation_0_to_1);
    int column_index = (int) feature1->image_position.X();
    assert((column_index >= 0) && (column_index < feature1->image->Width()));
    R2Point projection0 = feature1->image->UndistortedPosition(point0, column_index);
    if ((projection0.X() < 0) || (projection0.Y() < 0)) return max_squared_distance;
    R2Line projection1(feature1->image_position, feature1->image_direction);
    RNLength distance = R2Distance(projection0, projection1);
    sum += distance * distance;
    count++;
  }

  // Return average
  if (count == 0) return RN_INFINITY;
  RNScalar squared_distance = sum / count;
  if (squared_distance >= max_squared_distance) squared_distance = max_squared_distance - 1;
  return squared_distance;
}



RNScalar GSVPoseOptimization::
SquaredWorldDistance(const GSVFeatureCluster *cluster0, const GSVFeature *feature1, RNBoolean apply_pose_transformation) const
{
  // This function assumes that the scan_position of image features is uptodate

  // Get transformation
  R3Affine transformation0 = OptimizedTransformation(cluster0);
  R3Affine transformation1 = OptimizedTransformation(feature1);

  // Compute squared distance between feature scan positions (in sq. meters)
  RNScalar sum = 0;
  int count = 0;

  if (((cluster0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE))) {
    R3Point point0 = cluster0->scan_position;
    R3Point point1 = feature1->scan_position;
    if (apply_pose_transformation) point0.Transform(transformation0);
    if (apply_pose_transformation) point1.Transform(transformation1);
    sum += R3SquaredDistance(point0, point1);
    count++;
  }

  if (((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE))) {
    R3Line line0(cluster0->scan_position, cluster0->scan_direction);
    R3Point point1 = feature1->scan_position;
    if (apply_pose_transformation) line0.Transform(transformation0);
    if (apply_pose_transformation) point1.Transform(transformation1);
    RNScalar distance = R3Distance(line0, point1);
    sum += distance * distance;
    count++;
  }    

  if (((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE))) {
    R3Plane plane0(cluster0->scan_position, cluster0->scan_normal);
    R3Point point1 = feature1->scan_position;
    if (apply_pose_transformation) plane0.Transform(transformation0);
    if (apply_pose_transformation) point1.Transform(transformation1);
    RNScalar distance = R3Distance(plane0, point1);
    sum += distance * distance;
    count++;
  }    

  if (((cluster0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE))) {
    R3Point point0 = cluster0->scan_position;
    R3Line line1(feature1->scan_position, feature1->scan_direction);
    if (apply_pose_transformation) point0.Transform(transformation0);
    if (apply_pose_transformation) line1.Transform(transformation1);
    RNScalar distance = R3Distance(line1, point0);
    sum += distance * distance;
    count++;
  }      

  if (((cluster0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE))) {
    R3Point point0 = cluster0->scan_position;
    R3Plane plane1(feature1->scan_position, feature1->scan_normal);
    if (apply_pose_transformation) point0.Transform(transformation0);
    if (apply_pose_transformation) plane1.Transform(transformation1);
    RNScalar distance = R3Distance(plane1, point0);
    sum += distance * distance;
    count++;
  }      

  if (((cluster0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE))) {
    R3Point point0 = cluster0->scan_position;
    R3Point point1 = feature1->scan_position;
    R2Point point1A = feature1->image_position - feature1->image_scale * feature1->image_direction;
    R2Point point1B = feature1->image_position + feature1->image_scale * feature1->image_direction;
    R3Ray ray1A = feature1->image->RayThroughUndistortedPosition(point1A);
    R3Ray ray1B = feature1->image->RayThroughUndistortedPosition(point1B);
    R3Vector normal1 = ray1A.Vector() % ray1B.Vector(); normal1.Normalize();
    R3Plane plane1(point1, normal1);
    if (apply_pose_transformation) point0.Transform(transformation0);
    if (apply_pose_transformation) plane1.Transform(transformation1);
    RNLength distance = R3Distance(point0, plane1);
    sum += distance * distance;
    count++;
  }      

  // Return average
  if (count == 0) return RN_INFINITY;
  return sum / count;
}



RNScalar GSVPoseOptimization::
SquaredImageDistance(const GSVFeatureCluster *cluster0, const GSVFeature *feature1, RNBoolean apply_pose_transformation) const
{
  // This function assumes that the scan_position of image features is uptodate

  // Check feature types
  if (!feature1->image) return RN_INFINITY;
  if (cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) return RN_INFINITY;

  // Set maximum error (for projected features that lie outside viewport)
  RNScalar max_squared_distance = 1000 * 1000;

  // Get transformations
  R3Affine transformation_0_to_1 = R3identity_affine;
  R3Affine transformation_1_to_0 = R3identity_affine;
  if (feature1->image) {
    R3Affine transformation0 = OptimizedTransformation(cluster0);
    R3Affine transformation1 = OptimizedTransformation(feature1);
    transformation_0_to_1.InverseTransform(transformation1);
    transformation_0_to_1.Transform(transformation0);
    transformation_1_to_0.InverseTransform(transformation0);
    transformation_1_to_0.Transform(transformation1);
  }


  // Compute squared distance between feature image positions (in sq. pixels)
  RNScalar sum = 0;
  int count = 0;

  if (((cluster0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)  && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)   && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)  && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE))) {
    R3Point point0 = cluster0->scan_position;
    if (apply_pose_transformation) point0.Transform(transformation_0_to_1);
    int column_index = (int) feature1->image_position.X();
    assert((column_index >= 0) && (column_index < feature1->image->Width()));
    R2Point projection0 = feature1->image->UndistortedPosition(point0, column_index);
    if ((projection0.X() < 0) || (projection0.Y() < 0)) return max_squared_distance;
    R2Point projection1 = feature1->image_position;
    sum += R2SquaredDistance(projection0, projection1);
    count++;
  }

  if (((cluster0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)  && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)   && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)  && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE))) {
    R3Point point0 = cluster0->scan_position;
    if (apply_pose_transformation) point0.Transform(transformation_0_to_1);
    int column_index = (int) feature1->image_position.X();
    assert((column_index >= 0) && (column_index < feature1->image->Width()));
    R2Point projection0 = feature1->image->UndistortedPosition(point0, column_index);
    if ((projection0.X() < 0) || (projection0.Y() < 0)) return max_squared_distance;
    R2Line projection1(feature1->image_position, feature1->image_direction);
    RNLength distance = R2Distance(projection0, projection1);
    sum += distance * distance;
    count++;
  }

  // Return average
  if (count == 0) return RN_INFINITY;
  RNScalar squared_distance = sum / count;
  if (squared_distance >= max_squared_distance) squared_distance = max_squared_distance - 1;
  return squared_distance;
}



RNScalar GSVPoseOptimization::
SquaredWorldDistance(const GSVFeatureCluster *cluster0, const GSVFeatureCluster *cluster1, RNBoolean apply_pose_transformation) const
{
  // This function assumes that the scan_position of image features is uptodate

  // Get transformation
  R3Affine transformation0 = OptimizedTransformation(cluster0);
  R3Affine transformation1 = OptimizedTransformation(cluster1);

  // Compute squared distance between feature scan positions (in sq. meters)
  RNScalar sum = 0;
  int count = 0;

  if ((cluster0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) {
    R3Point point0 = cluster0->scan_position;
    R3Point point1 = cluster1->scan_position;
    if (apply_pose_transformation) point0.Transform(transformation0);
    if (apply_pose_transformation) point1.Transform(transformation1);
    sum += R3SquaredDistance(point0, point1);
    count++;
  }

  if (((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE))) {
    R3Line line0(cluster0->scan_position, cluster0->scan_direction);
    R3Point point1 = cluster1->scan_position;
    if (apply_pose_transformation) line0.Transform(transformation0);
    if (apply_pose_transformation) point1.Transform(transformation1);
    RNScalar distance = R3Distance(line0, point1);
    sum += distance * distance;
    count++;
  }    

  if (((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE))) {
    R3Plane plane0(cluster0->scan_position, cluster0->scan_normal);
    R3Point point1 = cluster1->scan_position;
    if (apply_pose_transformation) plane0.Transform(transformation0);
    if (apply_pose_transformation) point1.Transform(transformation1);
    RNScalar distance = R3Distance(plane0, point1);
    sum += distance * distance;
    count++;
  }    

  if (((cluster0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE))) {
    R3Point point0 = cluster0->scan_position;
    R3Line line1(cluster1->scan_position, cluster1->scan_direction);
    if (apply_pose_transformation) point0.Transform(transformation0);
    if (apply_pose_transformation) line1.Transform(transformation1);
    RNScalar distance = R3Distance(line1, point0);
    sum += distance * distance;
    count++;
  }      

  if (((cluster0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
      ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE))) {
    R3Point point0 = cluster0->scan_position;
    R3Plane plane1(cluster1->scan_position, cluster1->scan_normal);
    if (apply_pose_transformation) point0.Transform(transformation0);
    if (apply_pose_transformation) plane1.Transform(transformation1);
    RNScalar distance = R3Distance(plane1, point0);
    sum += distance * distance;
    count++;
  }      

  // Return average
  if (count == 0) return RN_INFINITY;
  return sum / count;
}



RNScalar GSVPoseOptimization::
SquaredDescriptorDistance(const GSVFeature *feature0, const GSVFeature *feature1) const
{
  // Return squared distance in descriptor space
  const GSVDescriptor& descriptor0 = feature0->descriptor;
  const GSVDescriptor& descriptor1 = feature1->descriptor;
  return descriptor0.SquaredDistance(descriptor1);
}



////////////////////////////////////////////////////////////////////////
// Pair creation functions
////////////////////////////////////////////////////////////////////////

static R3Point 
WorldPositionCallback(GSVFeature *feature, void *data)
{
  GSVPoseOptimization *optimization = (GSVPoseOptimization *) data;
  return optimization->WorldPosition(feature);
}



static RNArray<GSVFeature *> *
CreateFeatureAssociations(GSVPoseOptimization *optimization,
  RNBoolean mask_pairs, RNBoolean mask_correspondences, RNBoolean mask_clusters)
{
  // Initialize result
  if (optimization->NFeatures() == 0) return NULL;
  if (((optimization->NPairs() == 0) || !mask_pairs) && 
      ((optimization->NCorrespondences() == 0) || !mask_correspondences) && 
      ((optimization->NClusters() == 0) || !mask_clusters)) {
    return NULL;
  }

  // Allocate data
  RNArray<GSVFeature *> *feature_masks = new RNArray<GSVFeature *>[optimization->NFeatures()];
  if (!feature_masks) return NULL;

  if (mask_pairs) {
    for (int i = 0; i < optimization->NPairs(); i++) {
      GSVFeaturePair *pair = optimization->Pair(i);
      GSVFeature *feature0 = pair->Feature(0);
      GSVFeature *feature1 = pair->Feature(1);
      if (!feature_masks[feature0->index].FindEntry(feature1)) 
        feature_masks[feature0->index].Insert(feature1);
      if (!feature_masks[feature1->index].FindEntry(feature0)) 
        feature_masks[feature1->index].Insert(feature0);
    }
  }

  if (mask_correspondences) {
    for (int i = 0; i < optimization->NCorrespondences(); i++) {
      GSVFeatureCorrespondence *correspondence = optimization->Correspondence(i);
      GSVFeature *feature0 = correspondence->Feature(0);
      GSVFeature *feature1 = correspondence->Feature(1);
      if (!feature_masks[feature0->index].FindEntry(feature1)) 
        feature_masks[feature0->index].Insert(feature1);
      if (!feature_masks[feature1->index].FindEntry(feature0)) 
        feature_masks[feature1->index].Insert(feature0);
    }
  }

  if (mask_clusters) {
    for (int i = 0; i < optimization->NClusters(); i++) {
      GSVFeatureCluster *cluster = optimization->Cluster(i);
      for (int j = 0; j < cluster->NFeatures(); j++) {
        GSVFeature *feature0 = cluster->Feature(j);
        for (int k = j+1; k < cluster->NFeatures(); k++) {
          GSVFeature *feature1 = cluster->Feature(k);
          if (!feature_masks[feature0->index].FindEntry(feature1)) 
            feature_masks[feature0->index].Insert(feature1);
          if (!feature_masks[feature1->index].FindEntry(feature0)) 
            feature_masks[feature1->index].Insert(feature0);
        }
      }
    }
  }

  // Return feature mask
  return feature_masks;
}



int GSVPoseOptimization::
CreateAllPairs(RNBoolean check_duplicates)
{
  // Create masks of features not to be paired with each feature
  RNArray<GSVFeature *> *feature_masks = NULL;
  if (check_duplicates) feature_masks = CreateFeatureAssociations(this, TRUE, TRUE, TRUE);

  // Consider every pair of features
  for (int i0 = 0; i0 < NFeatures(); i0++) {
    GSVFeature *feature0 = Feature(i0);
    for (int i1 = 0; i1 < NFeatures(); i1++) {
      GSVFeature *feature1 = Feature(i1);
      if (feature_masks) {
        if (feature_masks[feature0->index].FindEntry(feature1)) continue;
        assert(!feature_masks[feature1->index].FindEntry(feature0));
        feature_masks[feature0->index].Insert(feature1);
        feature_masks[feature1->index].Insert(feature0);
      }
      GSVFeaturePair *pair = new GSVFeaturePair(feature0, feature1, score);
      InsertPair(pair);
    }
  }

  // Delete feature mask
  if (feature_masks) delete [] feature_masks;

  // Return success
  return 1;
}




int GSVPoseOptimization::
CreateFeasiblePairs(RNBoolean check_duplicates)
{
  // Create masks of features not to be paired with each feature
  RNArray<GSVFeature *> *feature_masks = NULL;
  if (check_duplicates) feature_masks = CreateFeatureAssociations(this, TRUE, TRUE, TRUE);

  // Consider every pair of features
  for (int i0 = 0; i0 < NFeatures(); i0++) {
    GSVFeature *feature0 = Feature(i0);
    for (int i1 = 0; i1 < NFeatures(); i1++) {
      GSVFeature *feature1 = Feature(i1);
      if (feature_masks) {
        if (feature_masks[feature0->index].FindEntry(feature1)) continue;
        assert(!feature_masks[feature1->index].FindEntry(feature0));
        feature_masks[feature0->index].Insert(feature1);
        feature_masks[feature1->index].Insert(feature0);
      }
      RNScalar score = PairScore(feature0, feature1);
      if (score < 0) continue;
      GSVFeaturePair *pair = new GSVFeaturePair(feature0, feature1, score);
      InsertPair(pair);
    }
  }

  // Delete feature pairs
  if (feature_masks) delete [] feature_masks;

  // Return success
  return 1;
}




int GSVPoseOptimization::
CreateCoplanarPairs(RNBoolean check_duplicates)
{
  // Parameters
  RNLength min_coplanar_pair_distance = 15;
  RNLength max_coplanar_pair_distance = 100;
  RNLength min_coplanar_pair_squared_distance = min_coplanar_pair_distance * min_coplanar_pair_distance;
  RNLength max_coplanar_pair_squared_distance = max_coplanar_pair_distance * max_coplanar_pair_distance;
  
  // Check features
  if (NFeatures() == 0) return 1;
  
  // Create masks of features not to be paired with each feature
  RNArray<GSVFeature *> *feature_masks = NULL;
  if (check_duplicates) feature_masks = CreateFeatureAssociations(this, TRUE, TRUE, TRUE);

  // Create pairs between features in same plane segment
  for (int i = 0; i < NFeatures(); i++) {
    GSVFeature *feature1 = Feature(i);
    if (feature1->feature_type != GSV_SCAN_PLANE_FEATURE_TYPE) continue;
    if (feature1->group == 0) continue;
    R3Point position1 = feature1->scan_position;
    for (int j = i+1; j < NFeatures(); j++) {
      GSVFeature *feature2 = Feature(j);
      if (feature2->feature_type != GSV_SCAN_PLANE_FEATURE_TYPE) continue;
      if (feature2->group == 0) continue;
      if (feature1->group != feature2->group) continue;
      R3Point position2 = feature2->scan_position;
      RNLength squared_distance = R3SquaredDistance(position1, position2);
      if (squared_distance < min_coplanar_pair_squared_distance) continue;
      if (squared_distance > max_coplanar_pair_squared_distance) continue;
      if (feature_masks) {
        if (feature_masks[feature1->index].FindEntry(feature2)) continue;
        assert(!feature_masks[feature2->index].FindEntry(feature1));
        feature_masks[feature1->index].Insert(feature2);
        feature_masks[feature2->index].Insert(feature1);
      }
      RNScalar score = feature1->score * feature2->score;
      GSVFeaturePair *pair = new GSVFeaturePair(feature1, feature2, score);
      InsertPair(pair);
    }
  }

  // Delete temporary data
  if (feature_masks) delete [] feature_masks;
 
  // Return success
  return 1;
}



int GSVPoseOptimization::
CreateClosePairs(RNBoolean check_duplicates)
{
  // Set parameters
  // mutually closest and ratio only consider other scan features within this cutoff
  RNScalar travel_distance_cutoff = 100; 

  ////////////////////

  // Create masks of features not to be paired with each feature
  RNArray<GSVFeature *> *feature_masks = NULL;
  if (check_duplicates) feature_masks = CreateFeatureAssociations(this, TRUE, TRUE, TRUE);

  // Create list of feature that make good candidates for every feature
  RNArray<GSVFeature *> *feature_candidates = new RNArray<GSVFeature *>[ NFeatures() ];
  if (max_pair_world_distance > 0) {
    // Create array of features (for kd-tree)
    RNArray<GSVFeature *> features;
    for (int i = 0; i < NFeatures(); i++) {
      GSVFeature *feature = Feature(i);
      if (feature->scan_position == R3unknown_point) continue;
      features.Insert(feature);
    }
    
    // Create kdtree of features
    R3Kdtree<GSVFeature *> kdtree(features, WorldPositionCallback, this);

    // Find candidate features 
    for (int index0 = 0; index0 < NFeatures(); index0++) {
      GSVFeature *feature0 = Feature(index0);
      R3Point position0 = WorldPosition(feature0);
      if (position0 == R3unknown_point) continue;
      RNArray<GSVFeature *> close_features;
      if (!kdtree.FindAll(position0, 0, max_pair_world_distance, close_features)) continue;
      for (int k = 0; k < close_features.NEntries(); k++) {
        GSVFeature *feature1 = close_features.Kth(k);
        if (feature1 == feature0) continue;
        int index1 = feature1->index;
        if (index1 <= index0) continue;
        assert(R3Distance(WorldPosition(feature0), WorldPosition(feature1)) <= max_pair_world_distance);
        if (feature_masks) {
          if (feature_masks[feature0->index].FindEntry(feature1)) continue;
          assert(!feature_masks[feature1->index].FindEntry(feature0));
          feature_masks[feature0->index].Insert(feature1);
          feature_masks[feature1->index].Insert(feature0);
        }
        if (PairScore(feature0, feature1) < 0) continue;
        feature_candidates[feature0->index].Insert(feature1);
        feature_candidates[feature1->index].Insert(feature0);
      }
    }
  }
  else {
    RNAbort("Not implemented");
    return 0;
  }

  ////////////////////

  // Create pairs from candidates
  for (int index0 = 0; index0 < NFeatures(); index0++) {
    GSVFeature *feature0 = Feature(index0);
    if (feature_candidates[index0].IsEmpty()) continue;
    R3Point position0 = WorldPosition(feature0);
    if (position0 == R3unknown_point) continue;
    GSVScanline *scanline0 = feature0->scanline;
    GSVScan *scan0 = (scanline0) ? scanline0->Scan() : NULL;
    GSVSegment *segment0 = (scan0) ? scan0->Segment() : NULL;
    RNLength travel0 = (scanline0) ? scanline0->TravelDistance() : 0.0;
    GSVDescriptor& descriptor0 = feature0->descriptor;
    GSVImage *image0 = feature0->image;

    // Consider all candidates for feature0
    for (int j = 0; j < feature_candidates[index0].NEntries(); j++) {
      GSVFeature *feature1 = feature_candidates[index0].Kth(j);
      int index1 = feature1->index;
      if (index1 <= index0) continue;
      
      // Compute score
      RNScalar score = PairScore(feature0, feature1);
      assert(score > 0);

      // Check relationships to other features
      if ((max_pair_world_distance_ratio > 0) || 
          (max_pair_descriptor_distance_ratio > 0)) {
        // Get stuff for feature1
        GSVScanline *scanline1 = feature1->scanline;
        GSVScan *scan1 = (scanline1) ? scanline1->Scan() : NULL;
        GSVSegment *segment1 = (scan1) ? scan1->Segment() : NULL;
        RNLength travel1 = (scanline1) ? scanline1->TravelDistance() : 0.0;
        GSVImage *image1 = feature1->image;
        R3Point position1 = WorldPosition(feature1);
        assert(position1 != R3unknown_point);
        GSVDescriptor& descriptor1 = feature1->descriptor;

        // Get distances between feature0 and feature1
        RNLength spatial_d01 = (max_pair_world_distance_ratio > 0.0) ? R3SquaredDistance(position0, position1) : 0.0;
        RNLength descriptor_d01 = (max_pair_descriptor_distance_ratio > 0.0) ? descriptor0.SquaredDistance(descriptor1) : 0.0;

        //////// 

        // Find distances from feature0 to closest other candidates
        RNScalar spatial_d02 = FLT_MAX;
        RNScalar descriptor_d02 = FLT_MAX;
        for (int k = 0; k < feature_candidates[index0].NEntries(); k++) {
          GSVFeature *feature2 = feature_candidates[index0].Kth(k);
          if (feature2 == feature1) continue;

          // Check if not in same region of same segment
          if (scanline1) {
            GSVScanline *scanline2 = feature2->scanline;
            if (scanline2) {
              GSVScan *scan2 = scanline2->Scan();
              GSVSegment *segment2 = scan2->Segment();
              if (segment1 != segment2) continue;
              if (travel_distance_cutoff > 0) {
                RNLength travel2 = scanline2->TravelDistance();
                if (fabs(travel1 - travel2) > travel_distance_cutoff) continue;
              }
            }
          }

          // Check if not in same image
          if (image1) {
            GSVImage *image2 = feature2->image;
            if (image2) {
              if (image1 != image2) continue;
            }
          }

          // Update closest spatial distance
          if (max_pair_world_distance_ratio > 0) {
            R3Point position2 = WorldPosition(feature2);
            RNLength spatial_d = R3SquaredDistance(position0, position2);
            if (spatial_d < spatial_d02) spatial_d02 = spatial_d; 
          }

          // Update closest descriptor distance
          if (max_pair_descriptor_distance_ratio > 0) {
            GSVDescriptor& descriptor2 = feature2->descriptor;
            RNLength descriptor_d = descriptor0.SquaredDistance(descriptor2);
            if (descriptor_d < descriptor_d02) descriptor_d02 = descriptor_d; 
          }
        }
          
        // Check relationships between spatial_d01 and spatial_d02
        if (max_pair_world_distance_ratio > 0) {
          if (spatial_d02 < FLT_MAX) {
            if (RNIsZero(spatial_d02)) continue;
            if (spatial_d01 / spatial_d02 > max_pair_world_distance_ratio) continue;
          }
        }

        // Check relationships between descriptor_d01 and descriptor_d02
        if (max_pair_descriptor_distance_ratio > 0) {
          if (descriptor_d02 < FLT_MAX) {
            if (RNIsZero(descriptor_d02)) continue;
            if (descriptor_d01 / descriptor_d02 > max_pair_descriptor_distance_ratio) continue;
          }
        }

        ////////////////

        // Find distances from feature1 to closest other candidates
        RNScalar spatial_d12 = FLT_MAX;
        RNScalar descriptor_d12 = FLT_MAX;
        for (int k = 0; k < feature_candidates[index1].NEntries(); k++) {
          GSVFeature *feature2 = feature_candidates[index1].Kth(k);
          if (feature2 == feature0) continue;

          // Check if not in same region of same segment
          if (scanline0) {
            GSVScanline *scanline2 = feature2->scanline;
            if (scanline2) {
              GSVScan *scan2 = scanline2->Scan();
              GSVSegment *segment2 = scan2->Segment();
              if (segment0 != segment2) continue;
              if (travel_distance_cutoff > 0) {
                RNLength travel2 = scanline2->TravelDistance();
                if (fabs(travel0 - travel2) > travel_distance_cutoff) continue;
              }
            }
          }

          // Check if not in same image
          if (image0) {
            GSVImage *image2 = feature2->image;
            if (image2) {
              if (image0 != image2) continue;
            }
          }

          // Update closest spatial distance
          if (max_pair_world_distance_ratio > 0) {
            R3Point position2 = WorldPosition(feature2);
            RNLength spatial_d = R3SquaredDistance(position1, position2);
            if (spatial_d < spatial_d12) spatial_d12 = spatial_d; 
          }

          // Update closest descriptor distance
          if (max_pair_descriptor_distance_ratio > 0) {
            GSVDescriptor& descriptor2 = feature2->descriptor;
            RNLength descriptor_d = descriptor1.SquaredDistance(descriptor2);
            if (descriptor_d < descriptor_d12) descriptor_d12 = descriptor_d; 
          }
        }
          
        // Check relationships between spatial_d01 and spatial_d12
        if (max_pair_world_distance_ratio > 0) {
          if (spatial_d12 < FLT_MAX) {
            if (RNIsZero(spatial_d12)) continue;
            if (spatial_d01 / spatial_d12 > max_pair_world_distance_ratio) continue;
          }
        }

        // Check relationships between descriptor_d01 and descriptor_d12
        if (max_pair_descriptor_distance_ratio > 0) {
          if (descriptor_d12 < FLT_MAX) {
            if (RNIsZero(descriptor_d12)) continue;
            if (descriptor_d01 / descriptor_d12 > max_pair_descriptor_distance_ratio) continue;
          }
        }
      }

      // Create pair
      GSVFeaturePair *pair = new GSVFeaturePair(feature0, feature1, score);
      InsertPair(pair);
    }
  }

  ////////////////////

  // Delete temporary memory
  if (feature_masks) delete [] feature_masks;
  if (feature_candidates) delete [] feature_candidates;

  // Return success
  return 1;
}
  


int GSVPoseOptimization::
CreatePixelPixelPlanePairs(RNBoolean check_duplicates)
{
  // Parameters
  RNLength max_pixel_pixel_mesh_distance = 0.5;
  RNScalar max_pixel_pixel_mesh_curvature = 1.0;

  // Check number of pairs
  if (NPairs() == 0) return 1;

  // Create segment pair index
  RNArray<GSVFeaturePair *> *segment_pairs = new RNArray<GSVFeaturePair *>[ scene->NSegments() ];
  for (int i = 0; i < NPairs(); i++) {
    GSVFeaturePair *pair = Pair(i);
    GSVFeature *feature0 = pair->Feature(0);
    if (feature0->feature_type != GSV_IMAGE_POINT_FEATURE_TYPE) continue;
    if (feature0->image_t <= 0.0) continue;
    GSVFeature *feature1 = pair->Feature(1);
    if (feature1->feature_type != GSV_IMAGE_POINT_FEATURE_TYPE) continue;
    if (feature1->image_t <= 0.0) continue;
    GSVImage *image0 = feature0->image;
    if (!image0) continue;
    GSVImage *image1 = feature1->image;
    if (!image1) continue;
    GSVSegment *segment0 = image0->Segment();
    segment_pairs[segment0->SceneIndex()].Insert(pair);
    GSVSegment *segment1 = image1->Segment();
    if (segment1 == segment0) continue;
    segment_pairs[segment1->SceneIndex()].Insert(pair);
  }

  // Create pixel-pixel-mesh pairs segment by segment
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      int segment_index = segment->SceneIndex();
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (scan->NScanlines() == 0) continue;
        if (ia == 0) continue;
        
        // Read the mesh
        GSVMesh *mesh = scan->Mesh();
        if (!mesh) continue;

        // Check each pair in segment
        for (int ic = 0; ic < segment_pairs[segment_index].NEntries(); ic++) {
          GSVFeaturePair *pair = segment_pairs[segment_index].Kth(ic);
          for (int k = 0; k < 2; k++) {
            GSVFeature *feature = pair->Feature(k);
            if (feature->image->Segment() != segment) continue;

            // Check if close to mesh
            R3MeshIntersection closest;
            RNScalar max_distance = max_pixel_pixel_mesh_distance;
            if (max_distance == RN_UNKNOWN) max_distance = FLT_MAX;
            if (mesh->ClosestPoint(feature->scan_position, &closest, 0, max_distance)) {
              GSVMeshVertex *vertex = NULL;
              if (closest.type == R3_MESH_VERTEX_TYPE) vertex = (GSVMeshVertex *) closest.vertex;
              else if (closest.type == R3_MESH_EDGE_TYPE) vertex = (GSVMeshVertex *) mesh->VertexOnEdge(closest.edge, 0);
              else if (closest.type == R3_MESH_FACE_TYPE) vertex = (GSVMeshVertex *) mesh->VertexOnFace(closest.face, 0);
              if (!mesh->IsVertexOnBoundary(vertex)) {
                if ((max_pixel_pixel_mesh_curvature == RN_UNKNOWN) || (fabs(mesh->VertexMeanCurvature(vertex)) < max_pixel_pixel_mesh_curvature)) {
                  GSVScanline *scanline = mesh->VertexScanline(vertex);
                  int point_index = mesh->VertexPointIndex(vertex);
                  R3Point scan_position = closest.point;
                  R3Vector scan_normal = mesh->FaceNormal(closest.face);
                  GSVFeature *feature2 = new GSVFeature(GSV_SCAN_PLANE_FEATURE_TYPE, scanline, point_index, 
                    scan_position, R3zero_vector, scan_normal, 1, 1);
                  InsertFeature(feature2);
                  GSVFeaturePair *pair2 = new GSVFeaturePair(feature, feature2, 1);
                  InsertPair(pair2);
                }
              }
            }
          }
        }

        // Delete mesh
        delete mesh;
      }
    }
  }

  // Delete segment pair index
  delete [] segment_pairs;

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Correspondence creation
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
CreateAllCorrespondences(RNBoolean check_duplicates)
{
  // Create masks of features not to be paired with each feature
  RNArray<GSVFeature *> *feature_masks = NULL;
  if (check_duplicates) feature_masks = CreateFeatureAssociations(this, FALSE, TRUE, TRUE);

  // Consider every correspondence of features
  for (int i = 0; i < NPairs(); i++) {
    GSVFeaturePair *pair = Pair(i);
    GSVFeature *feature0 = pair->Feature(0);
    GSVFeature *feature1 = pair->Feature(1);
    if (feature_masks) {
      if (feature_masks[feature0->index].FindEntry(feature1)) continue;
      assert(!feature_masks[feature1->index].FindEntry(feature0));
      feature_masks[feature0->index].Insert(feature1);
      feature_masks[feature1->index].Insert(feature0);
    }
    GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(feature0, feature1, pair->score);
    InsertCorrespondence(correspondence);
  }

  // Delete feature correspondence index
  if (feature_masks) delete [] feature_masks;

  // Return success
  return 1;
}




int GSVPoseOptimization::
CreateFeasibleCorrespondences(RNBoolean check_duplicates)
{
  // Create masks of features not to be paired with each feature
  RNArray<GSVFeature *> *feature_masks = NULL;
  if (check_duplicates) feature_masks = CreateFeatureAssociations(this, FALSE, TRUE, TRUE);

  // Consider every pair
  for (int i = 0; i < NPairs(); i++) {
    GSVFeaturePair *pair = Pair(i);
    GSVFeature *feature0 = pair->Feature(0);
    GSVFeature *feature1 = pair->Feature(1);
    if (feature_masks) {
      if (feature_masks[feature0->index].FindEntry(feature1)) continue;
      assert(!feature_masks[feature1->index].FindEntry(feature0));
      feature_masks[feature0->index].Insert(feature1);
      feature_masks[feature1->index].Insert(feature0);
    }
    RNScalar score = CorrespondenceScore(feature0, feature1);
    if (score <= 0) continue;
    GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(feature0, feature1, score);
    InsertCorrespondence(correspondence);
  }

  // Delete feature correspondence index
  if (feature_masks) delete [] feature_masks;

  // Return success
  return 1;
}



int GSVPoseOptimization::
CreateCloseCorrespondences(RNBoolean check_duplicates)
{
  // Check pairs
  if (NFeatures() == 0) return 1;
  if (NPairs() == 0) return 1;
  
  // Set parameters
  // mutually closest and ratio only consider other scan features within this cutoff
  RNScalar travel_distance_cutoff = 100; 

  ////////////////////

  // Create masks of features not to be paired with each feature
  RNArray<GSVFeature *> *feature_masks = NULL;
  if (check_duplicates) feature_masks = CreateFeatureAssociations(this, FALSE, TRUE, TRUE);

  // Create array of candidates for every feature
  RNArray<GSVFeature *> *feature_candidates = new RNArray<GSVFeature *>[ NFeatures() ];
  for (int i = 0; i < NPairs(); i++) {
    GSVFeaturePair *pair = Pair(i);
    GSVFeature *feature0 = pair->Feature(0);
    GSVFeature *feature1 = pair->Feature(1);
    if (feature_masks) {
      if (feature_masks[feature0->index].FindEntry(feature1)) continue;
      assert(!feature_masks[feature1->index].FindEntry(feature0));
      feature_masks[feature0->index].Insert(feature1);
      feature_masks[feature1->index].Insert(feature0);
    }
    feature_candidates[feature0->index].Insert(feature1);
    feature_candidates[feature1->index].Insert(feature0);
  }

  ////////////////////

  // Create correspondences from candidates
  for (int index0 = 0; index0 < NFeatures(); index0++) {
    GSVFeature *feature0 = Feature(index0);
    if (feature_candidates[index0].IsEmpty()) continue;
    R3Point position0 = WorldPosition(feature0);
    if (position0 == R3unknown_point) continue;
    GSVScanline *scanline0 = feature0->scanline;
    GSVScan *scan0 = (scanline0) ? scanline0->Scan() : NULL;
    GSVSegment *segment0 = (scan0) ? scan0->Segment() : NULL;
    RNLength travel0 = (scanline0) ? scanline0->TravelDistance() : 0.0;
    GSVDescriptor& descriptor0 = feature0->descriptor;
    GSVImage *image0 = feature0->image;

    // Consider all candidates for feature0
    for (int j = 0; j < feature_candidates[index0].NEntries(); j++) {
      GSVFeature *feature1 = feature_candidates[index0].Kth(j);
      int index1 = feature1->index;
      if (index1 <= index0) continue;

      // Compute score
      RNScalar score = CorrespondenceScore(feature0, feature1);
      if (score < 0) continue;

      // Check relationships to other features
      if ((max_correspondence_world_distance_ratio > 0) || 
          (max_correspondence_descriptor_distance_ratio > 0)) {
        // Get stuff for feature1
        GSVScanline *scanline1 = feature1->scanline;
        GSVScan *scan1 = (scanline1) ? scanline1->Scan() : NULL;
        GSVSegment *segment1 = (scan1) ? scan1->Segment() : NULL;
        RNLength travel1 = (scanline1) ? scanline1->TravelDistance() : 0.0;
        GSVImage *image1 = feature1->image;
        R3Point position1 = WorldPosition(feature1);
        assert(position1 != R3unknown_point);
        GSVDescriptor& descriptor1 = feature1->descriptor;

        // Get distances between feature0 and feature1
        RNLength spatial_d01 = (max_correspondence_world_distance_ratio > 0.0) ? R3SquaredDistance(position0, position1) : 0.0;
        RNLength descriptor_d01 = (max_correspondence_descriptor_distance_ratio > 0.0) ? descriptor0.SquaredDistance(descriptor1) : 0.0;

        //////// 

        // Find distances from feature0 to closest other candidates
        RNScalar spatial_d02 = FLT_MAX;
        RNScalar descriptor_d02 = FLT_MAX;
        for (int k = 0; k < feature_candidates[index0].NEntries(); k++) {
          GSVFeature *feature2 = feature_candidates[index0].Kth(k);
          if (feature2 == feature1) continue;

          // Check if not in same region of same segment
          if (scanline1) {
            GSVScanline *scanline2 = feature2->scanline;
            if (scanline2) {
              GSVScan *scan2 = scanline2->Scan();
              GSVSegment *segment2 = scan2->Segment();
              if (segment1 != segment2) continue;
              if (travel_distance_cutoff > 0) {
                RNLength travel2 = scanline2->TravelDistance();
                if (fabs(travel1 - travel2) > travel_distance_cutoff) continue;
              }
            }
          }

          // Check if not in same image
          if (image1) {
            GSVImage *image2 = feature2->image;
            if (image2) {
              if (image1 != image2) continue;
            }
          }

          // Update closest spatial distance
          if (max_correspondence_world_distance_ratio > 0) {
            R3Point position2 = WorldPosition(feature2);
            RNLength spatial_d = R3SquaredDistance(position0, position2);
            if (spatial_d < spatial_d02) spatial_d02 = spatial_d; 
          }

          // Update closest descriptor distance
          if (max_correspondence_descriptor_distance_ratio > 0) {
            GSVDescriptor& descriptor2 = feature2->descriptor;
            RNLength descriptor_d = descriptor0.SquaredDistance(descriptor2);
            if (descriptor_d < descriptor_d02) descriptor_d02 = descriptor_d; 
          }
        }
          
        // Check relationships between spatial_d01 and spatial_d02
        if (max_correspondence_world_distance_ratio > 0) {
          if (spatial_d02 < FLT_MAX) {
            if (RNIsZero(spatial_d02)) continue;
            if (spatial_d01 / spatial_d02 > max_correspondence_world_distance_ratio) continue;
          }
        }

        // Check relationships between descriptor_d01 and descriptor_d02
        if (max_correspondence_descriptor_distance_ratio > 0) {
          if (descriptor_d02 < FLT_MAX) {
            if (RNIsZero(descriptor_d02)) continue;
            if (descriptor_d01 / descriptor_d02 > max_correspondence_descriptor_distance_ratio) continue;
          }
        }

        ////////////////

        // Find distances from feature1 to closest other candidates
        RNScalar spatial_d12 = FLT_MAX;
        RNScalar descriptor_d12 = FLT_MAX;
        for (int k = 0; k < feature_candidates[index1].NEntries(); k++) {
          GSVFeature *feature2 = feature_candidates[index1].Kth(k);
          if (feature2 == feature0) continue;

          // Check if not in same region of same segment
          if (scanline0) {
            GSVScanline *scanline2 = feature2->scanline;
            if (scanline2) {
              GSVScan *scan2 = scanline2->Scan();
              GSVSegment *segment2 = scan2->Segment();
              if (segment0 != segment2) continue;
              if (travel_distance_cutoff > 0) {
                RNLength travel2 = scanline2->TravelDistance();
                if (fabs(travel0 - travel2) > travel_distance_cutoff) continue;
              }
            }
          }

          // Check if not in same image
          if (image0) {
            GSVImage *image2 = feature2->image;
            if (image2) {
              if (image0 != image2) continue;
            }
          }

          // Update closest spatial distance
          if (max_correspondence_world_distance_ratio > 0) {
            R3Point position2 = WorldPosition(feature2);
            RNLength spatial_d = R3SquaredDistance(position1, position2);
            if (spatial_d < spatial_d12) spatial_d12 = spatial_d; 
          }

          // Update closest descriptor distance
          if (max_correspondence_descriptor_distance_ratio > 0) {
            GSVDescriptor& descriptor2 = feature2->descriptor;
            RNLength descriptor_d = descriptor1.SquaredDistance(descriptor2);
            if (descriptor_d < descriptor_d12) descriptor_d12 = descriptor_d; 
          }
        }
          
        // Check relationships between spatial_d01 and spatial_d12
        if (max_correspondence_world_distance_ratio > 0) {
          if (spatial_d12 < FLT_MAX) {
            if (RNIsZero(spatial_d12)) continue;
            if (spatial_d01 / spatial_d12 > max_correspondence_world_distance_ratio) continue;
          }
        }

        // Check relationships between descriptor_d01 and descriptor_d12
        if (max_correspondence_descriptor_distance_ratio > 0) {
          if (descriptor_d12 < FLT_MAX) {
            if (RNIsZero(descriptor_d12)) continue;
            if (descriptor_d01 / descriptor_d12 > max_correspondence_descriptor_distance_ratio) continue;
          }
        }
      }

      // Create correspondence
      GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(feature0, feature1, score);
      InsertCorrespondence(correspondence);
    }
  }

  ////////////////////

  // Delete temporary memory
  if (feature_masks) delete [] feature_masks;
  if (feature_candidates) delete [] feature_candidates;

  ////////////////////

  // Return success
  return 1;
}
  



int GSVPoseOptimization::
CreateCoplanarCorrespondences(RNBoolean check_duplicates)
{
  // Parameters
  RNLength min_coplanar_correspondence_distance = 15;
  RNLength max_coplanar_correspondence_distance = 100;
  RNLength min_coplanar_correspondence_squared_distance = min_coplanar_correspondence_distance * min_coplanar_correspondence_distance;
  RNLength max_coplanar_correspondence_squared_distance = max_coplanar_correspondence_distance * max_coplanar_correspondence_distance;

  // Check features
  if (NFeatures() == 0) return 1;
  
  // Create masks of features not to be paired with each feature
  RNArray<GSVFeature *> *feature_masks = NULL;
  if (check_duplicates) feature_masks = CreateFeatureAssociations(this, FALSE, TRUE, TRUE);

  // Create correspondences between coplanar features
  for (int i = 0; i < NPairs(); i++) {
    GSVFeaturePair *pair = Pair(i);
    GSVFeature *feature0 = pair->Feature(0);
    GSVFeature *feature1 = pair->Feature(1);
    if (feature0->feature_type != GSV_SCAN_PLANE_FEATURE_TYPE) continue;
    if (feature1->feature_type != GSV_SCAN_PLANE_FEATURE_TYPE) continue;
    if (feature0->group == 0) continue;
    if (feature1->group == 0) continue;
    if (feature0->group != feature1->group) continue;
    R3Point position0 = feature0->scan_position;
    R3Point position1 = feature1->scan_position;
    RNLength squared_distance = R3SquaredDistance(position0, position1);
    if (squared_distance < min_coplanar_correspondence_squared_distance) continue;
    if (squared_distance > max_coplanar_correspondence_squared_distance) continue;
    if (feature_masks) {
      if (feature_masks[feature0->index].FindEntry(feature1)) continue;
      assert(!feature_masks[feature1->index].FindEntry(feature0));
      feature_masks[feature0->index].Insert(feature1);
      feature_masks[feature1->index].Insert(feature0);
    }
    GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(feature0, feature1, pair->score);
    InsertCorrespondence(correspondence);
  }

  // Delete temporary data
  if (feature_masks) delete [] feature_masks;
 
  // Return success
  return 1;
}



int GSVPoseOptimization::
CreatePixelPixelPlaneCorrespondences(RNBoolean check_duplicates)
{
  // Parameters
  RNLength max_pixel_pixel_mesh_distance = 0.1;
  RNScalar max_pixel_pixel_mesh_curvature = 0.1;

  // Check number of correspondences
  if (NCorrespondences() == 0) return 1;

  // Create segment correspondence index
  RNArray<GSVFeatureCorrespondence *> *segment_correspondences = new RNArray<GSVFeatureCorrespondence *>[ scene->NSegments() ];
  for (int i = 0; i < NCorrespondences(); i++) {
    GSVFeatureCorrespondence *correspondence = Correspondence(i);
    GSVFeature *feature0 = correspondence->Feature(0);
    if (feature0->feature_type != GSV_IMAGE_POINT_FEATURE_TYPE) continue;
    if (feature0->image_t <= 0.0) continue;
    GSVFeature *feature1 = correspondence->Feature(1);
    if (feature1->feature_type != GSV_IMAGE_POINT_FEATURE_TYPE) continue;
    if (feature1->image_t <= 0.0) continue;
    GSVImage *image0 = feature0->image;
    if (!image0) continue;
    GSVImage *image1 = feature1->image;
    if (!image1) continue;
    GSVSegment *segment0 = image0->Segment();
    segment_correspondences[segment0->SceneIndex()].Insert(correspondence);
    GSVSegment *segment1 = image1->Segment();
    if (segment1 == segment0) continue;
    segment_correspondences[segment1->SceneIndex()].Insert(correspondence);
  }

  // Create pixel-pixel-mesh correspondences segment by segment
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      int segment_index = segment->SceneIndex();
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (scan->NScanlines() == 0) continue;
        if (ia == 0) continue;
        
        // Read the mesh
        GSVMesh *mesh = scan->Mesh();
        if (!mesh) continue;

        // Check each correspondence in segment
        for (int ic = 0; ic < segment_correspondences[segment_index].NEntries(); ic++) {
          GSVFeatureCorrespondence *correspondence = segment_correspondences[segment_index].Kth(ic);
          for (int k = 0; k < 2; k++) {
            GSVFeature *feature = correspondence->Feature(k);
            if (feature->image->Segment() != segment) continue;

            // Check if close to mesh
            R3MeshIntersection closest;
            RNScalar max_distance = max_pixel_pixel_mesh_distance;
            if (max_distance == RN_UNKNOWN) max_distance = FLT_MAX;
            if (mesh->ClosestPoint(feature->scan_position, &closest, 0, max_distance)) {
              GSVMeshVertex *vertex = NULL;
              if (closest.type == R3_MESH_VERTEX_TYPE) vertex = (GSVMeshVertex *) closest.vertex;
              else if (closest.type == R3_MESH_EDGE_TYPE) vertex = (GSVMeshVertex *) mesh->VertexOnEdge(closest.edge, 0);
              else if (closest.type == R3_MESH_FACE_TYPE) vertex = (GSVMeshVertex *) mesh->VertexOnFace(closest.face, 0);
              if (!mesh->IsVertexOnBoundary(vertex)) {
                if ((max_pixel_pixel_mesh_curvature == RN_UNKNOWN) || (fabs(mesh->VertexMeanCurvature(vertex)) < max_pixel_pixel_mesh_curvature)) {
                  GSVScanline *scanline = mesh->VertexScanline(vertex);
                  int point_index = mesh->VertexPointIndex(vertex);
                  R3Point scan_position = closest.point;
                  R3Vector scan_normal = mesh->FaceNormal(closest.face);
                  GSVFeature *feature2 = new GSVFeature(GSV_SCAN_PLANE_FEATURE_TYPE, scanline, point_index, 
                    scan_position, R3zero_vector, scan_normal, 1, 1);
                  InsertFeature(feature2);
                  GSVFeatureCorrespondence *correspondence2 = new GSVFeatureCorrespondence(feature, feature2, 1);
                  InsertCorrespondence(correspondence2);
                }
              }
            }
          }
        }

        // Delete mesh
        delete mesh;
      }
    }
  }

  // Delete segment correspondence index
  delete [] segment_correspondences;

  // Return success
  return 1;
}



int GSVPoseOptimization::
CreateTransitiveClosureCorrespondences(RNBoolean check_duplicates)
{
  // Check 
  if (NFeatures() == 0) return 1;
  if (NCorrespondences() == 0) return 1;

  // Create list of feature correspondences
  RNArray<GSVFeatureCorrespondence *> *feature_correspondences = NULL;
  feature_correspondences = new RNArray<GSVFeatureCorrespondence *>[ NFeatures() ];
  for (int i = 0; i < NCorrespondences(); i++) {
    GSVFeatureCorrespondence *correspondence = Correspondence(i);
    GSVFeature *feature0 = correspondence->Feature(0);
    GSVFeature *feature1 = correspondence->Feature(1);
    feature_correspondences[feature0->index].Insert(correspondence);
    feature_correspondences[feature1->index].Insert(correspondence);
  }

  // Create transitive correspondences
  const int max_iterations = 64;
  for (int iteration = 0; iteration < max_iterations; iteration++) {
    RNBoolean done = TRUE;
    for (int i = 0; i < NFeatures(); i++) {
      GSVFeature *feature0 = Feature(i);
      for (int j = 0; j < feature_correspondences[i].NEntries(); j++) {
        GSVFeatureCorrespondence *correspondence1 = feature_correspondences[i].Kth(j); 
        GSVFeature *feature1 = correspondence1->Feature(0);
        if (feature1 == feature0) feature1 = correspondence1->Feature(1);
        for (int k = j+1; k < feature_correspondences[i].NEntries(); k++) {
          GSVFeatureCorrespondence *correspondence2 = feature_correspondences[i].Kth(k); 
          GSVFeature *feature2 = correspondence2->Feature(0);
          if (feature2 == feature0) feature2 = correspondence2->Feature(1);
          if (feature1 == feature2) continue;
          if (feature1->scanline && feature2->scanline && (feature1->scanline == feature2->scanline)) continue;
          if (feature1->image && feature2->image && (feature1->image == feature2->image)) continue;
          
          // Check for existing correspondence between feature1 and feature2
          if (check_duplicates) {
            RNBoolean found = FALSE;
            for (int s = 0; s < feature_correspondences[feature1->index].NEntries(); s++) {
              GSVFeatureCorrespondence *c = feature_correspondences[feature1->index].Kth(s); 
              if ((c->Feature(0) == feature1) && (c->Feature(1) == feature2)) { found = TRUE; break; }
              if ((c->Feature(1) == feature1) && (c->Feature(0) == feature2)) { found = TRUE; break; }
            }
            if (found) continue;
          }

          // Create new correspondence between feature1 and feature2
          RNScalar score = correspondence1->score * correspondence2->score;
          GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(feature1, feature2, score);
          feature_correspondences[feature1->index].Insert(correspondence);
          feature_correspondences[feature2->index].Insert(correspondence);
          InsertCorrespondence(correspondence);
          done = FALSE;
        }
      }
    }
    if (done) break;
  }

  // Delete arrays of feature correspondences
  if (feature_correspondences) delete [] feature_correspondences;
 
  // Return success
  return 1;
}



int GSVPoseOptimization::
CreateICPCorrespondences(RNBoolean check_duplicates)
{
  // Just checking
  if (icp_max_world_distance_end <= 0) icp_max_world_distance_end = icp_max_world_distance_start;
  if (icp_max_world_distance_end >= icp_max_world_distance_start) icp_max_world_distance_end = icp_max_world_distance_start;
  if (icp_max_image_distance_end <= 0) icp_max_image_distance_end = icp_max_image_distance_start;
  if (icp_max_image_distance_end >= icp_max_image_distance_start) icp_max_image_distance_end = icp_max_image_distance_start;

  // Initialize parameters
  const int max_iterations = 32;
  int saved_correspondences = NCorrespondences();
  RNScalar saved_max_correspondence_world_distance = max_correspondence_world_distance;
  RNScalar saved_max_correspondence_image_distance = max_correspondence_image_distance;
  RNScalar max_world_distance = icp_max_world_distance_start;
  RNScalar max_image_distance = icp_max_image_distance_start;

  // Iteratively create correspondences
  for (int i = 0; i < max_iterations; i++) {
    // Create correspondences
    TruncateCorrespondences(saved_correspondences);
    max_correspondence_world_distance = max_world_distance;
    max_correspondence_image_distance = max_image_distance;
    CreateCloseCorrespondences(check_duplicates);
  
    // Compute rigidity 
    // OPTIMIZATION MAKES NO CHANGE IF RIGIDITY IS TOO HIGH
    RNScalar rigidity = 0.1 * max_world_distance;

    // Solve for translations only
    Solve(FALSE, FALSE,  FALSE, FALSE,  TRUE, FALSE,  TRUE, FALSE,  TRUE,  TRUE, TRUE,  rigidity);

    // Solve for translations and rotations
    if (max_world_distance < 2) Solve(FALSE, FALSE,  FALSE, FALSE,  TRUE, TRUE,  TRUE, TRUE,  TRUE,  TRUE, TRUE,  rigidity);

    // Compute score
    printf("%d : %d %d %d %d : %g %g : %g %g\n", i+1, NFeatures(), NPairs(), NCorrespondences(), NClusters(),
      max_world_distance, max_image_distance, WorldRMSD(), ImageRMSD());

    // Check for termination
    if (((max_world_distance == 0) || (max_world_distance == icp_max_world_distance_end)) && 
        ((max_image_distance == 0) || (max_image_distance == icp_max_image_distance_end))) {
      break;
    }

    // Update max distance
    max_world_distance *= 0.5;
    max_image_distance *= 0.5;
    if (max_world_distance < icp_max_world_distance_end) max_world_distance = icp_max_world_distance_end;
    if (max_image_distance < icp_max_image_distance_end) max_image_distance = icp_max_image_distance_end;
  }    

  // Restore parameters
  max_correspondence_world_distance = saved_max_correspondence_world_distance;
  max_correspondence_image_distance = saved_max_correspondence_image_distance;

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Cluster creation
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
CreateAllClusters(RNBoolean check_duplicates)
{
  // Create all clusters from pairs
  int saved_ncorrespondences = NCorrespondences();
  if (!CreateAllCorrespondences(check_duplicates)) return 0;
  // if (!CreateCorrespondenceClusters(check_duplicates)) return 0;
  int status = CreateHierarchicalClusters(check_duplicates);
  TruncateCorrespondences(saved_ncorrespondences);

  // Return status
  return status;
}



int GSVPoseOptimization::
CreateCloseClusters(RNBoolean check_duplicates)
{
  // Create close clusters from pairs
  int saved_ncorrespondences = NCorrespondences();
  if (!CreateCloseCorrespondences(check_duplicates)) return 0;
  // if (!CreateCorrespondenceClusters(check_duplicates)) return 0;
  int status = CreateHierarchicalClusters(check_duplicates);
  TruncateCorrespondences(saved_ncorrespondences);

  // Return status
  return status;
}



int GSVPoseOptimization::
CreateCorrespondenceClusters(RNBoolean check_duplicates)
{
  // Create masks of features not to be paired with each feature
  RNArray<GSVFeature *> *feature_masks = NULL;
  if (check_duplicates) feature_masks = CreateFeatureAssociations(this, FALSE, FALSE, TRUE);

  // Consider every correspondence
  for (int i = 0; i < NCorrespondences(); i++) {
    GSVFeatureCorrespondence *correspondence = Correspondence(i);
    GSVFeature *feature0 = correspondence->Feature(0);
    GSVFeature *feature1 = correspondence->Feature(1);

    // Check if duplicate
    if (feature_masks) {
      if (feature_masks[feature0->index].FindEntry(feature1)) continue;
      assert(!feature_masks[feature1->index].FindEntry(feature0));
      feature_masks[feature0->index].Insert(feature1);
      feature_masks[feature1->index].Insert(feature0);
    }

    // Create cluster
    GSVFeatureCluster *cluster = new GSVFeatureCluster();
    cluster->InsertFeature(feature0);
    cluster->InsertFeature(feature1);
    cluster->Update(this);
    InsertCluster(cluster);
  }

  // Delete feature mask
  if (feature_masks) delete [] feature_masks;

  // Return success
  return 1;
}



int GSVPoseOptimization::
CreateICPClusters(RNBoolean check_duplicates)
{
  // Just checking
  if (icp_max_world_distance_end <= 0) icp_max_world_distance_end = icp_max_world_distance_start;
  if (icp_max_world_distance_end >= icp_max_world_distance_start) icp_max_world_distance_end = icp_max_world_distance_start;
  if (icp_max_image_distance_end <= 0) icp_max_image_distance_end = icp_max_image_distance_start;
  if (icp_max_image_distance_end >= icp_max_image_distance_start) icp_max_image_distance_end = icp_max_image_distance_start;

  // Initialize parameters
  const int max_iterations = 32;
  int saved_clusters = NClusters();
  RNScalar saved_max_correspondence_world_distance = max_correspondence_world_distance;
  RNScalar saved_max_correspondence_image_distance = max_correspondence_image_distance;
  RNScalar max_world_distance = icp_max_world_distance_start;
  RNScalar max_image_distance = icp_max_image_distance_start;

  // Iteratively create clusters
  for (int i = 0; i < max_iterations; i++) {
    // Create clusters
    TruncateClusters(saved_clusters);
    max_correspondence_world_distance = max_world_distance;
    max_correspondence_image_distance = max_image_distance;
    CreateCloseClusters(check_duplicates);
  
    // Compute rigidity 
    // OPTIMIZATION MAKES NO CHANGE IF RIGIDITY IS TOO HIGH
    RNScalar rigidity = 0.1 * max_world_distance;

    // Solve for translations only
    Solve(FALSE, FALSE,  FALSE, FALSE,  TRUE, FALSE,  TRUE, FALSE,  TRUE,  TRUE, TRUE,  rigidity);

    // Solve for translations and rotations
    if (max_world_distance < 2) Solve(FALSE, FALSE,  FALSE, FALSE,  TRUE, TRUE,  TRUE, TRUE,  TRUE,  TRUE, TRUE,  rigidity);

    // Compute score
    printf("%d : %d %d %d %d : %g %g : %g %g\n", i+1, NFeatures(), NPairs(), NCorrespondences(), NClusters(),
      max_world_distance, max_image_distance, WorldRMSD(), ImageRMSD());

    // Check for termination
    if (((max_world_distance == 0) || (max_world_distance == icp_max_world_distance_end)) && 
        ((max_image_distance == 0) || (max_image_distance == icp_max_image_distance_end))) {
      break;
    }

    // Update max distance
    max_world_distance *= 0.5;
    max_image_distance *= 0.5;
    if (max_world_distance < icp_max_world_distance_end) max_world_distance = icp_max_world_distance_end;
    if (max_image_distance < icp_max_image_distance_end) max_image_distance = icp_max_image_distance_end;
  }    

  // Restore parameters
  max_correspondence_world_distance = saved_max_correspondence_world_distance;
  max_correspondence_image_distance = saved_max_correspondence_image_distance;

  // Return success
  return 1;
}



int GSVPoseOptimization::
CreateHierarchicalClusters(RNBoolean check_duplicates)
{
  // Create hierarchical clusters from small clusters and correspondences

  // Check features
  if (NFeatures() == 0) return 1;
  if (NCorrespondences() == 0) return 1;
  if (NClusters() == 0) return 1;

  // Allocate temporary data
  RNArray<GSVFeatureCluster *> *cluster_neighbors = new RNArray<GSVFeatureCluster *>[ NClusters() ];
  GSVFeatureCluster **cluster_parents = new GSVFeatureCluster *[ NClusters() ];
  GSVFeatureCluster **feature_parents = new GSVFeatureCluster *[ NFeatures() ];
  for (int i = 0; i < NClusters(); i++) cluster_parents[i] = Cluster(i);
  for (int i = 0; i < NFeatures(); i++) feature_parents[i] = NULL;

  // Log cluster parents
  for (int i = 0; i < NClusters(); i++) {
    GSVFeatureCluster *cluster0 = Cluster(i);
    for (int j = 0; j < cluster0->NClusters(); j++) {
      GSVFeatureCluster *cluster1 = cluster0->Cluster(j);
      cluster_parents[cluster1->index] = cluster0;
    }
  }

  // Log feature parents
  for (int i = 0; i < NClusters(); i++) {
    GSVFeatureCluster *cluster0 = Cluster(i);
    for (int j = 0; j < cluster0->NFeatures(); j++) {
      GSVFeature *feature1 = cluster0->Feature(j);
      GSVFeatureCluster *parent = cluster0;
      while (cluster_parents[parent->index] != parent) 
        parent = cluster_parents[parent->index];
      feature_parents[feature1->index] = parent;
    }
  }

#if 0
  // Log cluster neighbors
  for (int i = 0; i < NClusters(); i++) {
    GSVFeatureCluster *cluster0 = Cluster(i);
    GSVFeatureCluster *parent0 = cluster_parents[cluster0->index];
    assert(cluster_parents[parent0->index] == parent0);
    for (int j = 0; j < cluster0->NFeatures(); j++) {
      GSVFeature *feature1 = cluster0->Feature(j);
      GSVFeatureCluster *parent1 = feature_parents[feature1->index];
      assert(cluster_parents[parent1->index] == parent1);
      if (parent0 == parent1) continue;
      if (cluster_neighbors[parent0->index].FindEntry(parent1)) continue;
      assert(!cluster_neighbors[parent1->index].FindEntry(parent0));
      cluster_neighbors[parent0->index].Insert(parent1);
      cluster_neighbors[parent1->index].Insert(parent0);
    }
    for (int j = 0; j < cluster0->NClusters(); j++) {
      GSVFeatureCluster *cluster1 = cluster0->Cluster(j);
      GSVFeatureCluster *parent1 = cluster_parents[cluster1->index];
      assert(cluster_parents[parent1->index] == parent1);
      if (parent0 == parent1) continue;
      if (cluster_neighbors[parent0->index].FindEntry(parent1)) continue;
      assert(!cluster_neighbors[parent1->index].FindEntry(parent0));
      cluster_neighbors[parent0->index].Insert(parent1);
      cluster_neighbors[parent1->index].Insert(parent0);
    }
  }

#else
  // Log cluster neighbors
  for (int i = 0; i < NCorrespondences(); i++) {
    GSVFeatureCorrespondence *correspondence = Correspondence(i);
    GSVFeature *feature0 = correspondence->Feature(0);
    GSVFeature *feature1 = correspondence->Feature(1);
    GSVFeatureCluster *cluster0 = feature_parents[feature0->index];
    GSVFeatureCluster *cluster1 = feature_parents[feature1->index];
    if (!cluster0) continue;
    if (!cluster1) continue;
    if (cluster0 == cluster1) continue;
    cluster_neighbors[cluster0->index].Insert(cluster1);
    cluster_neighbors[cluster1->index].Insert(cluster0);
  }
#endif

  // Create cluster for each connected component of neighbor graph
  int seed_index = 0;
  RNArray<GSVFeatureCluster *> created_clusters;
  while (seed_index < NClusters()) {
    GSVFeatureCluster *seed_cluster = Cluster(seed_index++);
    if (cluster_parents[seed_cluster->index] == NULL) continue;
    if (cluster_parents[seed_cluster->index] != seed_cluster) continue;
    if (cluster_neighbors[seed_cluster->index].IsEmpty()) continue;

    // Find connected component
    RNArray<GSVFeatureCluster *> connected_component;
    RNArray<GSVFeatureCluster *> stack;
    stack.Insert(seed_cluster);
    cluster_parents[seed_cluster->index] = NULL;
    while (!stack.IsEmpty()) {
      GSVFeatureCluster *cluster = stack.Tail();
      stack.RemoveTail();
      connected_component.Insert(cluster);
      for (int i = 0; i < cluster_neighbors[cluster->index].NEntries(); i++) {
        GSVFeatureCluster *neighbor_cluster = cluster_neighbors[cluster->index].Kth(i);
        if (cluster_parents[neighbor_cluster->index] == NULL) continue;
        if (cluster_parents[neighbor_cluster->index] != neighbor_cluster) continue;
        cluster_parents[neighbor_cluster->index] = NULL;
        stack.Insert(neighbor_cluster);
      }
    }

    // Create hierarchical cluster for connected component
    if (connected_component.NEntries() > 1) {
      // Create cluster
      GSVFeatureCluster *created_cluster = new GSVFeatureCluster();
      created_clusters.Insert(created_cluster);

      // Insert connected component as children of created cluster
      for (int i = 0; i < connected_component.NEntries(); i++) {
        GSVFeatureCluster *cluster = connected_component.Kth(i);
        created_cluster->InsertCluster(cluster);
      }
    }
  }

  // Insert clusters
  for (int i = 0; i < created_clusters.NEntries(); i++) {
    GSVFeatureCluster *cluster = created_clusters.Kth(i);
    cluster->Update(this);
    InsertCluster(cluster);
  }

  // Delete temporary data
  if (feature_parents) delete [] feature_parents;
  if (cluster_parents) delete [] cluster_parents;
  if (cluster_neighbors) delete [] cluster_neighbors;

  // Return success
  return 1;
}



int GSVPoseOptimization::
CreateMergedClusters(RNBoolean check_duplicates)
{
  // Create big clusters from small overlapping clusters and correspondences

  // Check features
  if (NFeatures() == 0) return 1;

  // Create arrays of features associated with each feature
  RNArray<GSVFeature *> *feature_associations = NULL;
  feature_associations = CreateFeatureAssociations(this, FALSE, TRUE, TRUE);

  // Remove existing clusters ???
  TruncateClusters(0);

  // Build merged clusters
  for (int i = 0; i < NFeatures(); i++) {
    GSVFeature *feature0 = Feature(i);
    if (!feature_associations[feature0->index].IsEmpty()) {
      // Create cluster
      GSVFeatureCluster *cluster = new GSVFeatureCluster();

      // Insert features
      cluster->InsertFeature(feature0);
      for (int j = 0; j < feature_associations[feature0->index].NEntries(); j++) {
        GSVFeature *feature1 = feature_associations[feature0->index].Kth(j);
        feature_associations[feature1->index].Empty();
        cluster->InsertFeature(feature1);
      }

      // Update cluster
      cluster->feature_type = GSV_NULL_FEATURE_TYPE;
      cluster->Update(this);
      InsertCluster(cluster);
    }
  }

  // Delete temporary memory
  if (feature_associations) delete [] feature_associations;

  // Return success
  return 1;
}



int GSVPoseOptimization::
CreateCoplanarClusters(RNBoolean check_duplicates)
{
  // Parameters
  const RNScalar min_feature_spacing = 5;
  const RNScalar min_segment_area = 10;
  // const RNScalar min_segment_connectivity = 0.75;
  const RNScalar min_segment_connectivity = 0.5;
  const RNScalar target_segment_area = 200;
  const int min_segment_nentries = 100;
  RNLength min_feature_spacing_squared = min_feature_spacing * min_feature_spacing;

  // Check scene
  if (!scene) return 0;

  // Get image directory name
  char laser_image_directory[4096];
  sprintf(laser_image_directory, "%s/laser_images", scene->CacheDataDirectoryName());

  // Create plane features and clusters
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (scan->NScanlines() == 0) continue;
        if (ia == 1) continue;

        // Read DA grids
        char image_name[4096];
        R2Grid da_scanline_index_grid;
        R2Grid da_point_index_grid;
        R2Grid da_position_x_grid;
        R2Grid da_position_y_grid;
        R2Grid da_position_z_grid;
        R2Grid da_plane_segment_A_grid;
        R2Grid da_plane_segment_B_grid;
        R2Grid da_plane_segment_C_grid;
        R2Grid da_plane_segment_D_grid;
        R2Grid da_plane_segment_id_grid;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_Scanline.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_scanline_index_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PointIndex.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_point_index_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionX.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_position_x_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionY.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_position_y_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_PositionZ.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_position_z_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_SmallPlaneA.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_plane_segment_A_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_SmallPlaneB.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_plane_segment_B_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_SmallPlaneC.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_plane_segment_C_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_SmallPlaneD.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_plane_segment_D_grid.Read(image_name)) return 0;
        sprintf(image_name, "%s/%s/%02d_%02d_DA_SmallPlaneId.grd", laser_image_directory, run->Name(), is, ia);
        if (!da_plane_segment_id_grid.Read(image_name)) return 0;

        // Check the grid resolution
        if (da_plane_segment_id_grid.YResolution() < 25) continue;

        // Mask scan points on GSV car
        for (int iy = 0; iy < 25; iy++) {
          for (int ix = 0; ix < da_plane_segment_id_grid.XResolution(); ix++) {
            da_scanline_index_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_point_index_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_position_x_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_position_y_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_position_z_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_plane_segment_A_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_plane_segment_B_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_plane_segment_C_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_plane_segment_D_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
            da_plane_segment_id_grid.SetGridValue(ix, iy, R2_GRID_UNKNOWN_VALUE);
          }
        }

        // Allocate temporary segment data 
        int nsegments = (int) (da_plane_segment_id_grid.Maximum() + 0.5) + 1;
        RNScalar *segment_connectivities = CreatePlaneSegmentConnectivities(nsegments, da_plane_segment_id_grid); 
        RNScalar *segment_areas = CreatePlaneSegmentAreas(nsegments, da_plane_segment_id_grid, 
          da_position_x_grid, da_position_y_grid, da_position_z_grid);
        RNArray<const RNScalar *> *segment_entries = new RNArray<const RNScalar *> [ nsegments ];
        for (int i = 0; i < da_plane_segment_id_grid.NEntries(); i++) {
          RNScalar segment_id_value = da_plane_segment_id_grid.GridValue(i);
          if (segment_id_value == R2_GRID_UNKNOWN_VALUE) continue;
          int segment_id = (int) (segment_id_value + 0.5);
          assert((segment_id >= 0) && (segment_id < nsegments));
          const RNScalar *segment_entry = &da_plane_segment_id_grid(i);
          segment_entries[segment_id].Insert(segment_entry);
        }

        // Randomize the order of the segment entries
        for (int segment_id = 0; segment_id < nsegments; segment_id++) {
          for (int k1 = 0; k1 < segment_entries[segment_id].NEntries(); k1++) {
            int k2 = (int) (RNRandomScalar() * segment_entries[segment_id].NEntries());
            segment_entries[segment_id].Swap(k1, k2);
          }
        }

        // Create features at points sampled from within plane segment
        for (int segment_id = 0; segment_id < nsegments; segment_id++) {
          // Check segment area
          RNScalar segment_area = segment_areas[segment_id];
          if (segment_area < min_segment_area) continue;

          // Check segment connectivity
          RNScalar segment_connectivity = segment_connectivities[segment_id];
          if (segment_connectivity < min_segment_connectivity) continue;
          
          // Check segment entries
          int segment_nentries = segment_entries[segment_id].NEntries();
          if (segment_nentries < min_segment_nentries) continue;

          // Create unique group number
          int group = (int) (RNRandomScalar() * INT_MAX);

          // printf("%6d %6d %9.3g %9.3g\n", segment_id, segment_entries[segment_id].NEntries(), segment_area, segment_connectivity);

          // Create features on segment with minimum spacing
          RNArray<GSVFeature *> created_features;
          for (int i = 0; i < segment_entries[segment_id].NEntries(); i++) {
            int grid_ix, grid_iy;
            const RNScalar *segment_entry = segment_entries[segment_id].Kth(i);
            int grid_index = segment_entry - da_plane_segment_id_grid.GridValues();
            da_plane_segment_id_grid.IndexToIndices(grid_index, grid_ix, grid_iy);

            // Determine the position
            RNScalar x = da_position_x_grid.GridValue(grid_index);
            RNScalar y = da_position_y_grid.GridValue(grid_index);
            RNScalar z = da_position_z_grid.GridValue(grid_index);
            if (x == R2_GRID_UNKNOWN_VALUE) continue;
            if (y == R2_GRID_UNKNOWN_VALUE) continue;
            if (z == R2_GRID_UNKNOWN_VALUE) continue;
            R3Point position(x, y, z);

            // Check if there is a previously created feature within minimum spacing
            RNBoolean well_spaced = TRUE;
            for (int j = 0; j < created_features.NEntries(); j++) {
              GSVFeature *created_feature = created_features.Kth(j);
              RNScalar distance_squared = R3SquaredDistance(position, created_feature->scan_position);
              if (distance_squared < min_feature_spacing_squared) { well_spaced = FALSE; break; }
            }

            // Create feature if not too close to previously created one
            if (well_spaced) {
              // Get feature information
              RNScalar scanline_index_value = da_scanline_index_grid.GridValue(grid_index);
              RNScalar point_index_value = da_point_index_grid.GridValue(grid_index);
              RNScalar a = da_plane_segment_A_grid.GridValue(grid_index);
              RNScalar b = da_plane_segment_B_grid.GridValue(grid_index);
              RNScalar c = da_plane_segment_C_grid.GridValue(grid_index);
              if (scanline_index_value == R2_GRID_UNKNOWN_VALUE) continue;
              if (point_index_value == R2_GRID_UNKNOWN_VALUE) continue;
              if (a == R2_GRID_UNKNOWN_VALUE) continue;
              if (b == R2_GRID_UNKNOWN_VALUE) continue;
              if (c == R2_GRID_UNKNOWN_VALUE) continue;
              int scanline_index = (int) (scanline_index_value + 0.5);
              int point_index = (int) (point_index_value + 0.5);
              GSVScanline *scanline = scan->Scanline(scanline_index);
              R3Vector normal(a, b, c);
              normal.Normalize();
              RNScalar score = (segment_area < target_segment_area) ? segment_area / target_segment_area : 1.0;
              GSVFeature *feature = new GSVFeature(GSV_SCAN_PLANE_FEATURE_TYPE, scanline, point_index,  
                position, R3zero_vector, normal, min_feature_spacing, score, group);
              created_features.Insert(feature);
              InsertFeature(feature);
            }
          }

          // Create cluster 
          if (created_features.NEntries() > 1) {
            R3Point centroid(0, 0, 0);
            for (int j = 0; j < created_features.NEntries(); j++) 
              centroid += created_features[j]->scan_position;
            centroid /= created_features.NEntries();
            R3Vector normal = created_features.Head()->scan_normal;
            GSVFeatureCluster *cluster = new GSVFeatureCluster(GSV_SCAN_PLANE_FEATURE_TYPE,
              centroid, R3zero_vector, normal, min_feature_spacing, 1.0);
            for (int j = 0; j < created_features.NEntries(); j++) {
              GSVFeature *feature = created_features.Kth(j);
              cluster->InsertFeature(feature);
            }
            InsertCluster(cluster);
          }
        }

        // Delete temporary segment data 
        delete [] segment_connectivities;
        delete [] segment_areas;
        delete [] segment_entries;
      }
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Scoring
////////////////////////////////////////////////////////////////////////

RNScalar GSVPoseOptimization::
WorldRMSD(void)
{
  // Initialize SSD
  int count = 0;
  RNScalar sum = 0;

  // Add squared residuals for correspondences
  for (int i = 0; i < NCorrespondences(); i++) {
    GSVFeatureCorrespondence *correspondence = Correspondence(i);
    GSVFeature *feature0 = correspondence->Feature(0);
    GSVFeature *feature1 = correspondence->Feature(1);
    if (feature0->image || feature1->image) continue;
    RNScalar dd = SquaredWorldDistance(feature0, feature1);
    if (dd == RN_INFINITY) continue;
    sum += dd;
    count++;
  }

  // Add squared residuals for clusters
  for (int i = 0; i < NClusters(); i++) {
    GSVFeatureCluster *cluster = Cluster(i);

    // Add residuals for cluster-feature
    for (int j = 0; j < cluster->NFeatures(); j++) {
      GSVFeature *feature = cluster->Feature(j);
      if (feature->image) continue;
      RNScalar dd = SquaredWorldDistance(cluster, feature);
      if (dd == RN_INFINITY) continue;
      sum += dd;
      count++;
    }

    // Add residuals for cluster-child
    for (int j = 0; j < cluster->NClusters(); j++) {
      GSVFeatureCluster *child = cluster->Cluster(j);
      RNScalar dd = SquaredWorldDistance(cluster, child);
      if (dd == RN_INFINITY) continue;
      sum += dd;
      count++;
    }
  }

  // Return RMSD
  return (count > 0) ? sqrt(sum/count) : 0;
}



RNScalar GSVPoseOptimization::
ImageRMSD(void)
{
  // Initialize SSD
  int count = 0;
  RNScalar sum = 0;

  // Add squared residuals for correspondences
  for (int i = 0; i < NCorrespondences(); i++) {
    GSVFeatureCorrespondence *correspondence = Correspondence(i);
    GSVFeature *feature0 = correspondence->Feature(0);
    GSVFeature *feature1 = correspondence->Feature(1);
    if (!feature0->image && !feature1->image) continue;
    RNScalar dd = SquaredImageDistance(feature0, feature1);
    if (dd == RN_INFINITY) continue;
    sum += dd;
    count++;
  }

  // Add squared residuals for clusters
  for (int i = 0; i < NClusters(); i++) {
    GSVFeatureCluster *cluster = Cluster(i);

    // Add residuals for cluster-child
    for (int j = 0; j < cluster->NFeatures(); j++) {
      GSVFeature *feature = cluster->Feature(j);
      if (!feature->image) continue;
      RNScalar dd = SquaredImageDistance(cluster, feature);
      if (dd == RN_INFINITY) continue;
      sum += dd;
      count++;
    }
  }

  // Return RMSD
  return (count > 0) ? sqrt(sum/count) : 0;
}



RNScalar GSVPoseOptimization::
Score(void)
{
  // Score is not implemented
  RNScalar score = 0.0;
      
  // Return score
  return score;
}



////////////////////////////////////////////////////////////////////////
// Input and Output
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
ReadFile(const char *filename)
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .fet)\n", filename);
    return 0;
  }

  // Read file of appropriate type
  if (!strncmp(extension, ".txt", 4)) {
    if (!ReadAsciiFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".fet", 4)) {
    if (!ReadBinaryFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".sft", 4)) {
    if (!ReadSiftFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".lin", 4)) {
    if (!ReadLineFile(filename)) return 0;
  }
  else {
    fprintf(stderr, "Unable to read file %s (unrecognized extension: %s)\n", filename, extension);
    return 0;
  }

  // Update pixel positions
  UpdatePixelPositionsFromRayIntersections(); 
  UpdateClusterPositionsFromFeaturePositions();  

  // Return success
  return 1;
}



int GSVPoseOptimization::
WriteFile(const char *filename, RNBoolean apply_pose_transformations) const
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .fet)\n", filename);
    return 0;
  }

  // Write file of appropriate type
  if (!strncmp(extension, ".txt", 4)) {
    if (!WriteAsciiFile(filename, apply_pose_transformations)) return 0;
  }
  else if (!strncmp(extension, ".fet", 4)) {
    if (!WriteBinaryFile(filename, apply_pose_transformations)) return 0;
  }
  else {
    fprintf(stderr, "Unable to read file %s (unrecognized extension: %s)\n", filename, extension);
    return 0;
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Ascii Input and Output
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
ReadAsciiFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }

  // Read file
  if (!ReadAscii(fp)) {
    fprintf(stderr, "Unable to read %s\n", filename);
    return 0;
  }

  // Close file
  fclose(fp);
  
  // Return success
  return 1;
}



int GSVPoseOptimization::
WriteAsciiFile(const char *filename, RNBoolean apply_pose_transformations) const
{
  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }

  // Write file
  if (!WriteAscii(fp, apply_pose_transformations)) {
    fprintf(stderr, "Unable to write %s\n", filename);
    return 0;
  }

  // Close file
  fclose(fp);
  
  // Return success
  return 1;
}



int GSVPoseOptimization::
ReadAscii(FILE *fp)
{
  // Check file
  if (!fp) fp = stdout;

  // Read file
  char keyword[1024];
  RNArray<GSVFeature *> read_features;
  RNArray<GSVFeatureCluster *> read_clusters;
  while (fscanf(fp, "%s", keyword) == (unsigned int) 1) {
    if (keyword[0] == '#') {
      char buffer[4096];
      fgets(buffer, 4096, fp);
    }
    else if (!strcmp(keyword, "F")) {
      // Parse feature type
      int feature_type;
      if (fscanf(fp, "%d", &feature_type) != (unsigned int) 1) {
        fprintf(stderr, "Invalid feature command\n");
        return 0;
      }

      // Parse feature info
      if (feature_type == GSV_SCAN_POINT_FEATURE_TYPE) {
        // Read scan info
        char run_name[1024];
        int segment_index, scan_index, scanline_index, point_index;
        double score, px, py, pz, dummy;
        if (fscanf(fp, "%s%d%d%d%d%lf%lf%lf%lf%lf%lf%lf", 
          run_name, &segment_index, &scan_index, &scanline_index, &point_index, &score,
          &px, &py, &pz, &dummy, &dummy, &dummy) != (unsigned int) 12) {
          fprintf(stderr, "Error parsing feature %d\n", NFeatures());
          return 0;
        }

        // Get scene info
        GSVRun *run = scene->Run(run_name);
        if (!run) { 
          fprintf(stderr, "Bad run %s\n", run_name);
          return 0;
        }
        if (segment_index >= run->NSegments()) { 
          fprintf(stderr, "Bad segment %d\n", segment_index);
          return 0;
        }
        GSVSegment *segment = run->Segment(segment_index);
        if (scan_index >= 0) {
          if (scan_index >= segment->NScans()) { 
            fprintf(stderr, "Bad scan %d\n", scan_index);
            return 0;
          }
        }
        GSVScan *scan = segment->Scan(scan_index);
        if (scanline_index >= scan->NScanlines()) {
          fprintf(stderr, "Bad scanline %d\n", scanline_index);
          return 0;
        }
        GSVScanline *scanline = scan->Scanline(scanline_index);
        R3Point position(px, py, pz);

        // Create feature
        GSVFeature *feature = new GSVFeature(feature_type, scanline, point_index, position, R3zero_vector, R3zero_vector, 1.0, score);
        read_features.Insert(feature);
        InsertFeature(feature);
      }
      else if (feature_type == GSV_SCAN_LINE_FEATURE_TYPE) {
        // Read scan info
        char run_name[1024];
        int segment_index, scan_index, scanline_index, point_index;
        double score, px, py, pz, dx, dy, dz, scale, dummy;
        if (fscanf(fp, "%s%d%d%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
          run_name, &segment_index, &scan_index, &scanline_index, &point_index, &score,
          &px, &py, &pz, &dx, &dy, &dz, &scale, &dummy, &dummy, &dummy, &dummy) != (unsigned int) 17) {
          fprintf(stderr, "Error parsing feature %d\n", NFeatures());
          return 0;
        }

        // Get scene info
        GSVRun *run = scene->Run(run_name);
        if (!run) { 
          fprintf(stderr, "Bad run %s\n", run_name);
          return 0;
        }
        if (segment_index >= run->NSegments()) { 
          fprintf(stderr, "Bad segment %d\n", segment_index);
          return 0;
        }
        GSVSegment *segment = run->Segment(segment_index);
        if (scan_index >= 0) {
          if (scan_index >= segment->NScans()) { 
            fprintf(stderr, "Bad scan %d\n", scan_index);
            return 0;
          }
        }
        GSVScan *scan = segment->Scan(scan_index);
        if (scanline_index >= scan->NScanlines()) {
          fprintf(stderr, "Bad scanline %d\n", scanline_index);
          return 0;
        }
        GSVScanline *scanline = scan->Scanline(scanline_index);
        R3Point position(px, py, pz);
        R3Vector direction(dx, dy, dz);

        // Create feature
        GSVFeature *feature = new GSVFeature(feature_type, scanline, point_index, position, direction, R3zero_vector, scale, score);
        read_features.Insert(feature);
        InsertFeature(feature);
      }
      else if (feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) {
        // Read scan info
        char run_name[1024];
        int segment_index, scan_index, scanline_index, point_index;
        double score, px, py, pz, nx, ny, nz;
        if (fscanf(fp, "%s%d%d%d%d%lf%lf%lf%lf%lf%lf%lf", 
          run_name, &segment_index, &scan_index, &scanline_index, &point_index, &score,
          &px, &py, &pz, &nx, &ny, &nz) != (unsigned int) 12) {
          fprintf(stderr, "Error parsing feature %d\n", NFeatures());
          return 0;
        }

        // Get scene info
        GSVRun *run = scene->Run(run_name);
        if (!run) { 
          fprintf(stderr, "Bad run %s\n", run_name);
          return 0;
        }
        if (segment_index >= run->NSegments()) { 
          fprintf(stderr, "Bad segment %d\n", segment_index);
          return 0;
        }
        GSVSegment *segment = run->Segment(segment_index);
        if (scan_index >= 0) {
          if (scan_index >= segment->NScans()) { 
            fprintf(stderr, "Bad scan %d\n", scan_index);
            return 0;
          }
        }
        GSVScan *scan = segment->Scan(scan_index);
        if (scanline_index >= scan->NScanlines()) {
          fprintf(stderr, "Bad scanline %d\n", scanline_index);
          return 0;
        }
        GSVScanline *scanline = scan->Scanline(scanline_index);
        R3Point position(px, py, pz);
        R3Vector normal(nx, ny, nz);

        // Create feature
        GSVFeature *feature = new GSVFeature(feature_type, scanline, point_index, position, R3zero_vector, normal, 1.0, score);
        read_features.Insert(feature);
        InsertFeature(feature);
      }
      else if (feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) {
        // Read image info
        char run_name[1024];
        int segment_index, panorama_index, image_index;
        double score, px, py, dummy;
        if (fscanf(fp, "%s%d%d%d%lf%lf%lf%lf%lf", 
          run_name, &segment_index, &panorama_index, &image_index, &score,
          &px, &py, &dummy, &dummy) != (unsigned int) 9) {
          fprintf(stderr, "Error parsing feature %d\n", read_features.NEntries());
          return 0;
        }

        // Get scene info
        GSVRun *run = scene->Run(run_name);
        if (!run) { 
          fprintf(stderr, "Bad run %s\n", run_name);
          return 0;
        }
        if (segment_index >= run->NSegments()) { 
          fprintf(stderr, "Bad segment %d\n", segment_index);
          return 0;
        }
        GSVSegment *segment = run->Segment(segment_index);
        if (panorama_index >= 0) {
          if (panorama_index >= segment->NPanoramas()) { 
            fprintf(stderr, "Bad panorama %d\n", panorama_index);
            return 0;
          }
        }
        GSVPanorama *panorama = segment->Panorama(panorama_index);
        if (image_index >= panorama->NImages()) {
          fprintf(stderr, "Bad image %d\n", image_index);
          return 0;
        }
        GSVImage *image = panorama->Image(image_index);

        // Create feature
        R2Point position(px, py);
        GSVFeature *feature = new GSVFeature(feature_type, image, position, R2zero_vector, 1.0, score);
        read_features.Insert(feature);
        InsertFeature(feature);
      }
      else if (feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) {
        // Read image info
        char run_name[1024];
        int segment_index, panorama_index, image_index;
        double score, px, py, dx, dy, scale, t, dummy;
        if (fscanf(fp, "%s%d%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
          run_name, &segment_index, &panorama_index, &image_index, &score,
          &px, &py, &dx, &dy, &scale, &t, &dummy, &dummy, &dummy, &dummy) != (unsigned int) 15) {
          fprintf(stderr, "Error parsing feature %d\n", read_features.NEntries());
          return 0;
        }

        // Get scene info
        GSVRun *run = scene->Run(run_name);
        if (!run) { 
          fprintf(stderr, "Bad run %s\n", run_name);
          return 0;
        }
        if (segment_index >= run->NSegments()) { 
          fprintf(stderr, "Bad segment %d\n", segment_index);
          return 0;
        }
        GSVSegment *segment = run->Segment(segment_index);
        if (panorama_index >= 0) {
          if (panorama_index >= segment->NPanoramas()) { 
            fprintf(stderr, "Bad panorama %d\n", panorama_index);
            return 0;
          }
        }
        GSVPanorama *panorama = segment->Panorama(panorama_index);
        if (image_index >= panorama->NImages()) {
          fprintf(stderr, "Bad image %d\n", image_index);
          return 0;
        }
        GSVImage *image = panorama->Image(image_index);
        R2Point position(px, py);
        R2Vector direction(dx, dy);

        // Create feature
        GSVFeature *feature = new GSVFeature(feature_type, image, position, direction, scale, score);
        feature->image_t = t;
        read_features.Insert(feature);
        InsertFeature(feature);
      }
      else {
        fprintf(stderr, "Unrecognized feature type %d\n", feature_type);
        return 0;
      }
    }
    else if (!strcmp(keyword, "D")) {
      // Read descriptor info
      int feature_index, descriptor_type, nvalues;
      if (fscanf(fp, "%d%d%d", &feature_index, &descriptor_type, &nvalues) != (unsigned int) 3) {
        fprintf(stderr, "Unable to parse descriptor\n");
        return 0;
      }

      // Check feature index
      if ((feature_index < 0) || (feature_index >= read_features.NEntries())) {
        fprintf(stderr, "Feature index %d is out of range\n", feature_index);
        return 0;
      }

      // Read descriptor
      GSVFeature *feature = read_features.Kth(feature_index);
      GSVDescriptor& descriptor = feature->descriptor;
      descriptor.descriptor_type = descriptor_type;
      descriptor.nvalues = nvalues;
      descriptor.values = new RNScalar [ nvalues ];
      for (int i = 0; i < nvalues; i++) {
        if (fscanf(fp, "%lf", &descriptor.values[i]) != (unsigned int) 1) {
          fprintf(stderr, "Unable to parse descriptor\n");
          return 0;
        }
      }
    }
    else if (!strcmp(keyword, "P")) {
      // Read pair
      int feature_index0, feature_index1;
      double score;
      if (fscanf(fp, "%d%d%lf", &feature_index0, &feature_index1, &score) != (unsigned int) 3) {
        fprintf(stderr, "Unable to parse pair %d\n", NPairs());
        return 0;
      }

      // Check feature indices
      if ((feature_index0 < 0) || (feature_index0 >= read_features.NEntries())) {
        fprintf(stderr, "Feature index %d is out of range\n", feature_index0);
        return 0;
      }
      if ((feature_index1 < 0) || (feature_index1 >= read_features.NEntries())) {
        fprintf(stderr, "Feature index %d is out of range\n", feature_index1);
        return 0;
      }

      // Create pair
      GSVFeature *feature0 = read_features.Kth(feature_index0);
      GSVFeature *feature1 = read_features.Kth(feature_index1);
      GSVFeaturePair *pair = new GSVFeaturePair(feature0, feature1, score);
      InsertPair(pair);
    }
    else if (!strcmp(keyword, "C")) {
      // Read correspondence
      int feature_index0, feature_index1;
      double score;
      if (fscanf(fp, "%d%d%lf", &feature_index0, &feature_index1, &score) != (unsigned int) 3) {
        fprintf(stderr, "Unable to parse correspondence %d\n", NCorrespondences());
        return 0;
      }

      // Check feature indices
      if ((feature_index0 < 0) || (feature_index0 >= read_features.NEntries())) {
        fprintf(stderr, "Feature index %d is out of range\n", feature_index0);
        return 0;
      }
      if ((feature_index1 < 0) || (feature_index1 >= read_features.NEntries())) {
        fprintf(stderr, "Feature index %d is out of range\n", feature_index1);
        return 0;
      }

      // Create correspondence
      GSVFeature *feature0 = read_features.Kth(feature_index0);
      GSVFeature *feature1 = read_features.Kth(feature_index1);
      GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence(feature0, feature1, score);
      InsertCorrespondence(correspondence);
    }
    else if (!strcmp(keyword, "G")) {
      // Read cluster
      int nfeatures, nchildren, feature_type, dummy;
      char separator[1024];
      double px, py, pz, dx, dy, dz, nx, ny, nz, tx, ty, tz, rx, ry, rz, scale, score;
      if (fscanf(fp, "%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%d%d%d%d%s", &nfeatures, &feature_type,
        &px, &py, &pz, &dx, &dy, &dz, &nx, &ny, &nz, &scale, &score, &tx, &ty, &tz, &rx, &ry, &rz, 
        &nchildren, &dummy, &dummy, &dummy, separator) != (unsigned int) 24) {
        fprintf(stderr, "Unable to parse cluster %d\n", NClusters());
        return 0;
      }

      // Check separator
      if (strcmp(separator, ":")) {
        fprintf(stderr, "Bad format for cluster %d\n", NClusters());
        return 0;
      }

      // Create cluster
      GSVFeatureCluster *cluster = new GSVFeatureCluster();
      cluster->feature_type = feature_type;
      cluster->scan_position = R3Point(px, py, pz);
      cluster->scan_direction = R3Vector(dx, dy, dz);
      cluster->scan_normal = R3Vector(nx, ny, nz);
      cluster->scan_scale = scale;
      cluster->score = score;
      cluster->translation = R3Vector(tx, ty, tz);
      cluster->rotation = R3Vector(rx, ry, rz);
      read_clusters.Insert(cluster);
      InsertCluster(cluster);

      // Read feature indices
      for (int j = 0; j < nfeatures; j++) {
        // Read index
        int feature_index;
        if (fscanf(fp, "%d", &feature_index) != (unsigned int) 1) {
          fprintf(stderr, "Unable to parse feature %d in cluster %d\n", j, NClusters());
          return 0;
        }
        
        // Check index
        if ((feature_index < 0) || (feature_index >= NFeatures())) {
          fprintf(stderr, "Invalid feature %d in cluster %d\n", feature_index, NClusters());
          return 0;
        }

        // Insert feature into cluster
        GSVFeature *feature = read_features.Kth(feature_index);
        cluster->InsertFeature(feature);
      }
    }
    else if (!strcmp(keyword, "GC")) {
      // Read cluster indices
      int cluster_index, child_index;
      if (fscanf(fp, "%d%d", &cluster_index, &child_index) != (unsigned int) 2) {
        fprintf(stderr, "Unable to parse cluster child\n");
        return 0;
      }

      // Check cluster index
      if ((cluster_index < 0) || (cluster_index >= read_clusters.NEntries())) {
        fprintf(stderr, "Invalid cluster index %d\n", cluster_index);
        return 0;
      }

      // Check child index
      if ((child_index < 0) || (child_index >= read_clusters.NEntries())) {
        fprintf(stderr, "Invalid child index %d\n", child_index);
        return 0;
      }

      // Insert child into cluster
      GSVFeatureCluster *cluster = read_clusters.Kth(cluster_index);
      GSVFeatureCluster *child = read_clusters.Kth(child_index);
      cluster->InsertCluster(child);
    }
    else if (!strcmp(keyword, "include")) {
      // Read filename
      char filename[2048];
      if (fscanf(fp, "%s", filename) != (unsigned int) 1) {
        fprintf(stderr, "Unable to parse include filename\n");
        return 0;
      }

      // Read included file
      if (!ReadFile(filename)) return 0;
    }
    else {
      // File format error
      fprintf(stderr, "Unrecognized keyword (%s) in correspondence file\n", keyword);
      return 0;
    }
  }

  // Return success
  return 1;
}



int GSVPoseOptimization::
WriteAscii(FILE *fp, RNBoolean apply_pose_transformations) const
{
  // Check file
  if (!fp) fp = stdout;

  // Write features
  for (int i = 0; i < NFeatures(); i++) {
    GSVFeature *feature = Feature(i);

    // Get scene info
    GSVScanline *scanline = feature->scanline;
    GSVImage *image = feature->image;
    GSVScan *scan = (scanline) ? scanline->Scan() : NULL;
    GSVPanorama *panorama = (image) ? image->Panorama() : NULL;
    assert(!panorama || !scan || (panorama->Segment() == scan->Segment()));
    GSVSegment *segment = (scan) ? scan->Segment() : NULL;
    GSVRun *run = (segment) ? segment->Run() : NULL;
    const char *run_name = (run) ? run->Name() : "None";
    int scanline_index = (scanline) ? scanline->ScanIndex() : -1;
    int scan_index = (scan) ? scan->SegmentIndex() : -1;
    int segment_index = (segment) ? segment->RunIndex() : -1;
    int image_index = (image) ? image->PanoramaIndex() : -1;
    int panorama_index = (panorama) ? panorama->SegmentIndex() : -1;

    // Write feature type
    fprintf(fp, "F %2d  ", feature->feature_type);

    // Write feature info
    if (feature->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) {
      // Write feature info
      R3Point scan_position = WorldPosition(feature, apply_pose_transformations); 
      fprintf(fp, "%30s %2d %2d %6d %3d   %8.6f   %12.6f %12.6f %12.6f  0 0 0\n", 
        run_name, segment_index, scan_index, scanline_index, feature->scan_point_index, feature->score,
        scan_position.X(), scan_position.Y(), scan_position.Z());
    }
    else if (feature->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) {
      // Write feature info
      R3Point scan_position = WorldPosition(feature, apply_pose_transformations); 
      R3Vector scan_direction = WorldDirection(feature, apply_pose_transformations); 
      fprintf(fp, "%30s %2d %2d %6d %3d   %8.6f   %12.6f %12.6f %12.6f  %8.6f %8.6f %8.6f  %g  0 0 0 0\n", 
        run_name, segment_index, scan_index, scanline_index, feature->scan_point_index, feature->score,
        scan_position.X(), scan_position.Y(), scan_position.Z(), 
        scan_direction.X(), scan_direction.Y(), scan_direction.Z(), 
        feature->scan_scale);
    }
    else if (feature->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) {
      // Write feature info
      R3Point scan_position = WorldPosition(feature, apply_pose_transformations); 
      R3Vector scan_normal = WorldNormal(feature, apply_pose_transformations); 
      fprintf(fp, "%30s %2d %2d %6d %3d   %8.6f   %12.6f %12.6f %12.6f  %8.6f %8.6f %8.6f\n", 
        run_name, segment_index, scan_index, scanline_index, feature->scan_point_index, feature->score,
        scan_position.X(), scan_position.Y(), scan_position.Z(), 
        scan_normal.X(), scan_normal.Y(), scan_normal.Z());
    }
    else if (feature->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) {
      // Write feature info
      fprintf(fp, "%30s %2d %6d %2d   %8.6f   %12.6f %12.6f  0  0\n", 
        run_name, segment_index,  panorama_index, image_index, feature->score,
        feature->image_position.X(), feature->image_position.Y());
    }
    else if (feature->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) {
      // Write feature info
      fprintf(fp, "%30s %2d %6d %2d   %8.6f   %12.6f %12.6f   %8.6f %8.6f  %8.6f %8.6f 0 0 0 0\n", 
        run_name, segment_index,  panorama_index, image_index, feature->score,
        feature->image_position.X(), feature->image_position.Y(), 
        feature->image_direction.X(), feature->image_direction.Y(),
        feature->image_scale, feature->image_t);
    }
  }

  ////////////////////////////////////////
  
  // Write descriptors
  for (int i = 0; i < NFeatures(); i++) {
    GSVFeature *feature = Feature(i);
    GSVDescriptor& descriptor = feature->descriptor;
    fprintf(fp, "D %6d %2d %3d  ", i, descriptor.descriptor_type, descriptor.nvalues);
    for (int i = 0; i < descriptor.nvalues; i++) {
      fprintf(fp, "%12.6f ", descriptor.values[i]);
    }
    fprintf(fp, "\n");
  }

  ////////////////////////////////////////
  
  // Write pairs
  for (int i = 0; i < NPairs(); i++) {
    GSVFeaturePair *pair = Pair(i);
    GSVFeature *feature0 = pair->Feature(0);
    GSVFeature *feature1 = pair->Feature(1);
    fprintf(fp, "P %d %d %g\n", feature0->index, feature1->index, pair->score);
  }
  
  ////////////////////////////////////////
  
  // Write correspondences
  for (int i = 0; i < NCorrespondences(); i++) {
    GSVFeatureCorrespondence *correspondence = Correspondence(i);
    GSVFeature *feature0 = correspondence->Feature(0);
    GSVFeature *feature1 = correspondence->Feature(1);
    fprintf(fp, "C %d %d %g\n", feature0->index, feature1->index, correspondence->score);
  }
  
  ////////////////////////////////////////
  
  // Write clusters
  for (int i = 0; i < NClusters(); i++) {
    GSVFeatureCluster *cluster = Cluster(i);
    fprintf(fp, "G  %d %d  %g %g %g  %g %g %g  %g %g %g  %g  %g  %g %g %g  %g %g %g  %d 0 0 0 : ", 
      cluster->NFeatures(), cluster->feature_type,
      cluster->scan_position.X(), cluster->scan_position.Y(), cluster->scan_position.Z(), 
      cluster->scan_direction.X(), cluster->scan_direction.Y(), cluster->scan_direction.Z(), 
      cluster->scan_normal.X(), cluster->scan_normal.Y(), cluster->scan_normal.Z(), 
      cluster->scan_scale, cluster->score,
      cluster->translation.X(), cluster->translation.Y(), cluster->translation.Z(), 
      cluster->rotation.X(), cluster->rotation.Y(), cluster->rotation.Z(), 
      cluster->NClusters());
    for (int j = 0; j < cluster->NFeatures(); j++) {
      GSVFeature *feature = cluster->Feature(j);
      fprintf(fp, "%d ", feature->index);
    }
    fprintf(fp, "\n");
  }
  
  // Write cluster hierarchy
  for (int i = 0; i < NClusters(); i++) {
    GSVFeatureCluster *cluster = Cluster(i);
    for (int j = 0; j < cluster->NClusters(); j++) {
      GSVFeatureCluster *child = cluster->Cluster(j);
      fprintf(fp, "GC %d %d\n", cluster->index, child->index);
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Binary Input and Output
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
ReadBinaryFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }

  // Read file
  if (!ReadBinary(fp)) {
    fprintf(stderr, "Unable to read %s\n", filename);
    return 0;
  }

  // Close file
  fclose(fp);
  
  // Return success
  return 1;
}



int GSVPoseOptimization::
WriteBinaryFile(const char *filename, RNBoolean apply_pose_transformations) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }

  // Write file
  if (!WriteBinary(fp, apply_pose_transformations)) {
    fprintf(stderr, "Unable to write %s\n", filename);
    return 0;
  }

  // Close file
  fclose(fp);
  
  // Return success
  return 1;
}



int GSVPoseOptimization::
ReadBinary(FILE *fp)
{
  // Check file
  if (!fp) fp = stdout;

  // Read endian test
  int endian_test;
  if (fread(&endian_test, sizeof(int), 1, fp) != (unsigned int) 1) {
    fprintf(stderr, "Unable to read correspondence file\n");
    return 0;
  }
  
  // Check endian test
  if (endian_test != 1) {
    fprintf(stderr, "Unable to read format of correspondence file\n");
    return 0;
  }

  // Read rest of header
  int dummy = 0;
  int major_version, minor_version;
  int nfeatures, npairs, ncorrespondences, nclusters;
  int read_cluster_children;
  fread(&major_version, sizeof(int), 1, fp);
  fread(&minor_version, sizeof(int), 1, fp);
  fread(&nfeatures, sizeof(int), 1, fp);
  fread(&npairs, sizeof(int), 1, fp);
  fread(&ncorrespondences, sizeof(int), 1, fp);
  fread(&nclusters, sizeof(int), 1, fp);
  fread(&read_cluster_children, sizeof(int), 1, fp);
  for (int i = 0; i < 56; i++) fread(&dummy, sizeof(int), 1, fp);

  // Read features
  RNArray<GSVFeature *> read_features;
  for (int i = 0; i < nfeatures; i++) {
    GSVFeature *feature = new GSVFeature();
    assert(feature);
    int scanline_index, image_index, descriptor_nvalues;
    fread(&feature->feature_type, sizeof(int), 1, fp);
    fread(&image_index, sizeof(int), 1, fp);
    fread(&scanline_index, sizeof(int), 1, fp);
    fread(&feature->scan_point_index, sizeof(int), 1, fp);
    fread(&feature->scan_position, sizeof(RNCoord), 3, fp);
    fread(&feature->scan_direction, sizeof(RNCoord), 3, fp);
    fread(&feature->scan_normal, sizeof(RNCoord), 3, fp);
    fread(&feature->scan_scale, sizeof(RNScalar), 1, fp);
    fread(&feature->image_position, sizeof(RNCoord), 2, fp);
    fread(&feature->image_direction, sizeof(RNCoord), 2, fp);
    fread(&feature->image_scale, sizeof(RNScalar), 1, fp);
    fread(&feature->image_t, sizeof(RNScalar), 1, fp);
    fread(&feature->score, sizeof(RNScalar), 1, fp);
    fread(&feature->group, sizeof(int), 1, fp);
    for (int j = 0; j < 20; j++) fread(&dummy, sizeof(int), 1, fp);
    fread(&feature->descriptor.descriptor_type, sizeof(int), 1, fp);
    fread(&descriptor_nvalues, sizeof(int), 1, fp);
    if (descriptor_nvalues > 0) {
      feature->descriptor.nvalues = descriptor_nvalues;
      feature->descriptor.values = new RNScalar [ feature->descriptor.nvalues ];
      fread(feature->descriptor.values, sizeof(RNScalar), feature->descriptor.nvalues, fp);
    }
    feature->scanline = (scanline_index >= 0) ? scanlines.Kth(scanline_index) : NULL;
    feature->image = (image_index >= 0) ? images.Kth(image_index) : NULL;
    read_features.Insert(feature);
    InsertFeature(feature);
  }

  // Read pairs
  for (int i = 0; i < npairs; i++) {
    GSVFeaturePair *pair = new GSVFeaturePair();
    RNScalar dummy[3];
    assert(pair);
    int feature_index0, feature_index1;
    fread(&feature_index0, sizeof(int), 1, fp);
    fread(&feature_index1, sizeof(int), 1, fp);
    fread(&pair->score, sizeof(RNScalar), 1, fp);
    fread(dummy, sizeof(RNScalar), 3, fp);
    for (int j = 0; j < 7; j++) fread(&dummy, sizeof(int), 1, fp);
    pair->features[0] = (feature_index0 >= 0) ? read_features.Kth(feature_index0) : NULL;
    pair->features[1] = (feature_index1 >= 0) ? read_features.Kth(feature_index1) : NULL;
    InsertPair(pair);
  }

  // Read correspondences
  for (int i = 0; i < ncorrespondences; i++) {
    GSVFeatureCorrespondence *correspondence = new GSVFeatureCorrespondence();
    assert(correspondence);
    int feature_index0, feature_index1;
    fread(&feature_index0, sizeof(int), 1, fp);
    fread(&feature_index1, sizeof(int), 1, fp);
    fread(&correspondence->score, sizeof(RNScalar), 1, fp);
    for (int j = 0; j < 13; j++) fread(&dummy, sizeof(int), 1, fp);
    correspondence->features[0] = (feature_index0 >= 0) ? read_features.Kth(feature_index0) : NULL;
    correspondence->features[1] = (feature_index1 >= 0) ? read_features.Kth(feature_index1) : NULL;
    InsertCorrespondence(correspondence);
  }

  // Read clusters
  RNArray<GSVFeatureCluster *> read_clusters;
  for (int i = 0; i < nclusters; i++) {
    int nfeatures, feature_index;
    GSVFeatureCluster *cluster = new GSVFeatureCluster();
    fread(&nfeatures, sizeof(int), 1, fp);
    fread(&cluster->feature_type, sizeof(int), 1, fp);
    fread(&cluster->scan_position, sizeof(RNScalar), 3, fp);
    fread(&cluster->scan_direction, sizeof(RNScalar), 3, fp);
    fread(&cluster->scan_normal, sizeof(RNScalar), 3, fp);
    fread(&cluster->scan_scale, sizeof(RNScalar), 1, fp);
    fread(&cluster->score, sizeof(RNScalar), 1, fp);
    fread(&cluster->translation, sizeof(RNScalar), 3, fp);
    fread(&cluster->rotation, sizeof(RNScalar), 3, fp);
    for (int j = 0; j < 16; j++) fread(&dummy, sizeof(int), 1, fp);
    for (int j = 0; j < nfeatures; j++) {
      fread(&feature_index, sizeof(int), 1, fp);
      GSVFeature *feature = read_features.Kth(feature_index);
      cluster->InsertFeature(feature);
    }
    read_clusters.Insert(cluster);
    InsertCluster(cluster);
  }

  // Read cluster children
  if (read_cluster_children) {
    int nchildren, child_index;
    for (int i = 0; i < read_clusters.NEntries(); i++) {
      GSVFeatureCluster *cluster = read_clusters.Kth(i);
      fread(&nchildren, sizeof(int), 1, fp);
      for (int j = 0; j < nchildren; j++) {
        fread(&child_index, sizeof(int), 1, fp);
        GSVFeatureCluster *child = read_clusters.Kth(child_index);
        cluster->InsertCluster(child);
      }
    }
  }

  // Return success
  return 1;
}



int GSVPoseOptimization::
WriteBinary(FILE *fp, RNBoolean apply_pose_transformations) const
{
  // Check file
  if (!fp) fp = stdout;

  // Convenient variables
  int dummy = 0;

  // Write endian test
  int endian_test = 1;
  if (fwrite(&endian_test, sizeof(int), 1, fp) != (unsigned int) 1) {
    fprintf(stderr, "Unable to write correspondence file\n");
    return 0;
  }

  // Write rest of header
  int major_version = 0;
  int minor_version = 1;
  int nfeatures = NFeatures();
  int npairs = NPairs();
  int ncorrespondences = NCorrespondences();
  int nclusters = NClusters();
  int write_cluster_children = 1;
  fwrite(&major_version, sizeof(int), 1, fp);
  fwrite(&minor_version, sizeof(int), 1, fp);
  fwrite(&nfeatures, sizeof(int), 1, fp);
  fwrite(&npairs, sizeof(int), 1, fp);
  fwrite(&ncorrespondences, sizeof(int), 1, fp);
  fwrite(&nclusters, sizeof(int), 1, fp);
  fwrite(&write_cluster_children, sizeof(int), 1, fp);
  for (int i = 0; i < 56; i++) fwrite(&dummy, sizeof(int), 1, fp);

  // Write features
  for (int i = 0; i < NFeatures(); i++) {
    GSVFeature *feature = Feature(i);
    GSVScanline *scanline = feature->scanline;
    ScanlineData *scanline_data = (scanline) ? (ScanlineData *) scanline->Data() : NULL;
    int scanline_index = (scanline_data) ? scanline_data->index : -1;
    R3Point scan_position = WorldPosition(feature, apply_pose_transformations); 
    R3Vector scan_direction = WorldDirection(feature, apply_pose_transformations); 
    R3Vector scan_normal = WorldNormal(feature, apply_pose_transformations); 
    GSVImage *image = feature->image;
    ImageData *image_data = (image) ? (ImageData *) image->Data() : NULL;
    int image_index = (image_data) ? image_data->index : -1;
    const GSVDescriptor& descriptor = feature->descriptor;
    fwrite(&feature->feature_type, sizeof(int), 1, fp);
    fwrite(&image_index, sizeof(int), 1, fp);
    fwrite(&scanline_index, sizeof(int), 1, fp);
    fwrite(&feature->scan_point_index, sizeof(int), 1, fp);
    fwrite(&scan_position, sizeof(RNCoord), 3, fp);
    fwrite(&scan_direction, sizeof(RNCoord), 3, fp);
    fwrite(&scan_normal, sizeof(RNCoord), 3, fp);
    fwrite(&feature->scan_scale, sizeof(RNScalar), 1, fp);
    fwrite(&feature->image_position, sizeof(RNCoord), 2, fp);
    fwrite(&feature->image_direction, sizeof(RNCoord), 2, fp);
    fwrite(&feature->image_scale, sizeof(RNScalar), 1, fp);
    fwrite(&feature->image_t, sizeof(RNScalar), 1, fp);
    fwrite(&feature->score, sizeof(RNScalar), 1, fp);
    fwrite(&feature->group, sizeof(int), 1, fp);
    for (int j = 0; j < 20; j++) fwrite(&dummy, sizeof(int), 1, fp);
    fwrite(&descriptor.descriptor_type, sizeof(int), 1, fp);
    fwrite(&descriptor.nvalues, sizeof(int), 1, fp);
    if (descriptor.nvalues > 0) {
      fwrite(descriptor.values, sizeof(RNScalar), descriptor.nvalues, fp);
    }
  }

  // Write pairs
  for (int i = 0; i < NPairs(); i++) {
    GSVFeaturePair *pair = Pair(i);
    const RNScalar dummy[3] = { 0.0, 0.0, 0.0 };
    int feature_index0 = pair->Feature(0)->index;
    int feature_index1 = pair->Feature(1)->index;
    fwrite(&feature_index0, sizeof(int), 1, fp);
    fwrite(&feature_index1, sizeof(int), 1, fp);
    fwrite(&pair->score, sizeof(RNScalar), 1, fp);
    fwrite(dummy, sizeof(RNScalar), 3, fp);
    for (int j = 0; j < 7; j++) fwrite(&dummy, sizeof(int), 1, fp);
  }

  // Write correspondences
  for (int i = 0; i < NCorrespondences(); i++) {
    GSVFeatureCorrespondence *correspondence = Correspondence(i);
    int feature_index0 = correspondence->Feature(0)->index;
    int feature_index1 = correspondence->Feature(1)->index;
    fwrite(&feature_index0, sizeof(int), 1, fp);
    fwrite(&feature_index1, sizeof(int), 1, fp);
    fwrite(&correspondence->score, sizeof(RNScalar), 1, fp);
    for (int j = 0; j < 13; j++) fwrite(&dummy, sizeof(int), 1, fp);
  }

  // Write clusters
  for (int i = 0; i < NClusters(); i++) {
    GSVFeatureCluster *cluster = Cluster(i);
    int nfeatures = cluster->NFeatures();
    fwrite(&nfeatures, sizeof(int), 1, fp);
    fwrite(&cluster->feature_type, sizeof(int), 1, fp);
    fwrite(&cluster->scan_position, sizeof(RNScalar), 3, fp);
    fwrite(&cluster->scan_direction, sizeof(RNScalar), 3, fp);
    fwrite(&cluster->scan_normal, sizeof(RNScalar), 3, fp);
    fwrite(&cluster->scan_scale, sizeof(RNScalar), 1, fp);
    fwrite(&cluster->score, sizeof(RNScalar), 1, fp);
    fwrite(&cluster->translation, sizeof(RNScalar), 3, fp);
    fwrite(&cluster->rotation, sizeof(RNScalar), 3, fp);
    for (int j = 0; j < 16; j++) fwrite(&dummy, sizeof(int), 1, fp);
    for (int j = 0; j < nfeatures; j++) {
      GSVFeature *feature = cluster->Feature(j);
      fwrite(&feature->index, sizeof(int), 1, fp);
    }
  }

  // Write cluster children
  if (write_cluster_children) {
    for (int i = 0; i < NClusters(); i++) {
      GSVFeatureCluster *cluster = Cluster(i);
      int nchildren = cluster->NClusters();
      fwrite(&nchildren, sizeof(int), 1, fp);
      for (int j = 0; j < nchildren; j++) {
        GSVFeatureCluster *child = cluster->Cluster(j);
        fwrite(&child->index, sizeof(int), 1, fp);
      }
    }
  }
    
  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Sift Feature Input 
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
ReadSiftFile(const char *filename)
{
  // Parse filename
  char run_buffer[4096];
  strncpy(run_buffer, filename, 4096);
  char *file_name_start = run_buffer;
  char *run_name_end = strrchr(run_buffer, '/');
  if (run_name_end) {
    file_name_start = run_name_end + 1;
    *run_name_end = '\0';
  }
  else {
    fprintf(stderr, "Unable to parse run directory from sift filename %s\n", filename); 
    return 0; 
  }

  // Get run name
  printf("%s\n", run_buffer);
  char *run_name = strrchr(run_buffer, '/');
  if (run_name) run_name++;
  else run_name = run_buffer;
  if (!run_name) {
    fprintf(stderr, "Bad run name in sift filename %s\n", filename); 
    return 0; 
  }

  // Parse segment, panorama, and image indicies from filename
  int segment_index, panorama_index, image_index;
  if (sscanf(file_name_start, "%d_%d_%d_Sift.sft",  
    &segment_index, &panorama_index, &image_index) != (unsigned int) 3) {
    fprintf(stderr, "Unable to parse image info from sift filename %s\n", filename);
    return 0;
  }

  // Get run
  GSVRun *run = scene->Run(run_name);
  if (!run) { 
    fprintf(stderr, "Bad run %s in sift filename %s\n", run_name, filename);
    return 0;
  }

  // Get segment index
  if ((segment_index < 0) || (segment_index >= run->NSegments())) { 
    fprintf(stderr, "Bad segment index %d in sift filename %s\n", segment_index, filename);
    return 0;
  }

  // Get segment
  GSVSegment *segment = run->Segment(segment_index);
  if (!segment) {
    fprintf(stderr, "Bad segment %d in sift filename %s\n", segment_index, filename);
    return 0;
  }

  // Get panorama index
  if ((panorama_index < 0) || (panorama_index >= segment->NPanoramas())) { 
    fprintf(stderr, "Bad panorama index %d in sift filename %s\n", panorama_index, filename);
    return 0;
  }

  // Get panorama
  GSVPanorama *panorama = segment->Panorama(panorama_index);
  if (!panorama) {
    fprintf(stderr, "Bad panorama %d in sift filename %s\n", panorama_index, filename);
    return 0;
  }

  // Get image index
  if ((image_index < 0) || (image_index >= panorama->NImages())) {
    fprintf(stderr, "Bad image index %d in sift filename %s\n", image_index, filename);
    return 0;
  }

  // Get image
  GSVImage *image = panorama->Image(image_index);
  if (!image) {
    fprintf(stderr, "Bad image %d in sift filename %s\n", image_index, filename);
    return 0;
  }

  // Read sift file
  return ReadSiftFile(image, filename);
}




int GSVPoseOptimization::
ReadSiftFile(GSVImage *image, const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }

  // Read header
  int nfeatures, descriptor_size;
  if (fscanf(fp, "%d%d", &nfeatures, &descriptor_size) != (unsigned int) 2) {
    fprintf(stderr, "Unable to read %s\n", filename);
    return 0;
  }

  // Read features
  double px, py, scale, orientation;
  for (int i = 0; i < nfeatures; i++) {
    // Read feature info
    if (fscanf(fp, "%lf%lf%lf%lf", &py, &px, &scale, &orientation) != (unsigned int) 4) {
      fprintf(stderr, "Unable to read feature %d from %s\n", i, filename);
      return 0;
    }

    // Apply scaling
    px /= sift_image_scale;
    py /= sift_image_scale;
    scale /= sift_image_scale;

    // Invert y coordinate
    py = image->Height() - py;

    // Create feature
    RNScalar score = scale;
    R2Point position(px, py);
    R2Vector direction(1, 0);
    direction.Rotate(orientation);
    GSVFeature *feature = new GSVFeature(GSV_IMAGE_POINT_FEATURE_TYPE, image, position, direction, scale, score);
    if (!feature) {
      fprintf(stderr, "Unable to create feature %d for %s\n", i, filename);
      return 0;
    }
      
    // Read feature descriptor
    if (descriptor_size > 0) {
      GSVDescriptor& descriptor = feature->descriptor;
      descriptor.descriptor_type = GSV_SIFT_DESCRIPTOR_TYPE;
      descriptor.nvalues = descriptor_size;
      descriptor.values = new RNScalar [ descriptor_size ];
      for (int j = 0; j < descriptor_size; j++) {
        if (fscanf(fp, "%lf", &descriptor.values[j]) != (unsigned int) 1) {
          fprintf(stderr, "Unable to read descriptor value %d for feature %d in %s\n", j, i, filename);
          return 0;
        }
      }
    }

    // Insert feature
    InsertFeature(feature);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Line Feature Input 
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
ReadLineFile(const char *filename)
{
  // Parse filename
  char run_buffer[4096];
  strncpy(run_buffer, filename, 4096);
  char *file_name_start = run_buffer;
  char *run_name_end = strrchr(run_buffer, '/');
  if (run_name_end) {
    file_name_start = run_name_end + 1;
    *run_name_end = '\0';
  }
  else {
    fprintf(stderr, "Unable to parse run directory from line filename %s\n", filename); 
    return 0; 
  }

  // Get run name
  printf("%s\n", run_buffer);
  char *run_name = strrchr(run_buffer, '/');
  if (run_name) run_name++;
  else run_name = run_buffer;
  if (!run_name) {
    fprintf(stderr, "Bad run name in line filename %s\n", filename); 
    return 0; 
  }

  // Parse segment, panorama, and image indicies from filename
  int segment_index, panorama_index, image_index;
  if (sscanf(file_name_start, "%d_%d_%d_Line.lin",  
    &segment_index, &panorama_index, &image_index) != (unsigned int) 3) {
    fprintf(stderr, "Unable to parse image info from line filename %s\n", filename);
    return 0;
  }

  // Get run
  GSVRun *run = scene->Run(run_name);
  if (!run) { 
    fprintf(stderr, "Bad run %s in line filename %s\n", run_name, filename);
    return 0;
  }

  // Get segment index
  if ((segment_index < 0) || (segment_index >= run->NSegments())) { 
    fprintf(stderr, "Bad segment index %d in line filename %s\n", segment_index, filename);
    return 0;
  }

  // Get segment
  GSVSegment *segment = run->Segment(segment_index);
  if (!segment) {
    fprintf(stderr, "Bad segment %d in line filename %s\n", segment_index, filename);
    return 0;
  }

  // Get panorama index
  if ((panorama_index < 0) || (panorama_index >= segment->NPanoramas())) { 
    fprintf(stderr, "Bad panorama index %d in line filename %s\n", panorama_index, filename);
    return 0;
  }

  // Get panorama
  GSVPanorama *panorama = segment->Panorama(panorama_index);
  if (!panorama) {
    fprintf(stderr, "Bad panorama %d in line filename %s\n", panorama_index, filename);
    return 0;
  }

  // Get image index
  if ((image_index < 0) || (image_index >= panorama->NImages())) {
    fprintf(stderr, "Bad image index %d in line filename %s\n", image_index, filename);
    return 0;
  }

  // Get image
  GSVImage *image = panorama->Image(image_index);
  if (!image) {
    fprintf(stderr, "Bad image %d in line filename %s\n", image_index, filename);
    return 0;
  }

  // Read line file
  return ReadLineFile(image, filename);
}




int GSVPoseOptimization::
ReadLineFile(GSVImage *image, const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }

  // Read header
  int nfeatures, descriptor_size;
  if (fscanf(fp, "%d%d", &nfeatures, &descriptor_size) != (unsigned int) 2) {
    fprintf(stderr, "Unable to read %s\n", filename);
    return 0;
  }

  // Read features
  double px1, py1, px2, py2, length, score;
  for (int i = 0; i < nfeatures; i++) {
    // Read feature info
    if (fscanf(fp, "%lf%lf%lf%lf%lf%lf", &px1, &py1, &px2, &py2, &length, &score) != (unsigned int) 6) {
      fprintf(stderr, "Unable to read feature %d from %s\n", i, filename);
      return 0;
    }

    // Read feature descriptor
    GSVDescriptor descriptor;
    if (descriptor_size > 0) {
      descriptor.descriptor_type = GSV_LINE_DESCRIPTOR_TYPE;
      descriptor.nvalues = descriptor_size;
      descriptor.values = new RNScalar [ descriptor_size ];
      for (int j = 0; j < descriptor_size; j++) {
        if (fscanf(fp, "%lf", &descriptor.values[j]) != (unsigned int) 1) {
          fprintf(stderr, "Unable to read descriptor value %d for feature %d in %s\n", j, i, filename);
          return 0;
        }
      }
    }

    // Get useful variables
    R2Point p1(px1, py1);
    R2Point p2(px2, py2);
    R2Vector direction = p2 - p1;
    direction.Normalize();
    RNScalar scale = length;

    // Create feature1
    GSVFeature *feature1 = new GSVFeature(GSV_IMAGE_LINE_FEATURE_TYPE, image, p1, direction, scale, score);
    if (!feature1) {
      fprintf(stderr, "Unable to create line feature %d for %s\n", i, filename);
      return 0;
    }
      
    // Create feature2
    GSVFeature *feature2 = new GSVFeature(GSV_IMAGE_LINE_FEATURE_TYPE, image, p2, -direction, scale, score);
    if (!feature2) {
      fprintf(stderr, "Unable to create line feature %d for %s\n", i, filename);
      return 0;
    }
      
    // Assign feature descriptor
    feature1->descriptor = descriptor;
    feature2->descriptor = descriptor;
      
    // Insert features
    InsertFeature(feature1);
    InsertFeature(feature2);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Transformations Input and Output
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
ReadTransformationsFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }

  // Read header
  int nlasers, ncameras, nvertices, dummy;
  if (fscanf(fp, "%d%d%d%d%d%d%d%d\n", &nlasers, &ncameras, &nvertices, 
    &dummy, &dummy, &dummy, &dummy, &dummy) != (unsigned int) 8) {
    fprintf(stderr, "Unable to read %s\n", filename);
    return 0;
  }

  // Read laser variables 
  for (int i = 0; i < lasers.NEntries(); i++) {
    GSVLaser *laser = lasers.Kth(i);
    LaserData *laser_data = (LaserData *) laser->Data();
    if (!laser_data) continue;

    // Read transformation
    double tx, ty, tz, rx, ry, rz;
    if (fscanf(fp, "%lf%lf%lf%lf%lf%lf", &tx, &ty, &tz, &rx, &ry, &rz) != (unsigned int) 6) {
      fprintf(stderr, "Unable to read transformation %d from %s\n", i, filename);
      return 0;
    }

    // Assign transformation
    laser_data->translation.Reset(tx, ty, tz);
    laser_data->rotation.Reset(rx, ry, rz);
  }

  // Read camera variables 
  for (int i = 0; i < cameras.NEntries(); i++) {
    GSVCamera *camera = cameras.Kth(i);
    CameraData *camera_data = (CameraData *) camera->Data();
    if (!camera_data) continue;

    // Read transformation
    double tx, ty, tz, rx, ry, rz;
    if (fscanf(fp, "%lf%lf%lf%lf%lf%lf", &tx, &ty, &tz, &rx, &ry, &rz) != (unsigned int) 6) {
      fprintf(stderr, "Unable to read transformation %d from %s\n", i, filename);
      return 0;
    }

    // Assign transformation
    camera_data->translation.Reset(tx, ty, tz);
    camera_data->rotation.Reset(rx, ry, rz);
  }

  // Read path vertex variables
  for (int i = 0; i < vertices.NEntries(); i++) {
    GSVPathVertex *vertex = vertices.Kth(i);

    // Read transformation
    double tx, ty, tz, rx, ry, rz;
    if (fscanf(fp, "%lf%lf%lf%lf%lf%lf", &tx, &ty, &tz, &rx, &ry, &rz) != (unsigned int) 6) {
      fprintf(stderr, "Unable to read transformation %d from %s\n", i, filename);
      return 0;
    }

    // Assign transformation
    vertex->translation.Reset(tx, ty, tz);
    vertex->rotation.Reset(rx, ry, rz);
  }

  // Close file
  fclose(fp);

  // Update pixel positions
  UpdatePixelPositionsFromRayIntersections(); 
  UpdateClusterPositionsFromFeaturePositions();  

  // Return success
  return 1;
}



int GSVPoseOptimization::
WriteTransformationsFile(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }

  // Write header
  fprintf(fp, "%d %d %d 0 0 0 0 0\n", lasers.NEntries(), cameras.NEntries(), vertices.NEntries());

  // Write laser variables to file
  for (int i = 0; i < lasers.NEntries(); i++) {
    GSVLaser *laser = lasers.Kth(i);
    LaserData *laser_data = (LaserData *) laser->Data();
    if (!laser_data) continue;
    const R3Vector& translation = laser_data->translation;
    const R3Vector& rotation = laser_data->rotation;
    fprintf(fp, "%15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n", 
            translation.X(), translation.Y(), translation.Z(),
            rotation.X(), rotation.Y(), rotation.Z());
  }

  // Write camera variables to file
  for (int i = 0; i < cameras.NEntries(); i++) {
    GSVCamera *camera = cameras.Kth(i);
    CameraData *camera_data = (CameraData *) camera->Data();
    if (!camera_data) continue;
    const R3Vector& translation = camera_data->translation;
    const R3Vector& rotation = camera_data->rotation;
    fprintf(fp, "%15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n", 
            translation.X(), translation.Y(), translation.Z(),
            rotation.X(), rotation.Y(), rotation.Z());
  }

  // Write path vertex variables to file
  for (int i = 0; i < vertices.NEntries(); i++) {
    GSVPathVertex *vertex = vertices.Kth(i);
    R3Vector translation = vertex->translation;
    R3Vector rotation = vertex->rotation;
    fprintf(fp, "%15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n", 
            translation.X(), translation.Y(), translation.Z(),
            rotation.X(), rotation.Y(), rotation.Z());
  }

  // Close file
  fclose(fp);
  
  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Inertias Input and Output
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
ReadInertiaFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }

  // Read header
  int nlasers, ncameras, nvertices, dummy;
  if (fscanf(fp, "%d%d%d%d%d%d%d%d\n", &nlasers, &ncameras, &nvertices, 
    &dummy, &dummy, &dummy, &dummy, &dummy) != (unsigned int) 8) {
    fprintf(stderr, "Unable to read %s\n", filename);
    return 0;
  }

  // Read laser variables 
  for (int i = 0; i < lasers.NEntries(); i++) {
    GSVLaser *laser = lasers.Kth(i);
    LaserData *laser_data = (LaserData *) laser->Data();
    if (!laser_data) continue;

    // Read inertias
    double tx, ty, tz, rx, ry, rz;
    if (fscanf(fp, "%lf%lf%lf%lf%lf%lf", &tx, &ty, &tz, &rx, &ry, &rz) != (unsigned int) 6) {
      fprintf(stderr, "Unable to read inertia %d from %s\n", i, filename);
      return 0;
    }

    // Assign inertias
    laser_data->inertia[TX] = tx;
    laser_data->inertia[TY] = ty;
    laser_data->inertia[TZ] = tz;
    laser_data->inertia[RX] = rx;
    laser_data->inertia[RY] = ry;
    laser_data->inertia[RZ] = rz;
  }

  // Read camera variables 
  for (int i = 0; i < cameras.NEntries(); i++) {
    GSVCamera *camera = cameras.Kth(i);
    CameraData *camera_data = (CameraData *) camera->Data();
    if (!camera_data) continue;

    // Read inertia
    double tx, ty, tz, rx, ry, rz;
    if (fscanf(fp, "%lf%lf%lf%lf%lf%lf", &tx, &ty, &tz, &rx, &ry, &rz) != (unsigned int) 6) {
      fprintf(stderr, "Unable to read inertia %d from %s\n", i, filename);
      return 0;
    }

    // Assign inertias
    camera_data->inertia[TX] = tx;
    camera_data->inertia[TY] = ty;
    camera_data->inertia[TZ] = tz;
    camera_data->inertia[RX] = rx;
    camera_data->inertia[RY] = ry;
    camera_data->inertia[RZ] = rz;
  }

  // Read path vertex variables
  for (int i = 0; i < vertices.NEntries(); i++) {
    GSVPathVertex *vertex = vertices.Kth(i);

    // Read inertia
    double tx, ty, tz, rx, ry, rz;
    if (fscanf(fp, "%lf%lf%lf%lf%lf%lf", &tx, &ty, &tz, &rx, &ry, &rz) != (unsigned int) 6) {
      fprintf(stderr, "Unable to read inertia %d from %s\n", i, filename);
      return 0;
    }

    // Assign inertia
    vertex->inertia[TX] = tx;
    vertex->inertia[TY] = ty;
    vertex->inertia[TZ] = tz;
    vertex->inertia[RX] = rx;
    vertex->inertia[RY] = ry;
    vertex->inertia[RZ] = rz;
  }

  // Close file
  fclose(fp);

  // Update cluster inertias
  UpdateClusterInertiasFromFeatureInertias();  

  // Return success
  return 1;
}



int GSVPoseOptimization::
WriteInertiaFile(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open %s\n", filename);
    return 0;
  }

  // Write header
  fprintf(fp, "%d %d %d 0 0 0 0 0\n", lasers.NEntries(), cameras.NEntries(), vertices.NEntries());

  // Write laser variables to file
  for (int i = 0; i < lasers.NEntries(); i++) {
    GSVLaser *laser = lasers.Kth(i);
    LaserData *laser_data = (LaserData *) laser->Data();
    if (!laser_data) continue;
    for (int j = 0; j < 6; j++) 
      fprintf(fp, "%15.9f ", laser_data->inertia[j]);
    fprintf(fp, "\n");
  }

  // Write camera variables to file
  for (int i = 0; i < cameras.NEntries(); i++) {
    GSVCamera *camera = cameras.Kth(i);
    CameraData *camera_data = (CameraData *) camera->Data();
    if (!camera_data) continue;
    for (int j = 0; j < 6; j++) 
      fprintf(fp, "%15.9f ", camera_data->inertia[j]);
    fprintf(fp, "\n");
  }

  // Write path vertex variables to file
  for (int i = 0; i < vertices.NEntries(); i++) {
    GSVPathVertex *vertex = vertices.Kth(i);
    for (int j = 0; j < 6; j++) 
      fprintf(fp, "%15.9f ", vertex->inertia[j]);
    fprintf(fp, "\n");
  }

  // Close file
  fclose(fp);
  
  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Options Input and Output
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
ReadOptionsFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open options file %s\n", filename);
    return 0;
  }

  // Read file
  int status = ReadOptions(fp);

  // Close file
  fclose(fp);

  // Return status
  return status;
}



int GSVPoseOptimization::
WriteOptionsFile(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open options file %s\n", filename);
    return 0;
  }

  // Write file
  int status = WriteOptions(fp);

  // Close file
  fclose(fp);

  // Return status
  return status;
}



int GSVPoseOptimization::
ParseOption(const char *keyword, const char *value) 
{
  // Check keyword
  if (!strcmp(keyword, "#")) return 1;
    
  // Pair creation options
  else if (!strcmp(keyword, "max_pair_world_distance"))
    max_pair_world_distance = atof(value);
  else if (!strcmp(keyword, "max_pair_world_distance_ratio"))
    max_pair_world_distance_ratio = atof(value);
  else if (!strcmp(keyword, "max_pair_image_distance"))
    max_pair_image_distance = atof(value);
  else if (!strcmp(keyword, "max_pair_world_direction_angle"))
    max_pair_world_direction_angle = atof(value);
  else if (!strcmp(keyword, "max_pair_world_normal_angle"))
    max_pair_world_normal_angle = atof(value);
  else if (!strcmp(keyword, "max_pair_image_distance_ratio"))
    max_pair_image_distance_ratio = atof(value);
  else if (!strcmp(keyword, "max_pair_image_direction_angle"))
    max_pair_image_direction_angle = atof(value);
  else if (!strcmp(keyword, "max_pair_spin_image_descriptor_distance"))
    max_pair_spin_image_descriptor_distance = atof(value);
  else if (!strcmp(keyword, "max_pair_shape_context_descriptor_distance"))
    max_pair_shape_context_descriptor_distance = atof(value);
  else if (!strcmp(keyword, "max_pair_sift_descriptor_distance"))
    max_pair_sift_descriptor_distance = atof(value);
  else if (!strcmp(keyword, "max_pair_line_descriptor_distance"))
    max_pair_line_descriptor_distance = atof(value);
  else if (!strcmp(keyword, "max_pair_descriptor_distance_ratio"))
    max_pair_descriptor_distance_ratio = atof(value);
  else if (!strcmp(keyword, "min_pair_path_parameter_difference"))
    min_pair_path_parameter_difference = atof(value);

  // Correspondence creation options
  else if (!strcmp(keyword, "max_correspondence_world_distance"))
    max_correspondence_world_distance = atof(value);
  else if (!strcmp(keyword, "max_correspondence_world_distance_ratio"))
    max_correspondence_world_distance_ratio = atof(value);
  else if (!strcmp(keyword, "max_correspondence_world_direction_angle"))
    max_correspondence_world_direction_angle = atof(value);
  else if (!strcmp(keyword, "max_correspondence_world_normal_angle"))
    max_correspondence_world_normal_angle = atof(value);
  else if (!strcmp(keyword, "max_correspondence_image_distance"))
    max_correspondence_image_distance = atof(value);
  else if (!strcmp(keyword, "max_correspondence_image_distance_ratio"))
    max_correspondence_image_distance_ratio = atof(value);
  else if (!strcmp(keyword, "max_correspondence_image_direction_angle"))
    max_correspondence_image_direction_angle = atof(value);
  else if (!strcmp(keyword, "max_correspondence_spin_image_descriptor_distance"))
    max_correspondence_spin_image_descriptor_distance = atof(value);
  else if (!strcmp(keyword, "max_correspondence_shape_context_descriptor_distance"))
    max_correspondence_shape_context_descriptor_distance = atof(value);
  else if (!strcmp(keyword, "max_correspondence_sift_descriptor_distance"))
    max_correspondence_sift_descriptor_distance = atof(value);
  else if (!strcmp(keyword, "max_correspondence_line_descriptor_distance"))
    max_correspondence_line_descriptor_distance = atof(value);
  else if (!strcmp(keyword, "max_correspondence_descriptor_distance_ratio"))
    max_correspondence_descriptor_distance_ratio = atof(value);
  else if (!strcmp(keyword, "min_correspondence_path_parameter_difference"))
    min_correspondence_path_parameter_difference = atof(value);
  else if (!strcmp(keyword, "icp_max_world_distance_start"))
    icp_max_world_distance_start = atof(value);
  else if (!strcmp(keyword, "icp_max_world_distance_end"))
    icp_max_world_distance_end = atof(value);
  else if (!strcmp(keyword, "icp_max_image_distance_start"))
    icp_max_image_distance_start = atof(value);
  else if (!strcmp(keyword, "icp_max_image_distance_end"))
    icp_max_image_distance_end = atof(value);

  // Optimization options
  else if (!strcmp(keyword, "expression_type"))
    expression_type = atoi(value);
  else if (!strcmp(keyword, "scan_residual_threshold"))
    scan_residual_threshold = atof(value);
  else if (!strcmp(keyword, "image_residual_threshold")) 
    image_residual_threshold = atof(value);
  else if (!strcmp(keyword, "dot_product_residual_threshold")) 
    dot_product_residual_threshold = atof(value);
  else if (!strcmp(keyword, "laser_inertia_tx_weight"))
    laser_inertia_weights[TX] = atof(value);
  else if (!strcmp(keyword, "laser_inertia_ty_weight"))
    laser_inertia_weights[TY] = atof(value);
  else if (!strcmp(keyword, "laser_inertia_tz_weight"))
    laser_inertia_weights[TZ] = atof(value);
  else if (!strcmp(keyword, "laser_inertia_rx_weight"))
    laser_inertia_weights[RX] = atof(value);
  else if (!strcmp(keyword, "laser_inertia_ry_weight"))
    laser_inertia_weights[RY] = atof(value);
  else if (!strcmp(keyword, "laser_inertia_rz_weight"))
    laser_inertia_weights[RZ] = atof(value);
  else if (!strcmp(keyword, "camera_intertia_tx_weight"))
    camera_inertia_weights[TX] = atof(value);
  else if (!strcmp(keyword, "camera_intertia_ty_weight"))
    camera_inertia_weights[TY] = atof(value);
  else if (!strcmp(keyword, "camera_intertia_tz_weight"))
    camera_inertia_weights[TZ] = atof(value);
  else if (!strcmp(keyword, "camera_intertia_rx_weight"))
    camera_inertia_weights[RX] = atof(value);
  else if (!strcmp(keyword, "camera_intertia_ry_weight"))
    camera_inertia_weights[RY] = atof(value);
  else if (!strcmp(keyword, "camera_intertia_rz_weight"))
    camera_inertia_weights[RZ] = atof(value);
  else if (!strcmp(keyword, "vertex_intertia_tx_weight"))
    vertex_inertia_weights[TX] = atof(value);
  else if (!strcmp(keyword, "vertex_intertia_ty_weight"))
    vertex_inertia_weights[TY] = atof(value);
  else if (!strcmp(keyword, "vertex_intertia_tz_weight"))
    vertex_inertia_weights[TZ] = atof(value);
  else if (!strcmp(keyword, "vertex_intertia_rx_weight"))
    vertex_inertia_weights[RX] = atof(value);
  else if (!strcmp(keyword, "vertex_intertia_ry_weight"))
    vertex_inertia_weights[RY] = atof(value);
  else if (!strcmp(keyword, "vertex_intertia_rz_weight"))
    vertex_inertia_weights[RZ] = atof(value);
  else if (!strcmp(keyword, "pixel_inertia_depth_weight"))
    pixel_feature_inertia_weights[T] = atof(value);
  else if (!strcmp(keyword, "cluster_inertia_tx_weight"))
    cluster_inertia_weights[TX] = atof(value);
  else if (!strcmp(keyword, "cluster_inertia_ty_weight"))
    cluster_inertia_weights[TY] = atof(value);
  else if (!strcmp(keyword, "cluster_inertia_tz_weight"))
    cluster_inertia_weights[TZ] = atof(value);
  else if (!strcmp(keyword, "cluster_inertia_rx_weight"))
    cluster_inertia_weights[RX] = atof(value);
  else if (!strcmp(keyword, "cluster_inertia_ry_weight"))
    cluster_inertia_weights[RY] = atof(value);
  else if (!strcmp(keyword, "cluster_inertia_rz_weight"))
    cluster_inertia_weights[RZ] = atof(value);
  else if (!strcmp(keyword, "path_rigidity_radius")) 
    path_rigidity_radius = atoi(value);
  else if (!strcmp(keyword, "path_rigidity_sigma"))
    path_rigidity_sigma = atof(value);      
  else if (!strcmp(keyword, "path_rigidity_weight"))
    path_rigidity_weight = atof(value);
  else if (!strcmp(keyword, "camera_camera_rigidity_weigh"))
    camera_camera_rigidity_weight = atof(value);
  else if (!strcmp(keyword, "laser_laser_rigidity_weight"))
    laser_laser_rigidity_weight = atof(value);
  else if (!strcmp(keyword, "camera_laser_rigidity_weight"))
    camera_laser_rigidity_weight = atof(value);
  else if (!strcmp(keyword, "scan_point_scan_point_correspondence_weight"))
    scan_point_scan_point_correspondence_weight = atof(value);
  else if (!strcmp(keyword, "scan_point_scan_line_correspondence_weight"))
    scan_point_scan_line_correspondence_weight = atof(value);
  else if (!strcmp(keyword, "scan_point_scan_plane_correspondence_weight"))
    scan_point_scan_plane_correspondence_weight = atof(value);
  else if (!strcmp(keyword, "scan_line_scan_line_correspondence_weight"))
    scan_line_scan_line_correspondence_weight = atof(value);
  else if (!strcmp(keyword, "scan_line_scan_plane_correspondence_weight"))
    scan_line_scan_plane_correspondence_weight = atof(value);
  else if (!strcmp(keyword, "scan_point_image_point_correspondence_weight"))
    scan_point_image_point_correspondence_weight = atof(value);
  else if (!strcmp(keyword, "scan_line_image_point_correspondence_weight"))
    scan_line_image_point_correspondence_weight = atof(value);
  else if (!strcmp(keyword, "scan_plane_image_point_correspondence_weight"))
    scan_plane_image_point_correspondence_weight = atof(value);
  else if (!strcmp(keyword, "image_point_image_point_correspondence_weight"))
    image_point_image_point_correspondence_weight = atof(value);
  else if (!strcmp(keyword, "image_point_image_line_correspondence_weight"))
    image_point_image_line_correspondence_weight = atof(value);
  else if (!strcmp(keyword, "image_line_image_line_correspondence_weight"))
    image_line_image_line_correspondence_weight = atof(value);
  else if (!strcmp(keyword, "scan_point_image_point_correspondence_reprojection_weight"))
    scan_point_image_point_correspondence_reprojection_weight = atof(value);
  else if (!strcmp(keyword, "scan_line_image_point_correspondence_reprojection_weight"))
    scan_line_image_point_correspondence_reprojection_weight = atof(value);
  else if (!strcmp(keyword, "scan_plane_image_point_correspondence_reprojection_weight"))
    scan_plane_image_point_correspondence_reprojection_weight = atof(value);
  else if (!strcmp(keyword, "image_point_image_point_correspondence_reprojection_weight"))
    image_point_image_point_correspondence_reprojection_weight = atof(value);
  else if (!strcmp(keyword, "image_point_image_line_correspondence_reprojection_weight"))
    image_point_image_line_correspondence_reprojection_weight = atof(value);
  else if (!strcmp(keyword, "image_line_image_line_correspondence_reprojection_weight"))
    image_line_image_line_correspondence_reprojection_weight = atof(value);
  else if (!strcmp(keyword, "scan_direction_scan_direction_correspondence_weight"))
    scan_direction_scan_direction_correspondence_weight = atof(value);

  // Unrecognized option
  else return 0;

  // Return success
  return 1;
}



int GSVPoseOptimization::
ReadOptions(FILE *fp)
{
  // Check file
  if (!fp) fp = stdout;

  // Parse file
  int status = 1;
  char buffer[4096];
  while (fgets(buffer, 4096, fp)) {
    char *bufferp = buffer;
    if (*bufferp == '#') continue;
    if (*bufferp == '\0') continue;
    char *keyword = strtok(bufferp, " \t\n=");
    if (!keyword) continue;
    char *value = strtok(NULL, " \t\n=");
    if (!value) continue;

    // Check keyword
    if (!strcmp(keyword, "#")) continue;
    
    // Parse keyword-value pair
    if (!ParseOption(keyword, value)) {
      fprintf(stderr, "Unrecognized keyword %s in options file\n", keyword);
      return 0;
    }
  }

  // Return whether everything was parsed
  return status;
}



int GSVPoseOptimization::
WriteOptions(FILE *fp) const
{
  // Check file
  if (!fp) fp = stdout;

  // Path parameters
  fprintf(fp, "# Path creation parameters\n");
  fprintf(fp, "max_path_vertex_spacing = %g\n", max_path_vertex_spacing);
  fprintf(fp, "\n");

  // Pair creation parameters
  fprintf(fp, "# Pair creation parameters\n");
  fprintf(fp, "max_pair_world_distance = %g\n", max_pair_world_distance);
  fprintf(fp, "max_pair_world_distance_ratio = %g\n", max_pair_world_distance_ratio);
  fprintf(fp, "max_pair_world_direction_angle = %g\n", max_pair_world_direction_angle);
  fprintf(fp, "max_pair_world_normal_angle = %g\n", max_pair_world_normal_angle);
  fprintf(fp, "max_pair_image_distance = %g\n", max_pair_image_distance);
  fprintf(fp, "max_pair_image_distance_ratio = %g\n", max_pair_image_distance_ratio);
  fprintf(fp, "max_pair_image_direction_angle = %g\n", max_pair_image_direction_angle);
  fprintf(fp, "max_pair_spin_image_descriptor_distance = %g\n", max_pair_spin_image_descriptor_distance);
  fprintf(fp, "max_pair_shape_context_descriptor_distance = %g\n", max_pair_shape_context_descriptor_distance);
  fprintf(fp, "max_pair_sift_descriptor_distance = %g\n", max_pair_sift_descriptor_distance);
  fprintf(fp, "max_pair_line_descriptor_distance = %g\n", max_pair_line_descriptor_distance);
  fprintf(fp, "max_pair_descriptor_distance_ratio = %g\n", max_pair_descriptor_distance_ratio);
  fprintf(fp, "min_pair_path_parameter_difference = %g\n", min_pair_path_parameter_difference);
  fprintf(fp, "\n");

  // Correspondence creation options
  fprintf(fp, "# Correspondence creation parameters\n");
  fprintf(fp, "max_correspondence_world_distance = %g\n", max_correspondence_world_distance);
  fprintf(fp, "max_correspondence_world_distance_ratio = %g\n", max_correspondence_world_distance_ratio);
  fprintf(fp, "max_correspondence_world_direction_angle = %g\n", max_correspondence_world_direction_angle);
  fprintf(fp, "max_correspondence_world_normal_angle = %g\n", max_correspondence_world_normal_angle);
  fprintf(fp, "max_correspondence_image_distance = %g\n", max_correspondence_image_distance);
  fprintf(fp, "max_correspondence_image_distance_ratio = %g\n", max_correspondence_image_distance_ratio);
  fprintf(fp, "max_correspondence_image_direction_angle = %g\n", max_correspondence_image_direction_angle);
  fprintf(fp, "max_correspondence_spin_image_descriptor_distance = %g\n", max_correspondence_spin_image_descriptor_distance);
  fprintf(fp, "max_correspondence_shape_context_descriptor_distance = %g\n", max_correspondence_shape_context_descriptor_distance);
  fprintf(fp, "max_correspondence_sift_descriptor_distance = %g\n", max_correspondence_sift_descriptor_distance);
  fprintf(fp, "max_correspondence_line_descriptor_distance = %g\n", max_correspondence_line_descriptor_distance);
  fprintf(fp, "max_correspondence_descriptor_distance_ratio = %g\n", max_correspondence_descriptor_distance_ratio);
  fprintf(fp, "min_correspondence_path_parameter_difference = %g\n", min_correspondence_path_parameter_difference);
  fprintf(fp, "icp_max_world_distance_start = %g\n", icp_max_world_distance_start);
  fprintf(fp, "icp_max_world_distance_end = %g\n", icp_max_world_distance_end);
  fprintf(fp, "icp_max_image_distance_start = %g\n", icp_max_image_distance_start);
  fprintf(fp, "icp_max_image_distance_end = %g\n", icp_max_image_distance_end);
  fprintf(fp, "\n");

  // Optimization options
  fprintf(fp, "expression_type = %d\n", expression_type);
  fprintf(fp, "scan_residual_threshold = %g\n", scan_residual_threshold);
  fprintf(fp, "image_residual_threshold = %g\n", image_residual_threshold);
  fprintf(fp, "dot_product_residual_threshold = %g\n", dot_product_residual_threshold);
  fprintf(fp, "laser_inertia_tx_weight = %g\n", laser_inertia_weights[TX]);
  fprintf(fp, "laser_inertia_ty_weight = %g\n", laser_inertia_weights[TY]);
  fprintf(fp, "laser_inertia_tz_weight = %g\n", laser_inertia_weights[TZ]);
  fprintf(fp, "laser_inertia_rx_weight = %g\n", laser_inertia_weights[RX]);
  fprintf(fp, "laser_inertia_ry_weight = %g\n", laser_inertia_weights[RY]);
  fprintf(fp, "laser_inertia_rz_weight = %g\n", laser_inertia_weights[RZ]);
  fprintf(fp, "camera_intertia_tx_weight = %g\n", camera_inertia_weights[TX]);
  fprintf(fp, "camera_intertia_ty_weight = %g\n", camera_inertia_weights[TY]);
  fprintf(fp, "camera_intertia_tz_weight = %g\n", camera_inertia_weights[TZ]);
  fprintf(fp, "camera_intertia_rx_weight = %g\n", camera_inertia_weights[RX]);
  fprintf(fp, "camera_intertia_ry_weight = %g\n", camera_inertia_weights[RY]);
  fprintf(fp, "camera_intertia_rz_weight = %g\n", camera_inertia_weights[RZ]);
  fprintf(fp, "vertex_intertia_tx_weight = %g\n", vertex_inertia_weights[TX]);
  fprintf(fp, "vertex_intertia_ty_weight = %g\n", vertex_inertia_weights[TY]);
  fprintf(fp, "vertex_intertia_tz_weight = %g\n", vertex_inertia_weights[TZ]);
  fprintf(fp, "vertex_intertia_rx_weight = %g\n", vertex_inertia_weights[RX]);
  fprintf(fp, "vertex_intertia_ry_weight = %g\n", vertex_inertia_weights[RY]);
  fprintf(fp, "vertex_intertia_rz_weight = %g\n", vertex_inertia_weights[RZ]);
  fprintf(fp, "pixel_inertia_depth_weight = %g\n", pixel_feature_inertia_weights[T]);
  fprintf(fp, "cluster_inertia_tx_weight = %g\n", cluster_inertia_weights[TX]);
  fprintf(fp, "cluster_inertia_ty_weight = %g\n", cluster_inertia_weights[TY]);
  fprintf(fp, "cluster_inertia_tz_weight = %g\n", cluster_inertia_weights[TZ]);
  fprintf(fp, "cluster_inertia_rx_weight = %g\n", cluster_inertia_weights[RX]);
  fprintf(fp, "cluster_inertia_ry_weight = %g\n", cluster_inertia_weights[RY]);
  fprintf(fp, "cluster_inertia_rz_weight = %g\n", cluster_inertia_weights[RZ]);
  fprintf(fp, "path_rigidity_radius = %d\n", path_rigidity_radius);
  fprintf(fp, "path_rigidity_sigma = %g\n", path_rigidity_sigma);      
  fprintf(fp, "path_rigidity_weight = %g\n", path_rigidity_weight);
  fprintf(fp, "camera_camera_rigidity_weight = %g\n", camera_camera_rigidity_weight);
  fprintf(fp, "laser_laser_rigidity_weight = %g\n", laser_laser_rigidity_weight);
  fprintf(fp, "camera_laser_rigidity_weight = %g\n", camera_laser_rigidity_weight);
  fprintf(fp, "scan_point_scan_point_correspondence_weight = %g\n", scan_point_scan_point_correspondence_weight);
  fprintf(fp, "scan_point_scan_line_correspondence_weight = %g\n", scan_point_scan_line_correspondence_weight);
  fprintf(fp, "scan_point_scan_plane_correspondence_weight = %g\n", scan_point_scan_plane_correspondence_weight);
  fprintf(fp, "scan_line_scan_line_correspondence_weight = %g\n", scan_line_scan_line_correspondence_weight);
  fprintf(fp, "scan_line_scan_plane_correspondence_weight = %g\n", scan_line_scan_plane_correspondence_weight);
  fprintf(fp, "scan_point_image_point_correspondence_weight = %g\n", scan_point_image_point_correspondence_weight);
  fprintf(fp, "scan_line_image_point_correspondence_weight = %g\n", scan_line_image_point_correspondence_weight);
  fprintf(fp, "scan_plane_image_point_correspondence_weight = %g\n", scan_plane_image_point_correspondence_weight);
  fprintf(fp, "image_point_image_point_correspondence_weight = %g\n", image_point_image_point_correspondence_weight);
  fprintf(fp, "image_point_image_line_correspondence_weight = %g\n", image_point_image_line_correspondence_weight);
  fprintf(fp, "image_line_image_line_correspondence_weight = %g\n", image_line_image_line_correspondence_weight);
  fprintf(fp, "scan_point_image_point_correspondence_reprojection_weight = %g\n", scan_point_image_point_correspondence_reprojection_weight);
  fprintf(fp, "scan_line_image_point_correspondence_reprojection_weight = %g\n", scan_line_image_point_correspondence_reprojection_weight);
  fprintf(fp, "scan_plane_image_point_correspondence_reprojection_weight = %g\n", scan_plane_image_point_correspondence_reprojection_weight);
  fprintf(fp, "image_point_image_point_correspondence_reprojection_weight = %g\n", image_point_image_point_correspondence_reprojection_weight);
  fprintf(fp, "image_point_image_line_correspondence_reprojection_weight = %g\n", image_point_image_line_correspondence_reprojection_weight);
  fprintf(fp, "image_line_image_line_correspondence_reprojection_weight = %g\n", image_line_image_line_correspondence_reprojection_weight);
  fprintf(fp, "scan_direction_scan_direction_correspondence_weight = %g\n", scan_direction_scan_direction_correspondence_weight);
  fprintf(fp, "\n");

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Update utility functions
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
ApplyOptimizedTransformationsToScene(void)
{
  // Check scene
  if (!scene) return 0;

  // Set new pose in every scanline
  for (int i = 0; i < scanlines.NEntries(); i++) {
    GSVScanline *scanline = scanlines.Kth(i);
    GSVPose pose = scanline->Pose();
    pose.Transform(OptimizedTransformation(scanline));
    scanline->SetPose(pose);
  }

  // Set new pose in every image
  for (int i = 0; i < images.NEntries(); i++) {
    GSVImage *image = images.Kth(i);
    GSVPose pose0 = image->Pose(0);
    GSVPose pose1 = image->Pose(image->Width()-1);
    R2Point point0(0, image->Height()/2);
    R2Point point1(image->Width()-1, image->Height()/2);
    pose0.Transform(OptimizedTransformation(image, point0));
    pose1.Transform(OptimizedTransformation(image, point1));
    image->SetPose(pose0, pose1);
  }

  // Clear optimization variables
  ClearOptimizationVariables();

  // Return success
  return 1;
}



int GSVPoseOptimization::
UpdatePixelPositionsFromCorrespondenceSolutions(void) 
{
  // Update scan_position from image_t
  for (int i = 0; i < NFeatures(); i++) {
    GSVFeature *feature = Feature(i);
    if (!feature->image) continue;
    if (feature->image_t <= 0) continue;
    R3Ray ray = feature->image->RayThroughUndistortedPosition(feature->image_position);
    feature->scan_position = ray.Point(feature->image_t);
    if (feature->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) {
      // ??? Should update scan_direction and scale ???
      R2Point pointA = feature->image_position - feature->image_scale * feature->image_direction;
      R2Point pointB = feature->image_position + feature->image_scale * feature->image_direction;
      R3Ray rayA = feature->image->RayThroughUndistortedPosition(pointA);
      R3Ray rayB = feature->image->RayThroughUndistortedPosition(pointB);
      feature->scan_normal = rayA.Vector() % rayB.Vector(); 
      feature->scan_normal.Normalize();
    }
  }

  // Return success
  return 1;
}




static int 
ComputeRayIntersection(GSVPoseOptimization *optimization,
  GSVFeature *feature0, GSVFeature *feature1, 
  R3Point& intersection_point0, R3Point& intersection_point1, 
  RNScalar& intersection_t0, RNScalar& intersection_t1)
{
  // Get transformations
  R3Affine transformation0 = optimization->OptimizedTransformation(feature0);
  R3Affine transformation1 = optimization->OptimizedTransformation(feature1);

  // Compute interesection
  if ((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) {
    R3Point point0 = feature0->scan_position;
    R3Ray ray1 = feature1->image->RayThroughUndistortedPosition(feature1->image_position);
    point0.Transform(transformation0);
    ray1.Transform(transformation1);
    RNScalar t = ray1.T(point0);
    if (RNIsPositive(t)) {
      intersection_point0 = point0;
      intersection_point1 = ray1.Point(t);
      intersection_t0 = 0;
      intersection_t1 = t;
      return 1;
    }
  }
  else if ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) {
    R3Ray ray0 = feature0->image->RayThroughUndistortedPosition(feature0->image_position);
    R3Point point1 = feature1->scan_position;
    ray0.Transform(transformation0);
    point1.Transform(transformation1);
    RNScalar t = ray0.T(point1);
    if (RNIsPositive(t)) {
      intersection_point0 = ray0.Point(t);
      intersection_point1 = point1;
      intersection_t0 = t;
      intersection_t1 = 0;
      return 1;
    }
  }
  else if ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) {
    R3Ray ray0 = feature0->image->RayThroughUndistortedPosition(feature0->image_position);
    R3Ray ray1 = feature1->image->RayThroughUndistortedPosition(feature1->image_position);
    ray0.Transform(transformation0);
    ray1.Transform(transformation1);
    const R3Vector& v0 = ray0.Vector();
    const R3Vector& v1 = ray1.Vector();
    RNScalar v0v0 = 1.0;  // v0.Dot(v0);
    RNScalar v1v1 = 1.0;  // v1.Dot(v1);
    RNScalar v0v1 = v0.Dot(v1);
    RNScalar denom = v0v1*v0v1 - v0v0*v1v1;
    if (RNIsNotZero(denom)) {
      R3Vector p0 = ray0.Start().Vector();
      R3Vector p1 = ray1.Start().Vector();
      RNScalar p0v0 = v0.Dot(p0);
      RNScalar p1v1 = v1.Dot(p1);
      RNScalar p0v1 = v1.Dot(p0);
      RNScalar p1v0 = v0.Dot(p1);
      RNScalar t0 = (v0v1*p1v1 + v1v1*p0v0 - v0v1*p0v1 - v1v1*p1v0) / denom;
      RNScalar t1 = (v0v1*p0v0 + v0v0*p1v1 - v0v1*p1v0 - v0v0*p0v1) / denom;
      if (RNIsPositive(t0) && RNIsPositive(t1)) {
        intersection_point0 = ray0.Point(t0);
        intersection_point1 = ray1.Point(t1); 
        intersection_t0 = t0;
        intersection_t1 = t1;
        return 1;
      }
    }
  }

  // No intersection returned
  return 0;
}



int GSVPoseOptimization::
UpdatePixelPositionsFromRayIntersections(void) 
{
  // Allocate data
  RNScalar *depths = new RNScalar [ NFeatures() ];
  RNScalar *weights = new RNScalar [ NFeatures() ];
  for (int i = 0; i < NFeatures(); i++) {
    depths[i] = 0.0;
    weights[i] = 0.0;
  }

  // Update scan positions for image features based on correspondences
  for (int i = 0; i < NCorrespondences(); i++) {
    GSVFeatureCorrespondence *correspondence = Correspondence(i);
    GSVFeature *feature0 = correspondence->Feature(0);
    GSVFeature *feature1 = correspondence->Feature(1);
    if (!feature0 || !feature1) continue;
    if ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) ||
        (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) {
      RNScalar t0 = 0, t1 = 0;
      R3Point intersection0, intersection1;
      if (ComputeRayIntersection(this, feature0, feature1, intersection0, intersection1, t0, t1)) {
        if (feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) {
          depths[feature0->index] += t0;
          weights[feature0->index] += 1.0;
        }
        if (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) {
          depths[feature1->index] += t1;
          weights[feature1->index] += 1.0;
        }
      }
    }
  }

  // Update scan positions for image features based on clusters
  for (int i = 0; i < NClusters(); i++) {
    GSVFeatureCluster *cluster = Cluster(i);
    R3Point cluster_position(0,0,0);
    RNScalar cluster_weight = 0.0;
    for (int j = 0; j < cluster->NFeatures(); j++) {
      GSVFeature *feature0 = cluster->Feature(j);
      for (int k = j+1; k < cluster->NFeatures(); k++) {
        GSVFeature *feature1 = cluster->Feature(k);
        if ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) ||
            (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) {
          RNScalar t0, t1;
          R3Point intersection0, intersection1;
          if (ComputeRayIntersection(this, feature0, feature1, intersection0, intersection1, t0, t1)) {
            cluster_position += intersection0;
            cluster_position += intersection1;
            cluster_weight += 2;
          }
        }
      }
    }
    if (cluster_weight > 0) {
      cluster_position /= cluster_weight;
      for (int j = 0; j < cluster->NFeatures(); j++) {
        GSVFeature *feature = cluster->Feature(j);
        if (feature->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) {
          R3Point p = cluster_position;
          p.InverseTransform(OptimizedTransformation(feature));
          R3Ray ray = feature->image->RayThroughUndistortedPosition(feature->image_position);
          depths[feature->index] += cluster_weight * ray.T(p);
          weights[feature->index] += cluster_weight;
        }
      }
    }
  }

  // Divide by sums to get means
  for (int i = 0; i < NFeatures(); i++) {
    GSVFeature *feature = Feature(i);
    if (weights[i] == 0.0) continue;
    feature->image_t = depths[i] / weights[i];
    R3Ray ray = feature->image->RayThroughUndistortedPosition(feature->image_position);
    feature->scan_position = ray.Point(feature->image_t);
  }

  // Delete data
  delete [] depths;
  delete [] weights;

  // Return success
  return 1;
}



int GSVPoseOptimization::
UpdatePixelPositionsFromMeshIntersections(void) 
{
  // Parameters
  RNLength max_distance = 1000;
  RNScalar max_normal_angle = 5.0*RN_PI/6.0;
  RNScalar max_curvature = 1;
  RNScalar min_cos_normal_angle = (max_normal_angle >= 0) ? cos(max_normal_angle) : -1;

  // Consider features scan-by-scan so that fetch mesh once for each scan
  for (int ir = 0; ir < scene->NRuns(); ir++) {
    GSVRun *run = scene->Run(ir);
    for (int is = 0; is < run->NSegments(); is++) {
      GSVSegment *segment = run->Segment(is);
      for (int ia = 0; ia < segment->NScans(); ia++) {
        GSVScan *scan = segment->Scan(ia);
        if (scan->NScanlines() == 0) continue;
        GSVMesh *mesh = scan->Mesh();
        if (!mesh) continue;
        printf("UPR %d %d %d\n", ir, is, ia);

        // Process all features in segment
        for (int k = 0; k < NFeatures(); k++) {
          GSVFeature *feature = Feature(k);
          GSVImage *image = feature->image;
          if (!image) continue;
          GSVPanorama *panorama = image->Panorama();
          if (!panorama) continue;
          if (panorama->Segment() != segment) continue;

          // Get ray through pixel
          R3Ray ray = feature->image->RayThroughUndistortedPosition(feature->image_position);
          R3Affine transformation = OptimizedTransformation(feature);
          ray.Transform(transformation);

          // Get ray intersection at midpoint
          R3MeshIntersection intersection;
          if (!mesh->Intersection(ray, &intersection, 0, max_distance)) continue;

          // Check if closest intersection
          if ((feature->image_t > 0) && (intersection.t > feature->image_t)) continue;
          
          // Get closest mesh vertex
          GSVMeshVertex *vertex = (GSVMeshVertex *) intersection.vertex;
          if (!vertex && intersection.edge) vertex = mesh->VertexOnEdge(intersection.edge);
          if (!vertex && intersection.face) vertex = mesh->VertexOnFace(intersection.face);
          if (!vertex) continue;
          
          // Check normal angle
          if (min_cos_normal_angle > -1) {
            RNScalar dot = -1.0 * ray.Vector().Dot(mesh->VertexNormal(vertex));
            if (dot < min_cos_normal_angle) continue;
          }
              
          // Check curvature
          if (max_curvature > 0) {
            RNScalar curvature = mesh->VertexCurvature(vertex);
            if (curvature > max_curvature) continue;
          }

          // Update feature info
          // GSVScanline *hit_scanline = mesh->VertexScanline(vertex);
          // int hit_point_index = mesh->VertexPointIndex(vertex);
          feature->scan_position = intersection.point;
          feature->scan_normal = mesh->VertexNormal(vertex);
          feature->image_t = intersection.t;

          // Determine scan direction
          if ((feature->image_direction.X() != RN_UNKNOWN) && (feature->image_direction.Y() != RN_UNKNOWN)) {
            // Reset feature info (require ray intersections at both endpoints)
            feature->image_t = 0;
            
            // Get ray intersection at first endpoint
            R2Point image_position1 = feature->image_position + feature->image_scale * feature->image_direction;
            R3Ray ray1 = feature->image->RayThroughUndistortedPosition(image_position1);
            ray1.Transform(transformation);
            R3MeshIntersection intersection1;
            if (!mesh->Intersection(ray1, &intersection1, 0, max_distance)) continue;
            GSVMeshVertex *vertex1 = (GSVMeshVertex *) intersection1.vertex;
            if (!vertex1 && intersection1.edge) vertex1 = mesh->VertexOnEdge(intersection1.edge);
            if (!vertex1 && intersection1.face) vertex1 = mesh->VertexOnFace(intersection1.face);
            if (!vertex1) continue;
            
            // Get ray intersection at first endpoint
            R2Point image_position2 = feature->image_position + feature->image_scale * feature->image_direction;
            R3Ray ray2 = feature->image->RayThroughUndistortedPosition(image_position2);
            ray2.Transform(transformation);
            R3MeshIntersection intersection2;
            if (!mesh->Intersection(ray2, &intersection2, 0, max_distance)) continue;
            GSVMeshVertex *vertex2 = (GSVMeshVertex *) intersection2.vertex;
            if (!vertex2 && intersection2.edge) vertex2 = mesh->VertexOnEdge(intersection2.edge);
            if (!vertex2 && intersection2.face) vertex2 = mesh->VertexOnFace(intersection2.face);
            if (!vertex2) continue;
            
            // Check if all three points are colinear
            // ???

            // Update feature info
            feature->scan_position = intersection.point;
            feature->scan_normal = mesh->VertexNormal(vertex);
            feature->scan_direction = 0.5 * (intersection2.point - intersection1.point);
            feature->scan_scale = feature->scan_direction.Length();
            if (RNIsNotZero(feature->scan_scale)) feature->scan_direction /= feature->scan_scale;
            feature->image_t = intersection.t;
          }
        }

        // Delete mesh
        delete mesh;
      }
    }
  }

  // Return success
  return 1;
}



int GSVPoseOptimization::
UpdateClusterPositionsFromFeaturePositions(void) 
{
  // Initialize marks
  int *marks = new int [ NClusters() ];
  for (int i = 0; i < NClusters(); i++) marks[i] = 0;

  // Create a stack with all clusters (sorted children last, clusters can be added multiple times)
  RNArray<GSVFeatureCluster *> stack;
  for (int i = 0; i < NClusters(); i++) {
    GSVFeatureCluster *cluster = Cluster(i);
    stack.Insert(cluster);
    for (int j = 0; j < cluster->NClusters(); j++) {
      GSVFeatureCluster *child = cluster->Cluster(j);
      stack.Insert(child);
    }
  }

  // Update scan positions for clusters
  while (!stack.IsEmpty()) {
    GSVFeatureCluster *cluster = stack.Tail();
    stack.RemoveTail();
    if (marks[cluster->index]) continue;
    marks[cluster->index] = 1;
    cluster->Update(this);
  }

  // Delete marks
  delete [] marks;

  // Return success
  return 1;
}



RNScalar
FeatureInertia(GSVFeature *feature, int variable_index)
{
  // Initialize result
  RNScalar inertia = 0;
  
  // Check feature type
  if (feature->image) {
    GSVImage *image = feature->image;
    ImageData *image_data = (ImageData *) image->Data();
    if (!image_data) return 0;
    GSVTapestry *tapestry = image->Tapestry();
    if (!tapestry) return 0;
    GSVSegment *segment = tapestry->Segment();
    if (!segment) return 0;
    SegmentData *segment_data = (SegmentData *) segment->Data();
    if (!segment_data) return 0;
    GSVPath *path = segment_data->path;
    if (!path) return 0;
    inertia = path->Inertia(image_data->path_parameter, variable_index);
  }
  else if (feature->scanline) {
    GSVScanline *scanline = feature->scanline;
    ScanlineData *scanline_data = (ScanlineData *) scanline->Data();
    if (!scanline_data) return 0;
    GSVScan *scan = scanline->Scan();
    if (!scan) return 0;
    GSVSegment *segment = scan->Segment();
    if (!segment) return 0;
    SegmentData *segment_data = (SegmentData *) segment->Data();
    if (!segment_data) return 0;
    GSVPath *path = segment_data->path;
    if (!path) return 0;
    inertia = path->Inertia(scanline_data->path_parameter, variable_index);
  }

  // Return inertia
  return inertia;
}



static void
UpdateClusterInertias(GSVFeatureCluster *cluster) 
{
  // Compute weighted average
  RNScalar total_inertia[6] = { 0.0 };
  RNScalar total_weight = 0;

  // Consider features
  for (int i = 0; i < cluster->NFeatures(); i++) {
    GSVFeature *feature = cluster->Feature(i);
    RNScalar weight = feature->score;
    total_weight += weight;
    for (int j = 0; j < 6; j++) {
      total_inertia[j] += weight * FeatureInertia(feature, j);
    }
  }

  // Consider subclusters
  for (int i = 0; i < cluster->NClusters(); i++) {
    GSVFeatureCluster *child = cluster->Cluster(i);
    RNScalar weight = child->score;
    total_weight += weight;
    for (int j = 0; j < 6; j++) {
      total_inertia[j] += weight * child->inertia[j];
    }
  }

  // Compute averages
  for (int j = 0; j < 6; j++) {
    cluster->inertia[j] = (total_weight > 0) ? total_inertia[j] / total_weight : RN_INFINITY;
  }
}



int GSVPoseOptimization::
UpdateClusterInertiasFromFeatureInertias(void) 
{
  // Initialize marks
  int *marks = new int [ NClusters() ];
  for (int i = 0; i < NClusters(); i++) marks[i] = 0;

  // Create a stack with all clusters (sorted children last, clusters can be added multiple times)
  RNArray<GSVFeatureCluster *> stack;
  for (int i = 0; i < NClusters(); i++) {
    GSVFeatureCluster *cluster = Cluster(i);
    stack.Insert(cluster);
    for (int j = 0; j < cluster->NClusters(); j++) {
      GSVFeatureCluster *child = cluster->Cluster(j);
      stack.Insert(child);
    }
  }

  // Update scan positions for clusters
  while (!stack.IsEmpty()) {
    GSVFeatureCluster *cluster = stack.Tail();
    stack.RemoveTail();
    if (marks[cluster->index]) continue;
    marks[cluster->index] = 1;
    UpdateClusterInertias(cluster);
  }

  // Delete marks
  delete [] marks;

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Transformation functions
////////////////////////////////////////////////////////////////////////

R3Affine GSVPoseOptimization::
OptimizedTransformation(const GSVFeature *feature) const
{
  // Return transformation to be applied to feature
  if (feature->image) return OptimizedTransformation(feature->image, feature->image_position);
  else if (feature->scanline) return OptimizedTransformation(feature->scanline);
  else return R3identity_affine;
}



R3Affine GSVPoseOptimization::
OptimizedTransformation(const GSVScanline *scanline) const
{
  // Get convenient variables
  ScanlineData *scanline_data = (ScanlineData *) scanline->Data();
  if (!scanline_data) return R3identity_affine;
  GSVScan *scan = scanline->Scan();
  if (!scan) return R3identity_affine;
  GSVLaser *laser = scanline->Scan()->Laser();
  if (!laser) return R3identity_affine;
  LaserData *laser_data = (LaserData *) laser->Data();
  if (!laser_data) return R3identity_affine;
  GSVSegment *segment = scan->Segment();
  if (!segment) return R3identity_affine;
  SegmentData *segment_data = (SegmentData *) segment->Data();
  if (!segment_data) return R3identity_affine;
  GSVPath *path = segment_data->path;
  if (!path) return R3identity_affine;

  // Get coordinate system 
  GSVPose pose = path->Pose(scanline_data->path_parameter);
  const R3Point& center = pose.Viewpoint();
  R3Vector right = pose.Right();
  R3Vector towards = pose.Towards();
  R3Vector up = pose.Up();

  // Get scanline viewpoint
  const GSVPose& scanline_pose = scanline->Pose();
  const R3Point& viewpoint = scanline_pose.Viewpoint();

  // Get transformation variables
  R3Vector center_translation = path->Translation(scanline_data->path_parameter);
  R3Vector center_rotation = path->Rotation(scanline_data->path_parameter);
  R3Vector viewpoint_translation = laser_data->translation;
  R3Vector viewpoint_rotation = laser_data->rotation;
  
  // Reparamerize viewpoint translation
  if (!viewpoint_translation.IsZero()) {
    R3Vector vt = viewpoint_translation;
    viewpoint_translation[0] = vt.X()*right.X() + vt.Y()*towards.X() + vt.Z()*up.X();
    viewpoint_translation[1] = vt.X()*right.Y() + vt.Y()*towards.Y() + vt.Z()*up.Y();
    viewpoint_translation[2] = vt.X()*right.Z() + vt.Y()*towards.Z() + vt.Z()*up.Z();
  }

  // Return transformation 
  R4Matrix matrix(R4identity_matrix);
  matrix.Translate(center_translation);
  if (!center_rotation.IsZero()) {
    matrix.Translate(center.Vector());
    matrix.Rotate(center_rotation);
    matrix.Translate(-(center.Vector()));
  }
  matrix.Translate(viewpoint_translation);
  if (!viewpoint_rotation.IsZero()) {
    matrix.Translate(viewpoint.Vector());
    matrix.Rotate(viewpoint_rotation);
    matrix.Translate(-(viewpoint.Vector()));
  }
  return R3Affine(matrix);
}



R3Affine GSVPoseOptimization::
OptimizedTransformation(const GSVImage *image, const R2Point& undistorted_position) const
{
  // Get convenient variables
  ImageData *image_data = (ImageData *) image->Data();
  if (!image_data) return R3identity_affine;
  GSVTapestry *tapestry = image->Tapestry();
  if (!tapestry) return R3identity_affine;
  GSVCamera *camera = image->Tapestry()->Camera();
  if (!camera) return R3identity_affine;
  CameraData *camera_data = (CameraData *) camera->Data();
  if (!camera_data) return R3identity_affine;
  GSVSegment *segment = tapestry->Segment();
  if (!segment) return R3identity_affine;
  SegmentData *segment_data = (SegmentData *) segment->Data();
  if (!segment_data) return R3identity_affine;
  GSVPath *path = segment_data->path;
  if (!path) return R3identity_affine;

  // Get coordinate system 
  GSVPose pose = path->Pose(image_data->path_parameter);
  const R3Point& center = pose.Viewpoint();
  const R3Vector right = pose.Right();
  const R3Vector towards = pose.Towards();
  const R3Vector up = pose.Up();

  // Get image viewpoint
  R3Ray ray = image->RayThroughUndistortedPosition(undistorted_position);
  const R3Point& viewpoint = ray.Start();

  // Get center translation
  R3Vector center_translation = R3zero_vector;
  R3Vector ct = path->Translation(image_data->path_parameter);
  if (!ct.IsZero()) {
    center_translation = ct;
    // center_translation[0] = ct.X()*right.X() + ct.Y()*towards.X() + ct.Z()*up.X();
    // center_translation[1] = ct.X()*right.Y() + ct.Y()*towards.Y() + ct.Z()*up.Y();
    // center_translation[2] = ct.X()*right.Z() + ct.Y()*towards.Z() + ct.Z()*up.Z();
  }

  // Get viewpoint translation
  R3Vector viewpoint_translation = R3zero_vector;
  R3Vector vt = camera_data->translation;
  if (!vt.IsZero()) {
    viewpoint_translation[0] = vt.X()*right.X() + vt.Y()*towards.X() + vt.Z()*up.X();
    viewpoint_translation[1] = vt.X()*right.Y() + vt.Y()*towards.Y() + vt.Z()*up.Y();
    viewpoint_translation[2] = vt.X()*right.Z() + vt.Y()*towards.Z() + vt.Z()*up.Z();
  }

  // Get center rotation
  R3Vector center_rotation = path->Rotation(image_data->path_parameter);

  // Get viewpoint rotation
  R3Vector viewpoint_rotation = camera_data->rotation;

  // Return transformation 
  R4Matrix matrix(R4identity_matrix);
  matrix.Translate(center_translation);
  if (!center_rotation.IsZero()) {
    matrix.Translate(center.Vector());
    matrix.Rotate(center_rotation);
    matrix.Translate(-(center.Vector()));
  }
  matrix.Translate(viewpoint_translation);
  if (!viewpoint_rotation.IsZero()) {
    matrix.Translate(viewpoint.Vector());
    matrix.Rotate(viewpoint_rotation);
    matrix.Translate(-(viewpoint.Vector()));
  }
  return R3Affine(matrix);
}



R3Affine GSVPoseOptimization::
OptimizedTransformation(const GSVFeatureCluster *cluster) const
{
  // Compute translation
  R3Vector translation = R3zero_vector;
  if (cluster->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) {
    translation = cluster->translation.X() * cluster->scan_normal;
  }
  else if (cluster->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) {
    RNDimension dim = cluster->scan_direction.MinDimension();
    R3Vector tmp = cluster->scan_direction % R3xyz_triad[dim];
    R3Vector axis1 = cluster->scan_direction % tmp; axis1.Normalize();
    R3Vector axis2 = cluster->scan_direction % axis1; axis2.Normalize();
    translation = cluster->translation.X()*axis1 + cluster->translation.Y()*axis2;
  }
  else {
    translation = cluster->translation;
  }

  // Return transformation 
  R4Matrix matrix(R4identity_matrix);
  matrix.Translate(translation);
  matrix.Translate(cluster->scan_position.Vector());
  matrix.Rotate(cluster->rotation);
  matrix.Translate(-(cluster->scan_position.Vector()));
  return R3Affine(matrix);
}



////////////////////////////////////////////////////////////////////////
// Optimization stuff
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
Solve(
  RNBoolean update_laser_translations, RNBoolean update_laser_rotations,
  RNBoolean update_camera_translations, RNBoolean update_camera_rotations,
  RNBoolean update_path_translations,  RNBoolean update_path_rotations,
  RNBoolean update_cluster_positions, RNBoolean update_cluster_rotations,
  RNBoolean update_pixel_depths, 
  RNBoolean include_correspondence_equations,
  RNBoolean include_cluster_equations,
  RNScalar rigidity, int solver)
{
  // Check consistency of options
  if (update_cluster_positions) include_cluster_equations = TRUE;
  if (update_cluster_rotations) include_cluster_equations = TRUE;
  if (include_cluster_equations) update_cluster_positions = TRUE;

  // Determine which laser variables to optimize
  laser_nv = 0;
  for (int i = 0; i < LASER_NV; i++) laser_v[i] = -1;
  if (update_laser_translations) {
    laser_v[TX] = laser_nv++;
    laser_v[TY] = laser_nv++;
    laser_v[TZ] = laser_nv++;
  }
  if (update_laser_rotations) {
    laser_v[RX] = laser_nv++;
    laser_v[RY] = laser_nv++;
    laser_v[RZ] = laser_nv++;
  }

  // Determine which camera variables to optimize
  camera_nv = 0;
  for (int i = 0; i < CAMERA_NV; i++) camera_v[i] = -1;
  if (update_camera_translations) {
    camera_v[TX] = camera_nv++;
    camera_v[TY] = camera_nv++;
    camera_v[TZ] = camera_nv++;
  }
  if (update_camera_rotations) {
    camera_v[RX] = camera_nv++;
    camera_v[RY] = camera_nv++;
    camera_v[RZ] = camera_nv++;
  }

  // Determine which path variables to optimize
  vertex_nv = 0;
  for (int i = 0; i < VERTEX_NV; i++) vertex_v[i] = -1;
  if (update_path_translations) {
    vertex_v[TX] = vertex_nv++;
    vertex_v[TY] = vertex_nv++;
    vertex_v[TZ] = vertex_nv++;
  }
  if (update_path_rotations) {
    vertex_v[RX] = vertex_nv++; 
    vertex_v[RY] = vertex_nv++; 
    vertex_v[RZ] = vertex_nv++;
  }

  // Determine which pixel feature variables to optimize
  pixel_feature_nv = 0;
  for (int i = 0; i < PIXEL_FEATURE_NV; i++) pixel_feature_v[i] = -1;
  if (update_pixel_depths) {
    pixel_feature_v[T] = pixel_feature_nv++;
  }

  // Determine which point cluster variables to optimize
  cluster_nv = 0;
  for (int i = 0; i < CLUSTER_NV; i++) cluster_v[i] = -1;
  if (update_cluster_positions) {
    cluster_v[TX] = cluster_nv++;
    cluster_v[TY] = cluster_nv++;
    cluster_v[TZ] = cluster_nv++;
  }
  if (update_cluster_rotations) {
    cluster_v[RX] = cluster_nv++;
    cluster_v[RY] = cluster_nv++;
    cluster_v[RZ] = cluster_nv++;
  }

  // Compute list of pixel features 
  pixel_features.Empty();
  assert(!pixel_feature_indices);
  if (NFeatures() > 0) {
    pixel_feature_indices = new int [ NFeatures() + 1 ];
    for (int i = 0; i < NFeatures(); i++) pixel_feature_indices[i] = -1;
    for (int i = 0; i < NCorrespondences(); i++) {
      GSVFeatureCorrespondence *correspondence = Correspondence(i);
      for (int j = 0; j < 2; j++) {
        GSVFeature *feature = correspondence->Feature(j);
        if (!feature->image) continue;
        if (pixel_feature_indices[feature->index] == -1) {
          pixel_feature_indices[feature->index] = pixel_features.NEntries();
          pixel_features.Insert(feature);
        }
      }
    }
    for (int i = 0; i < NClusters(); i++) {
      GSVFeatureCluster *cluster = Cluster(i);
      for (int j = 0; j < cluster->NFeatures(); j++) {
        GSVFeature *feature = cluster->Feature(j);
        if (!feature->image) continue;
        if (pixel_feature_indices[feature->index] == -1) {
          pixel_feature_indices[feature->index] = pixel_features.NEntries();
          pixel_features.Insert(feature);
        }
      }
    }
  }

  // Get total number of variables
  int n = 0;
  n += laser_nv * lasers.NEntries();
  n += camera_nv * cameras.NEntries();
  n += vertex_nv * vertices.NEntries();
  n += pixel_feature_nv * pixel_features.NEntries();
  n += cluster_nv * clusters.NEntries();

  // Check number of variables
  if (n == 0) {
    if (pixel_feature_indices) delete [] pixel_feature_indices;
    pixel_feature_indices = NULL;
    return 1;
  }

  // Create system of equations
  RNSystemOfEquations *equations = new RNSystemOfEquations(n);
  if (!equations) {
    fprintf(stderr, "Unable to allocate system of equations\n");
    if (pixel_feature_indices) delete pixel_feature_indices;
    pixel_feature_indices = NULL;
    return 0;
  }
  
  // Insert equations 
  printf("HEREA %g %g\n",  WorldRMSD(), ImageRMSD()); RNTime t; t.Read();
  AddInertiaEquations(equations);
  printf("HEREB %.3f\n", t.Elapsed()); t.Read();
  AddRigidityEquations(equations, rigidity);
  printf("HEREC %.3f\n", t.Elapsed()); t.Read();
  if (include_correspondence_equations) {
    AddCorrespondenceEquations(equations);
    printf("HERED %.3f\n", t.Elapsed()); t.Read();
  }
  if (include_cluster_equations) {
    AddClusterEquations(equations); 
    printf("HEREE %.3f\n", t.Elapsed()); t.Read();
  }

  // Allocate x variables
  assert(n == equations->NVariables());
  double *x = new double [ n ];
  for (int i = 0; i < n; i++) x[i] = 0;

  // Initialize optimization variables
  if (!InitializeOptimizationVariables(x)) {
    fprintf(stderr, "Unable to initialize optimization variables\n");
    if (pixel_feature_indices) delete pixel_feature_indices;
    pixel_feature_indices = NULL;
    delete equations;
    delete [] x;
    return 0;
  }

  // Print errors
  printf("Initial Errors:\n");
  PrintErrors(equations, x, rigidity, solver, TRUE);
  printf("\n");

  printf("HEREF %.3f\n", t.Elapsed()); t.Read();

  // Optimize
  if (!ExecuteOptimization(equations, x, solver)) {
    fprintf(stderr, "Unable to execute optimization\n");
    if (pixel_feature_indices) delete pixel_feature_indices;
    pixel_feature_indices = NULL;
    delete equations;
    delete [] x;
    return 0;
  }

  printf("HEREG %.3f\n", t.Elapsed()); t.Read();

  // Print final errors
  printf("Final Errors:\n");
  PrintErrors(equations, x, rigidity, solver, TRUE);
  printf("\n");

  // Extract optimization variables
  if (!ExtractOptimizationVariables(x)) {
    fprintf(stderr, "Unable to extract optimization variables\n");
    if (pixel_feature_indices) delete pixel_feature_indices;
    pixel_feature_indices = NULL;
    delete equations;
    delete [] x;
    return 0;
  }

  printf("HEREH %.3f\n", t.Elapsed()); t.Read();

  // Delete everything
  if (pixel_feature_indices) delete pixel_feature_indices;
  pixel_feature_indices = NULL;
  delete equations;
  delete [] x;

  printf("HEREI %g %g\n", WorldRMSD(), ImageRMSD());

  // Return success
  return 1;
}



int GSVPoseOptimization::
ClearOptimizationVariables(void)
{
  // Clear all laser transformations
  for (int i = 0; i < lasers.NEntries(); i++) {
    GSVLaser *laser = lasers.Kth(i);
    LaserData *laser_data = (LaserData *) laser->Data();
    laser_data->translation = R3zero_vector;
  }

  // Clear all camera transformations
  for (int i = 0; i < cameras.NEntries(); i++) {
    GSVCamera *camera = cameras.Kth(i);
    CameraData *camera_data = (CameraData *) camera->Data();
    camera_data->translation = R3zero_vector;
  }

  // Clear all spline vertex transformations
  for (int i = 0; i < vertices.NEntries(); i++) {
    GSVPathVertex *vertex = vertices.Kth(i);
    vertex->translation = R3zero_vector;
    vertex->rotation = R3zero_vector;
  }

  // Clear all pixel depths
  for (int i = 0; i < pixel_features.NEntries(); i++) {
    const double default_image_t = 10;
    GSVFeature *feature = pixel_features.Kth(i);
    R3Ray ray = feature->image->RayThroughUndistortedPosition(feature->image_position);
    feature->scan_position = ray.Point(default_image_t);
    feature->image_t = 0;
  }

  // Clear all feature cluster variables
  for (int i = 0; i < clusters.NEntries(); i++) {
    GSVFeatureCluster *cluster = clusters.Kth(i);
    cluster->translation = R3zero_vector;
    cluster->rotation = R3zero_vector;
    cluster->Update(this);
  }

  // Return success
  return 1;
}



int GSVPoseOptimization::
InitializeOptimizationVariables(double *x)
{
  // Initialize all laser variables
  if (laser_nv > 0) {
    for (int i = 0; i < lasers.NEntries(); i++) {
      for (int j = 0; j < LASER_NV; j++) {
        int v = LaserVariableIndex(i, j);
        if (v == -1) continue;
        x[v] = LaserVariableValue(i, j);
      }
    }
  }

  // Initialize all camera variables
  if (camera_nv > 0) {
    for (int i = 0; i < cameras.NEntries(); i++) {
      for (int j = 0; j < CAMERA_NV; j++) {
        int v = CameraVariableIndex(i, j);
        if (v == -1) continue;
        x[v] = CameraVariableValue(i, j);
      }
    }
  }

  // Initialize all vertex variables
  if (vertex_nv > 0) {
    for (int i = 0; i < vertices.NEntries(); i++) {
      for (int j = 0; j < VERTEX_NV; j++) {
        int v = VertexVariableIndex(i, j);
        if (v == -1) continue;
        x[v] = VertexVariableValue(i, j);
      }
    }
  }

  // Initialize all pixel feature variables
  if (pixel_feature_nv > 0) {
    for (int i = 0; i < pixel_features.NEntries(); i++) {
      for (int j = 0; j < PIXEL_FEATURE_NV; j++) {
        int v = PixelFeatureVariableIndex(i, j);
        if (v == -1) continue;
        x[v] = PixelFeatureVariableValue(i, j);
        if ((j == T) && (x[v] <= 0)) x[v] = 10;
      }
    }
  }

  // Initialize all cluster variables
  if (cluster_nv > 0) {
    for (int i = 0; i < clusters.NEntries(); i++) {
      for (int j = 0; j < CLUSTER_NV; j++) {
        int v = ClusterVariableIndex(i, j);
        if (v == -1) continue;
        x[v] = ClusterVariableValue(i, j);
      }
    }
  }

  // Return success
  return 1;
}



int GSVPoseOptimization::
ExtractOptimizationVariables(double *x)
{
  // Extract all laser variables
  if (laser_nv > 0) {
    printf("Lasers:\n");
    for (int i = 0; i < lasers.NEntries(); i++) {
      printf("L %d ", i);
      for (int j = 0; j < LASER_NV; j++) {
        SetLaserVariableValue(i, j, x);
        printf("%12.9f ", LaserVariableValue(i, j));
      }
      printf("\n");
    }
  }

  // Extract all camera variables
  if (camera_nv > 0) {
    // printf("Cameras:\n");
    for (int i = 0; i < cameras.NEntries(); i++) {
      printf("C %d ", i);
      for (int j = 0; j < CAMERA_NV; j++) {
        SetCameraVariableValue(i, j, x);
        printf("%12.9f ", CameraVariableValue(i, j));
      }
      printf("\n");
    }
  }

  // Extract all vertex variables
  if (vertex_nv > 0) {
    for (int i = 0; i < vertices.NEntries(); i++) {
      for (int j = 0; j < VERTEX_NV; j++) {
        SetVertexVariableValue(i, j, x);
      }
    }
  }

  // Extract all pixel feature variables
  if (pixel_feature_nv > 0) {
    for (int i = 0; i < pixel_features.NEntries(); i++) {
      for (int j = 0; j < PIXEL_FEATURE_NV; j++) {
        SetPixelFeatureVariableValue(i, j, x);
      }
    }
  }

 // Extract all cluster variables
  if (cluster_nv > 0) {
    for (int i = 0; i < clusters.NEntries(); i++) {
      // printf("G %d ", i);
      for (int j = 0; j < CLUSTER_NV; j++) {
        SetClusterVariableValue(i, j, x);
        // printf("%12.9f ", ClusterVariableValue(i, j));
      }
      // printf("\n");
    }
  }

  // Update pixel positions
  UpdatePixelPositionsFromRayIntersections(); 
  UpdateClusterPositionsFromFeaturePositions();  

  // Return success
  return 1;
}



int GSVPoseOptimization::
ExecuteOptimization(RNSystemOfEquations *equations, double *x, int solver)
{
  // Solve system of equations
  if (!equations->Minimize(x, solver)) {
    fprintf(stderr, "Unable to minimize system of equations\n");
    return 0;
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////
// Intertia equation stuff
////////////////////////////////////////////////////////

void GSVPoseOptimization::
AddInertiaEquations(RNSystemOfEquations *system)
{
  // Convenient variables
  RNScalar individual_inertia_scale = 1E3;

  // Insert equations for laser variables
  if (laser_nv > 0) {
    for (int i = 0; i < lasers.NEntries(); i++) {
      GSVLaser *laser = lasers.Kth(i);
      LaserData *laser_data = (LaserData *) laser->Data();
      if (!laser_data) continue;
      for (int j = 0; j < LASER_NV; j++) {
        if (LaserVariableIndex(i, j) == -1) continue;
        RNScalar weight = laser_inertia_weights[j] + individual_inertia_scale * laser_data->inertia[j];
        RNPolynomial *equation = new RNPolynomial();
        AddLaserVariableToEquation(equation, i, j, weight);
        *equation -= LaserVariableValue(i, j) * weight;
        system->InsertEquation(equation);
      }
    }
  }

  // Insert equations for camera variables
  if (camera_nv > 0) {
    for (int i = 0; i < cameras.NEntries(); i++) {
      GSVCamera *camera = cameras.Kth(i);
      CameraData *camera_data = (CameraData *) camera->Data();
      if (!camera_data) continue;
      for (int j = 0; j < CAMERA_NV; j++) {
        if (CameraVariableIndex(i, j) == -1) continue;
        RNScalar weight = camera_inertia_weights[j] + individual_inertia_scale * camera_data->inertia[j];
        RNPolynomial *equation = new RNPolynomial();
        AddCameraVariableToEquation(equation, i, j, weight);
        *equation -= CameraVariableValue(i, j) * weight;
        system->InsertEquation(equation);
      }
    }
  }

  // Insert equations for path variables
  if (vertex_nv > 0) {
    for (int i = 0; i < vertices.NEntries(); i++) {
      GSVPathVertex *vertex = vertices.Kth(i);
      for (int j = 0; j < VERTEX_NV; j++) {
        if (VertexVariableIndex(i, j) == -1) continue;
        RNScalar weight = vertex_inertia_weights[j] + individual_inertia_scale * vertex->inertia[j];
        RNPolynomial *equation = new RNPolynomial();
        AddPathVertexVariableToEquation(equation, i, j, weight);
        *equation -= VertexVariableValue(i, j) * weight;
        system->InsertEquation(equation);
      }
    }
  }

  // Insert equations for pixel feature variables
  if (pixel_feature_nv > 0) {
    for (int i = 0; i < pixel_features.NEntries(); i++) {
      for (int j = 0; j < PIXEL_FEATURE_NV; j++) {
        if (PixelFeatureVariableIndex(i, j) == -1) continue;
        RNScalar weight = pixel_feature_inertia_weights[j];
        RNPolynomial *equation = new RNPolynomial();
        AddPixelFeatureVariableToEquation(equation, i, j, weight);
        *equation -= PixelFeatureVariableValue(i, j) * weight;
        system->InsertEquation(equation);
      }
    }
  }

  // Insert equations for point cluster variables
  if (cluster_nv > 0) {
    for (int i = 0; i < clusters.NEntries(); i++) {
      GSVFeatureCluster *cluster = clusters.Kth(i);
      for (int j = 0; j < CLUSTER_NV; j++) {
        if (ClusterVariableIndex(i, j) == -1) continue;
        RNScalar weight = cluster_inertia_weights[j] + individual_inertia_scale * cluster->inertia[j];
        RNPolynomial *equation = new RNPolynomial();
        AddClusterVariableToEquation(equation, i, j, weight);
        *equation -= ClusterVariableValue(i, j) * weight;
        system->InsertEquation(equation);
      }
    }
  }
}



////////////////////////////////////////////////////////
// Rigidity equation stuff
////////////////////////////////////////////////////////

void GSVPoseOptimization::
AddRigidityEquations(RNSystemOfEquations *system, RNScalar rigidity)
{
#if 0
  // Insert equations to ensure smooth variations in path variables
  RNScalar path_variable_rigidity_weight = 100 * pow(0.5 + rigidity, 5);
  for (int is = 0; is < segments.NEntries(); is++) {
    GSVSegment *segment = segments.Kth(is);
    SegmentData *segment_data = (SegmentData *) segment->Data();
    if (!segment_data) continue;
    GSVPath *path = segment_data->path;
    if (!path) continue;
    for (int i = 1; i < path->NVertices(); i++) {
      GSVPathVertex *vertex0 = path->Vertex(i-1);
      GSVPathVertex *vertex1 = path->Vertex(i);
      for (int j = 0; j < VERTEX_NV; j++) {
        if (vertex_v[j] == -1) continue;
        RNPolynomial *equation = new RNPolynomial();
        AddPathVertexVariableToEquation(equation, vertex0->index, j, path_variable_rigidity_weight);
        AddPathVertexVariableToEquation(equation, vertex1->index, j, -path_variable_rigidity_weight);
        system->InsertEquation(equation);
      }
    }
  }
#endif

  // Insert equations to maintain the shape of the path
  RNScalar path_rigidity_sigma_squared = path_rigidity_sigma * path_rigidity_sigma; 
  RNScalar path_rigidity_scale_factor = 100 * pow(0.5 + rigidity, 5);
  for (int is = 0; is < segments.NEntries(); is++) {
    GSVSegment *segment = segments.Kth(is);
    SegmentData *segment_data = (SegmentData *) segment->Data();
    if (!segment_data) continue;
    GSVPath *path = segment_data->path;
    if (!path) continue;
    for (int i0 = 0; i0 < path->NVertices(); i0++) {
      GSVPathVertex *vertex0 = path->Vertex(i0);
      const R3Point& position0 = path->VertexPosition(i0);
      for (int k = 1; k <= path_rigidity_radius; k++) {
        int i1 = i0 + k;
        if (i1 >= path->NVertices()) continue;
        GSVPathVertex *vertex1 = path->Vertex(i1);
        const R3Point& position1 = path->VertexPosition(i1);
        RNScalar s = exp(-2 * k * k / (path_rigidity_sigma_squared));
        RNScalar w = s * path_rigidity_scale_factor * path_rigidity_weight;
        AddPathRigidityEquations(system, vertex0->index, vertex1->index, position0, position1, w);
        AddPathRigidityEquations(system, vertex1->index, vertex0->index, position1, position0, w);
      }
    }
  }

#if 0
  // Insert equations to maintain arrangment of cameras within same panorama
  for (int is = 0; is < segments.NEntries(); is++) {
    GSVSegment *segment = segments.Kth(is);
    for (int ip = 0; ip < segment->NPanoramas(); ip++) {
      GSVPanorama *panorama = segment->Panorama(ip);
      for (int ii = 0; ii < panorama->NImages(); ii++) {
        GSVImage *image0 = panorama->Image(ii);
        for (int jj = ii+1; jj < panorama->NImages(); jj++) {
          GSVImage *image1 = panorama->Image(jj);
          R3Point viewpoint0 = image0->Pose().Viewpoint();
          R3Point viewpoint1 = image1->Pose().Viewpoint();
          // AddPointPointCorrespondenceEquations(system, image0, viewpoint0, image1, viewpoint0, camera_camera_rigidity_weight);
          // AddPointPointCorrespondenceEquations(system, image0, viewpoint1, image1, viewpoint1, camera_camera_rigidity_weight);
          AddCameraRigidityEquations(system, image0, image1, camera_camera_rigidity_weight);
        }
      }
    }
  }

  // Insert equations to maintain arrangment of lasers at same timestamp
  for (int is = 0; is < segments.NEntries(); is++) {
    GSVSegment *segment = segments.Kth(is);
    if (segment->NScans() < 3) continue;
    GSVScan *scan0 = segment->Scan(0);
    GSVScan *scan2 = segment->Scan(2);
    for (int ie = 0; ie < scan0->NScanlines(); ie++) {
      GSVScanline *scanline0 = scan0->Scanline(ie);
      GSVScanline *scanline1 = scan2->FindScanlineBeforeTimestamp(scanline0->Timestamp());
      if (!scanline1) continue;
      R3Point viewpoint0 = scanline0->Pose().Viewpoint();
      R3Point viewpoint1 = scanline1->Pose().Viewpoint();
      AddPointPointCorrespondenceEquations(system, scanline0, viewpoint0, scanline1, viewpoint0, laser_laser_rigidity_weight);
      AddPointPointCorrespondenceEquations(system, scanline0, viewpoint1, scanline1, viewpoint1, laser_laser_rigidity_weight);
      // AddLaserRigidityEquations(system, scanline0, scanline1, laser_laser_rigidity_weight);
    }
  }
#endif
}



void GSVPoseOptimization::
AddPathRigidityEquations(RNSystemOfEquations *system,
  int vertex_index0, int vertex_index1, 
  const R3Point& v0, const R3Point& v1, RNScalar w)
{
  //        |  1  -rz  ry |
  // Rxyz ~ |  rz  1  -rx |
  //        | -ry  rx  1  |

  // (R0(d) + v0 + t0) - (v1 + t1) = 0
  // (R0(d) - d + t0 - t1 = 0

  // X:      d.X() + -rz0*d.Y() +  ry0*d.Z() + -d.X() + tx0 - s1tx = 0
  // Y:  rz0*d.X() +      d.Y() + -rx0*d.Z() + -d.Y() + ty0 - s1ty = 0
  // Z: -ry0*d.X() +  rx0*d.Y() +      d.Z() + -d.Z() + tz0 - s1tz = 0

  // X:      -rz0*d.Y() +   ry0*d.Z() + tx0 - s1tx = 0
  // Y:       rz0*d.X() +  -rx0*d.Z() + ty0 - s1ty = 0
  // Z:      -ry0*d.X() +   rx0*d.Y() + tz0 - s1tz = 0

  // Get convenient variables
  RNPolynomial *equation;
  R3Vector d = v1 - v0;

  // Align X coordinates
  equation = new RNPolynomial();
  AddPathVertexVariableToEquation(equation, vertex_index0, TX,  1.0,   1.0, TRUE); 
  AddPathVertexVariableToEquation(equation, vertex_index0, RY,  d.Z(), 1.0, TRUE); 
  AddPathVertexVariableToEquation(equation, vertex_index0, RZ, -d.Y(), 1.0, TRUE); 
  AddPathVertexVariableToEquation(equation, vertex_index1, TX, -1.0,   1.0, TRUE); 
  equation->Multiply(w);                              
  system->InsertEquation(equation);                 
                                                      
  // Align Y coordinates                              
  equation = new RNPolynomial();                      
  AddPathVertexVariableToEquation(equation, vertex_index0, TY,  1.0,   1.0, TRUE); 
  AddPathVertexVariableToEquation(equation, vertex_index0, RX, -d.Z(), 1.0, TRUE); 
  AddPathVertexVariableToEquation(equation, vertex_index0, RZ,  d.X(), 1.0, TRUE); 
  AddPathVertexVariableToEquation(equation, vertex_index1, TY, -1.0,   1.0, TRUE); 
  equation->Multiply(w);                              
  system->InsertEquation(equation);                 
                                                      
  // Align Z coordinates                              
  equation = new RNPolynomial();                      
  AddPathVertexVariableToEquation(equation, vertex_index0, TZ,  1.0,   1.0, TRUE);      
  AddPathVertexVariableToEquation(equation, vertex_index0, RX,  d.Y(), 1.0, TRUE);
  AddPathVertexVariableToEquation(equation, vertex_index0, RY, -d.X(), 1.0, TRUE);
  AddPathVertexVariableToEquation(equation, vertex_index1, TZ, -1.0,   1.0, TRUE);      
  equation->Multiply(w);
  system->InsertEquation(equation);
}



void GSVPoseOptimization::
AddLaserRigidityEquations(RNSystemOfEquations *system,
  GSVScanline *scanline0, GSVScanline *scanline1, RNScalar w)
{
  // Check weight
  if (w == 0) return;

  // Get viewpoints
  R3Point viewpoint0 = scanline0->Pose().Viewpoint();
  R3Point viewpoint1 = scanline1->Pose().Viewpoint();

  // Add equations for consistency of viewpoint0
  if (1) {
    // Get transformed coordinates of viewpoint0
    RNAlgebraic *px0, *py0, *pz0, *px1, *py1, *pz1;
    if (!ComputeTransformedPointCoordinates(scanline0, viewpoint0, px0, py0, pz0)) return;
    if (!ComputeTransformedPointCoordinates(scanline1, viewpoint0, px1, py1, pz1)) return;
    AddPointPointDistanceEquations(system, px0, py0, pz0, px1, py1, pz1, w);
  }

  // Add equations for consistency of viewpoint1
  if (1) {
    // Get transformed coordinates of viewpoint1
    RNAlgebraic *px0, *py0, *pz0, *px1, *py1, *pz1;
    if (!ComputeTransformedPointCoordinates(scanline0, viewpoint1, px0, py0, pz0)) return;
    if (!ComputeTransformedPointCoordinates(scanline1, viewpoint1, px1, py1, pz1)) return;
    AddPointPointDistanceEquations(system, px0, py0, pz0, px1, py1, pz1, w);
  }
}



void GSVPoseOptimization::
AddCameraRigidityEquations(RNSystemOfEquations *system,
  GSVImage *image0, GSVImage *image1, RNScalar w)
{
  // Check weight
  if (w == 0) return;

  // Get viewpoints
  R3Point viewpoint0 = image0->Pose().Viewpoint();
  R3Point viewpoint1 = image1->Pose().Viewpoint();

  // Add equations for consistency of viewpoint0
  if (1) {
    // Get transformed coordinates of viewpoint0
    RNAlgebraic *px0, *py0, *pz0, *px1, *py1, *pz1;
    if (!ComputeTransformedPointCoordinates(image0, viewpoint0, px0, py0, pz0)) return;
    if (!ComputeTransformedPointCoordinates(image1, viewpoint0, px1, py1, pz1)) return;
    AddPointPointDistanceEquations(system, px0, py0, pz0, px1, py1, pz1, w);
  }

  // Add equations for consistency of viewpoint1
  if (1) {
    // Get transformed coordinates of viewpoint1
    RNAlgebraic *px0, *py0, *pz0, *px1, *py1, *pz1;
    if (!ComputeTransformedPointCoordinates(image0, viewpoint1, px0, py0, pz0)) return;
    if (!ComputeTransformedPointCoordinates(image1, viewpoint1, px1, py1, pz1)) return;
    AddPointPointDistanceEquations(system, px0, py0, pz0, px1, py1, pz1, w);
  }
}



////////////////////////////////////////////////////////
// Correspondence equation stuff
////////////////////////////////////////////////////////

void GSVPoseOptimization::
AddCorrespondenceEquations(RNSystemOfEquations *system)
{
  // Add equations for every correspondence
  for (int i = 0; i < NCorrespondences(); i++) {
    GSVFeatureCorrespondence *correspondence = Correspondence(i);
    GSVFeature *feature0 = correspondence->Feature(0);
    GSVFeature *feature1 = correspondence->Feature(1);
    if (!feature0 || !feature1) continue;
    AddCorrespondenceEquations(system, feature0, feature1, correspondence->score);
  }
}



void GSVPoseOptimization::
AddCorrespondenceEquations(RNSystemOfEquations *system, 
  GSVFeature *feature0, GSVFeature *feature1, 
  RNScalar w)
{
  // Compute the equations weights
  if (w == 0) return;
  RNScalar w3 = w * scan_point_scan_point_correspondence_weight;
  if (feature0->image) w3 = w * scan_point_image_point_correspondence_weight;
  if (feature1->image) w3 = w * scan_point_image_point_correspondence_weight;
  RNScalar w2 = w * image_point_image_point_correspondence_reprojection_weight;
  RNScalar wd = w * scan_direction_scan_direction_correspondence_weight;

  // 3D distances
  if (w3 > 0) {
    // SCANPOINT <-> SCANPOINT
    // SCANPOINT <-> IMAGEPOINT
    // IMAGEPOINT <-> SCANPOINT
    // IMAGEPOINT <-> IMAGEPOINT
    if (((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE))) {
      RNAlgebraic *px0, *py0, *pz0, *px1, *py1, *pz1;
      if (!ComputeTransformedPointCoordinates(feature0, px0, py0, pz0)) return;
      if (!ComputeTransformedPointCoordinates(feature1, px1, py1, pz1)) return;
      AddPointPointDistanceEquations(system, px0, py0, pz0, px1, py1, pz1, w3);
    }

    // SCANLINE <- SCANPOINT
    // SCANLINE <- SCANLINE
    // SCANLINE <- IMAGEPOINT
    // SCANLINE <- IMAGELINE
    if (((feature0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE))) {
      RNAlgebraic *px0, *py0, *pz0, *dx0, *dy0, *dz0, *px1, *py1, *pz1;
      R3Line line0(feature0->scan_position, feature0->scan_direction);
      R3Point center0 = feature1->scan_position; center0.Project(line0);
      if (!ComputeTransformedPointCoordinates(feature0->scanline, center0, px0, py0, pz0)) return;
      if (!ComputeTransformedVectorCoordinates(feature0->scanline, feature0->scan_direction, dx0, dy0, dz0)) return;
      if (!ComputeTransformedPointCoordinates(feature1, px1, py1, pz1)) return;
      AddPointLineDistanceEquations(system, px1, py1, pz1, px0, py0, pz0, dx0, dy0, dz0, w3);
    }

    // SCANPOINT -> SCANLINE
    // SCANLINE -> SCANLINE
    // IMAGEPOINT -> SCANLINE
    // IMAGELINE -> SCANLINE
    if (((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE))) {
      RNAlgebraic *px0, *py0, *pz0, *px1, *py1, *pz1, *dx1, *dy1, *dz1;
      R3Line line1(feature1->scan_position, feature1->scan_direction);
      R3Point center1 = feature0->scan_position; center1.Project(line1);
      if (!ComputeTransformedPointCoordinates(feature1->scanline, center1, px1, py1, pz1)) return;
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, feature1->scan_direction, dx1, dy1, dz1)) return;
      if (!ComputeTransformedPointCoordinates(feature0, px0, py0, pz0)) return;
      AddPointLineDistanceEquations(system, px0, py0, pz0, px1, py1, pz1, dx1, dy1, dz1, w3);
    }

    // SCANPLANE <- SCANPOINT
    // SCANPLANE <- SCANLINE
    // SCANPLANE <- SCANLINE
    // SCANPLANE <- IMAGEPOINT
    // SCANPLANE <- IMAGELINE
    if (((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE))) {
      RNAlgebraic *px0, *py0, *pz0, *nx0, *ny0, *nz0, *px1, *py1, *pz1;
      R3Plane plane0(feature0->scan_position, feature0->scan_normal);
      R3Point center0 = feature1->scan_position; center0.Project(plane0);
      if (!ComputeTransformedPointCoordinates(feature0->scanline, center0, px0, py0, pz0)) return;
      if (!ComputeTransformedVectorCoordinates(feature0->scanline, feature0->scan_normal, nx0, ny0, nz0)) return;
      if (!ComputeTransformedPointCoordinates(feature1, px1, py1, pz1)) return;
      AddPointPlaneDistanceEquations(system, px1, py1, pz1, px0, py0, pz0, nx0, ny0, nz0, w3);
    }

    // SCANPOINT -> SCANPLANE
    // SCANLINE  -> SCANPLANE
    // SCANPLANE -> SCANPLANE
    // IMAGEPOINT -> SCANPLANE
    // IMAGELINE -> SCANPLANE
    if (((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)  && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)  && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)  && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)  && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE))) {
      RNAlgebraic *px0, *py0, *pz0, *px1, *py1, *pz1, *nx1, *ny1, *nz1;
      R3Plane plane1(feature1->scan_position, feature1->scan_normal);
      R3Point center1 = feature0->scan_position; center1.Project(plane1);
      if (!ComputeTransformedPointCoordinates(feature1->scanline, center1, px1, py1, pz1)) return;
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, feature1->scan_normal, nx1, ny1, nz1)) return;
      if (!ComputeTransformedPointCoordinates(feature0, px0, py0, pz0)) return;
      AddPointPlaneDistanceEquations(system, px0, py0, pz0, px1, py1, pz1, nx1, ny1, nz1, w3);
    }
  }

  // 2D distances
  if (w2 > 0) {
    // IMAGEPOINT <- SCANPOINT
    // IMAGEPOINT <- IMAGEPOINT
    if (((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE))) {
      RNAlgebraic *ix1, *iy1;
      if (!ComputeTransformedImageCoordinates(feature0, feature1, ix1, iy1)) return;
      AddPointPointDistanceEquations(system, ix1, iy1, feature0->image_position, w2);
    }

    // SCANPOINT -> IMAGEPOINT
    // IMAGEPOINT -> IMAGEPOINT 
    if (((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE))) {
      RNAlgebraic *ix0, *iy0;
      if (!ComputeTransformedImageCoordinates(feature1, feature0, ix0, iy0)) return;
      AddPointPointDistanceEquations(system, ix0, iy0, feature1->image_position, w2);
    }

    // IMAGELINE <- SCANPOINT
    // IMAGELINE <- SCANLINE
    // IMAGELINE <- IMAGEPOINT
    // IMAGELINE <- IMAGELINE
    if (((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE))) {
      RNAlgebraic *ix1, *iy1;
      if (!ComputeTransformedImageCoordinates(feature0, feature1, ix1, iy1)) return;
      AddPointLineDistanceEquations(system, ix1, iy1, feature0->image_position, feature0->image_direction, w2);
    }

    // SCANPOINT -> IMAGELINE
    // SCANLINE -> IMAGELINE
    // IMAGEPOINT -> IMAGELINE
    // IMAGELINE -> IMAGELINE
    if (((feature0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) ||
        ((feature0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE))) {
      RNAlgebraic *ix0, *iy0;
      if (!ComputeTransformedImageCoordinates(feature1, feature0, ix0, iy0)) return;
      AddPointLineDistanceEquations(system, ix0, iy0, feature1->image_position, feature1->image_direction, w2);
    }
  }

  // Direction compatibility
  if (wd > 0) {
    if ((feature0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) {
      // Keep vectors peripendicular to direction of feature1 perpendicular to direction of feature0
      RNAlgebraic *dx0, *dy0, *dz0, *nx1, *ny1, *nz1;
      RNDimension dim1 = feature1->scan_direction.MinDimension();
      R3Vector normal1A = feature1->scan_direction % R3xyz_triad[dim1]; normal1A.Normalize();
      R3Vector normal1B = feature1->scan_direction % normal1A; normal1B.Normalize();
      if (!ComputeTransformedVectorCoordinates(feature0->scanline, feature0->scan_direction, dx0, dy0, dz0)) return;
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, normal1A, nx1, ny1, nz1)) return;
      AddVectorVectorDotProductEquations(system, dx0, dy0, dz0, nx1, ny1, nz1, 0.25 * wd);
      if (!ComputeTransformedVectorCoordinates(feature0->scanline, feature0->scan_direction, dx0, dy0, dz0)) return;
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, normal1B, nx1, ny1, nz1)) return;
      AddVectorVectorDotProductEquations(system, dx0, dy0, dz0, nx1, ny1, nz1, 0.25 * wd);
      
      // Keep vectors peripendicular to direction of feature0 perpendicular to direction of feature1
      RNAlgebraic *dx1, *dy1, *dz1, *nx0, *ny0, *nz0;
      RNDimension dim0 = feature0->scan_direction.MinDimension();
      R3Vector normal0A = feature0->scan_direction % R3xyz_triad[dim0]; normal0A.Normalize();
      R3Vector normal0B = feature0->scan_direction % normal0A; normal0B.Normalize();
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, feature1->scan_direction, dx1, dy1, dz1)) return;
      if (!ComputeTransformedVectorCoordinates(feature0->scanline, normal0A, nx0, ny0, nz0)) return;
      AddVectorVectorDotProductEquations(system, dx1, dy1, dz1, nx0, ny0, nz0, 0.25 * wd);
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, feature1->scan_direction, dx1, dy1, dz1)) return;
      if (!ComputeTransformedVectorCoordinates(feature0->scanline, normal0B, nx0, ny0, nz0)) return;
      AddVectorVectorDotProductEquations(system, dx1, dy1, dz1, nx0, ny0, nz0, 0.25 * wd);
    }

    if ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) {
      // Keep vectors peripendicular to direction of normal1 perpendicular to direction of normal0
      RNAlgebraic *nx0, *ny0, *nz0, *dx1, *dy1, *dz1;
      RNDimension dim1 = feature1->scan_normal.MinDimension();
      R3Vector direction1A = feature1->scan_normal % R3xyz_triad[dim1]; direction1A.Normalize();
      R3Vector direction1B = feature1->scan_normal % direction1A; direction1B.Normalize();
      if (!ComputeTransformedVectorCoordinates(feature0->scanline, feature0->scan_normal, nx0, ny0, nz0)) return;
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, direction1A, dx1, dy1, dz1)) return;
      AddVectorVectorDotProductEquations(system, nx0, ny0, nz0, dx1, dy1, dz1, 0.25 * wd);
      if (!ComputeTransformedVectorCoordinates(feature0->scanline, feature0->scan_normal, nx0, ny0, nz0)) return;
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, direction1B, dx1, dy1, dz1)) return;
      AddVectorVectorDotProductEquations(system, nx0, ny0, nz0, dx1, dy1, dz1, 0.25 * wd);
      
      // Keep vectors peripendicular to direction of normal0 perpendicular to direction of normal1
      RNAlgebraic *nx1, *ny1, *nz1, *dx0, *dy0, *dz0;
      RNDimension dim0 = feature0->scan_normal.MinDimension();
      R3Vector direction0A = feature0->scan_normal % R3xyz_triad[dim0]; direction0A.Normalize();
      R3Vector direction0B = feature0->scan_normal % direction0A; direction0B.Normalize();
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, feature1->scan_normal, nx1, ny1, nz1)) return;
      if (!ComputeTransformedVectorCoordinates(feature0->scanline, direction0A, dx0, dy0, dz0)) return;
      AddVectorVectorDotProductEquations(system, nx1, ny1, nz1, dx0, dy0, dz0, 0.25 * wd);
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, feature1->scan_normal, nx1, ny1, nz1)) return;
      if (!ComputeTransformedVectorCoordinates(feature0->scanline, direction0B, dx0, dy0, dz0)) return;
      AddVectorVectorDotProductEquations(system, nx1, ny1, nz1, dx0, dy0, dz0, 0.25 * wd);
    }

    if ((feature0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) {
      // Keep direction of feature1 perpendicular to normal of feature0
      RNAlgebraic *nx0, *ny0, *nz0, *dx1, *dy1, *dz1;
      if (!ComputeTransformedVectorCoordinates(feature0->scanline, feature0->scan_normal, nx0, ny0, nz0)) return;
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, feature1->scan_direction, dx1, dy1, dz1)) return;
      AddVectorVectorDotProductEquations(system, nx0, ny0, nz0, dx1, dy1, dz1, wd);
    }

    if ((feature0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) {
      // Keep direction of feature0 perpendicular to normal of feature1
      RNAlgebraic *dx0, *dy0, *dz0, *nx1, *ny1, *nz1;
      if (!ComputeTransformedVectorCoordinates(feature0->scanline, feature0->scan_direction, dx0, dy0, dz0)) return;
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, feature1->scan_normal, nx1, ny1, nz1)) return;
      AddVectorVectorDotProductEquations(system, dx0, dy0, dz0, nx1, ny1, nz1, wd);
    }
  }
}



////////////////////////////////////////////////////////
// Cluster equation stuff
////////////////////////////////////////////////////////

void GSVPoseOptimization::
AddClusterEquations(RNSystemOfEquations *system)
{
  // Add equations for every cluster
  for (int i = 0; i < NClusters(); i++) {
    GSVFeatureCluster *cluster = Cluster(i);

    // Add equation to align with its features
    for (int j = 0; j < cluster->NFeatures(); j++) {
      GSVFeature *feature = cluster->Feature(j);
      AddClusterEquations(system, cluster, feature, feature->score);
    }

    // Add equation to align with its clusters
    for (int j = 0; j < cluster->NClusters(); j++) {
      GSVFeatureCluster *cluster2 = cluster->Cluster(j);
      AddClusterEquations(system, cluster, cluster2, cluster2->score);
    }
  }
}



void GSVPoseOptimization::
AddClusterEquations(RNSystemOfEquations *system, 
  GSVFeatureCluster *cluster0, GSVFeature *feature1, 
  RNScalar w)
{
  // Compute the equations weights
  if (w == 0) return;
  RNScalar w3 = w * scan_point_scan_point_correspondence_weight;
  if (feature1->image) w3 = w * scan_point_image_point_correspondence_weight;
  RNScalar w2 = w * image_point_image_point_correspondence_reprojection_weight;
  RNScalar wd = w * scan_direction_scan_direction_correspondence_weight;

  // 3D distances
  if (w3 > 0) {
    // SCANPOINT <-> SCANPOINT
    // SCANPOINT <-> IMAGEPOINT
    if (((cluster0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
        ((cluster0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE))) {
      RNAlgebraic *px0, *py0, *pz0, *px1, *py1, *pz1;
      if (!ComputeTransformedPointCoordinates(cluster0, cluster0->scan_position, px0, py0, pz0)) return;
      if (!ComputeTransformedPointCoordinates(feature1, px1, py1, pz1)) return;
      AddPointPointDistanceEquations(system, px0, py0, pz0, px1, py1, pz1, w3);
    }

    // SCANLINE <- SCANPOINT
    // SCANLINE <- SCANLINE
    // SCANLINE <- IMAGEPOINT
    // SCANLINE <- IMAGELINE
    if (((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
        ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
        ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
        ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE))) {
      RNAlgebraic *px0, *py0, *pz0, *dx0, *dy0, *dz0, *px1, *py1, *pz1;
      R3Line line0(cluster0->scan_position, cluster0->scan_direction);
      R3Point center0 = feature1->scan_position; center0.Project(line0);
      if (!ComputeTransformedPointCoordinates(cluster0, cluster0->scan_position, px0, py0, pz0)) return;
      if (!ComputeTransformedVectorCoordinates(cluster0, cluster0->scan_direction, dx0, dy0, dz0)) return;
      if (!ComputeTransformedPointCoordinates(feature1, px1, py1, pz1)) return;
      AddPointLineDistanceEquations(system, px1, py1, pz1, px0, py0, pz0, dx0, dy0, dz0, w3);
    }

    // SCANPOINT -> SCANLINE
    // SCANLINE -> SCANLINE
    if (((cluster0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
        ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE))) {
      RNAlgebraic *px0, *py0, *pz0, *px1, *py1, *pz1, *dx1, *dy1, *dz1;
      R3Line line1(feature1->scan_position, feature1->scan_direction);
      R3Point center1 = cluster0->scan_position; center1.Project(line1);
      if (!ComputeTransformedPointCoordinates(feature1->scanline, center1, px1, py1, pz1)) return;
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, feature1->scan_direction, dx1, dy1, dz1)) return;
      if (!ComputeTransformedPointCoordinates(cluster0, cluster0->scan_position, px0, py0, pz0)) return;
      AddPointLineDistanceEquations(system, px0, py0, pz0, px1, py1, pz1, dx1, dy1, dz1, w3);
    }

    // SCANPLANE <- SCANPOINT
    // SCANPLANE <- SCANLINE
    // SCANPLANE <- SCANLINE
    // SCANPLANE <- IMAGEPOINT
    // SCANPLANE <- IMAGELINE
    if (((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
        ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
        ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
        ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) ||
        ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE))) {
      RNAlgebraic *px0, *py0, *pz0, *nx0, *ny0, *nz0, *px1, *py1, *pz1;
      if (!ComputeTransformedPointCoordinates(cluster0, cluster0->scan_position, px0, py0, pz0)) return;
      if (!ComputeTransformedVectorCoordinates(cluster0, cluster0->scan_normal, nx0, ny0, nz0)) return;
      if (!ComputeTransformedPointCoordinates(feature1, px1, py1, pz1)) return;
      AddPointPlaneDistanceEquations(system, px1, py1, pz1, px0, py0, pz0, nx0, ny0, nz0, w3);
    }

    // SCANPOINT -> SCANPLANE
    // SCANLINE  -> SCANPLANE
    // SCANPLANE -> SCANPLANE
    if (((cluster0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
        ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)  && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
        ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)  && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE))) {
      RNAlgebraic *px0, *py0, *pz0, *px1, *py1, *pz1, *nx1, *ny1, *nz1;
      R3Plane plane1(feature1->scan_position, feature1->scan_normal);
      R3Point center1 = cluster0->scan_position; center1.Project(plane1);
      if (!ComputeTransformedPointCoordinates(feature1->scanline, center1, px1, py1, pz1)) return;
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, feature1->scan_normal, nx1, ny1, nz1)) return;
      if (!ComputeTransformedPointCoordinates(cluster0, cluster0->scan_position, px0, py0, pz0)) return;
      AddPointPlaneDistanceEquations(system, px0, py0, pz0, px1, py1, pz1, nx1, ny1, nz1, w3);
    }
  }

  // 2D distances
  if (w2 > 0) {
    // SCANPOINT -> IMAGEPOINT
    if ((cluster0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (feature1->feature_type == GSV_IMAGE_POINT_FEATURE_TYPE)) {
      RNAlgebraic *ix0, *iy0;
      if (!ComputeTransformedImageCoordinates(feature1, cluster0, ix0, iy0)) return;
      AddPointPointDistanceEquations(system, ix0, iy0, feature1->image_position, w2);
    }

    // SCANPOINT -> IMAGELINE
    // SCANLINE -> IMAGELINE
    if (((feature1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (cluster0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE)) ||
        ((feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (cluster0->feature_type == GSV_IMAGE_LINE_FEATURE_TYPE))) {
      RNAlgebraic *ix0, *iy0;
      if (!ComputeTransformedImageCoordinates(feature1, cluster0, ix0, iy0)) return;
      AddPointLineDistanceEquations(system, ix0, iy0, feature1->image_position, feature1->image_direction, w2);
    }
  }

  // Direction compatibility
  if (wd > 0) {
    if ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) {
      // Keep vectors peripendicular to direction of feature1 perpendicular to direction of cluster0
      RNAlgebraic *dx0, *dy0, *dz0, *nx1, *ny1, *nz1;
      RNDimension dim1 = feature1->scan_direction.MinDimension();
      R3Vector normal1A = feature1->scan_direction % R3xyz_triad[dim1]; normal1A.Normalize();
      R3Vector normal1B = feature1->scan_direction % normal1A; normal1B.Normalize();
      if (!ComputeTransformedVectorCoordinates(cluster0, cluster0->scan_direction, dx0, dy0, dz0)) return;
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, normal1A, nx1, ny1, nz1)) return;
      AddVectorVectorDotProductEquations(system, dx0, dy0, dz0, nx1, ny1, nz1, 0.25 * wd);
      if (!ComputeTransformedVectorCoordinates(cluster0, cluster0->scan_direction, dx0, dy0, dz0)) return;
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, normal1B, nx1, ny1, nz1)) return;
      AddVectorVectorDotProductEquations(system, dx0, dy0, dz0, nx1, ny1, nz1, 0.25 * wd);
      
      // Keep vectors peripendicular to direction of cluster0 perpendicular to direction of feature1
      RNAlgebraic *dx1, *dy1, *dz1, *nx0, *ny0, *nz0;
      RNDimension dim0 = cluster0->scan_direction.MinDimension();
      R3Vector normal0A = cluster0->scan_direction % R3xyz_triad[dim0]; normal0A.Normalize();
      R3Vector normal0B = cluster0->scan_direction % normal0A; normal0B.Normalize();
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, feature1->scan_direction, dx1, dy1, dz1)) return;
      if (!ComputeTransformedVectorCoordinates(cluster0, normal0A, nx0, ny0, nz0)) return;
      AddVectorVectorDotProductEquations(system, dx1, dy1, dz1, nx0, ny0, nz0, 0.25 * wd);
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, feature1->scan_direction, dx1, dy1, dz1)) return;
      if (!ComputeTransformedVectorCoordinates(cluster0, normal0B, nx0, ny0, nz0)) return;
      AddVectorVectorDotProductEquations(system, dx1, dy1, dz1, nx0, ny0, nz0, 0.25 * wd);
    }

    if ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) {
      // Keep vectors peripendicular to direction of normal1 perpendicular to direction of normal0
      RNAlgebraic *nx0, *ny0, *nz0, *dx1, *dy1, *dz1;
      RNDimension dim1 = feature1->scan_normal.MinDimension();
      R3Vector direction1A = feature1->scan_normal % R3xyz_triad[dim1]; direction1A.Normalize();
      R3Vector direction1B = feature1->scan_normal % direction1A; direction1B.Normalize();
      if (!ComputeTransformedVectorCoordinates(cluster0, cluster0->scan_normal, nx0, ny0, nz0)) return;
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, direction1A, dx1, dy1, dz1)) return;
      AddVectorVectorDotProductEquations(system, nx0, ny0, nz0, dx1, dy1, dz1, 0.25 * wd);
      if (!ComputeTransformedVectorCoordinates(cluster0, cluster0->scan_normal, nx0, ny0, nz0)) return;
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, direction1B, dx1, dy1, dz1)) return;
      AddVectorVectorDotProductEquations(system, nx0, ny0, nz0, dx1, dy1, dz1, 0.25 * wd);
      
      // Keep vectors peripendicular to direction of normal0 perpendicular to direction of normal1
      RNAlgebraic *nx1, *ny1, *nz1, *dx0, *dy0, *dz0;
      RNDimension dim0 = cluster0->scan_normal.MinDimension();
      R3Vector direction0A = cluster0->scan_normal % R3xyz_triad[dim0]; direction0A.Normalize();
      R3Vector direction0B = cluster0->scan_normal % direction0A; direction0B.Normalize();
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, feature1->scan_normal, nx1, ny1, nz1)) return;
      if (!ComputeTransformedVectorCoordinates(cluster0, direction0A, dx0, dy0, dz0)) return;
      AddVectorVectorDotProductEquations(system, nx1, ny1, nz1, dx0, dy0, dz0, 0.25 * wd);
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, feature1->scan_normal, nx1, ny1, nz1)) return;
      if (!ComputeTransformedVectorCoordinates(cluster0, direction0B, dx0, dy0, dz0)) return;
      AddVectorVectorDotProductEquations(system, nx1, ny1, nz1, dx0, dy0, dz0, 0.25 * wd);
    }

    if ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) {
      // Keep direction of feature1 perpendicular to normal of cluster0
      RNAlgebraic *nx0, *ny0, *nz0, *dx1, *dy1, *dz1;
      if (!ComputeTransformedVectorCoordinates(cluster0, cluster0->scan_normal, nx0, ny0, nz0)) return;
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, feature1->scan_direction, dx1, dy1, dz1)) return;
      AddVectorVectorDotProductEquations(system, nx0, ny0, nz0, dx1, dy1, dz1, wd);
    }

    if ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (feature1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) {
      // Keep direction of cluster0 perpendicular to normal of feature1
      RNAlgebraic *dx0, *dy0, *dz0, *nx1, *ny1, *nz1;
      if (!ComputeTransformedVectorCoordinates(cluster0, cluster0->scan_direction, dx0, dy0, dz0)) return;
      if (!ComputeTransformedVectorCoordinates(feature1->scanline, feature1->scan_normal, nx1, ny1, nz1)) return;
      AddVectorVectorDotProductEquations(system, dx0, dy0, dz0, nx1, ny1, nz1, wd);
    }
  }
}



void GSVPoseOptimization::
AddClusterEquations(RNSystemOfEquations *system, 
  GSVFeatureCluster *cluster0, GSVFeatureCluster *cluster1, 
  RNScalar w)
{
  // Compute the equations weights
  if (w == 0) return;
  RNScalar w3 = w * scan_point_scan_point_correspondence_weight;
  RNScalar wd = w * scan_direction_scan_direction_correspondence_weight;

  // 3D distances
  if (w3 > 0) {
    // SCANPOINT <-> SCANPOINT
    if ((cluster0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) {
      RNAlgebraic *px0, *py0, *pz0, *px1, *py1, *pz1;
      if (!ComputeTransformedPointCoordinates(cluster0, cluster0->scan_position, px0, py0, pz0)) return;
      if (!ComputeTransformedPointCoordinates(cluster1, cluster1->scan_position, px1, py1, pz1)) return;
      AddPointPointDistanceEquations(system, px0, py0, pz0, px1, py1, pz1, w3);
    }

    // SCANLINE <- SCANPOINT
    // SCANLINE <- SCANLINE
    if (((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
        ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE))) {
      RNAlgebraic *px0, *py0, *pz0, *dx0, *dy0, *dz0, *px1, *py1, *pz1;
      R3Line line0(cluster0->scan_position, cluster0->scan_direction);
      R3Point center0 = cluster1->scan_position; center0.Project(line0);
      if (!ComputeTransformedPointCoordinates(cluster0, cluster0->scan_position, px0, py0, pz0)) return;
      if (!ComputeTransformedVectorCoordinates(cluster0, cluster0->scan_direction, dx0, dy0, dz0)) return;
      if (!ComputeTransformedPointCoordinates(cluster1, cluster1->scan_position, px1, py1, pz1)) return;
      AddPointLineDistanceEquations(system, px1, py1, pz1, px0, py0, pz0, dx0, dy0, dz0, w3);
    }

    // SCANPOINT -> SCANLINE
    // SCANLINE -> SCANLINE
    if (((cluster0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
        ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE))) {
      RNAlgebraic *px0, *py0, *pz0, *px1, *py1, *pz1, *dx1, *dy1, *dz1;
      R3Line line1(cluster1->scan_position, cluster1->scan_direction);
      R3Point center1 = cluster0->scan_position; center1.Project(line1);
      if (!ComputeTransformedPointCoordinates(cluster1, cluster1->scan_position, px1, py1, pz1)) return;
      if (!ComputeTransformedVectorCoordinates(cluster1, cluster1->scan_direction, dx1, dy1, dz1)) return;
      if (!ComputeTransformedPointCoordinates(cluster0, cluster0->scan_position, px0, py0, pz0)) return;
      AddPointLineDistanceEquations(system, px0, py0, pz0, px1, py1, pz1, dx1, dy1, dz1, w3);
    }

    // SCANPLANE <- SCANPOINT
    // SCANPLANE <- SCANLINE
    // SCANPLANE <- SCANLINE
    if (((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_POINT_FEATURE_TYPE)) ||
        ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) ||
        ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE))) {
      RNAlgebraic *px0, *py0, *pz0, *nx0, *ny0, *nz0, *px1, *py1, *pz1;
      if (!ComputeTransformedPointCoordinates(cluster0, cluster0->scan_position, px0, py0, pz0)) return;
      if (!ComputeTransformedVectorCoordinates(cluster0, cluster0->scan_normal, nx0, ny0, nz0)) return;
      if (!ComputeTransformedPointCoordinates(cluster1, cluster1->scan_position, px1, py1, pz1)) return;
      AddPointPlaneDistanceEquations(system, px1, py1, pz1, px0, py0, pz0, nx0, ny0, nz0, w3);
    }

    // SCANPOINT -> SCANPLANE
    // SCANLINE  -> SCANPLANE
    // SCANPLANE -> SCANPLANE
    if (((cluster0->feature_type == GSV_SCAN_POINT_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
        ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)  && (cluster1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) ||
        ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)  && (cluster1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE))) {
      RNAlgebraic *px0, *py0, *pz0, *px1, *py1, *pz1, *nx1, *ny1, *nz1;
      if (!ComputeTransformedPointCoordinates(cluster1, cluster1->scan_position, px1, py1, pz1)) return;
      if (!ComputeTransformedVectorCoordinates(cluster1, cluster1->scan_normal, nx1, ny1, nz1)) return;
      if (!ComputeTransformedPointCoordinates(cluster0, cluster0->scan_position, px0, py0, pz0)) return;
      AddPointPlaneDistanceEquations(system, px0, py0, pz0, px1, py1, pz1, nx1, ny1, nz1, w3);
    }
  }

  // Direction compatibility
  if (wd > 0) {
    if ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) {
      // Keep vectors peripendicular to direction of cluster1 perpendicular to direction of cluster0
      RNAlgebraic *dx0, *dy0, *dz0, *nx1, *ny1, *nz1;
      RNDimension dim1 = cluster1->scan_direction.MinDimension();
      R3Vector normal1A = cluster1->scan_direction % R3xyz_triad[dim1]; normal1A.Normalize();
      R3Vector normal1B = cluster1->scan_direction % normal1A; normal1B.Normalize();
      if (!ComputeTransformedVectorCoordinates(cluster0, cluster0->scan_direction, dx0, dy0, dz0)) return;
      if (!ComputeTransformedVectorCoordinates(cluster1, normal1A, nx1, ny1, nz1)) return;
      AddVectorVectorDotProductEquations(system, dx0, dy0, dz0, nx1, ny1, nz1, 0.25 * wd);
      if (!ComputeTransformedVectorCoordinates(cluster0, cluster0->scan_direction, dx0, dy0, dz0)) return;
      if (!ComputeTransformedVectorCoordinates(cluster1, normal1B, nx1, ny1, nz1)) return;
      AddVectorVectorDotProductEquations(system, dx0, dy0, dz0, nx1, ny1, nz1, 0.25 * wd);
      
      // Keep vectors peripendicular to direction of cluster0 perpendicular to direction of cluster1
      RNAlgebraic *dx1, *dy1, *dz1, *nx0, *ny0, *nz0;
      RNDimension dim0 = cluster0->scan_direction.MinDimension();
      R3Vector normal0A = cluster0->scan_direction % R3xyz_triad[dim0]; normal0A.Normalize();
      R3Vector normal0B = cluster0->scan_direction % normal0A; normal0B.Normalize();
      if (!ComputeTransformedVectorCoordinates(cluster1, cluster1->scan_direction, dx1, dy1, dz1)) return;
      if (!ComputeTransformedVectorCoordinates(cluster0, normal0A, nx0, ny0, nz0)) return;
      AddVectorVectorDotProductEquations(system, dx1, dy1, dz1, nx0, ny0, nz0, 0.25 * wd);
      if (!ComputeTransformedVectorCoordinates(cluster1, cluster1->scan_direction, dx1, dy1, dz1)) return;
      if (!ComputeTransformedVectorCoordinates(cluster0, normal0B, nx0, ny0, nz0)) return;
      AddVectorVectorDotProductEquations(system, dx1, dy1, dz1, nx0, ny0, nz0, 0.25 * wd);
    }

    if ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) {
      // Keep vectors peripendicular to direction of normal1 perpendicular to direction of normal0
      RNAlgebraic *nx0, *ny0, *nz0, *dx1, *dy1, *dz1;
      RNDimension dim1 = cluster1->scan_normal.MinDimension();
      R3Vector direction1A = cluster1->scan_normal % R3xyz_triad[dim1]; direction1A.Normalize();
      R3Vector direction1B = cluster1->scan_normal % direction1A; direction1B.Normalize();
      if (!ComputeTransformedVectorCoordinates(cluster0, cluster0->scan_normal, nx0, ny0, nz0)) return;
      if (!ComputeTransformedVectorCoordinates(cluster1, direction1A, dx1, dy1, dz1)) return;
      AddVectorVectorDotProductEquations(system, nx0, ny0, nz0, dx1, dy1, dz1, 0.25 * wd);
      if (!ComputeTransformedVectorCoordinates(cluster0, cluster0->scan_normal, nx0, ny0, nz0)) return;
      if (!ComputeTransformedVectorCoordinates(cluster1, direction1B, dx1, dy1, dz1)) return;
      AddVectorVectorDotProductEquations(system, nx0, ny0, nz0, dx1, dy1, dz1, 0.25 * wd);
      
      // Keep vectors peripendicular to direction of normal0 perpendicular to direction of normal1
      RNAlgebraic *nx1, *ny1, *nz1, *dx0, *dy0, *dz0;
      RNDimension dim0 = cluster0->scan_normal.MinDimension();
      R3Vector direction0A = cluster0->scan_normal % R3xyz_triad[dim0]; direction0A.Normalize();
      R3Vector direction0B = cluster0->scan_normal % direction0A; direction0B.Normalize();
      if (!ComputeTransformedVectorCoordinates(cluster1, cluster1->scan_normal, nx1, ny1, nz1)) return;
      if (!ComputeTransformedVectorCoordinates(cluster0, direction0A, dx0, dy0, dz0)) return;
      AddVectorVectorDotProductEquations(system, nx1, ny1, nz1, dx0, dy0, dz0, 0.25 * wd);
      if (!ComputeTransformedVectorCoordinates(cluster1, cluster1->scan_normal, nx1, ny1, nz1)) return;
      if (!ComputeTransformedVectorCoordinates(cluster0, direction0B, dx0, dy0, dz0)) return;
      AddVectorVectorDotProductEquations(system, nx1, ny1, nz1, dx0, dy0, dz0, 0.25 * wd);
    }

    if ((cluster0->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_LINE_FEATURE_TYPE)) {
      // Keep direction of cluster1 perpendicular to normal of cluster0
      RNAlgebraic *nx0, *ny0, *nz0, *dx1, *dy1, *dz1;
      if (!ComputeTransformedVectorCoordinates(cluster0, cluster0->scan_normal, nx0, ny0, nz0)) return;
      if (!ComputeTransformedVectorCoordinates(cluster1, cluster1->scan_direction, dx1, dy1, dz1)) return;
      AddVectorVectorDotProductEquations(system, nx0, ny0, nz0, dx1, dy1, dz1, wd);
    }

    if ((cluster0->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) && (cluster1->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE)) {
      // Keep direction of cluster0 perpendicular to normal of cluster1
      RNAlgebraic *dx0, *dy0, *dz0, *nx1, *ny1, *nz1;
      if (!ComputeTransformedVectorCoordinates(cluster0, cluster0->scan_direction, dx0, dy0, dz0)) return;
      if (!ComputeTransformedVectorCoordinates(cluster1, cluster1->scan_normal, nx1, ny1, nz1)) return;
      AddVectorVectorDotProductEquations(system, dx0, dy0, dz0, nx1, ny1, nz1, wd);
    }
  }
}



////////////////////////////////////////////////////////
// Mid-level equation utilty stuff
////////////////////////////////////////////////////////

void GSVPoseOptimization::
AddPointPointDistanceEquations(RNSystemOfEquations *system, 
  RNAlgebraic *px0, RNAlgebraic *py0, RNAlgebraic *pz0, 
  RNAlgebraic *px1, RNAlgebraic *py1, RNAlgebraic *pz1,
  RNScalar w)
{
  // Check weight
  if (w == 0) return;

  // Add equation representing differences between point coordinates
  // *ex = w * (px1 - px0);
  // *ey = w * (py1 - py0);
  // *ez = w * (pz1 - pz0);
  RNAlgebraic *ex = new RNAlgebraic(RN_SUBTRACT_OPERATION, px1, px0);
  RNAlgebraic *ey = new RNAlgebraic(RN_SUBTRACT_OPERATION, py1, py0);
  RNAlgebraic *ez = new RNAlgebraic(RN_SUBTRACT_OPERATION, pz1, pz0);
  ex->Multiply(w);
  ey->Multiply(w);
  ez->Multiply(w);
  system->InsertEquation(ex, w * scan_residual_threshold);
  system->InsertEquation(ey, w * scan_residual_threshold);
  system->InsertEquation(ez, w * scan_residual_threshold);
}



void GSVPoseOptimization::
AddPointLineDistanceEquations(RNSystemOfEquations *system, 
  RNAlgebraic *px0, RNAlgebraic *py0, RNAlgebraic *pz0, 
  RNAlgebraic *px1, RNAlgebraic *py1, RNAlgebraic *pz1,
  RNAlgebraic *dx1, RNAlgebraic *dy1, RNAlgebraic *dz1,
  RNScalar w)
{
  // Check weight
  if (w == 0) return;

  // Compute vector from point0 to point1
  RNAlgebraic *vx = new RNAlgebraic(RN_SUBTRACT_OPERATION, px1, px0);
  RNAlgebraic *vy = new RNAlgebraic(RN_SUBTRACT_OPERATION, py1, py0);
  RNAlgebraic *vz = new RNAlgebraic(RN_SUBTRACT_OPERATION, pz1, pz0);
  
  // Add equations for distance between line1 and point0
  // || d % v || = 0;
  // ex = w * (dy1*vz - dz1*vy);
  // ey = w * (dz1*vx - dx1*vz);
  // ez = w * (dx1*vy - dy1*vx);
  RNAlgebraic *ex, *ey, *ez, *dxA, *dxB, *dyA, *dyB, *dzA, *dzB;
  dxA = new RNAlgebraic(RN_MULTIPLY_OPERATION, new RNAlgebraic(*dy1), new RNAlgebraic(*vz));
  dxB = new RNAlgebraic(RN_MULTIPLY_OPERATION, new RNAlgebraic(*dz1), new RNAlgebraic(*vy));
  ex = new RNAlgebraic(RN_SUBTRACT_OPERATION, dxA, dxB);
  ex = new RNAlgebraic(RN_MULTIPLY_OPERATION, w, ex);
  dyA = new RNAlgebraic(RN_MULTIPLY_OPERATION, dz1, new RNAlgebraic(*vx));
  dyB = new RNAlgebraic(RN_MULTIPLY_OPERATION, new RNAlgebraic(*dx1), vz);
  ey = new RNAlgebraic(RN_SUBTRACT_OPERATION, dyA, dyB);
  ey = new RNAlgebraic(RN_MULTIPLY_OPERATION, w, ey);
  dzA = new RNAlgebraic(RN_MULTIPLY_OPERATION, dx1, vy);
  dzB = new RNAlgebraic(RN_MULTIPLY_OPERATION, dy1, vx);
  ez = new RNAlgebraic(RN_SUBTRACT_OPERATION, dzA, dzB);
  ez = new RNAlgebraic(RN_MULTIPLY_OPERATION, w, ez);
  system->InsertEquation(ex, w * scan_residual_threshold);
  system->InsertEquation(ey, w * scan_residual_threshold);
  system->InsertEquation(ez, w * scan_residual_threshold);
}



void GSVPoseOptimization::
AddPointPlaneDistanceEquations(RNSystemOfEquations *system, 
  RNAlgebraic *px0, RNAlgebraic *py0, RNAlgebraic *pz0, 
  RNAlgebraic *px1, RNAlgebraic *py1, RNAlgebraic *pz1,
  RNAlgebraic *nx1, RNAlgebraic *ny1, RNAlgebraic *nz1,
  RNScalar w)
{
  // Check weight
  if (w == 0) return;

  // Add equation representing distance from point0 to plane1
  // *equation = w * ( (px0 - px1)*nx1 + (py0 - py1)*ny1 + (pz0 - pz1)*nz1 );
  RNAlgebraic *dx, *dy, *dz, *d;
  dx = new RNAlgebraic(RN_SUBTRACT_OPERATION, px0, px1);
  dx = new RNAlgebraic(RN_MULTIPLY_OPERATION, dx, nx1);
  dy = new RNAlgebraic(RN_SUBTRACT_OPERATION, py0, py1);
  dy = new RNAlgebraic(RN_MULTIPLY_OPERATION, dy, ny1);
  dz = new RNAlgebraic(RN_SUBTRACT_OPERATION, pz0, pz1);
  dz = new RNAlgebraic(RN_MULTIPLY_OPERATION, dz, nz1);
  d = new RNAlgebraic(RN_ADD_OPERATION, dx, dy);
  d = new RNAlgebraic(RN_ADD_OPERATION, d, dz);
  d = new RNAlgebraic(RN_MULTIPLY_OPERATION, d, w);
  system->InsertEquation(d, w * scan_residual_threshold);
}



void GSVPoseOptimization::
AddPointPointDistanceEquations(RNSystemOfEquations *system, 
  RNAlgebraic *ix0, RNAlgebraic *iy0, 
  const R2Point& p1, RNScalar w)
{
  // Check weight
  if (w == 0) return;

  // Add equations for 2D point-point distance
  RNAlgebraic *ex = new RNAlgebraic(RN_SUBTRACT_OPERATION, ix0, p1.X());
  RNAlgebraic *ey = new RNAlgebraic(RN_SUBTRACT_OPERATION, iy0, p1.Y());
  ex = new RNAlgebraic(RN_MULTIPLY_OPERATION, ex, w);
  ey = new RNAlgebraic(RN_MULTIPLY_OPERATION, ey, w);
  system->InsertEquation(ex, w * image_residual_threshold);
  system->InsertEquation(ey, w * image_residual_threshold);
}



void GSVPoseOptimization::
AddPointLineDistanceEquations(RNSystemOfEquations *system, 
  RNAlgebraic *ix0, RNAlgebraic *iy0, 
  const R2Point& p1, const R2Vector& d1, RNScalar w)
{
  // Check weight
  if (w == 0) return;

  // Compute vector normal to line1
  RNScalar nx1 = -d1.Y();
  RNScalar ny1 = d1.X();

  // Add equation representing distance from point0 to plane1
  // *equation = w * ( (px0 - px1)*nx1 + (py0 - py1)*ny1 );
  RNAlgebraic *dx, *dy, *d;
  dx = new RNAlgebraic(RN_SUBTRACT_OPERATION, ix0, p1.X());
  dx = new RNAlgebraic(RN_MULTIPLY_OPERATION, dx, nx1);
  dy = new RNAlgebraic(RN_SUBTRACT_OPERATION, iy0, p1.Y());
  dy = new RNAlgebraic(RN_MULTIPLY_OPERATION, dy, ny1);
  d = new RNAlgebraic(RN_ADD_OPERATION, dx, dy);
  d = new RNAlgebraic(RN_MULTIPLY_OPERATION, d, w);
  system->InsertEquation(d, w * image_residual_threshold);
}



void GSVPoseOptimization::
AddVectorVectorDotProductEquations(RNSystemOfEquations *system, 
  RNAlgebraic *dx0, RNAlgebraic *dy0, RNAlgebraic *dz0, 
  RNAlgebraic *dx1, RNAlgebraic *dy1, RNAlgebraic *dz1,
  RNScalar w)
{
  // Check weight
  if (w == 0) return;

  // Add equation representing differences between point coordinates
  // e = w * (dx0*dx1 + dy0*dy1 + dz0*dz1)
  RNAlgebraic *mx = new RNAlgebraic(RN_MULTIPLY_OPERATION, dx0, dx1);
  RNAlgebraic *my = new RNAlgebraic(RN_MULTIPLY_OPERATION, dy0, dy1);
  RNAlgebraic *mz = new RNAlgebraic(RN_MULTIPLY_OPERATION, dz0, dz1);
  RNAlgebraic *e = new RNAlgebraic(RN_ADD_OPERATION, mx, my);
  e = new RNAlgebraic(RN_ADD_OPERATION, e, mz);
  e->Multiply(w);
  system->InsertEquation(e, w * dot_product_residual_threshold);
}



////////////////////////////////////////////////////////////////////////
// Low level equation functions
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
ComputeTransformedPointCoordinates(GSVFeature *feature,
  RNAlgebraic *& px, RNAlgebraic *& py, RNAlgebraic *& pz)
{
  // Check feature type
  switch (feature->feature_type) {
  case GSV_SCAN_POINT_FEATURE_TYPE:
  case GSV_SCAN_LINE_FEATURE_TYPE:
  case GSV_SCAN_PLANE_FEATURE_TYPE:
    if (!ComputeTransformedPointCoordinates(feature->scanline, feature->scan_position, px, py, pz)) return 0;
    break;

  case GSV_IMAGE_POINT_FEATURE_TYPE:
  case GSV_IMAGE_LINE_FEATURE_TYPE: {
    int pixel_feature_index = pixel_feature_indices[feature->index];
    if (!ComputeTransformedPointCoordinates(feature->image, feature->image_position, pixel_feature_index, px, py, pz)) return 0;
    break; }

  default:
    fprintf(stderr, "Unrecognized feature type: %d\n", feature->feature_type);
    return 0;
  }

  // Return success
  return 1;
}



int GSVPoseOptimization::
ComputeTransformedPointCoordinates(GSVFeatureCluster *cluster, const R3Point& point,
  RNAlgebraic *& px, RNAlgebraic *& py, RNAlgebraic *& pz)
{
  // Get cluster variables
  const R3Point& center = cluster->scan_position;
  
  // Start with center
  px = new RNAlgebraic(center.X(), 0); 
  py = new RNAlgebraic(center.Y(), 0); 
  pz = new RNAlgebraic(center.Z(), 0); 
  
  // Check if point is not at center of rotation
  if (!R3Contains(center, point)) {
    // Get vector from center to point
    R3Vector d = point - center;

    // Get rotation
    RNPolynomial rx, ry, rz;
    AddClusterVariableToEquation(&rx, cluster->index, RX); 
    AddClusterVariableToEquation(&ry, cluster->index, RY); 
    AddClusterVariableToEquation(&rz, cluster->index, RZ); 

    // Get rotated vector from center to point
    // dx =     d.X() - rz*d.Y() + ry*d.Z();
    // dy =  rz*d.X()      d.Y() - rx*d.Z();
    // dz = -ry*d.X() + rx*d.Y()      d.Z();
    RNAlgebraic *dx, *dy, *dz;
    RNAlgebraic *dxX = new RNAlgebraic(1.0, 0); dxX->Multiply( d.X());
    RNAlgebraic *dxY = new RNAlgebraic(rz, 0);  dxY->Multiply(-d.Y());
    RNAlgebraic *dxZ = new RNAlgebraic(ry, 0);  dxZ->Multiply( d.Z());
    RNAlgebraic *dyX = new RNAlgebraic(rz, 0);  dyX->Multiply( d.X());
    RNAlgebraic *dyY = new RNAlgebraic(1.0, 0); dyY->Multiply( d.Y());
    RNAlgebraic *dyZ = new RNAlgebraic(rx, 0);  dyZ->Multiply(-d.Z());
    RNAlgebraic *dzX = new RNAlgebraic(ry, 0);  dzX->Multiply(-d.X());
    RNAlgebraic *dzY = new RNAlgebraic(rx, 0);  dzY->Multiply( d.Y());
    RNAlgebraic *dzZ = new RNAlgebraic(1.0, 0); dzZ->Multiply( d.Z());
    dx = new RNAlgebraic(                      dxX);
    dx = new RNAlgebraic(RN_ADD_OPERATION, dx, dxY);
    dx = new RNAlgebraic(RN_ADD_OPERATION, dx, dxZ);
    dy = new RNAlgebraic(                      dyX);
    dy = new RNAlgebraic(RN_ADD_OPERATION, dy, dyY);
    dy = new RNAlgebraic(RN_ADD_OPERATION, dy, dyZ);
    dz = new RNAlgebraic(                      dzX);
    dz = new RNAlgebraic(RN_ADD_OPERATION, dz, dzY);
    dz = new RNAlgebraic(RN_ADD_OPERATION, dz, dzZ);

    // Add rotated vector
    px = new RNAlgebraic(RN_ADD_OPERATION, px, dx);
    py = new RNAlgebraic(RN_ADD_OPERATION, py, dy);
    pz = new RNAlgebraic(RN_ADD_OPERATION, pz, dz);
  }
  
  // Get translation
  RNPolynomial *tx = new RNPolynomial();
  RNPolynomial *ty = new RNPolynomial();
  RNPolynomial *tz = new RNPolynomial();
  if (cluster->feature_type == GSV_SCAN_PLANE_FEATURE_TYPE) {
    // Translation is defined by one variable in original normal direction
    RNPolynomial t;
    AddClusterVariableToEquation(&t, cluster->index, TX);  
    *tx = t * cluster->scan_normal.X();
    *ty = t * cluster->scan_normal.Y();
    *tz = t * cluster->scan_normal.Z();
  }
  else if (cluster->feature_type == GSV_SCAN_LINE_FEATURE_TYPE) {
    RNPolynomial t1, t2;
    AddClusterVariableToEquation(&t1, cluster->index, TX);  
    AddClusterVariableToEquation(&t2, cluster->index, TY);  
    RNDimension dim = cluster->scan_direction.MinDimension();
    R3Vector tmp = cluster->scan_direction % R3xyz_triad[dim];
    R3Vector axis1 = cluster->scan_direction % tmp; axis1.Normalize();
    R3Vector axis2 = cluster->scan_direction % axis1; axis2.Normalize();
    *tx = t1*axis1.X() + t2*axis2.X();
    *ty = t1*axis1.Y() + t2*axis2.Y();
    *tz = t1*axis1.Z() + t2*axis2.Z();
  }
  else {
    // Translation is defined by three variables
    AddClusterVariableToEquation(tx, cluster->index, TX); 
    AddClusterVariableToEquation(ty, cluster->index, TY); 
    AddClusterVariableToEquation(tz, cluster->index, TZ); 
  }

  // Add translation
  px = new RNAlgebraic(RN_ADD_OPERATION, px, tx);
  py = new RNAlgebraic(RN_ADD_OPERATION, py, ty);
  pz = new RNAlgebraic(RN_ADD_OPERATION, pz, tz);

  // Return success
  return 1;
}



int GSVPoseOptimization::
ComputeTransformedPointCoordinates(GSVScanline *scanline, const R3Point& point,
  RNAlgebraic *& px, RNAlgebraic *& py, RNAlgebraic *& pz)
{
  // Compute transformed point
  //          |   1  -crz  cry |           |   1  -rz   ry |
  // p' = c + |  crz   1  -crx | (v - c) + |  rz    1  -rx | (p - v)
  //          | -cry  crx   1  |           | -ry   rx   1  |

  ////////////////////////////////  

  // Get laser and path
  if (!scanline) return 0;
  ScanlineData *scanline_data = (ScanlineData *) scanline->Data();
  if (!scanline_data) return 0;
  int scanline_index = scanline_data->index;
  if (scanline_index < 0) return 0;
  GSVScan *scan = scanline->Scan();
  if (!scan) return 0;
  GSVLaser *laser = scan->Laser();
  if (!laser) return 0;
  LaserData *laser_data = (LaserData *) laser->Data();
  if (!laser_data) return 0;
  int laser_index = laser_data->index;
  if (laser_index < 0) return 0;
  GSVSegment *segment = scan->Segment();
  if (!segment) return 0;
  SegmentData *segment_data = (SegmentData *) segment->Data();
  if (!segment_data) return 0;
  GSVPath *path = segment_data->path;
  if (!path) return 0;

  // Get center coordinate system 
  GSVPose path_pose = path->Pose(scanline_data->path_parameter);
  const R3Point& center = path_pose.Viewpoint();
  R3Vector right = path_pose.Right();
  R3Vector towards = path_pose.Towards();
  R3Vector up = path_pose.Up();

  // Get scanline viewpoint
  const GSVPose& scanline_pose = scanline->Pose();
  const R3Point& viewpoint = scanline_pose.Viewpoint();

  // Get vector from viewpoint to point
  R3Vector d = point - viewpoint;


  ////////////////////////////////  

  // Get center translation
  RNPolynomial *ctx = new RNPolynomial();
  RNPolynomial *cty = new RNPolynomial();
  RNPolynomial *ctz = new RNPolynomial();
  AddScanlineVariableToEquation(ctx, scanline_index, TX); 
  AddScanlineVariableToEquation(cty, scanline_index, TY); 
  AddScanlineVariableToEquation(ctz, scanline_index, TZ); 

  // Get center rotation 
  RNPolynomial crx, cry, crz;
  AddScanlineVariableToEquation(&crx, scanline_index, RX); 
  AddScanlineVariableToEquation(&cry, scanline_index, RY); 
  AddScanlineVariableToEquation(&crz, scanline_index, RZ); 

  // Get laser translation
  RNPolynomial vtx, vty, vtz;
  AddLaserVariableToEquation(&vtx, laser_index, TX); 
  AddLaserVariableToEquation(&vty, laser_index, TY); 
  AddLaserVariableToEquation(&vtz, laser_index, TZ); 

  // Get combined center and laser rotation
  RNPolynomial vrx, vry, vrz;
  AddScanlineVariableToEquation(&vrx, scanline_index, RX); 
  AddScanlineVariableToEquation(&vry, scanline_index, RY); 
  AddScanlineVariableToEquation(&vrz, scanline_index, RZ); 
  AddLaserVariableToEquation(&vrx, laser_index, RX, 1.0, 1.0, TRUE); 
  AddLaserVariableToEquation(&vry, laser_index, RY, 1.0, 1.0, TRUE); 
  AddLaserVariableToEquation(&vrz, laser_index, RZ, 1.0, 1.0, TRUE); 


  ////////////////////////////////  

  // Get center location
  // cx = center.X(); 
  // cy = center.Y(); 
  // cz = center.Z(); 
  RNAlgebraic *cx = new RNAlgebraic(center.X(), 0);  
  RNAlgebraic *cy = new RNAlgebraic(center.Y(), 0);  
  RNAlgebraic *cz = new RNAlgebraic(center.Z(), 0);  

  // Add center translation
  // cx += ctx;
  // cy += cty;
  // cz += ctz;
  cx = new RNAlgebraic(RN_ADD_OPERATION, cx, new RNAlgebraic(ctx));
  cy = new RNAlgebraic(RN_ADD_OPERATION, cy, new RNAlgebraic(cty));
  cz = new RNAlgebraic(RN_ADD_OPERATION, cz, new RNAlgebraic(ctz));
  
  ////////////////////////////////  

  // Get offset vector (from center to viewpoint)
  // ox = viewpoint - center
  RNAlgebraic *ox = new RNAlgebraic(viewpoint.X() - center.X(), 0);
  RNAlgebraic *oy = new RNAlgebraic(viewpoint.Y() - center.Y(), 0);
  RNAlgebraic *oz = new RNAlgebraic(viewpoint.Z() - center.Z(), 0);

  // Translate offset vector 
  // ox += vtx*right.X() + vty*towards.X() + vtz*up.X();
  // oy += vtx*right.Y() + vty*towards.Y() + vtz*up.Y();
  // oz += vtx*right.Z() + vty*towards.Z() + vtz*up.Z();
  RNAlgebraic *vtxX = new RNAlgebraic(vtx, 0);  vtxX->Multiply(right.X());   
  RNAlgebraic *vtxY = new RNAlgebraic(vty, 0);  vtxY->Multiply(towards.X()); 
  RNAlgebraic *vtxZ = new RNAlgebraic(vtz, 0);  vtxZ->Multiply(up.X());      
  RNAlgebraic *vtyX = new RNAlgebraic(vtx, 0);  vtyX->Multiply(right.Y());   
  RNAlgebraic *vtyY = new RNAlgebraic(vty, 0);  vtyY->Multiply(towards.Y()); 
  RNAlgebraic *vtyZ = new RNAlgebraic(vtz, 0);  vtyZ->Multiply(up.Y());      
  RNAlgebraic *vtzX = new RNAlgebraic(vtx, 0);  vtzX->Multiply(right.Z());   
  RNAlgebraic *vtzY = new RNAlgebraic(vty, 0);  vtzY->Multiply(towards.Z()); 
  RNAlgebraic *vtzZ = new RNAlgebraic(vtz, 0);  vtzZ->Multiply(up.Z());      
  ox = new RNAlgebraic(RN_ADD_OPERATION, ox, vtxX);
  ox = new RNAlgebraic(RN_ADD_OPERATION, ox, vtxY);
  ox = new RNAlgebraic(RN_ADD_OPERATION, ox, vtxZ);
  oy = new RNAlgebraic(RN_ADD_OPERATION, oy, vtyX);
  oy = new RNAlgebraic(RN_ADD_OPERATION, oy, vtyY);
  oy = new RNAlgebraic(RN_ADD_OPERATION, oy, vtyZ);
  oz = new RNAlgebraic(RN_ADD_OPERATION, oz, vtzX);
  oz = new RNAlgebraic(RN_ADD_OPERATION, oz, vtzY);
  oz = new RNAlgebraic(RN_ADD_OPERATION, oz, vtzZ);

  // Rotate offset vector 
  // ox =      ox - crz*oy + cry*oz;   
  // oy =  crz*ox +     oy - crx*oz;
  // oz = -cry*ox + crx*oy +     oz;
  RNAlgebraic *crxY = new RNAlgebraic(*oy); crxY->Multiply(new RNAlgebraic(crz, 0)); crxY->Negate();
  RNAlgebraic *crxZ = new RNAlgebraic(*oz); crxZ->Multiply(new RNAlgebraic(cry, 0)); 
  RNAlgebraic *cryX = new RNAlgebraic(*ox); cryX->Multiply(new RNAlgebraic(crz, 0)); 
  RNAlgebraic *cryZ = new RNAlgebraic(*oz); cryZ->Multiply(new RNAlgebraic(crx, 0)); cryZ->Negate();
  RNAlgebraic *crzX = new RNAlgebraic(*ox); crzX->Multiply(new RNAlgebraic(cry, 0)); crzX->Negate();
  RNAlgebraic *crzY = new RNAlgebraic(*oy); crzY->Multiply(new RNAlgebraic(crx, 0)); 
  ox = new RNAlgebraic(RN_ADD_OPERATION, ox, crxY);
  ox = new RNAlgebraic(RN_ADD_OPERATION, ox, crxZ);
  oy = new RNAlgebraic(RN_ADD_OPERATION, oy, cryX);
  oy = new RNAlgebraic(RN_ADD_OPERATION, oy, cryZ);
  oz = new RNAlgebraic(RN_ADD_OPERATION, oz, crzX);
  oz = new RNAlgebraic(RN_ADD_OPERATION, oz, crzY);

  // Get viewpoint position (add center position and offset vector)
  // vx = cx + ox
  // vy = cy + oy
  // vz = cz + oz
  RNAlgebraic *vx = new RNAlgebraic(RN_ADD_OPERATION, cx, ox);
  RNAlgebraic *vy = new RNAlgebraic(RN_ADD_OPERATION, cy, oy);
  RNAlgebraic *vz = new RNAlgebraic(RN_ADD_OPERATION, cz, oz);

  ////////////////////////////////  
 
  // Get rotated ray direction
  // dx =     d.X() - rz*d.Y() + ry*d.Z();
  // dy =  rz*d.X()      d.Y() - rx*d.Z();
  // dz = -ry*d.X() + rx*d.Y()      d.Z();
  RNAlgebraic *dx, *dy, *dz;
  RNAlgebraic *dxX = new RNAlgebraic(1.0, 0);  dxX->Multiply( d.X());
  RNAlgebraic *dxY = new RNAlgebraic(vrz, 0);  dxY->Multiply(-d.Y());
  RNAlgebraic *dxZ = new RNAlgebraic(vry, 0);  dxZ->Multiply( d.Z());
  RNAlgebraic *dyX = new RNAlgebraic(vrz, 0);  dyX->Multiply( d.X());
  RNAlgebraic *dyY = new RNAlgebraic(1.0, 0);  dyY->Multiply( d.Y());
  RNAlgebraic *dyZ = new RNAlgebraic(vrx, 0);  dyZ->Multiply(-d.Z());
  RNAlgebraic *dzX = new RNAlgebraic(vry, 0);  dzX->Multiply(-d.X());
  RNAlgebraic *dzY = new RNAlgebraic(vrx, 0);  dzY->Multiply( d.Y());
  RNAlgebraic *dzZ = new RNAlgebraic(1.0, 0);  dzZ->Multiply( d.Z());
  dx = new RNAlgebraic(                      dxX);
  dx = new RNAlgebraic(RN_ADD_OPERATION, dx, dxY);
  dx = new RNAlgebraic(RN_ADD_OPERATION, dx, dxZ);
  dy = new RNAlgebraic(                      dyX);
  dy = new RNAlgebraic(RN_ADD_OPERATION, dy, dyY);
  dy = new RNAlgebraic(RN_ADD_OPERATION, dy, dyZ);
  dz = new RNAlgebraic(                      dzX);
  dz = new RNAlgebraic(RN_ADD_OPERATION, dz, dzY);
  dz = new RNAlgebraic(RN_ADD_OPERATION, dz, dzZ);

  // Add viewpoint position (v) to rotated ray vector (d)
  // px = vx + dx
  // py = vy + dy
  // pz = vz + dz
  px = new RNAlgebraic(RN_ADD_OPERATION, vx, dx);
  py = new RNAlgebraic(RN_ADD_OPERATION, vy, dy);
  pz = new RNAlgebraic(RN_ADD_OPERATION, vz, dz);

  // Return success
  return 1;
}



int GSVPoseOptimization::
ComputeTransformedPointCoordinates(GSVImage *image, const R3Point& point,
  RNAlgebraic *& px, RNAlgebraic *& py, RNAlgebraic *& pz)
{
  // Compute transformed point
  //          |   1  -crz  cry |           |   1  -rz   ry |
  // p' = c + |  crz   1  -crx | (v - c) + |  rz    1  -rx | (p - v)
  //          | -cry  crx   1  |           | -ry   rx   1  |

  ////////////////////////////////  

  // Get camera and path
  if (!image) return 0;
  ImageData *image_data = (ImageData *) image->Data();
  if (!image_data) return 0;
  int image_index = image_data->index;
  if (image_index < 0) return 0;
  GSVTapestry *tapestry = image->Tapestry();
  if (!tapestry) return 0;
  GSVCamera *camera = tapestry->Camera();
  if (!camera) return 0;
  CameraData *camera_data = (CameraData *) camera->Data();
  if (!camera_data) return 0;
  int camera_index = camera_data->index;
  if (camera_index < 0) return 0;
  GSVSegment *segment = tapestry->Segment();
  if (!segment) return 0;
  SegmentData *segment_data = (SegmentData *) segment->Data();
  if (!segment_data) return 0;
  GSVPath *path = segment_data->path;
  if (!path) return 0;

  // Get center coordinate system 
  GSVPose path_pose = path->Pose(image_data->path_parameter);
  const R3Point& center = path_pose.Viewpoint();
  R3Vector right = path_pose.Right();
  R3Vector towards = path_pose.Towards();
  R3Vector up = path_pose.Up();

  // Get image viewpoint
  const GSVPose& image_pose = image->Pose();
  const R3Point& viewpoint = image_pose.Viewpoint();

  // Get vector from viewpoint to point
  R3Vector d = point - viewpoint;

  ////////////////////////////////  

  // Get center translation
  RNPolynomial *ctx = new RNPolynomial();
  RNPolynomial *cty = new RNPolynomial();
  RNPolynomial *ctz = new RNPolynomial();
  AddImageVariableToEquation(ctx, image_index, TX); 
  AddImageVariableToEquation(cty, image_index, TY); 
  AddImageVariableToEquation(ctz, image_index, TZ); 

  // Get center rotation 
  RNPolynomial crx, cry, crz;
  AddImageVariableToEquation(&crx, image_index, RX); 
  AddImageVariableToEquation(&cry, image_index, RY); 
  AddImageVariableToEquation(&crz, image_index, RZ); 

  // Get camera translation
  RNPolynomial vtx, vty, vtz;
  AddCameraVariableToEquation(&vtx, camera_index, TX); 
  AddCameraVariableToEquation(&vty, camera_index, TY); 
  AddCameraVariableToEquation(&vtz, camera_index, TZ); 

  // Get combined center and camera rotation
  RNPolynomial vrx, vry, vrz;
  AddImageVariableToEquation(&vrx, image_index, RX); 
  AddImageVariableToEquation(&vry, image_index, RY); 
  AddImageVariableToEquation(&vrz, image_index, RZ); 
  AddCameraVariableToEquation(&vrx, camera_index, RX, 1.0, 1.0, TRUE); 
  AddCameraVariableToEquation(&vry, camera_index, RY, 1.0, 1.0, TRUE); 
  AddCameraVariableToEquation(&vrz, camera_index, RZ, 1.0, 1.0, TRUE); 


  ////////////////////////////////  

  // Get center location
  // cx = center.X(); 
  // cy = center.Y(); 
  // cz = center.Z(); 
  RNAlgebraic *cx = new RNAlgebraic(center.X(), 0);  
  RNAlgebraic *cy = new RNAlgebraic(center.Y(), 0);  
  RNAlgebraic *cz = new RNAlgebraic(center.Z(), 0);  

  // Add center translation
  // cx += ctx;
  // cy += cty;
  // cz += ctz;
  cx = new RNAlgebraic(RN_ADD_OPERATION, cx, new RNAlgebraic(ctx));
  cy = new RNAlgebraic(RN_ADD_OPERATION, cy, new RNAlgebraic(cty));
  cz = new RNAlgebraic(RN_ADD_OPERATION, cz, new RNAlgebraic(ctz));
  
  ////////////////////////////////  

  // Get offset vector (from center to viewpoint)
  // ox = viewpoint - center
  RNAlgebraic *ox = new RNAlgebraic(viewpoint.X() - center.X(), 0);
  RNAlgebraic *oy = new RNAlgebraic(viewpoint.Y() - center.Y(), 0);
  RNAlgebraic *oz = new RNAlgebraic(viewpoint.Z() - center.Z(), 0);

  // Translate offset vector 
  // ox += vtx*right.X() + vty*towards.X() + vtz*up.X();
  // oy += vtx*right.Y() + vty*towards.Y() + vtz*up.Y();
  // oz += vtx*right.Z() + vty*towards.Z() + vtz*up.Z();
  RNAlgebraic *vtxX = new RNAlgebraic(vtx, 0);  vtxX->Multiply(right.X());   
  RNAlgebraic *vtxY = new RNAlgebraic(vty, 0);  vtxY->Multiply(towards.X()); 
  RNAlgebraic *vtxZ = new RNAlgebraic(vtz, 0);  vtxZ->Multiply(up.X());      
  RNAlgebraic *vtyX = new RNAlgebraic(vtx, 0);  vtyX->Multiply(right.Y());   
  RNAlgebraic *vtyY = new RNAlgebraic(vty, 0);  vtyY->Multiply(towards.Y()); 
  RNAlgebraic *vtyZ = new RNAlgebraic(vtz, 0);  vtyZ->Multiply(up.Y());      
  RNAlgebraic *vtzX = new RNAlgebraic(vtx, 0);  vtzX->Multiply(right.Z());   
  RNAlgebraic *vtzY = new RNAlgebraic(vty, 0);  vtzY->Multiply(towards.Z()); 
  RNAlgebraic *vtzZ = new RNAlgebraic(vtz, 0);  vtzZ->Multiply(up.Z());      
  ox = new RNAlgebraic(RN_ADD_OPERATION, ox, vtxX);
  ox = new RNAlgebraic(RN_ADD_OPERATION, ox, vtxY);
  ox = new RNAlgebraic(RN_ADD_OPERATION, ox, vtxZ);
  oy = new RNAlgebraic(RN_ADD_OPERATION, oy, vtyX);
  oy = new RNAlgebraic(RN_ADD_OPERATION, oy, vtyY);
  oy = new RNAlgebraic(RN_ADD_OPERATION, oy, vtyZ);
  oz = new RNAlgebraic(RN_ADD_OPERATION, oz, vtzX);
  oz = new RNAlgebraic(RN_ADD_OPERATION, oz, vtzY);
  oz = new RNAlgebraic(RN_ADD_OPERATION, oz, vtzZ);

  // Rotate offset vector 
  // ox =      ox - crz*oy + cry*oz;   
  // oy =  crz*ox +     oy - crx*oz;
  // oz = -cry*ox + crx*oy +     oz;
  RNAlgebraic *crxY = new RNAlgebraic(*oy); crxY->Multiply(new RNAlgebraic(crz, 0)); crxY->Negate();
  RNAlgebraic *crxZ = new RNAlgebraic(*oz); crxZ->Multiply(new RNAlgebraic(cry, 0)); 
  RNAlgebraic *cryX = new RNAlgebraic(*ox); cryX->Multiply(new RNAlgebraic(crz, 0)); 
  RNAlgebraic *cryZ = new RNAlgebraic(*oz); cryZ->Multiply(new RNAlgebraic(crx, 0)); cryZ->Negate();
  RNAlgebraic *crzX = new RNAlgebraic(*ox); crzX->Multiply(new RNAlgebraic(cry, 0)); crzX->Negate();
  RNAlgebraic *crzY = new RNAlgebraic(*oy); crzY->Multiply(new RNAlgebraic(crx, 0)); 
  ox = new RNAlgebraic(RN_ADD_OPERATION, ox, crxY);
  ox = new RNAlgebraic(RN_ADD_OPERATION, ox, crxZ);
  oy = new RNAlgebraic(RN_ADD_OPERATION, oy, cryX);
  oy = new RNAlgebraic(RN_ADD_OPERATION, oy, cryZ);
  oz = new RNAlgebraic(RN_ADD_OPERATION, oz, crzX);
  oz = new RNAlgebraic(RN_ADD_OPERATION, oz, crzY);

  // Get viewpoint position (add center position and offset vector)
  // vx = cx + ox
  // vy = cy + oy
  // vz = cz + oz
  RNAlgebraic *vx = new RNAlgebraic(RN_ADD_OPERATION, cx, ox);
  RNAlgebraic *vy = new RNAlgebraic(RN_ADD_OPERATION, cy, oy);
  RNAlgebraic *vz = new RNAlgebraic(RN_ADD_OPERATION, cz, oz);

  ////////////////////////////////  
 
  // Get rotated ray direction
  // dx =     d.X() - rz*d.Y() + ry*d.Z();
  // dy =  rz*d.X()      d.Y() - rx*d.Z();
  // dz = -ry*d.X() + rx*d.Y()      d.Z();
  RNAlgebraic *dx, *dy, *dz;
  RNAlgebraic *dxX = new RNAlgebraic(1.0, 0);  dxX->Multiply( d.X());
  RNAlgebraic *dxY = new RNAlgebraic(vrz, 0);  dxY->Multiply(-d.Y());
  RNAlgebraic *dxZ = new RNAlgebraic(vry, 0);  dxZ->Multiply( d.Z());
  RNAlgebraic *dyX = new RNAlgebraic(vrz, 0);  dyX->Multiply( d.X());
  RNAlgebraic *dyY = new RNAlgebraic(1.0, 0);  dyY->Multiply( d.Y());
  RNAlgebraic *dyZ = new RNAlgebraic(vrx, 0);  dyZ->Multiply(-d.Z());
  RNAlgebraic *dzX = new RNAlgebraic(vry, 0);  dzX->Multiply(-d.X());
  RNAlgebraic *dzY = new RNAlgebraic(vrx, 0);  dzY->Multiply( d.Y());
  RNAlgebraic *dzZ = new RNAlgebraic(1.0, 0);  dzZ->Multiply( d.Z());
  dx = new RNAlgebraic(                      dxX);
  dx = new RNAlgebraic(RN_ADD_OPERATION, dx, dxY);
  dx = new RNAlgebraic(RN_ADD_OPERATION, dx, dxZ);
  dy = new RNAlgebraic(                      dyX);
  dy = new RNAlgebraic(RN_ADD_OPERATION, dy, dyY);
  dy = new RNAlgebraic(RN_ADD_OPERATION, dy, dyZ);
  dz = new RNAlgebraic(                      dzX);
  dz = new RNAlgebraic(RN_ADD_OPERATION, dz, dzY);
  dz = new RNAlgebraic(RN_ADD_OPERATION, dz, dzZ);

  ////////////////////////////////  
 
  // Add viewpoint position (v) to rotated ray vector (d)
  // px = vx + dx
  // py = vy + dy
  // pz = vz + dz
  px = new RNAlgebraic(RN_ADD_OPERATION, vx, dx);
  py = new RNAlgebraic(RN_ADD_OPERATION, vy, dy);
  pz = new RNAlgebraic(RN_ADD_OPERATION, vz, dz);

  // Return success
  return 1;
}



int GSVPoseOptimization::
ComputeTransformedPointCoordinates(GSVImage *image, const R2Point& point, int pixel_feature_index,
  RNAlgebraic *& px, RNAlgebraic *& py, RNAlgebraic *& pz)
{
  // Compute transformed point
  // p' = v' + t * d'
  // 
  // Compute transformed point
  //          |   1  -crz  cry |        
  // v' = c + |  crz   1  -crx | (v - c)
  //          | -cry  crx   1  |        
  //
  //      |  1  -rz   ry |
  // d' = |  rz   1  -rx | d
  //      | -ry  rx   1  |
  //
  ////////////////////////////////  

  ////////////////////////////////  

  // Get camera and path
  if (!image) return 0;
  if (pixel_feature_index < 0) return 0;
  R3Ray ray = image->RayThroughUndistortedPosition(point);
  ImageData *image_data = (ImageData *) image->Data();
  if (!image_data) return 0;
  int image_index = image_data->index;
  if (image_index < 0) return 0;
  GSVTapestry *tapestry = image->Tapestry();
  if (!tapestry) return 0;
  GSVCamera *camera = tapestry->Camera();
  if (!camera) return 0;
  CameraData *camera_data = (CameraData *) camera->Data();
  if (!camera_data) return 0;
  int camera_index = camera_data->index;
  if (camera_index < 0) return 0;
  GSVSegment *segment = tapestry->Segment();
  if (!segment) return 0;
  SegmentData *segment_data = (SegmentData *) segment->Data();
  if (!segment_data) return 0;
  GSVPath *path = segment_data->path;
  if (!path) return 0;

  // Get center coordinate system 
  GSVPose path_pose = path->Pose(image_data->path_parameter);
  const R3Point& center = path_pose.Viewpoint();
  R3Vector right = path_pose.Right();
  R3Vector towards = path_pose.Towards();
  R3Vector up = path_pose.Up();
  const R3Vector& d = ray.Vector();

  ////////////////////////////////  

  // Get center translation
  RNPolynomial *ctx = new RNPolynomial();
  RNPolynomial *cty = new RNPolynomial();
  RNPolynomial *ctz = new RNPolynomial();
  AddImageVariableToEquation(ctx, image_index, TX); 
  AddImageVariableToEquation(cty, image_index, TY); 
  AddImageVariableToEquation(ctz, image_index, TZ); 

  // Get center rotation 
  RNPolynomial crx, cry, crz;
  AddImageVariableToEquation(&crx, image_index, RX); 
  AddImageVariableToEquation(&cry, image_index, RY); 
  AddImageVariableToEquation(&crz, image_index, RZ); 

  // Get camera translation
  RNPolynomial vtx, vty, vtz;
  AddCameraVariableToEquation(&vtx, camera_index, TX); 
  AddCameraVariableToEquation(&vty, camera_index, TY); 
  AddCameraVariableToEquation(&vtz, camera_index, TZ); 

  // Get combined center and camera rotation
  RNPolynomial vrx, vry, vrz;
  AddImageVariableToEquation(&vrx, image_index, RX); 
  AddImageVariableToEquation(&vry, image_index, RY); 
  AddImageVariableToEquation(&vrz, image_index, RZ); 
  AddCameraVariableToEquation(&vrx, camera_index, RX, 1.0, 1.0, TRUE); 
  AddCameraVariableToEquation(&vry, camera_index, RY, 1.0, 1.0, TRUE); 
  AddCameraVariableToEquation(&vrz, camera_index, RZ, 1.0, 1.0, TRUE); 

  // Get length of ray vector
  RNPolynomial t;
  AddPixelFeatureVariableToEquation(&t, pixel_feature_index, T);

  ////////////////////////////////  

  // Get center location
  // cx = center.X(); 
  // cy = center.Y(); 
  // cz = center.Z(); 
  RNAlgebraic *cx = new RNAlgebraic(center.X(), 0);  
  RNAlgebraic *cy = new RNAlgebraic(center.Y(), 0);  
  RNAlgebraic *cz = new RNAlgebraic(center.Z(), 0);  

  // Add center translation
  // cx += ctx;
  // cy += cty;
  // cz += ctz;
  cx = new RNAlgebraic(RN_ADD_OPERATION, cx, new RNAlgebraic(ctx));
  cy = new RNAlgebraic(RN_ADD_OPERATION, cy, new RNAlgebraic(cty));
  cz = new RNAlgebraic(RN_ADD_OPERATION, cz, new RNAlgebraic(ctz));
  
  ////////////////////////////////  

  // Get offset vector (from center to viewpoint)
  // ox = viewpoint - center
  RNAlgebraic *ox = new RNAlgebraic(ray.Start().X() - center.X(), 0);
  RNAlgebraic *oy = new RNAlgebraic(ray.Start().Y() - center.Y(), 0);
  RNAlgebraic *oz = new RNAlgebraic(ray.Start().Z() - center.Z(), 0);

  // Translate offset vector 
  // ox += vtx*right.X() + vty*towards.X() + vtz*up.X();
  // oy += vtx*right.Y() + vty*towards.Y() + vtz*up.Y();
  // oz += vtx*right.Z() + vty*towards.Z() + vtz*up.Z();
  RNAlgebraic *vtxX = new RNAlgebraic(vtx, 0);  vtxX->Multiply(right.X());   
  RNAlgebraic *vtxY = new RNAlgebraic(vty, 0);  vtxY->Multiply(towards.X()); 
  RNAlgebraic *vtxZ = new RNAlgebraic(vtz, 0);  vtxZ->Multiply(up.X());      
  RNAlgebraic *vtyX = new RNAlgebraic(vtx, 0);  vtyX->Multiply(right.Y());   
  RNAlgebraic *vtyY = new RNAlgebraic(vty, 0);  vtyY->Multiply(towards.Y()); 
  RNAlgebraic *vtyZ = new RNAlgebraic(vtz, 0);  vtyZ->Multiply(up.Y());      
  RNAlgebraic *vtzX = new RNAlgebraic(vtx, 0);  vtzX->Multiply(right.Z());   
  RNAlgebraic *vtzY = new RNAlgebraic(vty, 0);  vtzY->Multiply(towards.Z()); 
  RNAlgebraic *vtzZ = new RNAlgebraic(vtz, 0);  vtzZ->Multiply(up.Z());      
  ox = new RNAlgebraic(RN_ADD_OPERATION, ox, vtxX);
  ox = new RNAlgebraic(RN_ADD_OPERATION, ox, vtxY);
  ox = new RNAlgebraic(RN_ADD_OPERATION, ox, vtxZ);
  oy = new RNAlgebraic(RN_ADD_OPERATION, oy, vtyX);
  oy = new RNAlgebraic(RN_ADD_OPERATION, oy, vtyY);
  oy = new RNAlgebraic(RN_ADD_OPERATION, oy, vtyZ);
  oz = new RNAlgebraic(RN_ADD_OPERATION, oz, vtzX);
  oz = new RNAlgebraic(RN_ADD_OPERATION, oz, vtzY);
  oz = new RNAlgebraic(RN_ADD_OPERATION, oz, vtzZ);

  // Rotate offset vector 
  // ox =      ox - crz*oy + cry*oz;   
  // oy =  crz*ox +     oy - crx*oz;
  // oz = -cry*ox + crx*oy +     oz;
  RNAlgebraic *crxY = new RNAlgebraic(*oy); crxY->Multiply(new RNAlgebraic(crz, 0)); crxY->Negate();
  RNAlgebraic *crxZ = new RNAlgebraic(*oz); crxZ->Multiply(new RNAlgebraic(cry, 0)); 
  RNAlgebraic *cryX = new RNAlgebraic(*ox); cryX->Multiply(new RNAlgebraic(crz, 0)); 
  RNAlgebraic *cryZ = new RNAlgebraic(*oz); cryZ->Multiply(new RNAlgebraic(crx, 0)); cryZ->Negate();
  RNAlgebraic *crzX = new RNAlgebraic(*ox); crzX->Multiply(new RNAlgebraic(cry, 0)); crzX->Negate();
  RNAlgebraic *crzY = new RNAlgebraic(*oy); crzY->Multiply(new RNAlgebraic(crx, 0)); 
  ox = new RNAlgebraic(RN_ADD_OPERATION, ox, crxY);
  ox = new RNAlgebraic(RN_ADD_OPERATION, ox, crxZ);
  oy = new RNAlgebraic(RN_ADD_OPERATION, oy, cryX);
  oy = new RNAlgebraic(RN_ADD_OPERATION, oy, cryZ);
  oz = new RNAlgebraic(RN_ADD_OPERATION, oz, crzX);
  oz = new RNAlgebraic(RN_ADD_OPERATION, oz, crzY);
  
  // Get viewpoint position (add center position and offset vector)
  // vx = cx + ox
  // vy = cy + oy
  // vz = cz + oz
  RNAlgebraic *vx = new RNAlgebraic(RN_ADD_OPERATION, cx, ox);
  RNAlgebraic *vy = new RNAlgebraic(RN_ADD_OPERATION, cy, oy);
  RNAlgebraic *vz = new RNAlgebraic(RN_ADD_OPERATION, cz, oz);

  ////////////////////////////////  

  // Get rotated ray direction
  // dx =      d.X() - vrz*d.Y() + vry*d.Z();
  // dy =  vrz*d.X()       d.Y() - vrx*d.Z();
  // dz = -vry*d.X() + vrx*d.Y()       d.Z();
  RNAlgebraic *dx, *dy, *dz;
  RNAlgebraic *dxX = new RNAlgebraic(1.0, 0);  dxX->Multiply( d.X());
  RNAlgebraic *dxY = new RNAlgebraic(vrz, 0);  dxY->Multiply(-d.Y());
  RNAlgebraic *dxZ = new RNAlgebraic(vry, 0);  dxZ->Multiply( d.Z());
  RNAlgebraic *dyX = new RNAlgebraic(vrz, 0);  dyX->Multiply( d.X());
  RNAlgebraic *dyY = new RNAlgebraic(1.0, 0);  dyY->Multiply( d.Y());
  RNAlgebraic *dyZ = new RNAlgebraic(vrx, 0);  dyZ->Multiply(-d.Z());
  RNAlgebraic *dzX = new RNAlgebraic(vry, 0);  dzX->Multiply(-d.X());
  RNAlgebraic *dzY = new RNAlgebraic(vrx, 0);  dzY->Multiply( d.Y());
  RNAlgebraic *dzZ = new RNAlgebraic(1.0, 0);  dzZ->Multiply( d.Z());
  dx = new RNAlgebraic(                      dxX);
  dx = new RNAlgebraic(RN_ADD_OPERATION, dx, dxY);
  dx = new RNAlgebraic(RN_ADD_OPERATION, dx, dxZ);
  dy = new RNAlgebraic(                      dyX);
  dy = new RNAlgebraic(RN_ADD_OPERATION, dy, dyY);
  dy = new RNAlgebraic(RN_ADD_OPERATION, dy, dyZ);
  dz = new RNAlgebraic(                      dzX);
  dz = new RNAlgebraic(RN_ADD_OPERATION, dz, dzY);
  dz = new RNAlgebraic(RN_ADD_OPERATION, dz, dzZ);

  // Get scaled and rotated ray vector
  // dx *= t;
  // dy *= t;
  // dz *= t;
  dx = new RNAlgebraic(RN_MULTIPLY_OPERATION, dx, new RNAlgebraic(t, 0));
  dy = new RNAlgebraic(RN_MULTIPLY_OPERATION, dy, new RNAlgebraic(t, 0));
  dz = new RNAlgebraic(RN_MULTIPLY_OPERATION, dz, new RNAlgebraic(t, 0));

  ////////////////////////////////  

  // Add viewpoint position (v) to scaled and rotated ray vector (d)
  // px = vx + dx
  // py = vy + dy
  // pz = vz + dz
  px = new RNAlgebraic(RN_ADD_OPERATION, vx, dx);
  py = new RNAlgebraic(RN_ADD_OPERATION, vy, dy);
  pz = new RNAlgebraic(RN_ADD_OPERATION, vz, dz);

  // Return success
  return 1;
}



int GSVPoseOptimization::
ComputeTransformedVectorCoordinates(GSVFeatureCluster *cluster, const R3Vector& d,
  RNAlgebraic *& dx, RNAlgebraic *& dy, RNAlgebraic *& dz)
{
  // Get rotation
  RNPolynomial rx, ry, rz;
  AddClusterVariableToEquation(&rx, cluster->index, RX); 
  AddClusterVariableToEquation(&ry, cluster->index, RY); 
  AddClusterVariableToEquation(&rz, cluster->index, RZ); 

  // Compute rotated vector coordinates
  // dx =     d.X() - rz*d.Y() + ry*d.Z();
  // dy =  rz*d.X() +    d.Y() - rx*d.Z();
  // dz = -ry*d.X() + rx*d.Y() +    d.Z();
  RNAlgebraic *dxX = new RNAlgebraic(1.0, 0); dxX->Multiply( d.X()); 
  RNAlgebraic *dxY = new RNAlgebraic(rz, 0);  dxY->Multiply(-d.Y()); 
  RNAlgebraic *dxZ = new RNAlgebraic(ry, 0);  dxZ->Multiply( d.Z()); 
  RNAlgebraic *dyX = new RNAlgebraic(rz, 0);  dyX->Multiply( d.X()); 
  RNAlgebraic *dyY = new RNAlgebraic(1.0, 0); dyY->Multiply( d.Y()); 
  RNAlgebraic *dyZ = new RNAlgebraic(rx, 0);  dyZ->Multiply(-d.Z()); 
  RNAlgebraic *dzX = new RNAlgebraic(ry, 0);  dzX->Multiply(-d.X()); 
  RNAlgebraic *dzY = new RNAlgebraic(rx, 0);  dzY->Multiply( d.Y()); 
  RNAlgebraic *dzZ = new RNAlgebraic(1.0, 0); dzZ->Multiply( d.Z()); 
  dx = new RNAlgebraic(                      dxX);
  dx = new RNAlgebraic(RN_ADD_OPERATION, dx, dxY);
  dx = new RNAlgebraic(RN_ADD_OPERATION, dx, dxZ);
  dy = new RNAlgebraic(                      dyX);
  dy = new RNAlgebraic(RN_ADD_OPERATION, dy, dyY);
  dy = new RNAlgebraic(RN_ADD_OPERATION, dy, dyZ);
  dz = new RNAlgebraic(                      dzX);
  dz = new RNAlgebraic(RN_ADD_OPERATION, dz, dzY);
  dz = new RNAlgebraic(RN_ADD_OPERATION, dz, dzZ);

  // Return success
  return 1;
}



int GSVPoseOptimization::
ComputeTransformedVectorCoordinates(GSVScanline *scanline, const R3Vector& d,
  RNAlgebraic *& dx, RNAlgebraic *& dy, RNAlgebraic *& dz)
{
  // Return vector after (linearized) rotation
  //        |  1  -rz  ry |
  // d' =   |  rz  1  -rx | d
  //        | -ry  rx  1  |

  // Get useful variables
  if (!scanline) return 0;
  ScanlineData *scanline_data = (ScanlineData *) scanline->Data();
  if (!scanline_data) return 0;
  int scanline_index = scanline_data->index;
  if (scanline_index < 0) return 0;
  GSVLaser *laser = scanline->Scan()->Laser();
  if (!laser) return 0;
  LaserData *laser_data = (LaserData *) laser->Data();
  if (!laser_data) return 0;
  int laser_index = laser_data->index;
  if (laser_index < 0) return 0;
  
  ////////////////////////////////  

  // Get combined path and laser rotations
  RNPolynomial vrx, vry, vrz;
  AddScanlineVariableToEquation(&vrx, scanline_index, RX); 
  AddScanlineVariableToEquation(&vry, scanline_index, RY); 
  AddScanlineVariableToEquation(&vrz, scanline_index, RZ); 
  AddLaserVariableToEquation(&vrx, laser_index, RX, 1.0, 1.0, TRUE); 
  AddLaserVariableToEquation(&vry, laser_index, RY, 1.0, 1.0, TRUE); 
  AddLaserVariableToEquation(&vrz, laser_index, RZ, 1.0, 1.0, TRUE); 

  ////////////////////////////////  

  // Compute rotated vector coordinates
  // dx =      d.X() - vrz*d.Y() + vry*d.Z();
  // dy =  vrz*d.X() +     d.Y() - vrx*d.Z();
  // dz = -vry*d.X() + vrx*d.Y() +     d.Z();
  RNAlgebraic *dxX = new RNAlgebraic(1.0, 0); dxX->Multiply( d.X()); 
  RNAlgebraic *dxY = new RNAlgebraic(vrz, 0); dxY->Multiply(-d.Y()); 
  RNAlgebraic *dxZ = new RNAlgebraic(vry, 0); dxZ->Multiply( d.Z()); 
  RNAlgebraic *dyX = new RNAlgebraic(vrz, 0); dyX->Multiply( d.X()); 
  RNAlgebraic *dyY = new RNAlgebraic(1.0, 0); dyY->Multiply( d.Y()); 
  RNAlgebraic *dyZ = new RNAlgebraic(vrx, 0); dyZ->Multiply(-d.Z()); 
  RNAlgebraic *dzX = new RNAlgebraic(vry, 0); dzX->Multiply(-d.X()); 
  RNAlgebraic *dzY = new RNAlgebraic(vrx, 0); dzY->Multiply( d.Y()); 
  RNAlgebraic *dzZ = new RNAlgebraic(1.0, 0); dzZ->Multiply( d.Z()); 
  dx = new RNAlgebraic(                      dxX);
  dx = new RNAlgebraic(RN_ADD_OPERATION, dx, dxY);
  dx = new RNAlgebraic(RN_ADD_OPERATION, dx, dxZ);
  dy = new RNAlgebraic(                      dyX);
  dy = new RNAlgebraic(RN_ADD_OPERATION, dy, dyY);
  dy = new RNAlgebraic(RN_ADD_OPERATION, dy, dyZ);
  dz = new RNAlgebraic(                      dzX);
  dz = new RNAlgebraic(RN_ADD_OPERATION, dz, dzY);
  dz = new RNAlgebraic(RN_ADD_OPERATION, dz, dzZ);

  // Return success
  return 1;
}



int GSVPoseOptimization::
ComputeTransformedVectorCoordinates(GSVImage *image, const R3Vector& d,
  RNAlgebraic *& dx, RNAlgebraic *& dy, RNAlgebraic *& dz)
{
  // Return vector after (linearized) rotation
  //        |  1  -rz  ry |
  // d' =   |  rz  1  -rx | d
  //        | -ry  rx  1  |

  // Get useful variables
  if (!image) return 0;
  ImageData *image_data = (ImageData *) image->Data();
  if (!image_data) return 0;
  int image_index = image_data->index;
  if (image_index < 0) return 0;
  GSVCamera *camera = image->Tapestry()->Camera();
  if (!camera) return 0;
  CameraData *camera_data = (CameraData *) camera->Data();
  if (!camera_data) return 0;
  int camera_index = camera_data->index;
  if (camera_index < 0) return 0;

  ////////////////////////////////  

  // Get combined path and camera rotations
  RNPolynomial vrx, vry, vrz;
  AddImageVariableToEquation(&vrx, image_index, RX); 
  AddImageVariableToEquation(&vry, image_index, RY); 
  AddImageVariableToEquation(&vrz, image_index, RZ); 
  AddCameraVariableToEquation(&vrx, camera_index, RX, 1.0, 1.0, TRUE); 
  AddCameraVariableToEquation(&vry, camera_index, RY, 1.0, 1.0, TRUE); 
  AddCameraVariableToEquation(&vrz, camera_index, RZ, 1.0, 1.0, TRUE); 

  ////////////////////////////////  

  // Compute rotated vector coordinates
  // dx =      d.X() - vrz*d.Y() + vry*d.Z();
  // dy =  vrz*d.X() +     d.Y() - vrx*d.Z();
  // dz = -vry*d.X() + vrx*d.Y() +     d.Z();
  RNAlgebraic *dxX = new RNAlgebraic(1.0, 0); dxX->Multiply( d.X()); 
  RNAlgebraic *dxY = new RNAlgebraic(vrz, 0); dxY->Multiply(-d.Y()); 
  RNAlgebraic *dxZ = new RNAlgebraic(vry, 0); dxZ->Multiply( d.Z()); 
  RNAlgebraic *dyX = new RNAlgebraic(vrz, 0); dyX->Multiply( d.X()); 
  RNAlgebraic *dyY = new RNAlgebraic(1.0, 0); dyY->Multiply( d.Y()); 
  RNAlgebraic *dyZ = new RNAlgebraic(vrx, 0); dyZ->Multiply(-d.Z()); 
  RNAlgebraic *dzX = new RNAlgebraic(vry, 0); dzX->Multiply(-d.X()); 
  RNAlgebraic *dzY = new RNAlgebraic(vrx, 0); dzY->Multiply( d.Y()); 
  RNAlgebraic *dzZ = new RNAlgebraic(1.0, 0); dzZ->Multiply( d.Z()); 
  dx = new RNAlgebraic(                      dxX);
  dx = new RNAlgebraic(RN_ADD_OPERATION, dx, dxY);
  dx = new RNAlgebraic(RN_ADD_OPERATION, dx, dxZ);
  dy = new RNAlgebraic(                      dyX);
  dy = new RNAlgebraic(RN_ADD_OPERATION, dy, dyY);
  dy = new RNAlgebraic(RN_ADD_OPERATION, dy, dyZ);
  dz = new RNAlgebraic(                      dzX);
  dz = new RNAlgebraic(RN_ADD_OPERATION, dz, dzY);
  dz = new RNAlgebraic(RN_ADD_OPERATION, dz, dzZ);

  // Return success
  return 1;
}



int GSVPoseOptimization::
ComputeTransformedImageCoordinates(GSVImage *image, int column_index,
  RNAlgebraic *px, RNAlgebraic *py, RNAlgebraic *pz,
  RNAlgebraic *&ix, RNAlgebraic *&iy)
{
  // Check column index
  if (column_index < 0) column_index = 0;
  if (column_index >= image->Width()) column_index = image->Width()-1;

  // Get image field of view
  RNScalar xfov = 0.5 * image->Camera()->XFov();
  RNScalar yfov = 0.5 * image->Camera()->YFov();
  if (!RNIsPositive(xfov) || !RNIsPositive(yfov)) return 0;

  // Get image pose
  const GSVPose& pose = image->Pose(column_index);
  RNAlgebraic *vx, *vy, *vz, *rx, *ry, *rz, *ux, *uy, *uz, *tx, *ty, *tz; 
  if (!ComputeTransformedPointCoordinates(image, pose.Viewpoint(), vx, vy, vz)) return 0;
  if (!ComputeTransformedVectorCoordinates(image, pose.Right(), rx, ry, rz)) return 0;
  if (!ComputeTransformedVectorCoordinates(image, pose.Up(), ux, uy, uz)) return 0;
  if (!ComputeTransformedVectorCoordinates(image, pose.Towards(), tx, ty, tz)) return 0;  

  // Compute vector from viewpoint to point
  // dx = px - vx;
  // dy = py - vy;
  // dz = pz - vz;
  RNAlgebraic *dx = new RNAlgebraic(RN_SUBTRACT_OPERATION, px, vx);
  RNAlgebraic *dy = new RNAlgebraic(RN_SUBTRACT_OPERATION, py, vy);
  RNAlgebraic *dz = new RNAlgebraic(RN_SUBTRACT_OPERATION, pz, vz);
  
  // Compute camera coordinates
  // cx = rx*dx + ry*dy + rz*dz;
  // cy = ux*dx + uy*dy + uz*dz;
  // cz = tx*dx + ty*dy + tz*dz;
  RNAlgebraic *cx, *cy, *cz;
  RNAlgebraic *cxX = new RNAlgebraic(RN_MULTIPLY_OPERATION, rx, new RNAlgebraic(*dx)); 
  RNAlgebraic *cxY = new RNAlgebraic(RN_MULTIPLY_OPERATION, ry, new RNAlgebraic(*dy)); 
  RNAlgebraic *cxZ = new RNAlgebraic(RN_MULTIPLY_OPERATION, rz, new RNAlgebraic(*dz)); 
  RNAlgebraic *cyX = new RNAlgebraic(RN_MULTIPLY_OPERATION, ux, new RNAlgebraic(*dx)); 
  RNAlgebraic *cyY = new RNAlgebraic(RN_MULTIPLY_OPERATION, uy, new RNAlgebraic(*dy)); 
  RNAlgebraic *cyZ = new RNAlgebraic(RN_MULTIPLY_OPERATION, uz, new RNAlgebraic(*dz)); 
  RNAlgebraic *czX = new RNAlgebraic(RN_MULTIPLY_OPERATION, tx, dx); 
  RNAlgebraic *czY = new RNAlgebraic(RN_MULTIPLY_OPERATION, ty, dy); 
  RNAlgebraic *czZ = new RNAlgebraic(RN_MULTIPLY_OPERATION, tz, dz); 
  cx = new RNAlgebraic(RN_ADD_OPERATION, cxX, cxY); 
  cx = new RNAlgebraic(RN_ADD_OPERATION, cx,  cxZ); 
  cy = new RNAlgebraic(RN_ADD_OPERATION, cyX, cyY); 
  cy = new RNAlgebraic(RN_ADD_OPERATION, cy,  cyZ); 
  cz = new RNAlgebraic(RN_ADD_OPERATION, czX, czY); 
  cz = new RNAlgebraic(RN_ADD_OPERATION, cz,  czZ); 

  // Compute image coordinates
  // ix = (0.5*image->Width()/tan(xfov)) * cx / cz;
  // iy = (0.5*image->Height()/tan(yfov)) * cy / cz;
  ix = new RNAlgebraic(RN_MULTIPLY_OPERATION, 0.5*image->Width()/tan(xfov), cx);
  ix = new RNAlgebraic(RN_DIVIDE_OPERATION, ix, new RNAlgebraic(*cz));
  ix = new RNAlgebraic(RN_ADD_OPERATION, ix, 0.5*image->Width());
  iy = new RNAlgebraic(RN_MULTIPLY_OPERATION, 0.5*image->Height()/tan(yfov), cy);
  iy = new RNAlgebraic(RN_DIVIDE_OPERATION, iy, cz);
  iy = new RNAlgebraic(RN_ADD_OPERATION, iy, 0.5*image->Height());

  // Return success
  return 1;
}



int GSVPoseOptimization::
ComputeTransformedImageCoordinates(GSVFeature *image_feature, GSVFeature *feature, 
  RNAlgebraic *&ix, RNAlgebraic *&iy)
{
  // Return projected image position of feature
  if (image_feature->image == feature->image) {
    ix = new RNAlgebraic(feature->image_position.X(), 0);
    iy = new RNAlgebraic(feature->image_position.Y(), 0);
  }
  else {
    RNAlgebraic *px, *py, *pz;
    assert(image_feature->image);
    int column_index = (int) image_feature->image_position.X();
    if (!ComputeTransformedPointCoordinates(feature, px, py, pz)) return 0;
    if (!ComputeTransformedImageCoordinates(image_feature->image, column_index, px, py, pz, ix, iy)) return 0;
  }

#if 0
  // For debugging
  int n = 0;
  n += laser_nv * lasers.NEntries();
  n += camera_nv * cameras.NEntries();
  n += vertex_nv * vertices.NEntries();
  n += pixel_feature_nv * pixel_features.NEntries();
  n += cluster_nv * clusters.NEntries();
  double *x = new double [ n ];
  for (int i = 0; i < n; i++) x[i] = 0;
  InitializeOptimizationVariables(x);
  RNScalar ixa = ix->Evaluate(x);
  RNScalar iya = iy->Evaluate(x);
  delete [] x;
  R3Point p = WorldPosition(feature);
  p.InverseTransform(OptimizedTransformation(image_feature));
  int column_index = (int) image_feature->image_position.X();
  R2Point b = image_feature->image->UndistortedPosition(p, column_index);
  RNScalar ixb = b.X();
  RNScalar iyb = b.Y();

  printf("A %12.6f %12.6f\n", ixa, iya);
  printf("B %12.6f %12.6f\n", ixb, iyb);
#endif  

  // Return success
  return 1;
}



int GSVPoseOptimization::
ComputeTransformedImageCoordinates(GSVFeature *image_feature, GSVFeatureCluster *cluster, 
  RNAlgebraic *&ix, RNAlgebraic *&iy)
{
  RNAlgebraic *px, *py, *pz;
  assert(image_feature->image);
  int column_index = (int) image_feature->image_position.X();
  if (!ComputeTransformedPointCoordinates(cluster, cluster->scan_position, px, py, pz)) return 0;
  if (!ComputeTransformedImageCoordinates(image_feature->image, column_index, px, py, pz, ix, iy)) return 0;

#if 0
  // For debugging
  int n = 0;
  n += laser_nv * lasers.NEntries();
  n += camera_nv * cameras.NEntries();
  n += vertex_nv * vertices.NEntries();
  n += pixel_feature_nv * pixel_features.NEntries();
  n += cluster_nv * clusters.NEntries();
  double *x = new double [ n ];
  for (int i = 0; i < n; i++) x[i] = 0;
  InitializeOptimizationVariables(x);
  RNScalar ixa = ix->Evaluate(x);
  RNScalar iya = iy->Evaluate(x);
  delete [] x;
  printf("a %12.6f %12.6f %12.6f   %12.6f %12.6f %12.6f\n",
    cluster->scan_position.X(), cluster->scan_position.Y(), cluster->scan_position.Z(),
    cluster->translation.X(), cluster->translation.Y(), cluster->translation.Z());
  R3Point p = WorldPosition(cluster);
  printf("b %12.6f %12.6f %12.6f\n", p.X(), p.Y(), p.Z());
  p.InverseTransform(OptimizedTransformation(image_feature));
  printf("c %12.6f %12.6f %12.6f\n", p.X(), p.Y(), p.Z());
  R2Point b = image_feature->image->UndistortedPosition(p, column_index);
  RNScalar ixb = b.X();
  RNScalar iyb = b.Y();

  printf("A %12.6f %12.6f\n", ixa, iya);
  printf("B %12.6f %12.6f\n", ixb, iyb);
#endif  

  // Return success
  return 1;
}



int GSVPoseOptimization::
ComputeTransformedPixelDepth(GSVImage *image, int pixel_feature_index, RNAlgebraic *&t)
{
  // Get variable
  RNPolynomial *p = new RNPolynomial();
  AddPixelFeatureVariableToEquation(p, pixel_feature_index, T); 
  t = new RNAlgebraic(p);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Very low level equation functions
////////////////////////////////////////////////////////////////////////

void GSVPoseOptimization::
AddLaserVariableToEquation(RNPolynomial *equation, 
  int laser_index, int variable_index, 
  double a, double e, RNBoolean already_unique)
{
  // Check inputs
  if (a == 0.0) return;
  assert(equation);
  assert(laser_index >= 0);
  assert(variable_index >= 0);
  assert(variable_index < LASER_NV);

  // Add variable term to equation
  int index = LaserVariableIndex(laser_index, variable_index);
  if (index != -1) equation->AddTerm(a, index, e, already_unique);
  else equation->AddTerm(a * LaserVariableValue(laser_index, variable_index));
}



void GSVPoseOptimization::
AddCameraVariableToEquation(RNPolynomial *equation, 
  int camera_index, int variable_index, 
  double a, double e, RNBoolean already_unique)
{
  // Check inputs
  if (a == 0.0) return;
  assert(equation);
  assert(camera_index >= 0);
  assert(variable_index >= 0);
  assert(variable_index < CAMERA_NV);

  // Add variable term to equation
  int index = CameraVariableIndex(camera_index, variable_index);
  if (index != -1) equation->AddTerm(a, index, e, already_unique);
  else equation->AddTerm(a * CameraVariableValue(camera_index, variable_index));
}



void GSVPoseOptimization::
AddPathVertexVariableToEquation(RNPolynomial *equation, 
  int vertex_index, int variable_index, 
  double a, double e, RNBoolean already_unique)
{
  // Check inputs
  if (a == 0.0) return;
  assert(equation);
  assert(vertex_index >= 0);
  assert(vertex_index < vertices.NEntries());
  assert(variable_index >= 0);
  assert(variable_index < VERTEX_NV);

  // Check if variable is constrained
  int index = VertexVariableIndex(vertex_index, variable_index);
  if (index != -1) equation->AddTerm(a, index, e, already_unique);
  else equation->AddTerm(a * VertexVariableValue(vertex_index, variable_index));
}



void GSVPoseOptimization::
AddPathVertexVariableToEquation(RNPolynomial *equation, 
  GSVPath *path, RNScalar path_parameter, int variable_index, 
  double a, double e, RNBoolean already_unique)
{
  // Check inputs
  if (a == 0.0) return;
  assert(equation);
  assert(variable_index >= 0);
  assert(variable_index < VERTEX_NV);
  assert(e == 1.0);

  // Get useful variables
  R3CatmullRomSpline *spline = path->spline;
  if (!spline) return;
  RNScalar u = path_parameter;
  int iu = spline->VertexIndex(u);
  GSVPathVertex *vertex = path->Vertex(iu);
  int k1 = vertex->index;
  int k0 = k1 - 1;
  int k2 = k1 + 1;
  int k3 = k1 + 2;
  
  // Compute blending weights
  RNScalar t = u - iu;
  RNScalar b0 = spline->BlendingWeight(t, -1);
  RNScalar b1 = spline->BlendingWeight(t, 0);
  RNScalar b2 = spline->BlendingWeight(t, 1);
  RNScalar b3 = spline->BlendingWeight(t, 2);
  if (iu == 0) { b1 += b0; b0 = 0; k0 = k1; }
  else if (iu == spline->NVertices()-1) { b1 += b2; b1 += b3; b2 = 0; b3 = 0; k2 = k1; k3 = k1; }
  else if (iu == spline->NVertices()-2) { b2 += b3; b3 = 0; k3 = k2; }

  // Add terms to equation
  if (b0 != 0.0) AddPathVertexVariableToEquation(equation, k0, variable_index, a * b0, e, already_unique);
  if (b1 != 0.0) AddPathVertexVariableToEquation(equation, k1, variable_index, a * b1, e, already_unique);
  if (b2 != 0.0) AddPathVertexVariableToEquation(equation, k2, variable_index, a * b2, e, already_unique);
  if (b3 != 0.0) AddPathVertexVariableToEquation(equation, k3, variable_index, a * b3, e, already_unique);
}



void GSVPoseOptimization::
AddScanlineVariableToEquation(RNPolynomial *equation, 
  int scanline_index, int variable_index, 
  double a, double e, RNBoolean already_unique)
{
  // Check inputs
  if (a == 0.0) return;
  assert(equation);
  assert(scanline_index >= 0);
  assert(scanline_index < scanlines.NEntries());
  assert(variable_index >= 0);
  assert(variable_index < VERTEX_NV);
  assert(e == 1.0);

  // Get useful variables
  GSVScanline *scanline = scanlines.Kth(scanline_index);
  ScanlineData *scanline_data = (ScanlineData *) scanline->Data();
  if (!scanline_data) return;
  GSVScan *scan = scanline->Scan();
  if (!scan) return;
  GSVSegment *segment = scan->Segment();
  if (!segment) return;
  SegmentData *segment_data = (SegmentData *) segment->Data();
  if (!segment_data) return;
  GSVPath *path = segment_data->path;
  if (!path) return;
  R3CatmullRomSpline *spline = path->spline;
  if (!spline) return;
  RNScalar u = scanline_data->path_parameter;
  int iu = spline->VertexIndex(u);
  GSVPathVertex *vertex = path->Vertex(iu);
  int k1 = vertex->index;
  int k0 = k1 - 1;
  int k2 = k1 + 1;
  int k3 = k1 + 2;
  
  // Compute blending weights
  RNScalar t = u - iu;
  RNScalar b0 = spline->BlendingWeight(t, -1);
  RNScalar b1 = spline->BlendingWeight(t, 0);
  RNScalar b2 = spline->BlendingWeight(t, 1);
  RNScalar b3 = spline->BlendingWeight(t, 2);
  if (iu == 0) { b1 += b0; b0 = 0; k0 = k1; }
  else if (iu == spline->NVertices()-1) { b1 += b2; b1 += b3; b2 = 0; b3 = 0; k2 = k1; k3 = k1; }
  else if (iu == spline->NVertices()-2) { b2 += b3; b3 = 0; k3 = k2; }

  // Add terms to equation
  if (b0 != 0.0) AddPathVertexVariableToEquation(equation, k0, variable_index, a * b0, e, already_unique);
  if (b1 != 0.0) AddPathVertexVariableToEquation(equation, k1, variable_index, a * b1, e, already_unique);
  if (b2 != 0.0) AddPathVertexVariableToEquation(equation, k2, variable_index, a * b2, e, already_unique);
  if (b3 != 0.0) AddPathVertexVariableToEquation(equation, k3, variable_index, a * b3, e, already_unique);
}



void GSVPoseOptimization::
AddImageVariableToEquation(RNPolynomial *equation, 
  int image_index, int variable_index, 
  double a, double e, RNBoolean already_unique)
{
  // Check inputs
  if (a == 0.0) return;
  assert(equation);
  assert(image_index >= 0);
  assert(image_index < images.NEntries());
  assert(variable_index >= 0);
  assert(variable_index < VERTEX_NV);
  assert(e == 1.0);

  // Get useful variables
  GSVImage *image = images.Kth(image_index);
  ImageData *image_data = (ImageData *) image->Data();
  if (!image_data) return;
  GSVPanorama *panorama = image->Panorama();
  if (!panorama) return;
  GSVSegment *segment = panorama->Segment();
  if (!segment) return;
  SegmentData *segment_data = (SegmentData *) segment->Data();
  if (!segment_data) return;
  GSVPath *path = segment_data->path;
  if (!path) return;
  R3CatmullRomSpline *spline = path->spline;
  if (!spline) return;
  RNScalar u = image_data->path_parameter;
  int iu = spline->VertexIndex(u);
  GSVPathVertex *vertex = path->Vertex(iu);
  int k1 = vertex->index;
  int k0 = k1 - 1;
  int k2 = k1 + 1;
  int k3 = k1 + 2;
  
  // Compute blending weights
  RNScalar t = u - iu;
  RNScalar b0 = spline->BlendingWeight(t, -1);
  RNScalar b1 = spline->BlendingWeight(t, 0);
  RNScalar b2 = spline->BlendingWeight(t, 1);
  RNScalar b3 = spline->BlendingWeight(t, 2);
  if (iu == 0) { b1 += b0; b0 = 0; k0 = k1; }
  else if (iu == spline->NVertices()-1) { b1 += b2; b1 += b3; b2 = 0; b3 = 0; k2 = k1; k3 = k1; }
  else if (iu == spline->NVertices()-2) { b2 += b3; b3 = 0; k3 = k2; }

  // Add terms to equation
  if (b0 != 0.0) AddPathVertexVariableToEquation(equation, k0, variable_index, a * b0, e, already_unique);
  if (b1 != 0.0) AddPathVertexVariableToEquation(equation, k1, variable_index, a * b1, e, already_unique);
  if (b2 != 0.0) AddPathVertexVariableToEquation(equation, k2, variable_index, a * b2, e, already_unique);
  if (b3 != 0.0) AddPathVertexVariableToEquation(equation, k3, variable_index, a * b3, e, already_unique);
}



void GSVPoseOptimization::
AddPixelFeatureVariableToEquation(RNPolynomial *equation, 
  int pixel_feature_index, int variable_index, 
  double a, double e, RNBoolean already_unique)
{
  // Check inputs
  if (a == 0.0) return;
  assert(equation);
  assert(pixel_feature_index >= 0);
  assert(variable_index >= 0);
  assert(variable_index < PIXEL_FEATURE_NV);

  // Add variable term to equation
  int index = PixelFeatureVariableIndex(pixel_feature_index, variable_index);
  if (index != -1) equation->AddTerm(a, index, e, already_unique);
  else equation->AddTerm(a * PixelFeatureVariableValue(pixel_feature_index, variable_index));
}



void GSVPoseOptimization::
AddClusterVariableToEquation(RNPolynomial *equation, 
  int cluster_index, int variable_index, 
  double a, double e, RNBoolean already_unique)
{
  // Check inputs
  if (a == 0.0) return;
  assert(equation);
  assert(cluster_index >= 0);
  assert(variable_index >= 0);
  assert(variable_index < CLUSTER_NV);

  // Add variable term to equation
  int index = ClusterVariableIndex(cluster_index, variable_index);
  if (index != -1) equation->AddTerm(a, index, e, already_unique);
  else equation->AddTerm(a * ClusterVariableValue(cluster_index, variable_index));
}



////////////////////////////////////////////////////////////////////////
// Variable access utility functions
////////////////////////////////////////////////////////////////////////

int GSVPoseOptimization::
LaserVariableIndex(int laser_index, int variable_index)
{
  // Check if laser variable is active
  if (laser_nv == 0) return -1;
  if (laser_v[variable_index] == -1) return -1;

  // Check if variable is fixed
  GSVLaser *laser = lasers.Kth(laser_index);
  LaserData *laser_data = (LaserData *) laser->Data();
  if (!laser_data) return -1;
  if (laser_data->inertia[variable_index] >= 1.0) return -1;

  // Return index of laser variable
  return laser_index*laser_nv + laser_v[variable_index];
}



int GSVPoseOptimization::
CameraVariableIndex(int camera_index, int variable_index)
{
  // Check if camera variable is active
  if (camera_nv == 0) return -1;
  if (camera_v[variable_index] == -1) return -1;

  // Check if variable is fixed
  GSVCamera *camera = cameras.Kth(camera_index);
  CameraData *camera_data = (CameraData *) camera->Data();
  if (!camera_data) return -1;
  if (camera_data->inertia[variable_index] >= 1.0) return -1;

  // Return index of camera variable
  return lasers.NEntries()*laser_nv + 
         camera_index*camera_nv + camera_v[variable_index];
}



int GSVPoseOptimization::
VertexVariableIndex(int vertex_index, int variable_index)
{
  // Check if cluster variable is active
  if (vertex_nv == 0) return -1;
  if (vertex_v[variable_index] == -1) return -1;

  // Check if variable is fixed
  GSVPathVertex *vertex = vertices.Kth(vertex_index);
  if (vertex->inertia[variable_index] >= 1.0) return -1;

  // Return index of path vertex variable
  return lasers.NEntries()*laser_nv + 
         cameras.NEntries()*camera_nv + 
         vertex_index*vertex_nv + vertex_v[variable_index];
}



int GSVPoseOptimization::
PixelFeatureVariableIndex(int pixel_feature_index, int variable_index)
{
  // Check if pixel variable is active
  if (pixel_feature_nv == 0) return -1;
  if (pixel_feature_v[variable_index] == -1) return -1;

  // Return index of pixel variable
  return lasers.NEntries()*laser_nv + 
         cameras.NEntries()*camera_nv + 
         vertices.NEntries()*vertex_nv + 
         pixel_feature_index*pixel_feature_nv + pixel_feature_v[variable_index];
}



int GSVPoseOptimization::
ClusterVariableIndex(int cluster_index, int variable_index)
{
  // Check if cluster variable is active
  if (cluster_nv == 0) return -1;
  if (cluster_v[variable_index] == -1) return -1;

  // Check if variable is fixed
  GSVFeatureCluster *cluster = Cluster(cluster_index);
  if (cluster->inertia[variable_index] >= 1.0) return -1;

  // Return index of point cluster variable
  return lasers.NEntries()*laser_nv + 
         cameras.NEntries()*camera_nv + 
         vertices.NEntries()*vertex_nv + 
         pixel_features.NEntries()*pixel_feature_nv +
         cluster_index*cluster_nv + cluster_v[variable_index];
}



////////////////////////////////////////



RNScalar GSVPoseOptimization::
LaserVariableValue(int laser_index, int variable_index, RNScalar *x)
{
  // Return value from x
  int index = -1;
  if (x && (laser_v[variable_index] != -1)) 
    index = LaserVariableIndex(laser_index, variable_index);
  if (index != -1) return x[index];

  // Return stored value
  GSVLaser *laser = lasers.Kth(laser_index);
  LaserData *laser_data = (LaserData *) laser->Data();
  if (variable_index <= TZ) return laser_data->translation[variable_index];
  else if (variable_index <= RZ) return laser_data->rotation[variable_index-3];
  RNAbort("Invalid variable index");
  return 0.0;
}



RNScalar GSVPoseOptimization::
CameraVariableValue(int camera_index, int variable_index, RNScalar *x)
{
  // Return value from x
  int index = -1;
  if (x && (camera_v[variable_index] != -1)) 
    index = CameraVariableIndex(camera_index, variable_index);
  if (index != -1) return x[index];

  // Return stored value
  GSVCamera *camera = cameras.Kth(camera_index);
  CameraData *camera_data = (CameraData *) camera->Data();
  if (variable_index <= TZ) return camera_data->translation[variable_index];
  else if (variable_index <= RZ) return camera_data->rotation[variable_index-3];
  RNAbort("Invalid variable index");
  return 0.0;
}




RNScalar GSVPoseOptimization::
VertexVariableValue(int vertex_index, int variable_index, RNScalar *x)
{
  // Return value from x
  int index = -1;
  if (x && (vertex_v[variable_index] != -1)) 
    index = VertexVariableIndex(vertex_index, variable_index);
  if (index != -1) return x[index];

  // Return stored value
  GSVPathVertex *vertex = vertices.Kth(vertex_index);
  if (variable_index <= TZ) return vertex->translation[variable_index];
  else if (variable_index <= RZ) return vertex->rotation[variable_index-3];
  RNAbort("Invalid variable index");
  return 0.0;
}



RNScalar GSVPoseOptimization::
PixelFeatureVariableValue(int pixel_feature_index, int variable_index, RNScalar *x)
{
  // Return value from x
  int index = -1;
  if (x && (pixel_feature_v[variable_index] != -1)) 
    index = PixelFeatureVariableIndex(pixel_feature_index, variable_index);
  if (index != -1) return x[index];

  // Return stored value
  GSVFeature *feature = pixel_features.Kth(pixel_feature_index);
  if (variable_index == T) return feature->image_t;
  RNAbort("Invalid variable index");
  return 0.0;
}



RNScalar GSVPoseOptimization::
ClusterVariableValue(int cluster_index, int variable_index, RNScalar *x)
{
  // Return value from x
  int index = -1;
  if (x && (cluster_v[variable_index] != -1)) 
    index = ClusterVariableIndex(cluster_index, variable_index);
  if (index != -1) return x[index];

  // Return stored value
  GSVFeatureCluster *cluster = clusters.Kth(cluster_index);
  if (variable_index <= TZ) return cluster->translation[variable_index];
  else if (variable_index <= RZ) return cluster->rotation[variable_index-3];
  RNAbort("Invalid variable index");
  return 0.0;
}



////////////////////////////////////////



int GSVPoseOptimization::
SetLaserVariableValue(int laser_index, int variable_index, RNScalar *x)
{
  // Get value from x
  if (!x) return 0;
  if (laser_v[variable_index] == -1) return 0; 
  int index = LaserVariableIndex(laser_index, variable_index);
  if (index == -1) return 0;

  /// Set stored value
  GSVLaser *laser = lasers.Kth(laser_index);
  LaserData *laser_data = (LaserData *) laser->Data();
  if (variable_index <= TZ) laser_data->translation[variable_index] = x[index];
  else if (variable_index <= RZ) laser_data->rotation[variable_index-3] = x[index];
  else RNAbort("Invalid variable index");

  // Return success
  return 1;
}




int GSVPoseOptimization::
SetCameraVariableValue(int camera_index, int variable_index, RNScalar *x)
{
  // Get value from x
  if (!x) return 0;
  if (camera_v[variable_index] == -1) return 0; 
  int index = CameraVariableIndex(camera_index, variable_index);
  if (index == -1) return 0;

  /// Set stored value
  GSVCamera *camera = cameras.Kth(camera_index);
  CameraData *camera_data = (CameraData *) camera->Data();
  if (variable_index <= TZ) camera_data->translation[variable_index] = x[index];
  else if (variable_index <= RZ) camera_data->rotation[variable_index-3] = x[index];
  else RNAbort("Invalid variable index");

  // Return success
  return 1;
}




int GSVPoseOptimization::
SetVertexVariableValue(int vertex_index, int variable_index, RNScalar *x)
{
  // Get value from x
  if (!x) return 0;
  if (vertex_v[variable_index] == -1) return 0; 
  int index = VertexVariableIndex(vertex_index, variable_index);
  if (index == -1) return 0;

  /// Set stored value
  GSVPathVertex *vertex = vertices.Kth(vertex_index);
  if (variable_index <= TZ) vertex->translation[variable_index] = x[index];
  else if (variable_index <= RZ) vertex->rotation[variable_index-3] = x[index];
  else RNAbort("Invalid variable index");

  // Return success
  return 1;
}



int GSVPoseOptimization::
SetPixelFeatureVariableValue(int pixel_feature_index, int variable_index, RNScalar *x)
{
  // Get value from x
  if (!x) return 0;
  if (pixel_feature_v[variable_index] == -1) return 0; 
  int index = PixelFeatureVariableIndex(pixel_feature_index, variable_index);
  if (index == -1) return 0;

  // Set stored value
  GSVFeature *feature = pixel_features.Kth(pixel_feature_index);
  if (variable_index == T) feature->image_t = x[index];
  else RNAbort("Invalid variable index");

  // Return success
  return 1;
}



int GSVPoseOptimization::
SetClusterVariableValue(int cluster_index, int variable_index, RNScalar *x)
{
  // Get value from x
  if (!x) return 0;
  if (cluster_v[variable_index] == -1) return 0; 
  int index = ClusterVariableIndex(cluster_index, variable_index);
  if (index == -1) return 0;

  // Set stored value
  GSVFeatureCluster *cluster = clusters.Kth(cluster_index);
  if (variable_index <= TZ) cluster->translation[variable_index] = x[index];
  else if (variable_index <= RZ) cluster->rotation[variable_index-3] = x[index];
  else RNAbort("Invalid variable index");

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Debugging functions
////////////////////////////////////////////////////////////////////////

void GSVPoseOptimization::
PrintErrors(RNSystemOfEquations *equations, double *x, double rigidity, int solver, 
  int print_details, int print_values, int print_residuals, 
  int print_partial_derivatives, int print_equations)
{
  // Allocate memory
  int n = equations->NVariables();
  int m = equations->NEquations();
  double *y = new double [ m ];

  // Print parameters
  printf("Errors: ");
  printf("( %d %g | %d %d %d %d %d | %d %d %d %d | %d %d %d | %d %d %d )\n", 
    solver, rigidity, 
    laser_nv, camera_nv, vertex_nv, pixel_feature_nv, cluster_nv, 
    lasers.NEntries(), cameras.NEntries(), vertices.NEntries(), 
    pixel_features.NEntries(),
    NFeatures(), NClusters(), NCorrespondences(), 
    equations->NVariables(), equations->NEquations(), equations->NPartialDerivatives());

  // Print overall error
  RNScalar overall_error = 0;
  equations->EvaluateResiduals(x, y);
  int overall_m = equations->NEquations();
  for (int i = 0; i < overall_m; i++) overall_error += y[i] * y[i];
  RNScalar overall_rmsd = (overall_m > 0) ? sqrt(overall_error / overall_m) : 0;
  printf("  %20s = %18.6f %12.6f %6d\n", "Overall Error", overall_error, overall_rmsd, overall_m);
  printf("----------------------------------------------------------\n");
    
  if (print_details) {
    // Print inertia error
    RNScalar inertia_error = 0;
    RNSystemOfEquations inertia_equations(n);
    AddInertiaEquations(&inertia_equations);
    int inertia_m = inertia_equations.NEquations();
    if (inertia_m > 0) inertia_equations.EvaluateResiduals(x, y);
    for (int i = 0; i < inertia_m; i++) inertia_error += y[i] * y[i];
    RNScalar inertia_rmsd = (inertia_m > 0) ? sqrt(inertia_error / inertia_m) : 0;
    printf("  %20s = %18.6f %12.6f %6d\n", "Inertia Error", inertia_error, inertia_rmsd, inertia_m);

    // Print rigidity error
    RNScalar rigidity_error = 0;
    RNSystemOfEquations rigidity_equations(n);
    AddRigidityEquations(&rigidity_equations, rigidity);
    int rigidity_m = rigidity_equations.NEquations();
    if (rigidity_m > 0) rigidity_equations.EvaluateResiduals(x, y);
    for (int i = 0; i < rigidity_m; i++) rigidity_error += y[i] * y[i];
    RNScalar rigidity_rmsd = (rigidity_m > 0) ? sqrt(rigidity_error / rigidity_m) : 0;
    printf("  %20s = %18.6f %12.6f %6d\n", "Rigidity Error", rigidity_error, rigidity_rmsd, rigidity_m);

    // Print correspondence error
    RNScalar correspondence_error = 0;
    RNSystemOfEquations correspondence_equations(n);
    AddCorrespondenceEquations(&correspondence_equations);
    int correspondence_m = correspondence_equations.NEquations();
    if (correspondence_m > 0) correspondence_equations.EvaluateResiduals(x, y);
    for (int i = 0; i < correspondence_m; i++) correspondence_error += y[i] * y[i];
    RNScalar correspondence_rmsd = (correspondence_m > 0) ? sqrt(correspondence_error / correspondence_m) : 0;
    printf("  %20s = %18.6f %12.6f %6d\n", "Correspondence Error", correspondence_error, correspondence_rmsd, correspondence_m);

    // Print cluster error
    RNScalar cluster_error = 0;
    RNSystemOfEquations cluster_equations(n);
    AddClusterEquations(&cluster_equations);
    int cluster_m = cluster_equations.NEquations();
    if (cluster_m > 0) cluster_equations.EvaluateResiduals(x, y);
    for (int i = 0; i < cluster_m; i++) cluster_error += y[i] * y[i];
    RNScalar cluster_rmsd = (cluster_m > 0) ? sqrt(cluster_error / cluster_m) : 0;
    printf("  %20s = %18.6f %12.6f %6d\n", "Cluster Error", cluster_error, cluster_rmsd, cluster_m);
  }

  // Delete memory
  delete [] y;

  // Print values
  if (print_values) {
    printf("Values: %d\n", equations->NVariables());
    for (int i = 0; i < equations->NVariables(); i++) {
      printf("  %12d %12.6f ", i, x[i]);
      if ((i > 0) && ((i % 6) == 0)) printf("\n");
    }
    printf("\n");
  }

  // Print residuals
  if (print_residuals) {
    double sum = 0;
    int m = equations->NEquations();
    double *y = new double [ m ];
    equations->EvaluateResiduals(x, y);
    printf("Residuals: %d\n", m);
    for (int i = 0; i < m; i++) {
      printf("  %12d %12.6f\n", i, y[i]);
      sum += y[i]*y[i];
    }
    printf("%g\n", sum);
    delete [] y;
  }

  // Print partial derivatives
  if (print_partial_derivatives) {
    printf("Partial derivatives: %d\n", equations->NEquations() * equations->NVariables());
    for (int i = 0; i < equations->NEquations(); i++) {
      RNEquation *equation = equations->Equation(i);
      for (int j = 0; j < equations->NVariables(); j++) {
        RNScalar d = equation->PartialDerivative(x, j);
        if (d == 0.0) continue;
        printf("  %12d %12d %12.6f\n", i, j, d);
      }
    }
  }

  // Print equations
  if (print_equations) {
    printf("Equations:\n");
    for (int i = 0; i < equations->NEquations(); i++) {
      RNEquation *equation = equations->Equation(i);
      equation->Print();
    }
  }
}



