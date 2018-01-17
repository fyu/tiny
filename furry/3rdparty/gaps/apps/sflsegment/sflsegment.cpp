// Source file for the surfel segmentation program



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "R3Surfels/R3Surfels.h"



////////////////////////////////////////////////////////////////////////
// Program arguments
////////////////////////////////////////////////////////////////////////

static char *scene_name = NULL;
static char *database_name = NULL;
static char *parent_object_name = NULL;
static char *parent_node_name = NULL;
static char *source_node_name = NULL;
static char *mask_grid_name = NULL;
static char *support_grid_name = NULL;
static int cluster_hierarchically = FALSE;
static int cluster_mean_shift = FALSE;
static int min_cluster_points = 100;
static double min_cluster_spacing = 0.5;
static double min_cluster_diameter = 1;
static double max_cluster_diameter = 0;
static double max_cluster_offplane_distance = 2;
static double max_cluster_normal_angle = 0;
static int max_neighbors = 16;
static double max_neighbor_distance = 0.5;
static double max_neighbor_normal_angle = RN_PI / 4.0;
static double max_pair_centroid_distance = 20;
static double max_pair_offplane_distance = 2;
static double max_pair_normal_angle = RN_PI / 2.0;
static double min_pair_affinity = 1.0E-6;
static double min_pair_affinity_ratio = 0;
static int print_verbose = 0;
static int print_debug = 0;




////////////////////////////////////////////////////////////////////////
// Surfel scene I/O Functions
////////////////////////////////////////////////////////////////////////

static R3SurfelScene *
OpenScene(const char *scene_name, const char *database_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate scene
  R3SurfelScene *scene = new R3SurfelScene();
  if (!scene) {
    fprintf(stderr, "Unable to allocate scene\n");
    return NULL;
  }

  // Open scene files
  if (!scene->OpenFile(scene_name, database_name, "r+", "r+")) {
    delete scene;
    return NULL;
  }

  // Print statistics
  if (print_verbose) {
    printf("Opened scene ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Objects = %d\n", scene->NObjects());
    printf("  # Labels = %d\n", scene->NLabels());
    printf("  # Assignments = %d\n", scene->NLabelAssignments());
    printf("  # Features = %d\n", scene->NFeatures());
    printf("  # Nodes = %d\n", scene->Tree()->NNodes());
    printf("  # Blocks = %d\n", scene->Tree()->Database()->NBlocks());
    printf("  # Surfels = %d\n", scene->Tree()->Database()->NSurfels());
    fflush(stdout);
  }

  // Return scene
  return scene;
}



static int
CloseScene(R3SurfelScene *scene)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Print statistics
  if (print_verbose) {
    printf("Closing scene ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Objects = %d\n", scene->NObjects());
    printf("  # Labels = %d\n", scene->NLabels());
    printf("  # Assignments = %d\n", scene->NLabelAssignments());
    printf("  # Features = %d\n", scene->NFeatures());
    printf("  # Nodes = %d\n", scene->Tree()->NNodes());
    printf("  # Blocks = %d\n", scene->Tree()->Database()->NBlocks());
    printf("  # Surfels = %d\n", scene->Tree()->Database()->NSurfels());
    fflush(stdout);
  }

  // Close scene files
  if (!scene->CloseFile()) {
    delete scene;
    return 0;
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Mean shift segmentation functions
////////////////////////////////////////////////////////////////////////

struct MeanShiftCluster {
  R2Point centroid;
};

static R2Point
ComputeMeanShiftPosition(R3SurfelPointSet *pointset, const R2Point& start_position, RNLength radius)
{
  // Initialize result
  R3Point current_position(start_position.X(), start_position.Y(), 0);

  // Iteratively move towards centroid of neighborhood
  for (int iter = 0; iter < 10; iter++) {
    // Estimate ground plane
    R3SurfelCylinderConstraint cylinder_constraint(current_position, 2);
    R3SurfelPointSet *pointset1 = CreatePointSet(pointset, &cylinder_constraint);
    if (!pointset1) break;
    if (pointset1->NPoints() == 0) { delete pointset1; break; }
    R3Plane bottom_plane = FitSupportPlane(pointset1);
    R3Plane top_plane(bottom_plane[0], bottom_plane[1], bottom_plane[2], bottom_plane[3] - 2.5);
    delete pointset1;

    // Extract neighborhood pointset
    R3SurfelCylinderConstraint radius_constraint(current_position, radius);
    R3SurfelPlaneConstraint bottom_plane_constraint(bottom_plane, FALSE, FALSE, TRUE, 0.25);
    R3SurfelPlaneConstraint top_plane_constraint(top_plane, TRUE, FALSE, FALSE, 0.25);
    R3SurfelMultiConstraint multiconstraint;
    multiconstraint.InsertConstraint(&radius_constraint);
    multiconstraint.InsertConstraint(&bottom_plane_constraint);
    multiconstraint.InsertConstraint(&top_plane_constraint);
    R3SurfelPointSet *pointset2 = CreatePointSet(pointset, &multiconstraint);
    if (!pointset2) break;

    // Sum positions of points in neighborhood
    RNScalar total_weight = 0;
    R3Point total_position = R3zero_point;
    for (int i = 0; i < pointset2->NPoints(); i++) {
      const R3SurfelPoint *point = pointset2->Point(i);
      R3Point position = point->Position();
      total_position += position;
      total_weight += 1;
    }

    // Delete pointset
    delete pointset2;

    // Find centroid
    if (total_weight == 0) break;
    R3Point centroid = total_position / total_weight;
    if (R3SquaredDistance(centroid, current_position) < 0.0001) break;

    // Set current position and iterate
    current_position = centroid;
  }

  // Return mean shift centroid
  return R2Point(current_position.X(), current_position.Y());
}



static int
CreateMeanShiftClusters(R3SurfelPointSet *pointset, RNScalar spacing,
  RNArray<MeanShiftCluster *>& clusters)
{
  // Compute sampling variables
  R3Box bbox = pointset->BBox();
  int nxsamples = (int) (bbox.XLength() / spacing) + 1;
  int nysamples = (int) (bbox.YLength() / spacing) + 1;
  RNScalar xstep = bbox.XLength() / nxsamples;
  RNScalar ystep = bbox.YLength() / nysamples;

  // Create mask 
  R2Box mask_bbox(bbox[0][0], bbox[0][1], bbox[1][0], bbox[1][1]);
  R2Grid mask(2*nxsamples, 2*nysamples, mask_bbox);

  // Sample positions and mean shift each of them
  for (int ix = 0; ix < nxsamples; ix++) {
    double x = bbox.XMin() + (ix+0.5)*xstep;
    for (int iy = 0; iy < nysamples; iy++) {
      double y = bbox.YMin() + (iy+0.5)*ystep;

      // Mean shift
      R2Point seed(x, y);
      seed[0] += (0.5 - RNRandomScalar()) * 0.5 * xstep;
      seed[1] += (0.5 - RNRandomScalar()) * 0.5 * ystep;
      R2Point centroid = ComputeMeanShiftPosition(pointset, seed, spacing); 

      // Check/update mask
      RNScalar mask_value = mask.WorldValue(centroid);
      if (mask_value != 0) continue;
      mask.RasterizeWorldCircle(centroid, spacing, 1);

      // Create cluster
      MeanShiftCluster *cluster = new MeanShiftCluster();
      cluster->centroid = centroid;
      clusters.Insert(cluster);

      printf("%g %g\n", centroid.X(), centroid.Y());
    }
  }

  // Return success
  return 1;
}



static int
ClusterMeanShift(R3SurfelScene *scene, 
  const char *parent_object_name, 
  const char *parent_node_name, 
  const char *source_node_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_verbose) {
    printf("Creating mean shift clusters ...\n");
    fflush(stdout);
  }

  // Get convenient variables
  R3SurfelTree *tree = scene->Tree();
  if (!tree) return 0;
  R3SurfelDatabase *database = tree->Database();
  if (!database) return 0;

  // Find parent object
  R3SurfelObject *parent_object = scene->RootObject();
  if (parent_object_name) {
    parent_object = scene->FindObjectByName(parent_object_name);
    if (!parent_object) {
      fprintf(stderr, "Unable to find parent object with name %s\n", parent_object_name);
      return 0;
    }
  }

  // Find parent node
  R3SurfelNode *parent_node = tree->RootNode();
  if (parent_node_name) {
    parent_node = tree->FindNodeByName(parent_node_name);
    if (!parent_node) {
      fprintf(stderr, "Unable to find parent node with name %s\n", parent_node_name);
      return 0;
    }
  }

  // Find source node
  R3SurfelNode *source_node = tree->RootNode();
  if (source_node_name) {
    source_node = tree->FindNodeByName(source_node_name);
    if (!source_node) {
      fprintf(stderr, "Unable to find source node with name %s\n", source_node_name);
      return 0;
    }
  }

  // Initialize array of clusters
  RNArray<MeanShiftCluster *> clusters;

  // Split nodes
  tree->SplitNodes(source_node);

  // Create a set of blocks in leaf nodes to split
  RNArray<R3SurfelBlock *> blocks;
  RNArray<R3SurfelNode *> stack;
  stack.InsertTail(source_node);
  while (!stack.IsEmpty()) {
    R3SurfelNode *node = stack.Tail();
    stack.RemoveTail();
    if (node->NParts() == 0) {
      for (int i = 0; i < node->NBlocks(); i++) {
        blocks.Insert(node->Block(i));
      }
    }
    else {
      for (int i = 0; i < node->NParts(); i++) {
        stack.InsertTail(node->Part(i));
      }
    }
  }

  // Segment blocks
  for (int i = 0; i < blocks.NEntries(); i++) {
    R3SurfelBlock *block = blocks.Kth(i);
    R3SurfelNode *node = block->Node();
    if (!node) continue;
    if (node->NParts() > 0) continue;
    R3SurfelObject *object = node->Object();
    if (object) continue;
    
    // Print debug statement
    if (print_debug) {
      printf("  %6d/%6d : ", i, blocks.NEntries());
      fflush(stdout);
    }

    // Create pointset from block
    R3SurfelPointSet *pointset = new R3SurfelPointSet(block);

    // Print debug statement
    if (print_debug) {
      printf("%9d ", pointset->NPoints());
      fflush(stdout);
    }

    // Create clusters
    if (!CreateMeanShiftClusters(pointset, min_cluster_spacing, clusters)) {
      fprintf(stderr, "Error in clustering\n");
      delete pointset; 
      return 0; 
    }

    // Print debug statement
    if (print_debug) {
      printf("%6d\n", clusters.NEntries());
      fflush(stdout);
    }
  }

  // Print message
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Hierarchical clustering functions
////////////////////////////////////////////////////////////////////////

struct Cluster {
  Cluster *parent;
  RNArray<struct ClusterPair *> pairs;
  RNArray<R3SurfelPoint *> points;
  R3Point centroid;
  R3Plane plane;
  RNScalar variances[3];
  RNScalar affinity;
};

struct ClusterPair {
  Cluster *clusters[2];
  RNScalar affinity;
  ClusterPair **heapentry;
};


static R3SurfelObject *
CreateObject(R3SurfelScene *scene, 
  const RNArray<R3SurfelPoint *>& points, 
  R3SurfelNode *parent_node, R3SurfelObject *parent_object,
  const char *node_name, const char *object_name)
{
  // Get convenient variables
  R3SurfelTree *tree = scene->Tree();
  if (!tree) return NULL;
  R3SurfelDatabase *database = tree->Database();
  if (!database) return NULL;
  if (points.IsEmpty()) return NULL;

  // Compute bounding box
  R3Box bbox = R3null_box;
  for (int j = 0; j < points.NEntries(); j++) {
    R3SurfelPoint *point = points.Kth(j);
    bbox.Union(point->Position());
  }
      
  // Create array of surfels
  R3Point origin = bbox.Centroid();
  R3Surfel *surfels = new R3Surfel [ points.NEntries() ];
  for (int j = 0; j < points.NEntries(); j++) {
    R3SurfelPoint *point = points.Kth(j);
    R3SurfelBlock *block = point->Block();
    const R3Surfel *surfel = point->Surfel();
    RNCoord x = surfel->X() + block->Origin().X() - origin.X();
    RNCoord y = surfel->Y() + block->Origin().Y() - origin.Y();
    RNCoord z = surfel->Z() + block->Origin().Z() - origin.Z();
    surfels[j].SetCoords(x, y, z);
    surfels[j].SetColor(surfel->Color());
    surfels[j].SetAerial(surfel->IsAerial());
  }

  // Create block
  R3SurfelBlock *block = new R3SurfelBlock(surfels, points.NEntries(), origin);
  if (!block) {
    fprintf(stderr, "Unable to create node\n");
    delete [] surfels;
    return NULL;
  }
  
  // Delete surfels
  delete [] surfels;
      
  // Insert block
  block->UpdateProperties();
  database->InsertBlock(block);

  // Create node
  R3SurfelNode *node = NULL;
  if (block && parent_node) {
    node = new R3SurfelNode(node_name);
    if (!node) {
      fprintf(stderr, "Unable to create node\n");
      delete block;
    return NULL;
    }
    
    // Insert node
    node->InsertBlock(block);
    node->UpdateProperties();
    tree->InsertNode(node, parent_node);
  }

  // Create object
  R3SurfelObject *object = NULL;
  if (node && parent_object) {
    object = new R3SurfelObject(object_name);
    if (!object) {
      fprintf(stderr, "Unable to create object\n");
      delete block;
      delete node;
      return NULL;
    }

    // Insert object
    scene->InsertObject(object, parent_object);
    object->InsertNode(node);
    object->UpdateProperties();
  }

  // Release block
  database->ReleaseBlock(block);
    
  // Return object
  return object;
}



#if 0

static int
ComputeClusterShapeParameters(Cluster *cluster, 
  R3Point *result_centroid = NULL, R3Plane *result_plane = NULL,
  R3Triad *result_triad = NULL, RNScalar *result_variances = NULL)
{
  // Make array of sample positions
  int nsamples = 0;
  const int max_samples = 1024;
  R3Point samples[max_samples];
  RNScalar npoints = cluster->points.NEntries();
  int step = (int) (npoints / max_samples) + 1;
  int start = (int) (0.5 * step);
  for (int i = start; i < cluster->points.NEntries(); i += step) {
    assert(nsamples < max_samples);
    R3SurfelPoint *point = cluster->points.Kth(i);
    samples[nsamples] = point->Position();
    nsamples++;
  }

  // Compute centroid
  R3Point centroid = R3Centroid(nsamples, samples);
  if (result_centroid) *result_centroid = centroid;

  // Compute plane, triad, variances
  if (result_plane || result_triad || result_variances) {
    R3Triad triad = R3PrincipleAxes(centroid, nsamples, samples, NULL, result_variances);
    if (result_plane) *result_plane = R3Plane(centroid, triad[2]);
    if (result_triad) *result_triad = triad;
  }

  // Return success
  return 1;
}

#endif



static int
ComputeClusterPairShapeParameters(Cluster *cluster0, Cluster *cluster1,
  R3Point *result_centroid = NULL, R3Plane *result_plane = NULL,
  R3Triad *result_triad = NULL, RNScalar *result_variances = NULL)
{
  // Compute combined number of points
  int npoints0 = cluster0->points.NEntries();
  int npoints1 = cluster1->points.NEntries();
  RNScalar npoints = npoints0 + npoints1;
  RNScalar t0 = npoints0 / npoints;
  RNScalar t1 = npoints1 / npoints;

  // Compute combined centroid
  R3Point centroid = t0 * cluster0->centroid + t1 * cluster1->centroid;
  if (result_centroid) *result_centroid = centroid;

  // Compute combined plane, triad, and variances
  const int min_samples = 0;
  const int max_samples = 0;
  if ((npoints < min_samples) || (npoints > max_samples)) {
    // Too few/many points - average
    if (result_plane || result_triad) {
      R3Vector normal0 = cluster0->plane.Normal();
      R3Vector normal1 = cluster1->plane.Normal();
      if (normal0.Dot(normal1) < 0) normal1.Flip();
      R3Vector normal = t0 * normal0 + t1 * normal1;
      normal.Normalize();
      if (result_plane) {
        *result_plane = R3Plane(centroid, normal);
      }
      if (result_triad) {
        R3Vector yaxis = R3xyz_triad.Axis(normal.MinDimension());
        yaxis.Cross(normal);
        R3Vector xaxis = yaxis % normal;
        xaxis.Normalize();
        *result_triad = R3Triad(xaxis, yaxis, normal);
      }
    }
    if (result_variances) {
      result_variances[0] = 0.25 * R3SquaredDistance(cluster0->centroid, cluster1->centroid);
      result_variances[1] = 0;
      result_variances[2] = 0;
    }
  }
  else {
    // Make combined array of sample positions
    if (result_plane || result_triad || result_variances) {
      int nsamples = 0;
      R3Point samples[max_samples];
      for (int i = 0; i < cluster0->points.NEntries(); i++) {
        assert(nsamples < max_samples);
        R3SurfelPoint *point = cluster0->points.Kth(i);
        samples[nsamples] = point->Position();
        nsamples++;
      }
      for (int i = 0; i < cluster1->points.NEntries(); i++) {
        assert(nsamples < max_samples);
        R3SurfelPoint *point = cluster1->points.Kth(i);
        samples[nsamples] = point->Position();
        nsamples++;
      }
      
      // Compute shape parameters from sampled points
      R3Triad triad = R3PrincipleAxes(centroid, nsamples, samples, NULL, result_variances);
      if (result_plane) *result_plane = R3Plane(centroid, triad[2]);
      if (result_triad) *result_triad = triad;
    }
  }
  
  // Return success
  return 1;
}  




static RNScalar
ClusterAffinity(Cluster *cluster0, Cluster *cluster1)
{
  // Compute centroid distance
  RNScalar cluster_distance_factor = 1;
  if (max_pair_centroid_distance > 0) {
    RNLength cluster_distance = R3Distance(cluster0->centroid, cluster1->centroid);
    if (cluster_distance > max_pair_centroid_distance) return -FLT_MAX;
    cluster_distance_factor = 1.0 - cluster_distance / max_pair_centroid_distance;
  }

  // Compute offplane distances
  RNScalar offplane0_distance_factor = 1;
  RNScalar offplane1_distance_factor = 1;
  if (max_pair_offplane_distance > 0) {
    // Compute point0-plane1 distance
    RNLength offplane0_distance = R3Distance(cluster0->plane, cluster1->centroid);
    if (offplane0_distance > max_pair_offplane_distance) return -FLT_MAX;
    offplane0_distance_factor = 1.0 - offplane0_distance / max_pair_offplane_distance;

    // Compute point1-plane0 distance
    RNLength offplane1_distance = R3Distance(cluster1->plane, cluster0->centroid);
    if (offplane1_distance > max_pair_offplane_distance) return -FLT_MAX;
    offplane1_distance_factor = 1.0 - offplane1_distance / max_pair_offplane_distance;
  }

  // Compute normal angle
  RNScalar normal_angle_factor = 1;
  if (max_pair_normal_angle > 0) {
    RNScalar dot = fabs(cluster0->plane.Normal().Dot(cluster1->plane.Normal()));
    RNAngle normal_angle = (dot < 1) ? acos(dot) : 0;
    if (normal_angle > max_pair_normal_angle) return -FLT_MAX;
    normal_angle_factor = 1.0 - normal_angle / max_pair_normal_angle;
  }

  // Compute affinity
  RNScalar affinity = 1;
  affinity *= cluster_distance_factor;
  affinity *= offplane0_distance_factor;
  affinity *= offplane1_distance_factor;
  affinity *= normal_angle_factor;

  // Return affinity
  return affinity;
}



static int
ClusterHierarchically(R3SurfelScene *scene, R3SurfelPointGraph *graph, 
  R3SurfelObject *parent_object, R3SurfelNode *parent_node)
{
  // Get convenient variables
  R3SurfelTree *tree = scene->Tree();
  if (!tree) return 0;
  R3SurfelDatabase *database = tree->Database();
  if (!database) return 0;
  if (graph->NPoints() < min_cluster_points) return 1;
  if (!parent_object) parent_object = scene->RootObject();
  if (!parent_node) parent_node = tree->RootNode();

  ////////////////////////////////////////////////////////////////////////

  printf("HEREA\n");

  // Create clusters
  Cluster *clusters = new Cluster [ graph->NPoints() ];
  for (int i = 0; i < graph->NPoints(); i++) {
    R3SurfelPoint *point = graph->Point(i);
    R3Point centroid = point->Position();
    R3Vector normal = graph->PointNormal(i);
    R3Plane p(centroid, normal);
    Cluster *cluster = &clusters[i];
    cluster->parent = NULL;
    cluster->points.Insert(point);
    cluster->centroid = centroid;
    cluster->plane.Reset(centroid, normal);
    cluster->variances[0] = RN_EPSILON;
    cluster->variances[1] = RN_EPSILON;
    cluster->variances[2] = RN_EPSILON;
    cluster->affinity = RN_EPSILON;
  }

  ////////////////////////////////////////////////////////////////////////
  
  printf("HEREB %d\n", graph->NPoints());

  // Create cluster pairs 
  ClusterPair tmp;
  RNHeap<ClusterPair *> heap(&tmp, &tmp.affinity, &tmp.heapentry, FALSE);;
  for (int index0 = 0; index0 < graph->NPoints(); index0++) {
    Cluster *cluster0 = &clusters[index0];
    for (int j = 0; j < graph->NNeighbors(index0); j++) {
      R3SurfelPoint *point1 = graph->Neighbor(index0, j);
      int index1 = graph->PointIndex(point1);
      if (index1 <= index0) continue;
      Cluster *cluster1 = &clusters[index1];

      // Compute affinity
      RNScalar affinity = ClusterAffinity(cluster0, cluster1);
      if (affinity < min_pair_affinity) continue;

      // Create pair
      ClusterPair *pair = new ClusterPair();
      pair->clusters[0] = cluster0;
      pair->clusters[1] = cluster1;
      pair->affinity = affinity;
      pair->heapentry = NULL;
      heap.Push(pair);

      // Add pair to clusters
      cluster0->pairs.Insert(pair);
      cluster1->pairs.Insert(pair);
    }
  }

  ////////////////////////////////////////////////////////////////////////

  printf("HEREC %d\n", heap.NEntries());

  // Merge clusters hierarchically
  int merge_count = 0;
  while (!heap.IsEmpty()) {
    // Get pair
    ClusterPair *pair01 = heap.Pop();

    // Check if we are done
    if (pair01->affinity < min_pair_affinity) break;

    // Get clusters
    Cluster *cluster0 = pair01->clusters[0];
    Cluster *cluster1 = pair01->clusters[1];

    // if (print_debug) {
    //   static unsigned long count = 0;
    //   if ((count++ % 1000) == 0) {
    //     printf("  %15.12f : %9d %9d : %15d\n", pair01->affinity, 
    //            cluster0->points.NEntries(), cluster1->points.NEntries(), 
    //            heap.NEntries());
    //   }
    // }
    
    // Check if clusters have already been merged
    if (cluster0->parent || cluster1->parent) {
      // Find clusters
      Cluster *ancestor0 = cluster0;
      Cluster *ancestor1 = cluster1;
      while (ancestor0->parent) ancestor0 = ancestor0->parent;
      while (ancestor1->parent) ancestor1 = ancestor1->parent;
      if (ancestor0 != ancestor1) {
        // Find pair
        ClusterPair *pair = NULL;
        for (int j = 0; j < ancestor0->pairs.NEntries(); j++) {
          ClusterPair *tmp = ancestor0->pairs.Kth(j);
          if (tmp->clusters[0] == ancestor1) { pair = tmp; break; }
          if (tmp->clusters[1] == ancestor1) { pair = tmp; break; }
        }      
        
        // Create pair
        if (!pair) {
          RNScalar affinity = ClusterAffinity(ancestor0, ancestor1);
          if (affinity >= min_pair_affinity) {
            RNScalar denom = (ancestor0->affinity > ancestor1->affinity) ? ancestor0->affinity : ancestor1->affinity;
            RNScalar affinity_ratio = (denom > 0) ? affinity / denom : RN_INFINITY;
            if (affinity_ratio >= min_pair_affinity_ratio) {
              pair = new ClusterPair();
              pair->clusters[0] = ancestor0;
              pair->clusters[1] = ancestor1;
              pair->affinity = affinity;
              pair->heapentry = NULL;
              ancestor0->pairs.Insert(pair);
              ancestor1->pairs.Insert(pair);
              heap.Push(pair);
            }
          }
        }
      }
    }
    else {
      // Merge cluster1 into cluster0
      ComputeClusterPairShapeParameters(cluster0, cluster1, 
        &cluster0->centroid, &cluster0->plane, NULL, cluster0->variances);
      cluster0->points.Append(cluster1->points);
      cluster0->affinity = pair01->affinity;
      cluster1->points.Empty(TRUE);
      cluster1->parent = cluster0;
      cluster1->affinity = 0;

      // Update statistics
      merge_count++;
    }

    // Delete pair
    cluster0->pairs.Remove(pair01);
    cluster1->pairs.Remove(pair01);
    delete pair01;
  }

  ////////////////////////////////////////////////////////////////////////

  printf("HERED %d\n", heap.NEntries());

  // Create objects 
  int object_count = 0;
  RNArray<R3SurfelPoint *> unclustered_points;
  for (int i = 0; i < graph->NPoints(); i++) {
    Cluster *cluster = &clusters[i];
    if (cluster->parent) continue;
    if (cluster->points.NEntries() < min_cluster_points) {
      // Add points to list
      for (int j = 0; j < cluster->points.NEntries(); j++) {
        R3SurfelPoint *point = cluster->points.Kth(j);
        unclustered_points.Insert(point);
      }
    }
    else {
      // Insert points into object
      char node_name[1024];
      sprintf(node_name, "Cluster%d", object_count);
      if (CreateObject(scene, cluster->points, parent_node, parent_object, node_name, node_name)) {
        object_count++;
      }
    }
  }

  // Insert unclustered points
  if (unclustered_points.NEntries() > 0) {
    char node_name[1024];
    sprintf(node_name, "Unclustered");
    CreateObject(scene, unclustered_points, parent_node, NULL, node_name, NULL);
  }
  
  printf("HEREE %d\n", object_count);

  // Delete clusters
  delete [] clusters;

  // Return success
  return 1;
}



static int
ClusterHierarchically(R3SurfelScene *scene, 
  const char *parent_object_name, 
  const char *parent_node_name, 
  const char *source_node_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_verbose) {
    printf("Creating hierarchical clustering segments ...\n");
    fflush(stdout);
  }

  // Get convenient variables
  R3SurfelTree *tree = scene->Tree();
  if (!tree) return 0;
  R3SurfelDatabase *database = tree->Database();
  if (!database) return 0;

  // Find parent object
  R3SurfelObject *parent_object = scene->RootObject();
  if (parent_object_name) {
    parent_object = scene->FindObjectByName(parent_object_name);
    if (!parent_object) {
      fprintf(stderr, "Unable to find parent object with name %s\n", parent_object_name);
      return 0;
    }
  }

  // Find parent node
  R3SurfelNode *parent_node = tree->RootNode();
  if (parent_node_name) {
    parent_node = tree->FindNodeByName(parent_node_name);
    if (!parent_node) {
      fprintf(stderr, "Unable to find parent node with name %s\n", parent_node_name);
      return 0;
    }
  }

  // Find source node
  R3SurfelNode *source_node = tree->RootNode();
  if (source_node_name) {
    source_node = tree->FindNodeByName(source_node_name);
    if (!source_node) {
      fprintf(stderr, "Unable to find source node with name %s\n", source_node_name);
      return 0;
    }
  }

  // Split nodes
  tree->SplitNodes(source_node);

  // Create a set of blocks in leaf nodes to split
  RNArray<R3SurfelBlock *> blocks;
  RNArray<R3SurfelNode *> stack;
  stack.InsertTail(source_node);
  while (!stack.IsEmpty()) {
    R3SurfelNode *node = stack.Tail();
    stack.RemoveTail();
    if (node->NParts() == 0) {
      for (int i = 0; i < node->NBlocks(); i++) {
        blocks.Insert(node->Block(i));
      }
    }
    else {
      for (int i = 0; i < node->NParts(); i++) {
        stack.InsertTail(node->Part(i));
      }
    }
  }

  // Segment blocks 
  for (int i = 0; i < blocks.NEntries(); i++) {
    R3SurfelBlock *block = blocks.Kth(i);
    R3SurfelNode *node = block->Node();
    if (!node) continue;
    if (node->NParts() > 0) continue;
    R3SurfelObject *object = node->Object();
    if (object) continue;
    
    // Print debug statement
    if (print_debug) {
      printf("  %6d/%6d : ", i, blocks.NEntries());
      fflush(stdout);
    }

    // Create graph from block
    R3SurfelPointSet *pointset = new R3SurfelPointSet(block);
    R3SurfelPointGraph *graph = new R3SurfelPointGraph(*pointset, max_neighbors, max_neighbor_distance);
    delete pointset;

    // Print debug statement
    if (print_debug) {
      printf("%s %9d ", node->Name(), graph->NPoints());
      fflush(stdout);
    }

    // Create objects
    if (!ClusterHierarchically(scene, graph, parent_object, node)) {
      fprintf(stderr, "Error in clustering\n");
      delete graph; 
      return 0; 
    }

    // Print debug statement
    if (print_debug) {
      printf("%6d\n", node->NParts());
      fflush(stdout);
    }

    // Delete graph
    delete graph;

    // Delete block
    node->RemoveBlock(block);
    database->RemoveBlock(block);
    delete block;
  }

  // Print message
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Support surface segmentation functions
////////////////////////////////////////////////////////////////////////

static int
SegmentSupportGrid(R3SurfelScene *scene, 
  const char *parent_object_name, 
  const char *source_node_name,
  const char *parent_node_name,
  const char *support_grid_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int nobjects = 0;
  if (print_verbose) {
    printf("Creating support grid object ...\n");
    fflush(stdout);
  }

  // Get tree
  R3SurfelTree *tree = scene->Tree();
  if (!tree) {
    fprintf(stderr, "Scene has no tree\n");
    return 0;
  }    

  // Find parent object
  R3SurfelObject *parent_object = scene->RootObject();
  if (parent_object_name) {
    parent_object = scene->FindObjectByName(parent_object_name);
    if (!parent_object) {
      fprintf(stderr, "Unable to find parent object with name %s\n", parent_object_name);
      return 0;
    }
  }

  // Find source node
  R3SurfelNode *source_node = tree->RootNode();
  if (source_node_name) {
    source_node = tree->FindNodeByName(source_node_name);
    if (!source_node) {
      fprintf(stderr, "Unable to find source node with name %s\n", source_node_name);
      return 0;
    }
  }

  // Read support grid
  R2Grid *support_grid = new R2Grid();
  if (!support_grid->ReadFile(support_grid_name)) return 0;

  // Create constraint
  R3SurfelMultiConstraint constraint;
  R3SurfelObjectConstraint non_object_constraint(NULL, TRUE);
  R3SurfelOverheadGridConstraint support_constraint(support_grid, R3_SURFEL_CONSTRAINT_EQUAL, 
   R3_SURFEL_CONSTRAINT_Z, R3_SURFEL_CONSTRAINT_VALUE, 0, 0, 0.25);
  constraint.InsertConstraint(&non_object_constraint);
  constraint.InsertConstraint(&support_constraint);

  // Split source node
  RNArray<R3SurfelNode *> support_nodes, not_support_nodes;
  tree->SplitLeafNodes(source_node, constraint, &support_nodes, &not_support_nodes);
  if (support_nodes.NEntries() > 0) {
    // Create object name from grid name
    char object_name[1024];
    const char *grid_namep = strstr(support_grid_name, "/");
    if (grid_namep) grid_namep++;
    else grid_namep = support_grid_name;
    strcpy(object_name, grid_namep);
    char *object_namep = strrchr(object_name, '.');
    if (object_namep) *object_namep = '\0';

    // Create object
    R3SurfelObject *object = new R3SurfelObject(object_name);

    // Set names of nodes and insert them into object
    for (int i = 0; i < support_nodes.NEntries(); i++) {
      R3SurfelNode *node = support_nodes.Kth(i);
      if (node->Object()) continue;
      char node_name[1024];
      sprintf(node_name, "Support_%d", i);
      node->SetName(node_name);
      object->InsertNode(node);
    }

    // Set names of other nodes
    for (int i = 0; i < not_support_nodes.NEntries(); i++) {
      R3SurfelNode *node = not_support_nodes.Kth(i);
      char node_name[1024];
      sprintf(node_name, "Not_Support_%d", i);
      node->SetName(node_name);
    }

    // Insert object into scene
    scene->InsertObject(object, parent_object);

    // Update properties
    object->UpdateProperties();

    // Update statistics
    nobjects++;
  }

  // Print message
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Objects = %d\n", nobjects);
    printf("  # Nodes = %d\n", support_nodes.NEntries());
    fflush(stdout);
  }

  // Delete grid
  if (support_grid) delete support_grid;

  // Return success
  return 1;
}



static int
SegmentMaskGrid(R3SurfelScene *scene, 
  const char *parent_object_name, 
  const char *source_node_name,
  const char *parent_node_name,
  const char *mask_grid_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int nobjects = 0;
  if (print_verbose) {
    printf("Creating mask grid object ...\n");
    fflush(stdout);
  }

  // Get tree
  R3SurfelTree *tree = scene->Tree();
  if (!tree) {
    fprintf(stderr, "Scene has no tree\n");
    return 0;
  }    

  // Find parent object
  R3SurfelObject *parent_object = scene->RootObject();
  if (parent_object_name) {
    parent_object = scene->FindObjectByName(parent_object_name);
    if (!parent_object) {
      fprintf(stderr, "Unable to find parent object with name %s\n", parent_object_name);
      return 0;
    }
  }

  // Find source node
  R3SurfelNode *source_node = tree->RootNode();
  if (source_node_name) {
    source_node = tree->FindNodeByName(source_node_name);
    if (!source_node) {
      fprintf(stderr, "Unable to find source node with name %s\n", source_node_name);
      return 0;
    }
  }

  // Read mask grid
  R2Grid *mask_grid = new R2Grid();
  if (!mask_grid->ReadFile(mask_grid_name)) return 0;

  // Create constraint
  R3SurfelMultiConstraint constraint;
  R3SurfelObjectConstraint non_object_constraint(NULL, TRUE);
  R3SurfelOverheadGridConstraint mask_constraint(mask_grid, R3_SURFEL_CONSTRAINT_GREATER, 
   R3_SURFEL_CONSTRAINT_OPERAND, R3_SURFEL_CONSTRAINT_VALUE, 0.5, 0, 0);
  constraint.InsertConstraint(&non_object_constraint);
  constraint.InsertConstraint(&mask_constraint);

  // Split source node
  RNArray<R3SurfelNode *> mask_nodes, not_mask_nodes;
  tree->SplitLeafNodes(source_node, constraint, &mask_nodes, &not_mask_nodes);
  if (mask_nodes.NEntries() > 0) {
    // Create object name from grid name
    char object_name[1024];
    const char *grid_namep = strstr(mask_grid_name, "/");
    if (grid_namep) grid_namep++;
    else grid_namep = mask_grid_name;
    strcpy(object_name, grid_namep);
    char *object_namep = strrchr(object_name, '.');
    if (object_namep) *object_namep = '\0';

    // Create object
    R3SurfelObject *object = new R3SurfelObject(object_name);

    // Insert nodes into object
    for (int i = 0; i < mask_nodes.NEntries(); i++) {
      R3SurfelNode *node = mask_nodes.Kth(i);
      if (node->Object()) continue;
      char node_name[1024];
      sprintf(node_name, "Mask_%d", i);
      node->SetName(node_name);
      object->InsertNode(node);
    }

    // Set names of other nodes
    for (int i = 0; i < not_mask_nodes.NEntries(); i++) {
      R3SurfelNode *node = not_mask_nodes.Kth(i);
      char node_name[1024];
      sprintf(node_name, "Not_Mask_%d", i);
      node->SetName(node_name);
    }

    // Insert object into scene
    scene->InsertObject(object, parent_object);

    // Update properties
    object->UpdateProperties();

    // Update statistics
    nobjects++;
  }

  // Print message
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Objects = %d\n", nobjects);
    printf("  # Nodes = %d\n", mask_nodes.NEntries());
    fflush(stdout);
  }

  // Delete grid
  if (mask_grid) delete mask_grid;

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Argument Parsing Functions
////////////////////////////////////////////////////////////////////////

static int 
ParseArgs(int argc, char **argv)
{
  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) print_verbose = 1;
      else if (!strcmp(*argv, "-debug")) print_debug = 1;
      else if (!strcmp(*argv, "-parent_object")) { argc--; argv++; parent_object_name = *argv; }
      else if (!strcmp(*argv, "-parent_node")) { argc--; argv++; parent_node_name = *argv; }
      else if (!strcmp(*argv, "-source_node")) { argc--; argv++; source_node_name = *argv; }
      else if (!strcmp(*argv, "-mask_grid")) { argc--; argv++; mask_grid_name = *argv; }
      else if (!strcmp(*argv, "-support_grid")) { argc--; argv++; support_grid_name = *argv; }
      else if (!strcmp(*argv, "-cluster_hierarchically")) cluster_hierarchically = 1;
      else if (!strcmp(*argv, "-cluster_mean_shift")) cluster_mean_shift = 1;
      else if (!strcmp(*argv, "-min_cluster_points")) { argc--; argv++; min_cluster_points = atoi(*argv); }
      else if (!strcmp(*argv, "-min_cluster_diameter")) { argc--; argv++; min_cluster_diameter = atof(*argv); }
      else if (!strcmp(*argv, "-max_cluster_diameter")) { argc--; argv++; max_cluster_diameter = atof(*argv); }
      else if (!strcmp(*argv, "-max_cluster_offplane_distance")) { argc--; argv++; max_cluster_offplane_distance = atof(*argv); }
      else if (!strcmp(*argv, "-max_cluster_normal_angle")) { argc--; argv++; max_cluster_normal_angle = atof(*argv); }
      else if (!strcmp(*argv, "-max_neighbors")) { argc--; argv++; max_neighbors = atoi(*argv); }
      else if (!strcmp(*argv, "-max_neighbor_distance")) { argc--; argv++; max_neighbor_distance = atof(*argv); }
      else if (!strcmp(*argv, "-max_neighbor_normal_angle")) { argc--; argv++; max_neighbor_normal_angle = atof(*argv); }
      else if (!strcmp(*argv, "-max_pair_centroid_distance")) { argc--; argv++; max_pair_centroid_distance = atof(*argv); }
      else if (!strcmp(*argv, "-max_pair_offplane_distance")) { argc--; argv++; max_pair_offplane_distance = atof(*argv); }
      else if (!strcmp(*argv, "-max_pair_normal_angle")) { argc--; argv++; max_pair_normal_angle = atof(*argv); }
      else if (!strcmp(*argv, "-min_pair_affinity")) { argc--; argv++; min_pair_affinity = atof(*argv); }
      else if (!strcmp(*argv, "-min_pair_affinity_ratio")) { argc--; argv++; min_pair_affinity_ratio = atof(*argv); }
      else if (!strcmp(*argv, "-large_planes")) {
        min_cluster_points = 10000;
        min_cluster_diameter = 4;
        max_cluster_diameter = 0;
        max_cluster_offplane_distance = 1;
        max_pair_centroid_distance = 50;
        max_pair_offplane_distance = 1;
        max_pair_normal_angle = RN_PI / 4.0;
        max_neighbor_distance = 0.5;
        max_neighbor_normal_angle = RN_PI / 8.0;
        min_pair_affinity = 1.0E-6;
      }
      else if (!strcmp(*argv, "-outdoor_objects") || !strcmp(*argv, "-small_objects")) {
        min_cluster_points = 10;
        min_cluster_diameter = 1;
        max_cluster_diameter = 20;
        max_cluster_offplane_distance = 0;
        max_pair_centroid_distance = 10;
        max_pair_offplane_distance = 0;
        max_pair_normal_angle = 0;
        max_neighbor_distance = 0.5;
        max_neighbor_normal_angle = 0;
        min_pair_affinity_ratio = 0.75;
        min_pair_affinity = 1.0E-6;
      }
      else if (!strcmp(*argv, "-indoor_planes")) {
        min_cluster_points = 100;
        min_cluster_diameter = 0.25;
        max_cluster_diameter = 0;
        max_cluster_offplane_distance = 0.1;
        max_pair_centroid_distance = 4;
        max_pair_offplane_distance = 0.1;
        max_pair_normal_angle = RN_PI / 3.0;
        max_neighbor_distance = 0.1;
        max_neighbor_normal_angle = RN_PI / 4.0;
        min_pair_affinity = 1.0E-6;
      }
      else if (!strcmp(*argv, "-indoor_superpixels")) {
        min_cluster_points = 10;
        min_cluster_diameter = 0.1;
        max_cluster_diameter = 0;
        max_cluster_offplane_distance = 0.1;
        max_pair_centroid_distance = 4;
        max_pair_offplane_distance = 0.1;
        max_pair_normal_angle = RN_PI / 3.0;
        max_neighbor_distance = 0.1;
        max_neighbor_normal_angle = RN_PI / 4.0;
        min_pair_affinity = 1.0E-6;
      }
      else if (!strcmp(*argv, "-indoor_objects")) {
        min_cluster_points = 10;
        min_cluster_diameter = 1;
        max_cluster_diameter = 10;
        max_cluster_offplane_distance = 0;
        max_pair_centroid_distance = 4;
        max_pair_offplane_distance = 0;
        max_pair_normal_angle = 0;
        max_neighbor_distance = 0.25;
        max_neighbor_normal_angle = 0;
        min_pair_affinity_ratio = 0.75;
        min_pair_affinity = 1.0E-6;
      }
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
    else {
      if (!scene_name) scene_name = *argv;
      else if (!database_name) database_name = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check filenames
  if (!scene_name || !database_name) {
    fprintf(stderr, "Usage: sflsegment scenefile databasefile [options]\n");
    return FALSE;
  }

  // Check for default settings
  if (!support_grid_name && !mask_grid_name && !cluster_mean_shift) cluster_hierarchically = 1;

  // Return OK status 
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // Parse program arguments
  if (!ParseArgs(argc, argv)) exit(-1);

  // Open scene
  R3SurfelScene *scene = OpenScene(scene_name, database_name);
  if (!scene) exit(-1);

  // Split out surfels selected by mask
  if (mask_grid_name) {
    if (!SegmentMaskGrid(scene, parent_object_name, parent_node_name, source_node_name, mask_grid_name)) return 0;
  }

  // Split out surfels selected by overhead grids
  if (support_grid_name) {
    if (!SegmentSupportGrid(scene, parent_object_name, parent_node_name, source_node_name, support_grid_name)) return 0;
  }

  // Cluster hierarchically
  if (cluster_hierarchically) {
    if (!ClusterHierarchically(scene, parent_object_name, parent_node_name, source_node_name)) return 0;
  }

  // Cluster with mean shift
  if (cluster_mean_shift) {
    if (!ClusterMeanShift(scene, parent_object_name, parent_node_name, source_node_name)) return 0;
  }

  // Close scene
  if (!CloseScene(scene)) exit(-1);

  // Return success 
  return 0;
}


