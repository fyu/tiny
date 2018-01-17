// Program to align models



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include <R3Shapes/R3Shapes.h>



////////////////////////////////////////////////////////////////////////
// Program arguments
////////////////////////////////////////////////////////////////////////

static char *query_mesh_name = NULL;
static char *reference_mesh_name = NULL;
static char *correspondence_name = NULL;
static char *query_point_name = NULL;
static char *reference_point_name = NULL;
static char *query_arff_name = NULL;
static char *reference_arff_name = NULL;
static char *query_distance_name = NULL;
static char *reference_distance_name = NULL;
static int num_points = 64;
static int distance_type = 0; // 0=dijkstra, 1=euclidean
static int search_type = 2; // 0=branch and bound, 1=priority driven, 2=beam
static int num_correspondences = 8;
static int max_beam_search_candidates = 32*1024;
static RNScalar max_span_length_distance = 0.5;
static RNScalar span_length_scale = 1;
static int print_verbose = 0;
static int print_debug = 0;



////////////////////////////////////////////////////////////////////////
// Type definitions
////////////////////////////////////////////////////////////////////////

struct Point {
  int model_id;
  int vertex_id;
  R3Point position;
  R3MeshVertex *vertex;
  float *distances;
  int ndistances;
  float *features;
  int nfeatures;
};

struct Model {
  char name[256];
  R3Mesh *mesh;
  RNArray<Point *> points;
};

struct Match {
  int npoints;
  static const int num_points = 64;
  Point *query_points[num_points];
  Point *reference_points[num_points];
  RNScalar feature_distance;
  RNScalar total_distance;
};



////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////

struct VertexData {
  R3MeshVertex *vertex;
  double distance;
  VertexData **heappointer;
};


static float *
ComputeDijkstraDistances(R3Mesh *mesh, const RNArray<R3MeshVertex *>& source_vertices)
{
  // Compute distances from source_vertices to all other vertices

  // Allocate temporary vertex data
  VertexData *vertex_data = new VertexData [ mesh->NVertices() ];
  if (!vertex_data) {
    fprintf(stderr, "Unable to allocate vertex data\n");
    return 0;
  }

  // Initialize vertex data
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(i);
    VertexData *data = &vertex_data[i];
    data->vertex = vertex;
    data->distance = FLT_MAX;
    data->heappointer = NULL;
  }

  // Initialize heap
  VertexData tmp;
  RNHeap<VertexData *> heap(&tmp, &(tmp.distance), &(tmp.heappointer));
  for (int i = 0; i < source_vertices.NEntries(); i++) {
    R3MeshVertex *source_vertex = source_vertices.Kth(i);
    VertexData *data = &vertex_data[mesh->VertexID(source_vertex)];
    data->distance = 0;
    heap.Push(data);
  }

  // Iteratively visit vertices
  while (!heap.IsEmpty()) {
    VertexData *data = heap.Pop();
    R3MeshVertex *vertex = data->vertex;
    for (int j = 0; j < mesh->VertexValence(vertex); j++) {
      R3MeshEdge *edge = mesh->EdgeOnVertex(vertex, j);
      R3MeshVertex *neighbor_vertex = mesh->VertexAcrossEdge(edge, vertex);
      VertexData *neighbor_data = &vertex_data[mesh->VertexID(neighbor_vertex)];
      RNScalar old_distance = neighbor_data->distance;
      RNScalar new_distance = mesh->EdgeLength(edge) + data->distance;
      if (new_distance < old_distance) {
        neighbor_data->distance = new_distance;
        if (old_distance < FLT_MAX) heap.Update(neighbor_data);
        else heap.Push(neighbor_data);
      }
    }
  }

  // Copy distance values into output array
  int ndistances = mesh->NVertices();
  float *distances = new float [ ndistances ];
  if (!distances) { fprintf(stderr, "Unable to allocate distances\n"); return 0; }
  for (int i = 0; i < ndistances; i++) distances[i] = vertex_data[i].distance;

  // Delete vertex data
  delete [] vertex_data;

  // Return distances
  return distances;
}



static float *
ComputeDijkstraDistances(R3Mesh *mesh, R3MeshVertex *source_vertex)
{
  // Compute distances from source_vertex to all other vertices
  RNArray<R3MeshVertex *> source_vertices;
  source_vertices.Insert(source_vertex);
  return ComputeDijkstraDistances(mesh, source_vertices);
}



static float *
ComputeEuclideanDistances(R3Mesh *mesh, const RNArray<R3MeshVertex *>& source_vertices)
{
  // Allocate memory for distances
  float *distances = new float [ mesh->NVertices() ];
  if (!distances) {
    fprintf(stderr, "Unable to allocate distances\n");
    return NULL;
  }

  // Compute distances
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(i);
    const R3Point& position = mesh->VertexPosition(vertex);

    // Compute distance to closest source vertex
    float distance = FLT_MAX;
    for (int j = 0; j < source_vertices.NEntries(); j++) {
      R3MeshVertex *source_vertex = source_vertices.Kth(j);
      const R3Point& source_position = mesh->VertexPosition(source_vertex);
      float d = (float) R3Distance(position, source_position);
      if (d < distance) distance = d;
    }

    // Assign distance
    distances[i] = distance;
  }

  // Return distances
  return distances;
}



static float *
ComputeDistances(R3Mesh *mesh, const RNArray<R3MeshVertex *>& source_vertices)
{
  // Compute distances from source_vertex to all other vertices
  if (distance_type == 0) return ComputeDijkstraDistances(mesh, source_vertices);
  else return ComputeEuclideanDistances(mesh, source_vertices);
}



static float *
ComputeDistances(R3Mesh *mesh, R3MeshVertex *source_vertex)
{
  // Compute distances from source_vertex to all other vertices
  RNArray<R3MeshVertex *> source_vertices;
  source_vertices.Insert(source_vertex);
  return ComputeDistances(mesh, source_vertices);
}



static int
CreatePoints(Model *model, int num_points) 
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Get convenient variables
  R3Mesh *mesh = model->mesh;

  // Select random starting vertex
  RNArray<R3MeshVertex *> selected_vertices;
  int i = (int) (RNRandomScalar() * mesh->NVertices());
  R3MeshVertex *vertex = mesh->Vertex(i);
  selected_vertices.Insert(vertex);

  // Iteratively select furthest vertex
  for (int i = 0; i < num_points; i++) {

    // Compute distances from selected vertices to all vertices
    float *distances = ComputeDistances(mesh, selected_vertices);

    // Find furthest vertex
    float furthest_distance = 0;
    R3MeshVertex *furthest_vertex = NULL;
    for (int j = 0; j < mesh->NVertices(); j++) {
      if (distances[j] > furthest_distance) {
        furthest_distance = distances[j];
        furthest_vertex = mesh->Vertex(j);
      }
    }

    // Check if found furthest vertex
    if (!furthest_vertex) break;

    // Add furthest vertex to selected set
    if (i == 0) selected_vertices.Truncate(0);
    selected_vertices.Insert(furthest_vertex);

    // Create point
    Point *point = new Point();
    point->model_id = model->points.NEntries();
    point->vertex_id = mesh->VertexID(furthest_vertex);
    point->position = mesh->VertexPosition(furthest_vertex);
    point->vertex = furthest_vertex;
    point->distances = NULL;
    point->ndistances = 0;
    point->features = NULL;
    point->nfeatures = 0;
    model->points.Insert(point);

    // Delete distances
    delete [] distances;
  }

  // Print statistics
  if (print_verbose) {
    printf("Created points ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Points = %d\n", model->points.NEntries());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
CreateEuclideanDistances(Model *model)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int ndistances = 0;

  // Allocate memory for distances
  int npoints = model->points.NEntries();
  for (int i = 0; i < npoints; i++) {
    Point *point = model->points.Kth(i);
    point->distances = new float [ npoints ];
    for (int j = 0; j < npoints; j++) point->distances[j] = 0;
  }

  // Compute distance normalization factor
  float radius = (float) model->mesh->AverageRadius();
  float normalization = (radius > 0) ? 1.0f / radius : 1.0f;

  // Compute distances
  for (int i = 0; i < model->points.NEntries(); i++) {
    Point *point1 = model->points.Kth(i);
    for (int j = i+1; j < model->points.NEntries(); j++) {
      Point *point2 = model->points.Kth(j);

      // Compute distance 
      float distance = (float) R3Distance(point1->position, point2->position);

      // Normalize and assign distance
      distance *= normalization;
      model->points[i]->distances[j] = distance;
      model->points[j]->distances[i] = distance;
      ndistances++;
    }
  }

  // Print statistics
  if (print_verbose) {
    printf("Created Euclidean distances ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Points = %d\n", model->points.NEntries());
    printf("  # Distances = %d\n", ndistances);
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
CreateDijkstraDistances(Model *model)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int ndistances = 0;

  // Allocate memory for distances
  int npoints = model->points.NEntries();
  for (int i = 0; i < npoints; i++) {
    Point *point = model->points.Kth(i);
    point->distances = new float [ npoints ];
    for (int j = 0; j < npoints; j++) point->distances[j] = 0;
  }

  // Compute distance normalization factor
  float radius = (float) model->mesh->AverageRadius();
  float normalization = (radius > 0) ? 1.0f / radius : 1.0f;

  // Compute distances
  for (int i = 0; i < model->points.NEntries(); i++) {
    Point *point1 = model->points.Kth(i);

    // Compute dijkstra distances to all vertices
    float *distances = ComputeDijkstraDistances(model->mesh, point1->vertex);

    // Assign distances to points
    for (int j = 0; j < model->points.NEntries(); j++) {
      Point *point2 = model->points.Kth(j);

      // Assign distance
      float distance = distances[point2->vertex_id];

      // Normalize and assign distance
      distance *= normalization;
      model->points[i]->distances[j] = distance;
      model->points[j]->distances[i] = distance;
      ndistances++;
    }

    // Delete dijkstra distances from point1
    delete [] distances;
  }

  // Print statistics
  if (print_verbose) {
    printf("Created dijkstra distances ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Points = %d\n", model->points.NEntries());
    printf("  # Distances = %d\n", ndistances);
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
CreateDistances(Model *model)
{
  // Compute pairwise distances between points
  if (distance_type == 0) return CreateDijkstraDistances(model);
  else return CreateEuclideanDistances(model);
}



static int
CreateDistanceHistogramDescriptors(Model *model)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate memory for features
  int nfeatures = 32;
  for (int i = 0; i < model->points.NEntries(); i++) {
    Point *point = model->points.Kth(i);
    point->features = new float [ nfeatures ];
    point->nfeatures = nfeatures;
    for (int j = 0; j < nfeatures; j++) point->features[j] = 0;
  }

  // Compute distance normalization factor
  float radius = (float) model->mesh->AverageRadius();
  float normalization = (radius > 0) ? 1.0f / radius : 1.0f;
  normalization *= nfeatures / 3.0;
  float vote = 1.0 / model->mesh->NVertices();

  // Compute distances
  for (int i = 0; i < model->points.NEntries(); i++) {
    Point *point = model->points.Kth(i);

    // Compute distances to all vertices
    float *distances = ComputeDistances(model->mesh, point->vertex);

    // Fill histogram of distances
    for (int j = 0; j < model->mesh->NVertices(); j++) {
      float bin = normalization * distances[j];
      int bin1 = (int) bin;
      int bin2 = bin1 + 1;
      float t = bin - bin1;
      if (bin1 >= nfeatures) bin1 = nfeatures-1;
      if (bin2 >= nfeatures) bin2 = nfeatures-1;
      point->features[bin1] += (1-t) * vote;
      point->features[bin2] += t * vote;
    }

    // Delete distances from point1
    delete [] distances;
  }

  // Print statistics
  if (print_verbose) {
    printf("Created histogram of distances ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Points = %d\n", model->points.NEntries());
    printf("  # Features = %d\n", nfeatures);
    printf("  Distance type = %d\n", distance_type);
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int 
CompareMatches(const void *m1, const void *m2)
{
  const Match *match1 = *((const Match **) m1);
  const Match *match2 = *((const Match **) m2);
  RNScalar delta = match1->total_distance - match2->total_distance;
  if (delta < 0) return -1;
  else if (delta > 0) return 1;
  else return 0;
}



////////////////////////////////////////////////////////////////////////
// Input functions
////////////////////////////////////////////////////////////////////////

static int
ReadMesh(Model *model, const char *mesh_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate mesh
  R3Mesh *mesh = new R3Mesh();
  if (!mesh) {
    fprintf(stderr, "Unable to allocate mesh\n");
    return 0;
  }

  // Read mesh from file
  if (!mesh->ReadFile(mesh_name)) {
    delete mesh;
    return 0;
  }

  // Associate mesh with model
  model->mesh = mesh;

  // Print statistics
  if (print_verbose) {
    printf("Read mesh from %s ...\n", mesh_name);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Faces = %d\n", mesh->NFaces());
    printf("  # Edges = %d\n", mesh->NEdges());
    printf("  # Vertices = %d\n", mesh->NVertices());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
ReadPoints(Model *model, const char *point_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Open file
  FILE *fp = fopen(point_name, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open points file: %s\n", point_name);
    return 0;
  }

  // Read points from file (vertexid is first column of each row)
  char buffer[2048];
  while (fgets(buffer, 2048, fp)) {
    int vertex_index;
    if (sscanf(buffer, "%d", &vertex_index) == 1) {
      // Check vertex
      if ((vertex_index < 0) || (vertex_index >= model->mesh->NVertices())) {
        fprintf(stderr, "Invalid vertex index %d in %s\n", vertex_index, point_name);
        return 0;
      }

      // Find vertex in mesh
      R3MeshVertex *vertex = model->mesh->Vertex(vertex_index);
      if (!vertex) {
        fprintf(stderr, "Unable to find vertex %d in %s\n", vertex_index, point_name);
        return 0;
      }

      // Create point
      Point *point = new Point();
      point->model_id = model->points.NEntries();
      point->vertex_id = model->mesh->VertexID(vertex);
      point->position = model->mesh->VertexPosition(vertex);
      point->vertex = vertex;
      point->distances = NULL;
      point->ndistances = 0;
      point->features = NULL;
      point->nfeatures = 0;

      // Insert point into model
      model->points.Insert(point);
    }
  }

  // Close file
  fclose(fp);

  // Print statistics
  if (print_verbose) {
    printf("Read points from %s ...\n", point_name);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Points = %d\n", model->points.NEntries());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
ReadDescriptors(Model *model, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Read properties from file
  R3MeshPropertySet properties(model->mesh);
  if (!properties.Read(filename)) return 0;

  // Allocate memory for features
  for (int i = 0; i < model->points.NEntries(); i++) {
    Point *point = model->points.Kth(i);
    point->features = new float [ properties.NProperties() ];
    point->nfeatures = properties.NProperties();
    for (int j = 0; j < properties.NProperties(); j++) {
      R3MeshProperty *property = properties.Property(j);
      point->features[j] = property->VertexValue(point->vertex_id) ;
    }
  }

  // Print statistics
  if (print_verbose) {
    printf("Read arff file %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Points = %d\n", model->points.NEntries());
    printf("  # Features = %d\n", properties.NProperties());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
ReadDistances(Model *model, const char *distance_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int ndistances = 0;

  // Open file
  FILE *fp = fopen(distance_name, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open distance file: %s\n", distance_name);
    return 0;
  }

  // Allocate memory for distances
  int npoints = model->points.NEntries();
  for (int i = 0; i < npoints; i++) {
    Point *point = model->points.Kth(i);
    point->distances = new float [ npoints ];
    for (int j = 0; j < npoints; j++) point->distances[j] = 0;
  }

  // Compute distance normalization factor
  float radius = (float) model->mesh->AverageRadius();
  float normalization = (radius > 0) ? 1.0f / radius : 1.0f;

  // Read distances
  for (int i = 0; i < npoints; i++) {
    for (int j = i+1; j < npoints; j++) {
      float distance;
      if (fscanf(fp, "%f", &distance) != 1) {
        fprintf(stderr, "Error reading distances file %s\n", distance_name);
        return 0;
      }

      // Normalize and assign distance
      distance *= normalization;
      model->points[i]->distances[j] = distance;
      model->points[j]->distances[i] = distance;
      ndistances++;
    }
  }

  // Close file
  fclose(fp);

  // Print statistics
  if (print_verbose) {
    printf("Read distances %s ...\n", distance_name);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Points = %d\n", npoints);
    printf("  # Distances = %d\n", ndistances);
    fflush(stdout);
  }

  // Return success
  return 1;
}



static Model *
ReadModel(const char *mesh_name, const char *point_name, const char *arff_name, const char *distance_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate model
  Model *model = new Model();
  if (!model) {
    fprintf(stderr, "Unable to allocate model\n");
    return NULL;
  }

  // Copy name
  strcpy(model->name, mesh_name);

  // Read mesh
  if (!ReadMesh(model, mesh_name)) return NULL;

  // Read/create points
  if (point_name) { if (!ReadPoints(model, point_name)) return NULL; }
  else { if (!CreatePoints(model, num_points)) return NULL; }

  // Read/create shape descriptor
  if (arff_name) { if (!ReadDescriptors(model, arff_name)) return NULL; }
  else { if (!CreateDistanceHistogramDescriptors(model)) return NULL; }

  // Read/create distances
  if (distance_name) { if (!ReadDistances(model, distance_name)) return NULL; }
  else { if (!CreateDistances(model)) return NULL; }

  // Return model
  return model;
}


  
////////////////////////////////////////////////////////////////////////
// File output functions
////////////////////////////////////////////////////////////////////////

static int
WriteMatch(Match *match, Model *query_model, Model *reference_model, const char *correspondence_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Open file
  FILE *fp = fopen(correspondence_name, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open query correspondence file: %s\n", correspondence_name);
    return 0;
  }

  // Write corresponding points (2  vertex_id x y z  vertex_id x y z)
  for (int i = 0; i < match->npoints; i++) {
    Point *query_point = match->query_points[i];
    Point *reference_point = match->reference_points[i];
    fprintf(fp, "%6d %6d\n", query_point->vertex_id, reference_point->vertex_id);
  }

  // Close file
  fclose(fp);

  // Print statistics
  if (print_verbose) {
    printf("Wrote match to %s ...\n", correspondence_name);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Points = %d\n", match->npoints);
    printf("  Total distance = %g\n", match->total_distance);
    printf("  Feature distance = %g\n", match->feature_distance);
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Pairwise match creation
////////////////////////////////////////////////////////////////////////

static RNArray<Match *> *
CreatePairs(Model *query_model, Model *reference_model)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate array of pairs to return
  RNArray<Match *> *pairs = new RNArray<Match *>();
  if (!pairs) {
    fprintf(stderr, "Unable to allocate pair array\n");
    return NULL;
  }

  // Allocate temporary arrays of pairs (for sorting and ranking)
  RNArray<Match *> *query_pairs = new RNArray<Match *> [ query_model->points.NEntries() ];
  RNArray<Match *> *reference_pairs = new RNArray<Match *> [ reference_model->points.NEntries() ];

  // Create pairs between all pairs of points
  for (int i = 0; i < query_model->points.NEntries(); i++) {
    Point *query_point = query_model->points.Kth(i);
    float *query_features = query_point->features;
    int nfeatures = query_point->nfeatures;

    for (int j = 0; j < reference_model->points.NEntries(); j++) {
      Point *reference_point = reference_model->points.Kth(j);
      if (reference_point == query_point) continue;
      float *reference_features = reference_point->features;
      if (reference_point->nfeatures != nfeatures) continue;

      // Compute distance
      float distance = 0;
      for (int k = 0; k < nfeatures; k++) {
        float delta = query_features[k] - reference_features[k];
        distance += delta * delta;
      }

      // Create match
      Match *match = new Match();
      match->query_points[0] = query_point;
      match->reference_points[0] = reference_point;
      match->feature_distance = distance;
      match->total_distance = distance;
      match->npoints = 1;

      // Add match to arrays 
      query_pairs[i].Insert(match);
      reference_pairs[j].Insert(match);
    }
  }

  // Determine maximum rank of any pair
  int max_query_rank = query_model->points.NEntries() / 4; 
  int max_reference_rank = reference_model->points.NEntries() / 4; 

  // Sort/rank pairs for each query point
  for (int i = 0; i < query_model->points.NEntries(); i++) {
    query_pairs[i].Sort(CompareMatches);
    for (int j = 0; j < query_pairs[i].NEntries(); j++) {
      Match *query_pair = query_pairs[i].Kth(j);
      if (j < max_query_rank) {
        query_pair->feature_distance = j;
        query_pair->total_distance = j;
      }
      else {
        query_pair->feature_distance = RN_INFINITY;
        query_pair->total_distance = RN_INFINITY;
      }
    }
  }

  // Sort/rank pairs for each reference point
  for (int i = 0; i < reference_model->points.NEntries(); i++) {
    reference_pairs[i].Sort(CompareMatches);
    for (int j = 0; j < reference_pairs[i].NEntries(); j++) {
      Match *reference_pair = reference_pairs[i].Kth(j);
      if (j < max_reference_rank) {
        reference_pair->feature_distance += j;
        reference_pair->total_distance += j;
      }
      else {
        reference_pair->feature_distance = RN_INFINITY;
        reference_pair->total_distance = RN_INFINITY;
      }
    }
  }

  // Make single list of pairs with good ranks
  for (int i = 0; i < query_model->points.NEntries(); i++) {
    for (int j = 0; j < query_pairs[i].NEntries(); j++) {
      Match *pair = query_pairs[i].Kth(j);
      if (pair->total_distance >= RN_INFINITY) delete pair;
      else pairs->Insert(pair);
    }
  }

  // Delete temporary arrays of pairs
  delete [] query_pairs;
  delete [] reference_pairs;

  // Print messages
  if (print_verbose) {
    printf("Created pairs ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Pairs = %d\n", pairs->NEntries());
    fflush(stdout); 
  }

  // Return array of pairs
  return pairs;
}



////////////////////////////////////////////////////////////////////////
// Combinatorial search for best match using priority driven search
////////////////////////////////////////////////////////////////////////

static Match *
FindBestMatchWithPriorityDrivenSearch(Model *query_model, Model *reference_model, int num_correspondences, int max_heap_size = 0)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int npops = 0;

  // Initialize best match distance
  Match *best_match = new Match();
  RNScalar best_distance = RN_INFINITY;

  // Compute pairs
  RNArray<Match *> *pairs = CreatePairs(query_model, reference_model);
  if (!pairs) {
    fprintf(stderr, "Unable to create pairs\n");
    return NULL;
  }

  // Intialize heap with pairs
  Match tmp; int value_offset = (unsigned char *) &(tmp.total_distance) - (unsigned char *) &tmp;
  RNHeap<Match *> heap(value_offset);
  for (int i = 0; i < pairs->NEntries(); i++) {
    Match *pair = pairs->Kth(i);
    heap.Push(pair);
  }

  // Iteratively grow matches until find one with num_correspondences
  while (TRUE) {
    // Match *match = heap.Pop();
    Match *match = heap.Pop();
    if (!match) break;
    npops++;

    // Print debug message
    if (print_debug) {
      printf("%d %d %g %g\n", npops, match->npoints, match->feature_distance, match->total_distance);
    }

    // Check if we found a complete match
    if (match->npoints == num_correspondences) {
      *best_match = *match;
      break;
    }

    // Consider all plausible extensions to match
    for (int i = 0; i < pairs->NEntries(); i++) {
      Match *pair = pairs->Kth(i);
      Point *query_point = pair->query_points[0];
      Point *reference_point = pair->reference_points[0];

      // Check if order is preserved among query points (avoid different permutations of same points)
      if (query_point->model_id <= match->query_points[match->npoints-1]->model_id) continue;

      // Check if points of pair already appear in match
      RNBoolean ok = TRUE;
      for (int j = 0; j < match->npoints; j++) {
        if ((query_point == match->query_points[j]) || (reference_point == match->reference_points[j])) {
          ok = FALSE;
          break;
        }
      }
      if (!ok) continue;

      // Start from original distance
      RNScalar total_distance = match->total_distance;
      RNScalar feature_distance = match->feature_distance;

      // Add distance due to extra pair
      total_distance += pair->total_distance;
      feature_distance += pair->feature_distance;

      // Check if as good as best complete match seen so far
      if (total_distance >= best_distance) continue;

      // Add distance due to span length differences
      RNBoolean plausible = TRUE;
      for (int j = 0; j < match->npoints; j++) {
        RNLength query_span_length = query_point->distances[match->query_points[j]->model_id];
        RNLength reference_span_length = reference_point->distances[match->reference_points[j]->model_id];
        RNScalar greatest_span_length = (query_span_length > reference_span_length) ? query_span_length : reference_span_length;
        RNScalar span_length_distance = fabs(query_span_length - reference_span_length);
        if (greatest_span_length > 0) span_length_distance /= greatest_span_length;
        if (span_length_distance > max_span_length_distance) { plausible = FALSE; break; }
        total_distance += 4 * span_length_scale * span_length_distance;
      }
      if (!plausible) continue;

      // Check if as good as best complete match seen so far
      if (total_distance >= best_distance) continue;

      // Compute new match
      Match *extended_match = new Match(*match);
      extended_match->total_distance = total_distance;
      extended_match->feature_distance = feature_distance;
      extended_match->query_points[extended_match->npoints] = pair->query_points[0];
      extended_match->reference_points[extended_match->npoints] = pair->reference_points[0];
      extended_match->npoints++;

      // Update best match seen so far
      if ((extended_match->npoints > best_match->npoints) ||
          ((extended_match->npoints == best_match->npoints) && 
           (extended_match->total_distance <= best_match->total_distance))) {
        // Update best match
        *best_match = *extended_match;

        // Update best distance
        if (best_match->npoints == num_correspondences) {
          best_distance = best_match->total_distance;
        }
      }

      // Push match onto heap
      heap.Push(extended_match);

      // Check if heap is too big
      if (max_heap_size > 0) {
        if (heap.NEntries() > max_heap_size) {
          int target_size = max_heap_size/2;
          heap.Sort(target_size);
          for (int j = target_size; j < heap.NEntries(); j++) {
            Match *m = heap.Kth(j);
            if (m->npoints > 1) delete m;
          }
          heap.Truncate(target_size, FALSE);
        }
      }
    }

    // Delete match
    if (match->npoints > 1) delete match;
  }

  // Print messages
  if (print_verbose) {
    if (best_match) {
      printf("Found match ...\n");
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      printf("  # Pops = %d\n", npops);
      printf("  # Correspondences = %d\n", best_match->npoints);
      printf("  Feature distance = %g\n", best_match->feature_distance);
      printf("  Total distance = %g\n", best_match->total_distance);
      printf("  Error = %g\n", 1000 * (num_correspondences - best_match->npoints) + best_match->total_distance);
    }
    else {
      printf("Unable to find match ...\n");
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      printf("  # Pops = %d\n", npops);
    }
    fflush(stdout); 
  }

  // Return best match
  return best_match;
}



////////////////////////////////////////////////////////////////////////
// Combinatorial search for best match using beam search
////////////////////////////////////////////////////////////////////////

static Match *
FindBestMatchWithBeamSearch(Model *query_model, Model *reference_model, int num_correspondences)
{
  // Use priority-driven search, but limit heap to max_beam_search_candidates
  return FindBestMatchWithPriorityDrivenSearch(query_model, reference_model, num_correspondences, max_beam_search_candidates);
}



////////////////////////////////////////////////////////////////////////
// Combinatorial search for best match using branch and bound
////////////////////////////////////////////////////////////////////////

void FindBestMatchWithBranchAndBoundSearch(Model *query_model, Model *reference_model, int num_correspondences, 
  const RNArray<Match *>& pairs, const Match& match, Match& best_match, RNScalar& best_distance, int& npops)
{
  // Update statistics
  npops++;

  // Print debug message
  if (print_debug) {
    printf("%d %d %g %g\n", npops, match.npoints, match.feature_distance, match.total_distance);
  }

  // Consider all plausible extensions to match
  if ((match.npoints < num_correspondences) && (match.total_distance < best_distance)) {
    for (int i = 0; i < pairs.NEntries(); i++) {
      Match *pair = pairs.Kth(i);
      Point *query_point = pair->query_points[0];
      Point *reference_point = pair->reference_points[0];

      // Check if order is preserved among query points (avoid different permutations of same points)
      if (query_point->model_id <= match.query_points[match.npoints-1]->model_id) continue;

      // Check if points of pair already appear in match
      RNBoolean ok = TRUE;
      for (int j = 0; j < match.npoints; j++) {
        if ((query_point == match.query_points[j]) || (reference_point == match.reference_points[j])) {
          ok = FALSE;
          break;
        }
      }
      if (!ok) continue;

      // Start from original distance
      RNScalar total_distance = match.total_distance;
      RNScalar feature_distance = match.feature_distance;

      // Add distance due to extra pair
      total_distance += pair->total_distance;
      feature_distance += pair->feature_distance;

      // Check if as good as best complete match seen so far
      if (total_distance >= best_distance) continue;

      // Add distance due to span length differences
      RNBoolean plausible = TRUE;
      for (int j = 0; j < match.npoints; j++) {
        RNLength query_span_length = query_point->distances[match.query_points[j]->model_id];
        RNLength reference_span_length = reference_point->distances[match.reference_points[j]->model_id];
        RNScalar greatest_span_length = (query_span_length > reference_span_length) ? query_span_length : reference_span_length;
        RNScalar span_length_distance = fabs(query_span_length - reference_span_length);
        if (greatest_span_length > 0) span_length_distance /= greatest_span_length;
        if (span_length_distance > max_span_length_distance) { plausible = FALSE; break; }
        total_distance += 4 * span_length_scale * span_length_distance;
      }
      if (!plausible) continue;

      // Check if as good as best complete match seen so far
      if (total_distance >= best_distance) continue;

      // Compute new match
      Match extended_match(match);
      extended_match.total_distance = total_distance;
      extended_match.feature_distance = feature_distance;
      extended_match.query_points[extended_match.npoints] = pair->query_points[0];
      extended_match.reference_points[extended_match.npoints] = pair->reference_points[0];
      extended_match.npoints++;

      // Update best match seen so far
      if ((extended_match.npoints > best_match.npoints) ||
          ((extended_match.npoints == best_match.npoints) && 
           (extended_match.total_distance <= best_match.total_distance))) {
        // Update best match
        best_match = extended_match;

        // Update best distance
        if (best_match.npoints == num_correspondences) {
          best_distance = best_match.total_distance;
        }
      }

      // Recurse
      FindBestMatchWithBranchAndBoundSearch(query_model, reference_model, num_correspondences, 
        pairs, extended_match, best_match, best_distance, npops);
    }
  }
}


static Match *
FindBestMatchWithBranchAndBoundSearch(Model *query_model, Model *reference_model, int num_correspondences)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Initialize state
  Match *best_match = new Match();
  RNScalar best_distance = RN_INFINITY;
  int npops = 0;

  // Compute pairs
  RNArray<Match *> *pairs = CreatePairs(query_model, reference_model);
  if (!pairs) {
    fprintf(stderr, "Unable to create pairs\n");
    return NULL;
  }

  // Recursively build matches starting with all pairs
  for (int i = 0; i < pairs->NEntries(); i++) {
    Match *pair = pairs->Kth(i);
    FindBestMatchWithBranchAndBoundSearch(query_model, reference_model, num_correspondences, 
        *pairs, *pair, *best_match, best_distance, npops);
  }

  // Print messages
  if (print_verbose) {
    if (best_match) {
      printf("Found match ...\n");
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      printf("  # Pops = %d\n", npops);
      printf("  # Correspondences = %d\n", best_match->npoints);
      printf("  Feature distance = %g\n", best_match->feature_distance);
      printf("  Total distance = %g\n", best_match->total_distance);
      printf("  Error = %g\n", 1000 * (num_correspondences - best_match->npoints) + best_match->total_distance);
    }
    else {
      printf("Unable to find match ...\n");
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      printf("  # Pops = %d\n", npops);
    }
    fflush(stdout); 
  }

  // Return best match
  return best_match;
}



////////////////////////////////////////////////////////////////////////
// Combinatorial search for best match
////////////////////////////////////////////////////////////////////////

static Match *
FindBestMatch(Model *query_model, Model *reference_model, int num_correspondences, int search_type)
{
  // Initialize result
  Match *match = NULL;

  // Compute match
  if (search_type == 0) {
    match = FindBestMatchWithBranchAndBoundSearch(query_model, reference_model, num_correspondences);
  }
  else if (search_type == 1) {
    match = FindBestMatchWithPriorityDrivenSearch(query_model, reference_model, num_correspondences);
  }
  else if (search_type == 2) {
    match = FindBestMatchWithBeamSearch(query_model, reference_model, num_correspondences);
  }
  else {
    fprintf(stderr, "Unrecognized search algorithm: %d\n", search_type);
  }

  // Return match
  return match;
}



////////////////////////////////////////////////////////////////////////
// Program argument parsing
////////////////////////////////////////////////////////////////////////

int
ParseArgs (int argc, char **argv)
{
  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) { print_verbose = 1; }
      else if (!strcmp(*argv, "-debug")) { print_debug = 1; }
      else if (!strcmp(*argv, "-branch_and_bound_search")) search_type = 0; 
      else if (!strcmp(*argv, "-priority_driven_search")) search_type = 1; 
      else if (!strcmp(*argv, "-search_type")) { ++argv; --argc; search_type = atoi(*argv); }
      else if (!strcmp(*argv, "-dijkstra_distance")) distance_type = 0; 
      else if (!strcmp(*argv, "-euclidean_distance")) distance_type = 1; 
      else if (!strcmp(*argv, "-distance_type")) { ++argv; --argc; distance_type = atoi(*argv); }
      else if (!strcmp(*argv, "-num_points")) { ++argv; --argc; num_points = atoi(*argv); }
      else if (!strcmp(*argv, "-num_correspondences")) { ++argv; --argc; num_correspondences = atoi(*argv); }
      else if (!strcmp(*argv, "-n")) { ++argv; --argc; num_correspondences = atoi(*argv); }
      else if (!strcmp(*argv, "-max_beam_search_candidates")) { ++argv; --argc; max_beam_search_candidates = atoi(*argv); }
      else if (!strcmp(*argv, "-max_span_length_distance")) { ++argv; --argc; max_span_length_distance = atof(*argv); }
      else if (!strcmp(*argv, "-span_length_scale")) { ++argv; --argc; span_length_scale = atof(*argv); }
      else if (!strcmp(*argv, "-points")) { 
        ++argv; --argc; query_point_name = *argv; 
        ++argv; --argc; reference_point_name = *argv; 
      }
      else if (!strcmp(*argv, "-arffs")) { 
        ++argv; --argc; query_arff_name = *argv; 
        ++argv; --argc; reference_arff_name = *argv; 
      }
      else if (!strcmp(*argv, "-distances")) { 
        ++argv; --argc; query_distance_name = *argv; 
        ++argv; --argc; reference_distance_name = *argv; 
      }
      else { 
        fprintf(stderr, "Invalid program argument: %s", *argv); 
        return 0; 
      }
      argv++; argc--;
    }
    else {
      if (!query_mesh_name) query_mesh_name = *argv;
      else if (!reference_mesh_name) reference_mesh_name = *argv;
      else if (!correspondence_name) correspondence_name = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); return 0; }
      argv++; argc--;
    }
  }

  // Check required filename arguments
  if (!query_mesh_name || !reference_mesh_name || !correspondence_name) {
    fprintf(stderr, "Usage: msh2cor query_mesh reference_mesh correspondence_file\n"
                    "       [-n <num_correspondences>] (default=8)\n"
                    "       [-points query_point_file reference_point_file]\n"
                    "       [-arffs query_arff_file reference_arff_file]\n"
                    "       [-distances query_distance_file reference_distance_file]"
                    "       [options]\n");
    return 0;
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program 
////////////////////////////////////////////////////////////////////////

int
main(int argc, char **argv)
{
  // Parse args
  if(!ParseArgs(argc, argv)) exit(-1);

  // Read first model
  Model *query_model = ReadModel(query_mesh_name, query_point_name, query_arff_name, query_distance_name);
  if (!query_model) exit(-1);

  // Read second model
  Model *reference_model = query_model;
  if ((query_mesh_name && reference_mesh_name && (strcmp(query_mesh_name, reference_mesh_name))) ||
      (query_point_name && reference_point_name && (strcmp(query_point_name, reference_point_name))) ||
      (query_arff_name && reference_arff_name && (strcmp(query_arff_name, reference_arff_name))) ||
      (query_distance_name && reference_distance_name && (strcmp(query_distance_name, reference_distance_name)))){
    reference_model = ReadModel(reference_mesh_name, reference_point_name, reference_arff_name, reference_distance_name);
    if (!reference_model) exit(-1);
  }

  // Compute best match
  Match *match = FindBestMatch(query_model, reference_model, num_correspondences, search_type);
  if (!match) exit(-1);

  // Write correspondences
  if (!WriteMatch(match, query_model, reference_model, correspondence_name)) exit(-1);

  // Return success
  return 0;
}







