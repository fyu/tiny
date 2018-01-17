// Source file for the program to compute dense correspondences from sparse ones



////////////////////////////////////////////////////////////////////////
// Include files 
////////////////////////////////////////////////////////////////////////

#include "R3Shapes/R3Shapes.h"



////////////////////////////////////////////////////////////////////////
// Program arguments
////////////////////////////////////////////////////////////////////////

static char *input_mesh1_name = NULL;
static char *input_mesh2_name = NULL;
static char *input_correspondence_name = NULL;
static char *output_map_name = NULL;
static double max_delta = 0;
static double min_sparse_correspondence_spacing = 0.1; 
static int max_sparse_correspondences = 128;
static int add_swapped_correspondences = 0;
static int print_verbose = 0;
static int print_debug = 0;


////////////////////////////////////////////////////////////////////////
// Type definitions
////////////////////////////////////////////////////////////////////////

struct SparseVertexCorrespondence {
  // vertices[0] and vertices[1] have same number of entries
  // For every i, vertices[0][i] and vertices[1][i] are corresponding
  // vertices in mesh[0] and mesh[1], respectively
  R3Mesh *mesh[2];
  RNArray<R3MeshVertex *> vertices[2];
  SparseVertexCorrespondence(R3Mesh *mesh0, R3Mesh *mesh1) { 
    mesh[0] = mesh0; 
    mesh[1] = mesh1; 
  };
};

struct DenseVertexCorrespondence {
  // vertices has mesh[0]->NVertices() entries
  // each vertices[i] is pointer to vertex in mesh[1]
  // corresponding to mesh[0]->Vertex(i) in mesh[0]
  // if vertices[i] is NULL, then mesh[0]->Vertex(i) is outlier
  R3Mesh *mesh[2];
  R3MeshVertex **vertices; 
  DenseVertexCorrespondence(R3Mesh *mesh0, R3Mesh *mesh1) { 
    mesh[0] = mesh0; 
    mesh[1] = mesh1; 
    vertices = new R3MeshVertex * [ mesh0->NVertices() ];
    for (int i = 0; i < mesh0->NVertices(); i++) vertices[i] = NULL;
  };
};



////////////////////////////////////////////////////////////////////////
// Input/Output
////////////////////////////////////////////////////////////////////////

static R3Mesh *
ReadMesh(char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate mesh
  R3Mesh *mesh = new R3Mesh();
  assert(mesh);

  // Read mesh from file
  if (!mesh->ReadFile(filename)) {
    fprintf(stderr, "Unable to read mesh %s\n", filename);
    return NULL;
  }

  // Normalize the mesh
  R3Point centroid = mesh->Centroid();
  RNScalar radius = sqrt( mesh->Area() );
  RNScalar scale = (radius > 0) ? 1.0 / radius : 1;
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(i);
    R3Point position = mesh->VertexPosition(vertex);
    position -= centroid.Vector();
    position *= scale;
    mesh->SetVertexPosition(vertex, position);
  }

  // Print statistics
  if (print_verbose) {
    printf("Read mesh from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Faces = %d\n", mesh->NFaces());
    printf("  # Edges = %d\n", mesh->NEdges());
    printf("  # Vertices = %d\n", mesh->NVertices());
    fflush(stdout);
  }

  // Return success
  return mesh;
}



static SparseVertexCorrespondence *
ReadSparseVertexCorrespondence(R3Mesh *mesh0, R3Mesh *mesh1, char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Parse filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .cor)\n", filename);
    return 0;
  }

  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open correspondence file %s\n", filename);
    return NULL;
  }

  // Create correspondence
  SparseVertexCorrespondence *correspondence = new SparseVertexCorrespondence(mesh0, mesh1);
  if (!correspondence) {
    fprintf(stderr, "Unable to allocate correspondence for %s\n", filename); 
    return NULL; 
  }

  // Check filename extension
  if (!strcmp(extension, ".cor")) {
    // Read correspondences
    int id0, id1;
    while (fscanf(fp, "%d%d", &id0, &id1) == (unsigned int) 2) {
      R3MeshVertex *vertex0 = mesh0->Vertex(id0);
      R3MeshVertex *vertex1 = mesh1->Vertex(id1);
      correspondence->vertices[0].Insert(vertex0);
      correspondence->vertices[1].Insert(vertex1);
      if (add_swapped_correspondences && (vertex0 != vertex1)) {
        correspondence->vertices[0].Insert(vertex1);
        correspondence->vertices[1].Insert(vertex0);
      }
    }
  }
  else if (!strcmp(extension, ".vk")) {
    // Read header
    double fdummy;
    char format[256], dummy[256];
    int nmeshes, ncorrespondences, idummy;
    if (fscanf(fp, "%s%s%s%s%s%d%s%s%s%d%s%d%d%d%d%d", dummy, format,  dummy, dummy,  dummy, &nmeshes,
               dummy, dummy,  dummy, &idummy,   dummy, &idummy, &idummy, &idummy, &idummy, &ncorrespondences) != (unsigned int) 16) {
      fprintf(stderr, "Unable to read %s\n", filename);
      return NULL;
    }

    // Read correspondences
    int id0, id1;
    int count = 0;
    while (fscanf(fp, "%d%lg%d%lg", &id0, &fdummy, &id1, &fdummy) == (unsigned int) 4) {
      R3MeshVertex *vertex0 = mesh0->Vertex(id0);
      R3MeshVertex *vertex1 = mesh1->Vertex(id1);
      correspondence->vertices[0].Insert(vertex0);
      correspondence->vertices[1].Insert(vertex1);
      if (add_swapped_correspondences && (vertex0 != vertex1)) {
        correspondence->vertices[0].Insert(vertex1);
        correspondence->vertices[1].Insert(vertex0);
      }
      count++;
    }

    // Check number of correspondences
    if (count != ncorrespondences) {
      fprintf(stderr, "Mismatching number of correspondences in %s\n", filename);
      return NULL;
    }
  }
  else {
    fprintf(stderr, "Unrecognized correspondence file extension: %s\n", extension);
    return NULL;
  }

  // Close file
  fclose(fp);

  // Print statistics
  if (print_verbose) {
    printf("Read sparse correspondences from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Correspondences = %d\n", correspondence->vertices[0].NEntries());
    fflush(stdout);
  }

  // Return correspondence
  return correspondence;
}



static int
WriteDenseVertexCorrespondence(DenseVertexCorrespondence *map, char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open map file %s\n", filename);
    return 0;
  }

  // Write correspondences
  for (int i = 0; i < map->mesh[0]->NVertices(); i++) {
    if (map->vertices[i] == NULL) fprintf(fp, "-1\n");
    else fprintf(fp, "%d\n", map->mesh[1]->VertexID(map->vertices[i]));
  }

  // Close file
  fclose(fp);

  // Print statistics
  if (print_verbose) {
    printf("Wrote dense correspondences to %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Correspondences = %d\n", map->mesh[0]->NVertices());
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Distance utility functions
////////////////////////////////////////////////////////////////////////

static RNLength ***
CreateDistances(SparseVertexCorrespondence *correspondence)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int count = 0;

  // Create distances
  RNLength ***distances = new RNLength ** [ 2 ];
  for (int i = 0; i < 2; i++) {
    R3Mesh *mesh = correspondence->mesh[i];
    distances[i] = new RNLength * [ mesh->NVertices() ];
    for (int j = 0; j < mesh->NVertices(); j++) distances[i][j] = NULL;
    for (int j = 0; j < correspondence->vertices[i].NEntries(); j++) {
      R3MeshVertex *vertex = correspondence->vertices[i].Kth(j);
      distances[i][mesh->VertexID(vertex)] = mesh->DijkstraDistances(vertex);
      count++;
    }
  }

  // Print statistics
  if (print_verbose) {
    printf("Computed distances ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Distances = %d\n", count);
    fflush(stdout);
  }

  // Return distances
  return distances;
}



static int
DeleteDistances(SparseVertexCorrespondence *correspondence, RNLength ***distances)
{
  // Delete distances
  for (int i = 0; i < 2; i++) {
    R3Mesh *mesh = correspondence->mesh[i];
    for (int j = 0; j < mesh->NVertices(); j++) {
      delete [] distances[i][j];
    }
    delete distances[i];
  }
  delete [] distances;

  // Return success
  return 1;
}



static RNLength 
ComputeDistanceDifference(SparseVertexCorrespondence *correspondence,
  R3MeshVertex *vertex0, R3MeshVertex *vertex1, 
  RNLength ***distances)
{
  // Initialize distance
  RNLength delta = 0;

  // Get convenient variables
  R3Mesh *mesh0 = correspondence->mesh[0];
  R3Mesh *mesh1 = correspondence->mesh[1];

  // Compute distance in feature space
  for (int k = 0; k < correspondence->vertices[0].NEntries(); k++) {
    int k0 = mesh0->VertexID(correspondence->vertices[0].Kth(k));
    int k1 = mesh1->VertexID(correspondence->vertices[1].Kth(k));
    double d0 = distances[0][k0][mesh0->VertexID(vertex0)];
    double d1 = distances[1][k1][mesh1->VertexID(vertex1)];
    double d = fabs(d1 - d0);
    if (d0 > 0) delta += d / d0;
    else if (d > 0) delta += RN_INFINITY; 
  }

  // Return feature space distance
  return delta;
}



////////////////////////////////////////////////////////////////////////
// Sparse Correspondence Expansion
////////////////////////////////////////////////////////////////////////

static RNArray<R3MeshVertex *> *
CreateCandidates(R3Mesh *mesh, RNArray<R3MeshVertex *>& correspondence_vertices, RNScalar *importances, int max_candidates)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate array of candidates
  RNArray<R3MeshVertex *> *candidates = new RNArray<R3MeshVertex *>();
  if (!candidates) {
    fprintf(stderr, "Unable to allocate candidates\n");
    return NULL;
  }

  // Create array of values (importance * distance)
  RNLength *dists = mesh->DijkstraDistances(correspondence_vertices);
  RNScalar *values = new RNScalar [ mesh->NVertices() ];
  for (int i = 0; i < mesh->NVertices(); i++) {
    values[i] = dists[i] * importances[i];
  }

  // Fill array of candidates
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(i);
    if (dists[i] <= 0) continue;
    if (importances[i] < 0) continue;
    candidates->Insert(vertex);
  }

  // Sort the best candidates (insertion sort)
  if (candidates->NEntries() > max_candidates) {
    for (int i = 0; i < max_candidates; i++) {
      for (int j = i+1; j < candidates->NEntries(); j++) {
        RNScalar value_i = values[mesh->VertexID(candidates->Kth(i))];
        RNScalar value_j = values[mesh->VertexID(candidates->Kth(j))];
        if (value_i < value_j) candidates->Swap(i, j);
      }
    }
  }

  // Truncate array of candidates
  candidates->Truncate(max_candidates);

  // Delete temporary data
  delete [] dists;
  delete [] values;

  // Print statistics
  if (print_debug) {
    printf("  Created candidates ...\n");
    printf("    Time = %.2f seconds\n", start_time.Elapsed());
    printf("    # Candidates = %d\n", candidates->NEntries());
    fflush(stdout);
  }

  // Return candidates
  return candidates;
}



static int 
ExpandSparseCorrespondence(SparseVertexCorrespondence *correspondence, RNLength ***distances, RNScalar **importances, int max_candidates, RNLength min_spacing, RNScalar max_delta)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  R3mesh_mark++;
  int count = 0;

  // Get convenient variables
  R3Mesh *mesh0 = correspondence->mesh[0];
  R3Mesh *mesh1 = correspondence->mesh[1];

  // Create candidates
  RNArray<R3MeshVertex *> *candidates0 = CreateCandidates(mesh0, correspondence->vertices[0], importances[0], max_candidates);
  if (!candidates0 || candidates0->IsEmpty()) return 0; 
  RNArray<R3MeshVertex *> *candidates1 = CreateCandidates(mesh1, correspondence->vertices[1], importances[1], max_candidates);
  if (!candidates1 || candidates1->IsEmpty()) return 0; 

  // Create struct to store new correspondences
  SparseVertexCorrespondence mutually_closest_candidates(mesh0, mesh1);

  // Consider each candidate in candidates0
  for (int i = 0; i < candidates0->NEntries(); i++) {
    R3MeshVertex *candidate0 = candidates0->Kth(i);
    if (mesh0->VertexMark(candidate0) == R3mesh_mark) continue;

    // Find closest candidate in candidates1 
    R3MeshVertex *candidate1 = NULL;
    RNScalar distance1 = (max_delta > 0) ? correspondence->vertices[0].NEntries() * max_delta : FLT_MAX;
    for (int j = 0; j < candidates1->NEntries(); j++) {
      R3MeshVertex *vertex1 = candidates1->Kth(j);
      if (mesh1->VertexMark(vertex1) == R3mesh_mark) continue;
      RNLength distance = ComputeDistanceDifference(correspondence, candidate0, vertex1, distances);
      if (distance < distance1) {
        candidate1 = vertex1;
        distance1 = distance;
      }
    }

    // Check if found candidate1
    if (!candidate1) continue;

    // Check if there is a candidate0 closer to candidate1
    for (int j = 0; j < candidates0->NEntries(); j++) {
      R3MeshVertex *vertex0 = candidates0->Kth(j);
      if (mesh0->VertexMark(candidate0) == R3mesh_mark) continue;
      RNLength distance = ComputeDistanceDifference(correspondence, vertex0, candidate1, distances);
      if (distance < distance1) {
        candidate0 = NULL;
        break;
      }
    }

    // Check if mutually closest 
    if (!candidate0) break;

    printf("HERE2\n"); fflush(stdout);

    // Found a mutually closest correspondence
    mutually_closest_candidates.vertices[0].Insert(candidate0);
    mutually_closest_candidates.vertices[1].Insert(candidate1);
    mesh0->SetVertexMark(candidate0, R3mesh_mark);
    mesh1->SetVertexMark(candidate1, R3mesh_mark);
    count++;

    // Print debug message
    if (print_debug) {
      printf("Found correspondence %d %d %g\n", mesh0->VertexID(candidate0), mesh1->VertexID(candidate1), distance1);
      fflush(stdout);
    }
  }

  // Add mutually closest candidates to sparse correspondence
  for (int i = 0; i < count; i++) {
    correspondence->vertices[0].Insert(mutually_closest_candidates.vertices[0].Kth(i));
    correspondence->vertices[1].Insert(mutually_closest_candidates.vertices[1].Kth(i));
  }

  // Print statistics
  if (print_verbose) {
    printf("  Expanded correspondences ...\n");
    printf("    Time = %.2f seconds\n", start_time.Elapsed());
    printf("    # Correspondences = %d\n", correspondence->vertices[0].NEntries());
    printf("    # Created = %d\n", count);
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int 
ExpandSparseCorrespondence(SparseVertexCorrespondence *correspondence, RNLength ***distances, int max_correspondences, RNLength min_spacing, RNScalar max_delta)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int niterations = 0;

  // Get convenient variables
  R3Mesh *mesh0 = correspondence->mesh[0];
  R3Mesh *mesh1 = correspondence->mesh[1];

  // Check/adjust maximum number of correspondences
  if (max_correspondences > mesh0->NVertices()) max_correspondences = mesh0->NVertices();
  if (max_correspondences > mesh1->NVertices()) max_correspondences = mesh1->NVertices();

  // Create importances
  RNScalar max_curvature = 10;
  RNLength **importances = new RNScalar * [ 2 ];
  for (int i = 0; i < 2; i++) {
    R3Mesh *mesh = correspondence->mesh[i];
    R3MeshProperty property(mesh);
    for (int j = 0; j < mesh->NVertices(); j++) { 
      R3MeshVertex *vertex = mesh->Vertex(j);
      RNScalar value = fabs(mesh->VertexMeanCurvature(vertex));
      if (value > max_curvature) value = max_curvature;
      property.SetVertexValue(j, value);
    }
    // property.Blur(0.2);
    importances[i] = new RNScalar [ mesh->NVertices() ];
    for (int j = 0; j < mesh->NVertices(); j++) { 
      importances[i][j] = property.VertexValue(j);
    }
  }

  if (print_debug) {
    printf("Created importances ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    fflush(stdout);
  }

  // Expand set of correspondences coarse to fine
  while (correspondence->vertices[0].NEntries() < max_correspondences) {
    int ncorrespondences = correspondence->vertices[0].NEntries();
    int max_candidates = max_correspondences - ncorrespondences;
    if (max_candidates > ncorrespondences) max_candidates = ncorrespondences;
    if (max_candidates == 0) break;

    // Create correspondences
    ExpandSparseCorrespondence(correspondence, distances, importances, max_candidates, min_spacing, max_delta);

    // Update distances
    RNBoolean updated = FALSE;
    for (int i = 0; i < 2; i++) {
      R3Mesh *mesh = correspondence->mesh[i];
      for (int j = 0; j < correspondence->vertices[i].NEntries(); j++) {
        R3MeshVertex *vertex = correspondence->vertices[i].Kth(j);
        if (distances[i][mesh->VertexID(vertex)] != NULL) continue;
        distances[i][mesh->VertexID(vertex)] = mesh->DijkstraDistances(vertex);
        updated = TRUE;
      }
    }
 
    // Check if updated anything
    if (!updated) break;

    // Update statistics
    niterations++;
  }

  // Delete importances
  for (int i = 0; i < 2; i++) {
    delete importances[i];
  }
  delete [] importances;

  // Print statistics
  if (print_verbose) {
    printf("Expanded sparse correspondence ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Correspondences = %d\n", correspondence->vertices[0].NEntries());
    printf("  # Iterations = %d\n", niterations);
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Dense Correspondence Creation
////////////////////////////////////////////////////////////////////////

static DenseVertexCorrespondence * 
CreateDenseVertexCorrespondence(SparseVertexCorrespondence *correspondence, RNLength ***distances, RNScalar max_delta)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Get convenient variables
  int nsparse = correspondence->vertices[0].NEntries();
  R3Mesh *mesh0 = correspondence->mesh[0];
  R3Mesh *mesh1 = correspondence->mesh[1];

  // Create map
  DenseVertexCorrespondence *map = new DenseVertexCorrespondence(mesh0, mesh1);
  if (!map) {
    fprintf(stderr, "Unable to allocate map\n");
    return NULL;
  }

  // Compute map
  // For each vertex in mesh0, find closest vertex in mesh1 in feature space
  int ndense = 0;
  for (int i = 0; i < mesh0->NVertices(); i++) {
    R3MeshVertex *vertex0 = mesh0->Vertex(i);

    // Compute closest vertex on mesh1 in feature space
    double min_delta = (max_delta > 0) ? nsparse * max_delta : FLT_MAX;
    for (int j = 0; j < mesh1->NVertices(); j++) {
      R3MeshVertex *vertex1 = mesh1->Vertex(j);
      double delta = ComputeDistanceDifference(correspondence, vertex0, vertex1, distances);
      if (delta < min_delta) {
        map->vertices[i] = vertex1;
        min_delta = delta;
      }

      // Update statistics
      if (map->vertices[i]) ndense++;
    }
  }

  // Print statistics
  if (print_verbose) {
    printf("Created map ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Mapped = %d\n", ndense);
    fflush(stdout);
  }

  // Return map
  return map;
}



////////////////////////////////////////////////////////////////////////
// Program argument parsing
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
      else if (!strcmp(*argv, "-add_swapped_correspondences")) add_swapped_correspondences = 1;
      else if (!strcmp(*argv, "-max_delta")) { argc--; argv++; max_delta = atof(*argv); }
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
    else {
      if (!input_mesh1_name) input_mesh1_name = *argv;
      else if (!input_mesh2_name) input_mesh2_name = *argv;
      else if (!input_correspondence_name) input_correspondence_name = *argv;
      else if (!output_map_name) output_map_name = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check input filename
  if (!input_mesh1_name || !input_mesh2_name || !input_correspondence_name || !output_map_name) {
    fprintf(stderr, "Usage: mshcorr input_mesh1 input_mesh2 input_correspondence output_map [options]\n");
    return 0;
  }

  // Return OK status 
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int 
main(int argc, char **argv)
{
  // Check number of arguments
  if (!ParseArgs(argc, argv)) exit(1);

  // Read meshes 
  R3Mesh *meshes[2] = { NULL, NULL };
  meshes[0] = ReadMesh(input_mesh1_name);
  if (!meshes[0]) exit(-1);
  meshes[1] = ReadMesh(input_mesh2_name);
  if (!meshes[1]) exit(-1);

  // Read sparse correspondence
  SparseVertexCorrespondence *correspondence = ReadSparseVertexCorrespondence(meshes[0], meshes[1], input_correspondence_name);
  if (!correspondence) exit(-1);

  // Create distances
  RNLength ***distances = CreateDistances(correspondence);
  if (!distances) exit(-1);

  // Expand sparse correspondence
  int status = ExpandSparseCorrespondence(correspondence, distances, max_sparse_correspondences, min_sparse_correspondence_spacing, max_delta);
  if (!status) exit(-1);

  // Create dense correspondence
  DenseVertexCorrespondence *map = CreateDenseVertexCorrespondence(correspondence, distances, max_delta);
  if (!map) exit(-1);

  // Write dense correspondence
  if (!WriteDenseVertexCorrespondence(map, output_map_name)) exit(-1);

  // Delete distances
  if (!DeleteDistances(correspondence, distances)) exit(-1);

  // Return success 
  return 0;
}

















