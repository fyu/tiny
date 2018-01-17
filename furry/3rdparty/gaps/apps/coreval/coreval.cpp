// Source file to print information about a mesh correspondence



// Include files 

#include "R3Shapes/R3Shapes.h"



// Program arguments

static char *mesh1_name = NULL;
static char *mesh2_name = NULL;
static char *predicted_correspondence_name = NULL;
static char *truth_correspondence_name = NULL;
static double normalized_distance_threshold = sqrt(1.0/(20.0*RN_PI));
static int print_verbose = 0;
static int print_debug = 0;



// Type definitions

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
    fprintf(stderr, "Unable to open vertex file %s\n", filename);
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



static DenseVertexCorrespondence *
ReadDenseVertexCorrespondence(R3Mesh *mesh0, R3Mesh *mesh1, char *filename)
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
    fprintf(stderr, "Unable to open vertex file %s\n", filename);
    return NULL;
  }

  // Create correspondence
  DenseVertexCorrespondence *correspondence = new DenseVertexCorrespondence(mesh0, mesh1);
  if (!correspondence) {
    fprintf(stderr, "Unable to allocate correspondence for %s\n", filename); 
    return NULL; 
  }

  // Check filename extension
  if (!strcmp(extension, ".map")) {
    // Read correspondences
    for (int i = 0; i < mesh0->NVertices(); i++) {
      int id1;
      if (fscanf(fp, "%d", &id1) != (unsigned int) 1) {
        fprintf(stderr, "Unable to read correspondence %d from %s\n", i, filename);
        return NULL;
      }
      R3MeshVertex *vertex1 = mesh1->Vertex(id1);
      correspondence->vertices[i] = vertex1;
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
    printf("Read dense correspondences from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Correspondences = %d\n", mesh0->NVertices());
    fflush(stdout);
  }

  // Return correspondence
  return correspondence;
}



static double
EvaluateSparseVertexCorrespondence(SparseVertexCorrespondence *predicted, DenseVertexCorrespondence *truth)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Get useful variables
  R3Mesh *mesh0 = predicted->mesh[0];
  R3Mesh *mesh1 = predicted->mesh[1];
  double normalization = sqrt(mesh1->Area());
  if (normalization <= 0) return RN_INFINITY;
  else normalization = 1.0 / normalization;

  // Consider each vertex of sparse correspondence
  int total_count = 0;
  int correct_count = 0;
  double normalized_sum = 0;
  for (int i = 0; i < predicted->vertices[0].NEntries(); i++) {
    R3MeshVertex *vertex0 = predicted->vertices[0].Kth(i);
    R3MeshVertex *predicted_vertex1 = predicted->vertices[1].Kth(i);
    R3MeshVertex *truth_vertex1 = truth->vertices[ mesh0->VertexID(vertex0) ];
    RNScalar distance = (truth_vertex1) ? mesh1->DijkstraDistance(predicted_vertex1, truth_vertex1) : RN_INFINITY;
    RNScalar normalized_distance = (truth_vertex1) ? normalization * distance : 10;
    normalized_sum += normalized_distance * normalized_distance;
    if (normalized_distance < normalized_distance_threshold) correct_count++;
    if (print_debug) printf("C %d %g\n", i, normalized_distance);
    total_count++;
  }

  // Compute average distance normalized by mesh size
  double normalized_rmsd = (total_count > 0) ? sqrt(normalized_sum / total_count) : 0;
  double normalized_classification_rate = (total_count > 0) ? (double) correct_count / (double) total_count : -1;

  // Print statistics
  if (print_verbose) {
    printf("Evaluated sparse correspondences ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Correspondences = %d\n", total_count);
    printf("  # Correct = %d\n", correct_count);
    printf("  Normalization factor = %g\n", normalization);
    printf("  Normalized RMSD = %g\n", normalized_rmsd);
    printf("  Normalized classification rate = %g\n", normalized_classification_rate);
    printf("  Normalized distance threshold = %g\n", normalized_distance_threshold);
    fflush(stdout);
  }

  // Return evaluation metric
  return normalized_rmsd;
}



static double
EvaluateDenseVertexCorrespondence(DenseVertexCorrespondence *predicted, DenseVertexCorrespondence *truth)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Get useful variables
  R3Mesh *mesh0 = predicted->mesh[0];
  R3Mesh *mesh1 = predicted->mesh[1];
  double normalization = sqrt(mesh1->Area());
  if (normalization <= 0) return RN_INFINITY;
  else normalization = 1.0 / normalization;

  // Consider each vertex of dense correspondence
  int total_count = 0;
  int correct_count = 0;
  double normalized_sum = 0;
  for (int i = 0; i < mesh0->NVertices(); i++) {
    R3MeshVertex *predicted_vertex1 = predicted->vertices[i];
    if (!predicted_vertex1) continue;
    R3MeshVertex *truth_vertex1 = truth->vertices[i];
    RNScalar distance = (truth_vertex1) ? mesh1->DijkstraDistance(predicted_vertex1, truth_vertex1) : RN_INFINITY;
    RNScalar normalized_distance = (truth_vertex1) ? normalization * distance : 10;
    normalized_sum += normalized_distance * normalized_distance;
    if (normalized_distance < normalized_distance_threshold) correct_count++;
    if (print_debug) printf("C %d %g\n", i, normalized_distance);
    total_count++;
  }

  // Compute average distance normalized by mesh size
  double normalized_rmsd = (total_count > 0) ? sqrt(normalized_sum) / total_count : 0;
  double normalized_classification_rate = (total_count > 0) ? (double) correct_count / (double) total_count : -1;

  // Print statistics
  if (print_verbose) {
    printf("Evaluated dense correspondences ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Correspondences = %d\n", total_count);
    printf("  # Correct = %d\n", correct_count);
    printf("  Normalization factor = %g\n", normalization);
    printf("  Normalized RMSD = %g\n", normalized_rmsd);
    printf("  Normalized classification rate = %g\n", normalized_classification_rate);
    printf("  Normalized distance threshold = %g\n", normalized_distance_threshold);
    fflush(stdout);
  }

  // Return evaluation metric
  return normalized_rmsd;
}



static double
EvaluateDenseVertexCorrespondence(DenseVertexCorrespondence *predicted, SparseVertexCorrespondence *truth)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Get useful variables
  R3Mesh *mesh0 = predicted->mesh[0];
  R3Mesh *mesh1 = predicted->mesh[1];
  double normalization = sqrt(mesh1->Area());
  if (normalization <= 0) return RN_INFINITY;
  else normalization = 1.0 / normalization;

  // Consider each vertex of sparse truth correspondence
  int total_count = 0;
  int correct_count = 0;
  double normalized_sum = 0;
  for (int i = 0; i < truth->vertices[0].NEntries(); i++) {
    R3MeshVertex *truth_vertex0 = truth->vertices[0].Kth(i);
    R3MeshVertex *truth_vertex1 = truth->vertices[1].Kth(i);
    R3MeshVertex *predicted_vertex1 = predicted->vertices[ mesh0->VertexID(truth_vertex0) ];
    RNScalar distance = (predicted_vertex1) ? mesh1->DijkstraDistance(truth_vertex1, predicted_vertex1) : RN_INFINITY;
    RNScalar normalized_distance = (predicted_vertex1) ? normalization * distance : 10;
    normalized_sum += normalized_distance * normalized_distance;
    if (normalized_distance < normalized_distance_threshold) correct_count++;
    if (print_debug) printf("C %d %g\n", i, normalized_distance);
    total_count++;
  }

  // Compute average distance normalized by mesh size
  double normalized_rmsd = (total_count > 0) ? sqrt(normalized_sum) / total_count : 0;
  double normalized_classification_rate = (total_count > 0) ? (double) correct_count / (double) total_count : -1;

  // Print statistics
  if (print_verbose) {
    printf("Evaluated dense correspondences ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Correspondences = %d\n", total_count);
    printf("  # Correct = %d\n", correct_count);
    printf("  Normalization factor = %g\n", normalization);
    printf("  Normalized RMSD = %g\n", normalized_rmsd);
    printf("  Normalized classification rate = %g\n", normalized_classification_rate);
    printf("  Normalized distance threshold = %g\n", normalized_distance_threshold);
    fflush(stdout);
  }

  // Return evaluation metric
  return normalized_rmsd;
}



static int 
ParseArgs(int argc, char **argv)
{
  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) print_verbose = 1;
      else if (!strcmp(*argv, "-debug")) print_debug = 1;
      else if (!strcmp(*argv, "-normalized_distance_threshold")) { argv++; argc--; normalized_distance_threshold = atof(*argv); }
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
    else {
      if (!mesh1_name) mesh1_name = *argv;
      else if (!mesh2_name) mesh2_name = *argv;
      else if (!predicted_correspondence_name) predicted_correspondence_name = *argv;
      else if (!truth_correspondence_name) truth_correspondence_name = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check input filename
  if (!mesh1_name || !mesh2_name || !predicted_correspondence_name || !truth_correspondence_name) {
    fprintf(stderr, "Usage: coreval mesh1 mesh2 predicted_correspondence truth_dense_correspondence [options]\n");
    return 0;
  }

  // Return OK status 
  return 1;
}



int 
main(int argc, char **argv)
{
  // Check number of arguments
  if (!ParseArgs(argc, argv)) exit(1);

  // Read mesh1 
  R3Mesh *mesh1 = ReadMesh(mesh1_name);
  if (!mesh1) exit(-1);

  // Read mesh2
  R3Mesh *mesh2 = ReadMesh(mesh2_name);
  if (!mesh2) exit(-1);

  // Check predicted correspondence type 
  RNScalar value = 0;
  if (strstr(predicted_correspondence_name, ".map")) {
    // Read predicted DENSE correspondence
    DenseVertexCorrespondence *predicted_correspondence = ReadDenseVertexCorrespondence(mesh1, mesh2, predicted_correspondence_name);
    if (!predicted_correspondence) exit(-1);

    if (strstr(truth_correspondence_name, ".map")) {
      // Read truth dense correspondence
      DenseVertexCorrespondence *truth_correspondence = ReadDenseVertexCorrespondence(mesh1, mesh2, truth_correspondence_name);
      if (!truth_correspondence) exit(-1);

      // Evaluate predicted correspondences with respect to truth correspondence
      value = EvaluateDenseVertexCorrespondence(predicted_correspondence, truth_correspondence);
    }
    else {
      // Read truth sparse correspondence
      SparseVertexCorrespondence *truth_correspondence = ReadSparseVertexCorrespondence(mesh1, mesh2, truth_correspondence_name);
      if (!truth_correspondence) exit(-1);

      // Evaluate predicted correspondences with respect to truth correspondence
      value = EvaluateDenseVertexCorrespondence(predicted_correspondence, truth_correspondence);
    }
  }
  else {
    // Read predicted SPARSE correspondence
    SparseVertexCorrespondence *predicted_correspondence = ReadSparseVertexCorrespondence(mesh1, mesh2, predicted_correspondence_name);
    if (!predicted_correspondence) exit(-1);

    if (strstr(truth_correspondence_name, ".map")) {
      // Read truth dense correspondence
      DenseVertexCorrespondence *truth_correspondence = ReadDenseVertexCorrespondence(mesh1, mesh2, truth_correspondence_name);
      if (!truth_correspondence) exit(-1);

      // Evaluate predicted correspondences with respect to truth correspondence
      value = EvaluateSparseVertexCorrespondence(predicted_correspondence, truth_correspondence);
    }
    else {
      fprintf(stderr, "Cannot evaluate sparse correspondences with sparse truth file\n");
      return 0;
    }
  }

  // Print value
  if (!print_verbose) printf("%g\n", value);

  // Return success 
  return 0;
}

















