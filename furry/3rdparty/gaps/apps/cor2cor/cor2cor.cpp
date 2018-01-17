// Source file for the program to process sparse correspondences 



// Include files 

#include "R3Shapes/R3Shapes.h"



////////////////////////////////////////////////////////////////////////
// Program arguments
////////////////////////////////////////////////////////////////////////

static char *input_mesh1_name = NULL;
static char *input_mesh2_name = NULL;
static char *input_correspondence_name = NULL;
static char *output_correspondence_name = NULL;
static int correspondences_trace_symmetry_axis = 0;
static int interpolation_min_correspondences = 0;
static int interpolation_max_correspondences = 0;
static RNLength interpolation_min_spacing = 0;
static RNLength interpolation_max_spacing = 0;
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
  int nvertices;
  SparseVertexCorrespondence(R3Mesh *mesh0, R3Mesh *mesh1) { 
    mesh[0] = mesh0; 
    mesh[1] = mesh1; 
    nvertices = 0;
  };
};



////////////////////////////////////////////////////////////////////////
// Input/output
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
      correspondence->nvertices++;
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
    while (fscanf(fp, "%d%lg%d%lg", &id0, &fdummy, &id1, &fdummy) == (unsigned int) 4) {
      R3MeshVertex *vertex0 = mesh0->Vertex(id0);
      R3MeshVertex *vertex1 = mesh1->Vertex(id1);
      correspondence->vertices[0].Insert(vertex0);
      correspondence->vertices[1].Insert(vertex1);
      correspondence->nvertices++;
    }

    // Check number of correspondences
    if (correspondence->nvertices != ncorrespondences) {
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
    printf("  # Correspondences = %d\n", correspondence->nvertices);
    fflush(stdout);
  }

  // Return correspondence
  return correspondence;
}



static int
WriteSparseVertexCorrespondence(SparseVertexCorrespondence *correspondence, char *filename)
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

  // Write corresponding vertex IDs
  for (int i = 0; i < correspondence->nvertices; i++) {
    R3MeshVertex *v0 = correspondence->vertices[0].Kth(i);
    R3MeshVertex *v1 = correspondence->vertices[1].Kth(i);
    int id0 = correspondence->mesh[0]->VertexID(v0);
    int id1 = correspondence->mesh[1]->VertexID(v1);
    fprintf(fp, "%6d %6d\n", id0, id1);
  }

  // Close file
  fclose(fp);

  // Print statistics
  if (print_verbose) {
    printf("Wrote correspondences to %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Correspondences = %d\n", correspondence->nvertices);
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Interpolation
////////////////////////////////////////////////////////////////////////

static RNLength 
ComputePropertyDifference0(R3MeshPropertySet **properties, R3MeshVertex *vertex0, R3MeshVertex *vertex1)
{
  // Initialize property
  RNLength delta = 0;
  
  // Compute property in feature space
  for (int i = 0; i < properties[0]->NProperties(); i++) {
    R3MeshProperty *property0 = properties[0]->Property(i);
    R3MeshProperty *property1 = properties[1]->Property(i);
    double value0 = property0->VertexValue(vertex0);
    double value1 = property1->VertexValue(vertex1);
    double d = fabs(value1 - value0);
    if (value0 > 0) delta += d / value0;
    // else if (d > 0) delta += RN_INFINITY; 
  }

  // Return L1 difference of properties
  return delta;
}



static RNLength
ComputePropertyDifference1(R3MeshPropertySet **properties, R3MeshVertex *vertex0, R3MeshVertex *vertex1)
{
  // Initialize property
  RNLength delta = 0;

  // Compute property in feature space
  for (int i = 0; i < properties[0]->NProperties(); i++) {
    R3MeshProperty *property0 = properties[0]->Property(i);
    R3MeshProperty *property1 = properties[1]->Property(i);
    double value0 = property0->VertexValue(vertex0);
    double value1 = property1->VertexValue(vertex1);
    double d = fabs(value1 - value0);
    if (value1 > 0) delta += d / value1;
    // else if (d > 0) delta += RN_INFINITY;
  }

  // Return L1 difference of properties
  return delta;
}



static RNScalar *
CreateSides(R3Mesh *mesh, RNArray<R3MeshVertex *>& correspondence_vertices)
{
  // Allocate/initialize sides
  RNScalar *sides = new RNScalar [ mesh->NVertices() ];
  for (int i = 0; i < mesh->NVertices(); i++) sides[i] = 0;

  // Compute axis
  R3mesh_mark++;
  RNArray<R3MeshEdge *> axis_edges;
  RNArray<R3MeshVertex *> axis_vertices;
  for (int i = 0; i < correspondence_vertices.NEntries(); i++) {
    RNArray<R3MeshEdge *> shortest_path_edges;
    R3MeshVertex *vertex0 = correspondence_vertices.Kth(i);
    R3MeshVertex *vertex1 = correspondence_vertices.Kth((i+1) % correspondence_vertices.NEntries());
    mesh->DijkstraDistance(vertex1, vertex0, &shortest_path_edges);
    R3MeshVertex *vertex = vertex0;
    for (int j = 0; j < shortest_path_edges.NEntries(); j++) {
      R3MeshEdge *edge = shortest_path_edges.Kth(j);
      assert(mesh->IsVertexOnEdge(vertex, edge));
      axis_edges.Insert(edge);
      axis_vertices.Insert(vertex);
      vertex = mesh->VertexAcrossEdge(edge, vertex);
      mesh->SetVertexMark(vertex, R3mesh_mark);
    }
  }

  // Find seed vertex on ccw side
  R3MeshVertex *ccw_seed = NULL;
  for (int j = 0; j < axis_edges.NEntries(); j++) {
    R3MeshVertex *vertex = axis_vertices.Kth(j);
    R3MeshEdge *edge = axis_edges.Kth(j);
    assert(mesh->IsVertexOnEdge(vertex, edge));
    R3MeshEdge *ccw_edge = mesh->EdgeOnVertex(vertex, edge, RN_CCW);
    R3MeshVertex *ccw_vertex = mesh->VertexAcrossEdge(ccw_edge, vertex);
    if (mesh->VertexMark(ccw_vertex) == R3mesh_mark) continue;
    ccw_seed = ccw_vertex;
    break;
  }

  // Check if found seed
  if (!ccw_seed) return sides;

  // Compute distances from axis 
  RNLength *distances = mesh->DijkstraDistances(axis_vertices);

  // Flood fill ccw side from seed
  RNArray<R3MeshVertex *> stack;
  stack.Insert(ccw_seed);
  mesh->SetVertexMark(ccw_seed, R3mesh_mark);
  while (!stack.IsEmpty()) {
    R3MeshVertex *vertex = stack.Tail();
    stack.RemoveTail();
    sides[mesh->VertexID(vertex)] = distances[mesh->VertexID(vertex)];
    for (int i = 0; i < mesh->VertexValence(vertex); i++) {
      R3MeshEdge *edge = mesh->EdgeOnVertex(vertex, i);
      R3MeshVertex *neighbor = mesh->VertexAcrossEdge(edge, vertex);
      if (mesh->VertexMark(neighbor) != R3mesh_mark) {
        mesh->SetVertexMark(neighbor, R3mesh_mark);
        stack.Insert(neighbor);
      }
    }
  }

  // Flood fill cw side from seed
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(i);
    if (mesh->VertexMark(vertex) == R3mesh_mark) continue;
    sides[mesh->VertexID(vertex)] = -distances[mesh->VertexID(vertex)];
  }

  // Delete distances
  delete [] distances;

  // Return sides
  return sides;
}



static int
CreateProperties(R3MeshPropertySet **properties, SparseVertexCorrespondence *correspondence, int previous_ncorrespondences)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_debug) {
    printf("  Updating properties ...\n");
    fflush(stdout);
  }

  // Create property set
  for (int i = 0; i < 2; i++) {
    R3Mesh *mesh = correspondence->mesh[i];

    // Compute distance properties
    for (int j = previous_ncorrespondences; j < correspondence->nvertices; j++) {
      R3MeshVertex *vertex = correspondence->vertices[i].Kth(j);
      char name[256];
      sprintf(name, "d%d\n", mesh->VertexID(vertex));
      RNLength *distances = mesh->DijkstraDistances(vertex);
      R3MeshProperty *property = new R3MeshProperty(mesh, name, distances);
      properties[i]->Insert(property);
      delete [] distances;
    }

    // Compute side properties
    if (correspondences_trace_symmetry_axis) {
      if (previous_ncorrespondences == 0) {
        RNScalar *sides = CreateSides(mesh, correspondence->vertices[i]);
        R3MeshProperty *property = new R3MeshProperty(mesh, "Side", sides);
        property->Multiply(1);
        properties[i]->Insert(property);
        delete [] sides;
      }
    }
  }

  // Debug information
  if (print_debug) {
    static int counter = 1;
    char buffer[256];
    sprintf(buffer, "properties0%d.arff", counter);
    properties[0]->Write(buffer);
    sprintf(buffer, "properties1%d.arff", counter);
    properties[1]->Write(buffer);
    counter++;
  }

  // Print messhage
  if (print_debug) {
    printf("    Time = %.2f seconds\n", start_time.Elapsed());
    printf("    # Properties = %d\n", properties[0]->NProperties());
    fflush(stdout);
  }

  // Return success
  return 1;
}


static R3MeshProperty *
CreateImportances(R3Mesh *mesh, RNArray<R3MeshVertex *>& correspondence_vertices)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_debug) {
    printf("    Creating importances ...\n");
    fflush(stdout);
  }

  // Allocate importances
  R3MeshProperty *importances = new R3MeshProperty(mesh, "Importance");
  if (!importances) {
    fprintf(stderr, "Unable to allocate importances\n");
    return NULL;
  }

  // Compute curvatures
  R3MeshProperty curvatures(mesh);
  for (int i = 0; i < mesh->NVertices(); i++) { 
    R3MeshVertex *vertex = mesh->Vertex(i);
    RNScalar value = mesh->VertexGaussCurvature(vertex);
    curvatures.SetVertexValue(i, value);
  }

  // Normalize curvatures
  RNScalar min_curvature = curvatures.Minimum();
  RNScalar max_curvature = curvatures.Maximum();
  if (min_curvature < 0) min_curvature = 0;
  if (max_curvature > 1000) max_curvature = 1000;
  if (max_curvature <= min_curvature) max_curvature = min_curvature + 1;
  for (int i = 0; i < mesh->NVertices(); i++) { 
    RNScalar value = curvatures.VertexValue(i);
    if (value < min_curvature) value = min_curvature;
    if (value > max_curvature) value = max_curvature;
    value = (value - min_curvature) / (max_curvature - min_curvature);
    curvatures.SetVertexValue(i, value);
  }

  // Blur curvatures
  RNScalar spacing = sqrt( 1.0 / (correspondence_vertices.NEntries() * RN_PI));
  RNScalar sigma = 0.2 * spacing;
  curvatures.Blur(sigma);

  // Compute importances
  for (int i = 0; i < mesh->NVertices(); i++) { 
    RNScalar value = curvatures.VertexValue(i);
    if (curvatures.IsLocalMaximum(i)) value *= 2;
    importances->SetVertexValue(i, value);
  }

  // Debugging
  if (print_debug) {
    static int iteration = 1;
    char buffer[256];
    sprintf(buffer, "importance%d.val", iteration++);
    importances->Write(buffer);
  }

  // Print messhage
  if (print_debug) {
    printf("      Time = %.2f seconds\n", start_time.Elapsed());
    printf("      # Importances = %d\n", mesh->NVertices());
    printf("      # Correspondences = %d\n", correspondence_vertices.NEntries());
    printf("      Spacing = %g\n", spacing);
    printf("      Sigma = %g\n", sigma);
    fflush(stdout);
  }

  // Return importances
  return importances;
}



static RNArray<R3MeshVertex *> *
CreateCandidates(R3Mesh *mesh, RNArray<R3MeshVertex *>& correspondence_vertices, int max_candidates)
{
  // Create importances
  R3MeshProperty *importances = CreateImportances(mesh, correspondence_vertices);
  if (!importances) return NULL;

  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_debug) {
    printf("    Creating candidates ...\n");
    fflush(stdout);
  }

 // Allocate array of candidates
  RNArray<R3MeshVertex *> *candidates = new RNArray<R3MeshVertex *>();
  if (!candidates) {
    fprintf(stderr, "Unable to allocate candidates\n");
    return NULL;
  }

  // Find candidates one at a time
  RNArray<R3MeshVertex *> current_vertices(correspondence_vertices);
  for (int i = 0; i < max_candidates; i++) {
    // Compute distances from current vertices
    RNLength *dists = mesh->DijkstraDistances(current_vertices);

    // Find candidate vertex with highest value (importance * distance)
    R3MeshVertex *candidate = NULL;
    RNScalar best_value = -FLT_MAX;
    for (int i = 0; i < mesh->NVertices(); i++) {
      if (dists[i] == 0) continue;
      RNScalar value = (1 + importances->VertexValue(i)) * dists[i];
      if (value > best_value) {
        candidate = mesh->Vertex(i);
        best_value = value;
      }
    }

    // Delete distances from current vertices
    delete [] dists;

    // Check if found candidate vertex
    if (!candidate) break;

    // Add candidate
    candidates->Insert(candidate);
    current_vertices.Insert(candidate);
  }

  // Debugging
  if (print_debug) {
    static int iteration = 1;
    char buffer[256];
    sprintf(buffer, "candidate%d.pid", iteration++);
    FILE *fp = fopen(buffer, "w");
    for (int i = 0; i < candidates->NEntries(); i++) 
      fprintf(fp, "%d\n", mesh->VertexID(candidates->Kth(i)));
    fclose(fp);
  }

  // Print statistics
  if (print_debug) {
    printf("      Time = %.2f seconds\n", start_time.Elapsed());
    printf("      # Candidates = %d\n", candidates->NEntries());
    fflush(stdout);
  }

  // Delete importances
  delete importances;

  // Return candidates
  return candidates;
}



static int 
InterpolateCorrespondence(SparseVertexCorrespondence *correspondence, 
  int min_correspondences, int max_correspondences, 
  RNLength min_spacing, RNLength max_spacing, 
  R3MeshPropertySet **properties, int max_candidates)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  R3mesh_mark++;
  int count = 0;
  if (print_debug) {
    static int iteration = 1;
    printf("  Iteration %d ...\n", iteration++);
    fflush(stdout);
  }

  // Get convenient variables
  R3Mesh *mesh0 = correspondence->mesh[0];
  R3Mesh *mesh1 = correspondence->mesh[1];

  // Create candidates
  RNArray<R3MeshVertex *> *candidates0 = CreateCandidates(mesh0, correspondence->vertices[0], max_candidates);
  if (!candidates0 || candidates0->IsEmpty()) return 0; 
  RNArray<R3MeshVertex *> *candidates1 = CreateCandidates(mesh1, correspondence->vertices[1], max_candidates);
  if (!candidates1 || candidates1->IsEmpty()) return 0; 

  // Create struct to store new correspondences
  SparseVertexCorrespondence mutually_closest_candidates(mesh0, mesh1);

  // Consider each candidate in candidates0
  for (int i = 0; i < candidates0->NEntries(); i++) {
    R3MeshVertex *candidate0 = candidates0->Kth(i);
    if (mesh0->VertexMark(candidate0) == R3mesh_mark) continue;

    // Find closest candidate in candidates1 
    RNScalar distance1 = FLT_MAX;
    R3MeshVertex *candidate1 = NULL;
    for (int j = 0; j < candidates1->NEntries(); j++) {
      R3MeshVertex *vertex1 = candidates1->Kth(j);
      if (mesh1->VertexMark(vertex1) == R3mesh_mark) continue;
      RNLength distance = ComputePropertyDifference0(properties, candidate0, vertex1);
      if (distance < distance1) {
        candidate1 = vertex1;
        distance1 = distance;
      }
    }

    // Check if found candidate1
    if (!candidate1) continue;

    // Check if there is a candidate0 closer to candidate1
    RNScalar distance2 = FLT_MAX;
    R3MeshVertex *candidate2 = NULL;
    for (int j = 0; j < candidates0->NEntries(); j++) {
      R3MeshVertex *vertex0 = candidates0->Kth(j);
      if (mesh0->VertexMark(candidate0) == R3mesh_mark) continue;
      RNLength distance = ComputePropertyDifference1(properties, vertex0, candidate1);
      if (distance < distance2) {
        candidate2 = vertex0;
        distance2 = distance;
      }
    }

    // Check if mutually closest 
    if (candidate0 != candidate2) {
      if (print_debug) {
        RNScalar side0 = properties[0]->Property(properties[0]->NProperties()-1)->VertexValue(candidate0);
        RNScalar side1 = properties[1]->Property(properties[1]->NProperties()-1)->VertexValue(candidate1);
        RNScalar side_fraction = (distance1 > 0) ? fabs(side0 - side1) / distance1 : 0;
        printf("    M %d %d : %g (%g) : %d %g\n", mesh0->VertexID(candidate0), mesh1->VertexID(candidate1), distance1, side_fraction,
              mesh0->VertexID(candidate2), distance2);
        fflush(stdout);
      }
      continue;
    }

    // Found a mutually closest correspondence
    mutually_closest_candidates.vertices[0].Insert(candidate0);
    mutually_closest_candidates.vertices[1].Insert(candidate1);
    mutually_closest_candidates.nvertices++;
    mesh0->SetVertexMark(candidate0, R3mesh_mark);
    mesh1->SetVertexMark(candidate1, R3mesh_mark);
    count++;

    // Print debug message
    if (print_debug) {
      RNScalar side0 = properties[0]->Property(properties[0]->NProperties()-1)->VertexValue(candidate0);
      RNScalar side1 = properties[1]->Property(properties[1]->NProperties()-1)->VertexValue(candidate1);
      RNScalar side_fraction = (distance1 > 0) ? fabs(side0 - side1) / distance1 : 0;
      printf("    C %d %d : %g (%g)\n", mesh0->VertexID(candidate0), mesh1->VertexID(candidate1), distance1, side_fraction);
      fflush(stdout);
    }
  }

  // Add mutually closest candidates to sparse correspondence
  for (int i = 0; i < count; i++) {
    correspondence->vertices[0].Insert(mutually_closest_candidates.vertices[0].Kth(i));
    correspondence->vertices[1].Insert(mutually_closest_candidates.vertices[1].Kth(i));
    correspondence->nvertices++;
  }

  // Print statistics
  if (print_debug) {
    printf("    Time = %.2f seconds\n", start_time.Elapsed());
    printf("    # Correspondences = %d\n", correspondence->nvertices);
    printf("    # New Correspondences = %d\n", count);
    printf("    # Candidates = %d %d\n", candidates0->NEntries(), candidates1->NEntries());
    fflush(stdout);
  }

  // Delete the candidates
  delete candidates0;
  delete candidates1;

  // Return success
  return 1;
}



static int 
InterpolateCorrespondence(SparseVertexCorrespondence *correspondence, 
  int min_correspondences, int max_correspondences, 
  RNLength min_spacing, RNLength max_spacing)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int niterations = 0;
  if (print_verbose) {
    printf("Interpolating correspondences ...\n");
    fflush(stdout);
  }

  // Get convenient variables
  R3Mesh *mesh0 = correspondence->mesh[0];
  R3Mesh *mesh1 = correspondence->mesh[1];

  // Check/adjust number of correspondences
  if (min_correspondences > mesh0->NVertices()) min_correspondences = mesh0->NVertices();
  if (min_correspondences > mesh1->NVertices()) min_correspondences = mesh1->NVertices();
  if (min_correspondences <= correspondence->nvertices) return 1;

  // Create properties
  R3MeshPropertySet **properties = new R3MeshPropertySet *[ 2 ];
  properties[0] = new R3MeshPropertySet(mesh0);
  properties[1] = new R3MeshPropertySet(mesh1);
  CreateProperties(properties, correspondence, 0);

  // Interpolate correspondences a bit at a time (coarse to fine)
  while (correspondence->nvertices < min_correspondences) {
    // Remember old number of correspondences
    int ncorrespondences = correspondence->nvertices;

    // Create correspondences
    if (!InterpolateCorrespondence(correspondence, 
       min_correspondences, max_correspondences, 
       min_spacing, max_spacing, 
       properties, 2 * ncorrespondences)) return 0;

    // Check if created any new correspondences
    if (correspondence->nvertices == ncorrespondences) break;

    // Update properties
    if (correspondence->nvertices < min_correspondences) {
      CreateProperties(properties, correspondence, ncorrespondences);
    }

    // Update statistics
    niterations++;
  }

  // Print statistics
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Iterations = %d\n", niterations);
    printf("  # Correspondences = %d\n", correspondence->nvertices);
    printf("  Min Correspondences = %d\n", min_correspondences);
    printf("  Min Spacing = %g\n", min_spacing);
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Parse program arguments
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
      else if (!strcmp(*argv, "-symmetry_axis")) correspondences_trace_symmetry_axis = 1; 
      else if (!strcmp(*argv, "-interpolate")) { argc--; argv++; interpolation_min_correspondences = atoi(*argv); }
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
    else {
      if (!input_mesh1_name) input_mesh1_name = *argv;
      else if (!input_mesh2_name) input_mesh2_name = *argv;
      else if (!input_correspondence_name) input_correspondence_name = *argv;
      else if (!output_correspondence_name) output_correspondence_name = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check input filename
  if (!input_mesh1_name || !input_mesh2_name || !input_correspondence_name || !output_correspondence_name) {
    fprintf(stderr, "Usage: mshcorr input_mesh1 input_mesh2 input_correspondence output_correspondence [options]\n");
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

  // Interpolate correspondences
  if (interpolation_min_correspondences) {
    if (!InterpolateCorrespondence(correspondence, 
      interpolation_min_correspondences, interpolation_max_correspondences, 
      interpolation_min_spacing, interpolation_max_spacing)) {
      exit(-1);
    }
  }

  // Write sparse correspondence
  int status = WriteSparseVertexCorrespondence(correspondence, output_correspondence_name);
  if (!status) exit(-1);

  // Return success 
  return 0;
}




