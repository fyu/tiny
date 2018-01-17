// Source file for mesh processing utilities



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV/GSV.h"
#include "process.h"



////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////

extern int print_debug;



////////////////////////////////////////////////////////////////////////
// Functions
////////////////////////////////////////////////////////////////////////

static RNScalar 
EdgeLengthCallback(GSVMeshEdge *edge, void *data)
{
  // Return value associated with edge for heap sorting
  const GSVMesh *mesh = (GSVMesh *) data;
  return mesh->EdgeLength(edge);
}



int SplitEdges(GSVMesh *mesh, RNLength max_edge_length, int selected_segment)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int count = 0;

  // Split long edges 
  if (max_edge_length <= 0) return 1;

  // Create priority queue for long edges between two selected faces
  RNHeap<GSVMeshEdge *> heap(EdgeLengthCallback, NULL, mesh, FALSE);
  for (int i = 0; i < mesh->NEdges(); i++) {
    GSVMeshEdge *edge = mesh->Edge(i);
    GSVMeshFace *face0 = mesh->FaceOnEdge(edge, 0);
    GSVMeshFace *face1 = mesh->FaceOnEdge(edge, 1);
    if (!face0 || !face1) continue;
    if ((selected_segment >= 0) && (mesh->FaceSegment(face0) != selected_segment)) continue;
    if ((selected_segment >= 0) && (mesh->FaceSegment(face1) != selected_segment)) continue;
    RNLength length = mesh->EdgeLength(edge);
    if (length <= max_edge_length) continue;
    heap.Push(edge);
  }

  // Subdivide edges in longest to shortest order
  while (!heap.IsEmpty()) {
    GSVMeshEdge *edge = heap.Pop();
    RNLength old_length = mesh->EdgeLength(edge);
    GSVMeshVertex *v0 = mesh->VertexOnEdge(edge, 0);
    GSVMeshVertex *v1 = mesh->VertexOnEdge(edge, 1);
    GSVMeshFace *face0 = mesh->FaceOnEdge(edge, 0);
    GSVMeshFace *face1 = mesh->FaceOnEdge(edge, 1);
    if (!face0 || !face1) continue;
    GSVMeshVertex *vA = mesh->VertexAcrossFace(face0, edge);
    GSVMeshVertex *vB = mesh->VertexAcrossFace(face0, edge);
    R3Point midpoint = 0.5 * (mesh->VertexPosition(v0) + mesh->VertexPosition(v1));
    if (R3Contains(midpoint, mesh->VertexPosition(vA))) continue;
    if (R3Contains(midpoint, mesh->VertexPosition(vB))) continue;
    GSVMeshVertex *vertex = mesh->SplitEdge(edge, midpoint);
    for (int i = 0; i < mesh->VertexValence(vertex); i++) {
      GSVMeshEdge *edge = mesh->EdgeOnVertex(vertex, i);
      GSVMeshFace *face0 = mesh->FaceOnEdge(edge, 0);
      GSVMeshFace *face1 = mesh->FaceOnEdge(edge, 1);
      if (!face0 || !face1) continue;
      if ((selected_segment >= 0) && (mesh->FaceSegment(face0) != selected_segment)) continue;
      if ((selected_segment >= 0) && (mesh->FaceSegment(face1) != selected_segment)) continue;
      RNLength length = mesh->EdgeLength(edge);
      if (length <= max_edge_length) continue;
      if (RNIsGreaterOrEqual(length, old_length)) continue;
      heap.Push(edge);
    }
    count++;
  }

  // Print message
  if (print_debug) {
    printf("    Split edges\n");
    printf("      Time = %.2f seconds\n", start_time.Elapsed());
    printf("      Maximum Edge Length = %g\n", max_edge_length);
    printf("      # Splits = %d\n", count);
    fflush(stdout);
  }

  // Return success
  return 1;
}



int CollapseEdges(GSVMesh *mesh, RNLength min_edge_length, int selected_segment)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int count = 0;

  printf("%d\n", selected_segment);

  // Collapse short edges
  if (min_edge_length <= 0) return 1;

  // Create priority queue for short edges
  R3mesh_mark++;
  RNHeap<GSVMeshEdge *> heap(EdgeLengthCallback, NULL, mesh, TRUE);
  for (int i = 0; i < mesh->NEdges(); i++) {
    GSVMeshEdge *edge = mesh->Edge(i);
    GSVMeshFace *face0 = mesh->FaceOnEdge(edge, 0);
    GSVMeshFace *face1 = mesh->FaceOnEdge(edge, 1);
    if (!face0 || !face1) continue;
    if ((selected_segment >= 0) && (mesh->FaceSegment(face0) != selected_segment)) continue;
    if ((selected_segment >= 0) && (mesh->FaceSegment(face1) != selected_segment)) continue;
    RNLength length = mesh->EdgeLength(edge);
    if (length >= min_edge_length) continue;
    mesh->SetEdgeMark(edge, R3mesh_mark);
    heap.Push(edge);
  }

  // Collapse edges in shortest to longest order
  while (!heap.IsEmpty()) {
    GSVMeshEdge *edge = heap.Pop();
    assert(mesh->EdgeMark(edge) == R3mesh_mark);
    assert(mesh->IsEdgeOnMesh(edge));

    // Check edge length
    if (mesh->EdgeLength(edge) > min_edge_length) break;

    // Remove edges to be deleted from queue (e01, e11)
    GSVMeshVertex *v1 = mesh->VertexOnEdge(edge, 1);
    GSVMeshFace *f0 = mesh->FaceOnEdge(edge, 0);  
    GSVMeshFace *f1 = mesh->FaceOnEdge(edge, 1);  
    GSVMeshEdge *e01 = (f0) ? mesh->EdgeAcrossVertex(v1, edge, f0) : NULL;
    GSVMeshEdge *e11 = (f1) ? mesh->EdgeAcrossVertex(v1, edge, f1) : NULL;
    if (e01) { heap.Remove(e01); mesh->SetEdgeMark(e01, 0); }
    if (e11) { heap.Remove(e11); mesh->SetEdgeMark(e11, 0); }

    // Collapse edge
    GSVMeshVertex *vertex = mesh->CollapseEdge(edge);
    if (!vertex) continue;

    // Update adjacent edges
    for (int i = 0; i < mesh->VertexValence(vertex); i++) {
      GSVMeshEdge *edge = mesh->EdgeOnVertex(vertex, i);
      if (mesh->EdgeMark(edge) == R3mesh_mark) {
        heap.Update(edge);
      }
      else {
        RNLength length = mesh->EdgeLength(edge);
        if (length < min_edge_length) {
          mesh->SetEdgeMark(edge, R3mesh_mark);
          heap.Push(edge);
        }
      }
    }

    // Update statistics
    count++;
  }

  // Print message
  if (print_debug) {
    printf("    Collapsed edges\n");
    printf("      Time = %.2f seconds\n", start_time.Elapsed());
    printf("      Minimum Edge Length = %g\n", min_edge_length);
    printf("      # Collapses = %d\n", count);
    fflush(stdout);
  }

  // Return success
  return 1;
}



int SmoothVertices(GSVMesh *mesh, RNScalar sigma, int max_iterations, int selected_segment)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int count = 0;

  // Check inputs
  if (sigma == 0) return 1;
  if (max_iterations == 0) return 1;

  // Create array of vertices to smooth
  R3mesh_mark++;
  RNArray<GSVMeshVertex *> vertices;
  for (int i = 0; i < mesh->NVertices(); i++) {
    RNBoolean member = TRUE;
    GSVMeshVertex *vertex = mesh->Vertex(i);
    for (int j = 0; j < mesh->VertexValence(vertex); j++) {
      GSVMeshEdge *edge = mesh->EdgeOnVertex(vertex, j);
      GSVMeshFace *face0 = mesh->FaceOnEdge(edge, 0);
      GSVMeshFace *face1 = mesh->FaceOnEdge(edge, 1);
      if (!face0 || !face1) { member = FALSE; break; }
      if ((selected_segment >= 0) && (mesh->FaceSegment(face0) != selected_segment)) { member = FALSE; break; }
      if ((selected_segment >= 0) && (mesh->FaceSegment(face1) != selected_segment)) { member = FALSE; break; }
    }
    if (!member) continue;
    mesh->SetVertexMark(vertex, R3mesh_mark);
    vertices.Insert(vertex);
  }

  // Iteratively smooth vertices
  // This is not right because of order dependency (should average old positions)
  // but it is faster and probably good enough for our purposes here
  for (int iter = 0; iter < max_iterations; iter++) {
    for (int i = 0; i < vertices.NEntries(); i++) {
      GSVMeshVertex *vertex = vertices.Kth(i);
      assert(mesh->VertexMark(vertex) == R3mesh_mark);

      // Compute vertex position
      RNScalar weight = 1;
      R3Point position = mesh->VertexPosition(vertex);
      RNScalar exp_factor = 1.0 / (-2.0 * sigma * sigma);
      for (int j = 0; j < mesh->VertexValence(vertex); j++) {
        GSVMeshEdge *edge = mesh->EdgeOnVertex(vertex, j);
        GSVMeshVertex *neighbor_vertex = mesh->VertexAcrossEdge(edge, vertex);
        const R3Point& neighbor_position = mesh->VertexPosition(neighbor_vertex);
        RNLength length = mesh->EdgeLength(edge);
        RNScalar w = exp(exp_factor * length * length);
        GSVMeshFace *face0 = mesh->FaceOnEdge(edge, 0);
        GSVMeshFace *face1 = mesh->FaceOnEdge(edge, 1);
        if (!face0 || !face1) continue;
        if (selected_segment >= 0) {
          if ((mesh->FaceSegment(face0) != selected_segment) || (mesh->FaceSegment(face1) != selected_segment)) {
            w *= 10;
          }
        }
        position += w * neighbor_position;
        weight += w;
      }
      
      // Update vertex position
      R3Point smooth_position = position / weight;
      mesh->SetVertexPosition(vertex, smooth_position);

      // Update statistics
      count++;
    }
  }

  // Print message
  if (print_debug) {
    printf("    Smoothed vertices\n");
    printf("      Time = %.2f seconds\n", start_time.Elapsed());
    printf("      Sigma = %g\n", sigma);
    printf("      # Moves = %d\n", count);
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Hole filling
////////////////////////////////////////////////////////////////////////

struct BoundaryVertexData {
  GSVMeshVertex *vertex;
  GSVMeshEdge *edge;
  BoundaryVertexData *prev;
  BoundaryVertexData *next;
  BoundaryVertexData **heappointer;
  RNScalar error;
};



static RNScalar
BoundaryVertexError(GSVMesh *mesh, GSVMeshVertex *vertex, 
  GSVMeshVertex *prev_vertex, GSVMeshVertex *next_vertex,
  RNLength max_edge_length_squared)
{
  // Check for repeated vertex (figure 8) -- would create degenerate face
  if (vertex == prev_vertex) return RN_INFINITY;
  if (vertex == next_vertex) return RN_INFINITY;
  if (prev_vertex == next_vertex) return RN_INFINITY;

  // Check if edge already exists (across back side) -- would create fin
  if (mesh->EdgeBetweenVertices(prev_vertex, next_vertex)) return RN_INFINITY;

  // Check if any vertex is on left/right/top/bottom boundary
  for (int i = 0; i < 3; i++) {
    int boundary_type = mesh->VertexBoundaryType(vertex);
    if (i == 1) boundary_type = mesh->VertexBoundaryType(prev_vertex);
    else if (i == 2) boundary_type = mesh->VertexBoundaryType(next_vertex);
    if (boundary_type == GSV_MESH_LEFT_BOUNDARY) return RN_INFINITY;
    if (boundary_type == GSV_MESH_RIGHT_BOUNDARY) return RN_INFINITY;
    if (boundary_type == GSV_MESH_TOP_BOUNDARY) return RN_INFINITY;
    if (boundary_type == GSV_MESH_BOTTOM_BOUNDARY) return RN_INFINITY;
  }

  // Check if edge would exceed max edge length
  if (max_edge_length_squared > 0) {
    const R3Point& prev_position = mesh->VertexPosition(prev_vertex);
    const R3Point& next_position = mesh->VertexPosition(next_vertex);
    if (R3SquaredDistance(prev_position, next_position) > max_edge_length_squared) return RN_INFINITY;
  }

  // Initialize error
  RNScalar error = 0;

  // Add difference between depths
  int prev_depth = mesh->VertexLaserDepth(prev_vertex);
  int next_depth = mesh->VertexLaserDepth(next_vertex);
  error += fabs(prev_depth - next_depth);

  // Add angle between edges
  const R3Vector& normal = mesh->VertexNormal(vertex);
  const R3Point& position = mesh->VertexPosition(vertex);
  const R3Point& prev_position = mesh->VertexPosition(prev_vertex);
  const R3Point& next_position = mesh->VertexPosition(next_vertex);
  R3Vector prev_vector = prev_position - position;
  R3Vector next_vector = next_position - position;
  RNAngle interior_angle = R3InteriorAngle(prev_vector, next_vector);
  R3Vector cross = next_vector % prev_vector;
  if (normal.Dot(cross) < 0) interior_angle = RN_TWO_PI - interior_angle;
  if (RNIsEqual(interior_angle, RN_TWO_PI)) return RN_INFINITY;
  // error += interior_angle;

#if 0
  // Add angle between normals
  RNAngle normal_angle = 0;
  int boundary_type = mesh->VertexBoundaryType(vertex);
  if (boundary_type != GSV_MESH_SILHOUETTE_BOUNDARY) {
    const R3Vector& prev_normal = mesh->VertexNormal(prev_vertex);
    const R3Vector& next_normal = mesh->VertexNormal(next_vertex);
    normal_angle += R3InteriorAngle(new_normal, normal);
    normal_angle += R3InteriorAngle(new_normal, prev_normal);
    normal_angle += R3InteriorAngle(new_normal, next_normal);
  }
  error += normal_angle / 3;

  // Add distance between points
  // const R3Point& prev_position = mesh->VertexPosition(prev_vertex);
  // const R3Point& next_position = mesh->VertexPosition(next_vertex);
  // error += R3Distance(prev_position, next_position);
#endif

  // Favor edges between vertices of the same type
  int prev_boundary_type = mesh->VertexBoundaryType(prev_vertex);
  int next_boundary_type = mesh->VertexBoundaryType(next_vertex);
  if (prev_boundary_type != next_boundary_type) error *= 1.5;

  // Return error
  return error;
}



static void
FillHoleBoundary(GSVMesh *mesh, BoundaryVertexData *head, RNLength max_edge_length_squared)
{
  // Initialize heap
  BoundaryVertexData tmp;
  RNHeap<BoundaryVertexData *> heap(&tmp, &(tmp.error), &(tmp.heappointer));
  BoundaryVertexData *vdata = head;
  if (vdata) do {
    heap.Push(vdata);
    vdata = vdata->next;
  } while (vdata != head);

  // Iteratively select boundary vertex to form apex of new triangle
  while (heap.NEntries() > 3) {
    // Pop vertex off heap
    BoundaryVertexData *vdata = heap.Pop();
    BoundaryVertexData *prev_vdata = vdata->prev;
    BoundaryVertexData *next_vdata = vdata->next;
    GSVMeshVertex *vertex = vdata->vertex;
    // GSVMeshEdge *edge = vdata->edge;
    RNScalar error = vdata->error;
    delete vdata;

    // Check error
    if (error == RN_INFINITY) {
      prev_vdata->next = next_vdata;
      next_vdata->prev = prev_vdata;
      break;
    }
    
    // Get vertices
    GSVMeshVertex *prev_vertex = prev_vdata->vertex;
    GSVMeshVertex *next_vertex = next_vdata->vertex;

    // Create face
    GSVMeshFace *face = mesh->CreateFace(prev_vertex, vertex, next_vertex);
    if (!face) { /* ERROR */ return; }

    // Update edge 
    prev_vdata->edge = mesh->EdgeAcrossFace(face, vertex);

    // Remove vdata from linked list
    prev_vdata->next = next_vdata;
    next_vdata->prev = prev_vdata;

    // Update error of prev_vdata and next_vdata
    for (int i = 0; i < 2; i++) {
      BoundaryVertexData *vdata = (i == 0) ? prev_vdata : next_vdata;
      if (!vdata) continue;
      GSVMeshVertex *prev_vertex = vdata->prev->vertex;
      GSVMeshVertex *next_vertex = vdata->next->vertex;
      vdata->error = BoundaryVertexError(mesh, vdata->vertex, prev_vertex, next_vertex, max_edge_length_squared);
      heap.Update(vdata);
    }
  }

  // Create last triangle
  if (heap.NEntries() == 3) {
    BoundaryVertexData *vdata = heap.Pop();
    GSVMeshVertex *vertex = vdata->vertex;
    GSVMeshVertex *prev_vertex = vdata->prev->vertex;
    GSVMeshVertex *next_vertex = vdata->next->vertex;
    mesh->CreateFace(prev_vertex, vertex, next_vertex);
    delete vdata;
  }

  // Delete remaining data (other two vertices of last triangle)
  while (!heap.IsEmpty()) {
    vdata = heap.Pop(); 
    delete vdata;
  }
}



static void 
FillHole(GSVMesh *mesh, GSVMeshEdge *seed_edge, RNLength max_edge_length_squared)
{
  // Check if seed edge is on boundary
  GSVMeshFace *seed_face = mesh->FaceOnEdge(seed_edge);
  if (!seed_face) return; 
  if (mesh->FaceAcrossEdge(seed_edge, seed_face)) return;
  GSVMeshVertex *seed_vertex = mesh->VertexOnEdge(seed_edge, seed_face, RN_CCW);

  // Create counterclockwise linked list of boundary vertices
  BoundaryVertexData *head = NULL;
  BoundaryVertexData *tail = NULL;
  GSVMeshEdge *edge = seed_edge;
  GSVMeshVertex *vertex = seed_vertex;
  GSVMeshFace *face = seed_face;
  do {
    // Create boundary vertex 
    BoundaryVertexData *vdata = new BoundaryVertexData();
    vdata->vertex = vertex;
    vdata->edge = edge;
    vdata->prev = NULL;
    vdata->next = NULL;
    vdata->heappointer = NULL;
    vdata->error = 0;

    // Add boundary vertex to list
    if (!head) head = vdata;
    if (tail) tail->next = vdata;
    vdata->prev = tail;
    tail = vdata;

    // Update vdata error 
    if (vdata->prev) {
      GSVMeshVertex *prev_vertex = vdata->prev->vertex;
      GSVMeshVertex *next_vertex = mesh->VertexAcrossEdge(vdata->edge, vdata->vertex);
      vdata->error = BoundaryVertexError(mesh, vertex, prev_vertex, next_vertex, max_edge_length_squared);
    }

    // Mark edge (so don't revisit in later hole)
    mesh->SetEdgeMark(edge, R3mesh_mark);

    // Find next vertex
    vertex = mesh->VertexAcrossEdge(edge, vertex);

    // Find next edge
    GSVMeshEdge *e = edge;
    while (e) {
      GSVMeshFace *f = mesh->FaceOnVertex(vertex, e, RN_CCW);
      if (!f) { assert(e != edge); edge = e; break; }
      e = mesh->EdgeOnFace(f, vertex, RN_CW);
      if (!e || (e == edge)) return; 
    } 

    // Find next face
    face = mesh->FaceOnVertex(vertex, edge, RN_CW);
  } while (vertex != seed_vertex);

  // Complete loop
  if (head && tail) {
    // Add link between head and tail
    head->prev = tail;
    tail->next = head;

    // Update head error 
    if (head) {
      BoundaryVertexData *vdata = head;
      GSVMeshVertex *prev_vertex = vdata->prev->vertex;
      GSVMeshVertex *next_vertex = vdata->next->vertex;
      vdata->error = BoundaryVertexError(mesh, vdata->vertex, prev_vertex, next_vertex, max_edge_length_squared);
    }
  }

  // Fill hole 
  FillHoleBoundary(mesh, head, max_edge_length_squared);
}



int FillHoles(GSVMesh *mesh, RNLength max_edge_length)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Iteratively find boundary edge and fill hole
  R3mesh_mark++;
  for (int iter = 0; iter < 4; iter++) {
    for (int i = 0; i < mesh->NEdges(); i++) {
      GSVMeshEdge *edge = mesh->Edge(i);
      if (mesh->EdgeMark(edge) == R3mesh_mark) continue;
      if (!mesh->IsEdgeOnBoundary(edge)) continue;
      FillHole(mesh, edge, max_edge_length * max_edge_length);
    }
  }

  // Print debug statistics
  if (print_debug) {
    printf("    Filled holes ...\n");
    printf("      Time = %.2f seconds\n", start_time.Elapsed());
    printf("      # Faces = %d\n", mesh->NFaces());
    printf("      # Edges = %d\n", mesh->NEdges());
    printf("      # Vertices = %d\n", mesh->NVertices());
    fflush(stdout);
  }

  // Return success
  return 1;
}





