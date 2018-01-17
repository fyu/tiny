/* Source file for the R3 scene class */



/* Include files */

#include "R3Graphics.h"




/* Member functions */

int 
R3InitScene()
{
  /* Return success */
  return TRUE;
}



void 
R3StopScene()
{
}



R3Scene::
R3Scene(void)
  : nodes(),
    lights(),
    ambient(0, 0, 0),
    background(0, 0, 0)
{
  // Create root node
  root = new R3SceneNode(this);

  // Initialize viewer
  viewer.SetCamera(R3default_camera);
  viewer.SetViewport(R2default_viewport);
}



R3Scene::
~R3Scene(void)
{
  // Delete everything
  // ???
}



void R3Scene::
InsertNode(R3SceneNode *node) 
{
  // Insert node
  assert(!node->scene);
  assert(node->scene_index == -1);
  node->scene = this;
  node->scene_index = nodes.NEntries();
  nodes.Insert(node);
}



void R3Scene::
RemoveNode(R3SceneNode *node) 
{
  // Remove node
  assert(node->scene == this);
  assert(node->scene_index >= 0);
  RNArrayEntry *entry = nodes.KthEntry(node->scene_index);
  assert(nodes.EntryContents(entry) == node);
  R3SceneNode *tail = nodes.Tail();
  nodes.EntryContents(entry) = tail;
  tail->scene_index = node->scene_index;
  nodes.RemoveTail();
  node->scene_index = -1;
  node->scene = NULL;
}



void R3Scene::
InsertLight(R3Light *light) 
{
  // Insert light
  lights.Insert(light);
}



void R3Scene::
RemoveLight(R3Light *light) 
{
  // Remove light
  lights.Remove(light);
}



void R3Scene::
SetCamera(const R3Camera& camera) 
{
  // Remember camera
  viewer.SetCamera(camera);
}



void R3Scene::
SetViewport(const R2Viewport& viewport) 
{
  // Remember viewport
  viewer.SetViewport(viewport);
}



void R3Scene::
SetViewer(const R3Viewer& viewer) 
{
  // Remember viewer
  this->viewer = viewer;
}



RNBoolean R3Scene::
Intersects(const R3Ray& ray, 
  R3SceneNode **hit_node, R3SceneElement **hit_element, R3Shape **hit_shape,
  R3Point *hit_point, R3Vector *hit_normal, RNScalar *hit_t) const
{
  // Intersect with root node
  if (!root->Intersects(ray, hit_node, hit_element, hit_shape, hit_point, hit_normal, hit_t)) return FALSE;

  // Normalize normal vector
  if (hit_normal) hit_normal->Normalize();

  // Return hit
  return TRUE;
}



void R3Scene::
Draw(const R3DrawFlags draw_flags, RNBoolean set_camera, RNBoolean set_lights) const
{
  // Draw null material
  R3null_material.Draw();

  // Set camera
  if (set_camera) {
    viewer.Camera().Load(); 
  }

  // Set lights
  if (set_lights) {
    for (int i = 0; i < lights.NEntries(); i++) {
      lights.Kth(i)->Draw(i);
    }
  }

  // Draw nodes
  root->Draw(draw_flags);

  // Draw null material
  R3null_material.Draw();
}



////////////////////////////////////////////////////////////////////////
// I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3Scene::
ReadFile(const char *filename)
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .txt)\n", filename);
    return 0;
  }

  // Read file of appropriate type
  if (!strncmp(extension, ".scn", 4)) {
    if (!ReadPrincetonFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".ssc", 4)) {
    if (!ReadParseFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".txt", 4)) {
    if (!ReadSupportHierarchyFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".hier", 5)) {
    if (!ReadGrammarHierarchyFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".obj", 4)) {
    if (!ReadObjFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".off", 4)) {
    if (!ReadMeshFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".ply", 4)) {
    if (!ReadMeshFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".rct", 4)) {
    if (!ReadRectangleFile(filename)) return 0;
  }
  else {
    fprintf(stderr, "Unable to read file %s (unrecognized extension: %s)\n", filename, extension);
    return 0;
  }

  // Provide default camera
  if (Camera() == R3default_camera) {
    double scene_radius = BBox().DiagonalRadius();
    R3Point scene_center = BBox().Centroid();
    R3Vector towards = R3Vector(0, 0, -1);
    R3Vector up = R3Vector(0, 1, 0);
    R3Point eye = scene_center - 3 * scene_radius * towards;
    R3Camera camera(eye, towards, up, 0.25, 0.25, 0.01 * scene_radius, 100 * scene_radius);
    SetCamera(camera);
  }

  // Provide default lights
  if (NLights() == 0) {
    // Create first directional light
    RNRgb color1(1,1,1);
    R3Vector direction1(-3,-4,-5);
    direction1.Normalize();
    R3DirectionalLight *light1 = new R3DirectionalLight(direction1, color1);
    InsertLight(light1);

    // Create second directional light
    RNRgb color2(0.5, 0.5, 0.5);
    R3Vector direction2(3,2,3);
    direction2.Normalize();
    R3DirectionalLight *light2 = new R3DirectionalLight(direction2, color2);
    InsertLight(light2);
  }

  // Return success
  return 1;
}



int R3Scene::
WriteFile(const char *filename) const
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .txt)\n", filename);
    return 0;
  }

  // Write file of appropriate type
  if (!strncmp(extension, ".txt", 4)) {
    if (!WriteSupportHierarchyFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".scn", 4)) {
    if (!WritePrincetonFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".ssc", 4)) {
    if (!WriteParseFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".obj", 4)) {
    if (!WriteObjFile(filename)) return 0;
  }
  else {
    fprintf(stderr, "Unable to write file %s (unrecognized extension: %s)\n", filename, extension);
    return 0;
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// OBJ FILE I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

static R3Material *
FindMaterial(const RNArray<R3Material *>& materials, const char *name)
{
  // Search for material
  for (int i = 0; i < materials.NEntries(); i++) {
    R3Material *material = materials.Kth(i);
    if (!strcmp(name, material->Name())) return material;
  }

  // Material not found
  return NULL;
}



static int
InsertSceneElement(R3Scene *scene, R3SceneNode *node, R3Material *material, 
  const RNArray<R3TriangleVertex *>& verts, const RNArray<R3Triangle *>& tris)
{
  // Find element
  R3SceneElement *element = NULL;
  for (int i = 0; i < node->NElements(); i++) {
    R3SceneElement *g = node->Element(i);
    if (g->Material() == material) { element = g; break; }
  }

  // Create element
  if (!element) {
    element = new R3SceneElement(material);
  }

  // Mark vertices used by tris
  for (int i = 0; i < verts.NEntries(); i++) {
    R3TriangleVertex *vertex = verts.Kth(i);
    vertex->SetMark(0);
  }

  // Create array of verts used by tris
  RNArray<R3TriangleVertex *> tri_verts;
  for (int i = 0; i < tris.NEntries(); i++) {
    R3Triangle *triangle = tris.Kth(i);
    for (int j = 0; j < 3; j++) {
      R3TriangleVertex *vertex = triangle->Vertex(j);
      if (vertex->Mark() == 0) {
        tri_verts.Insert(vertex);
        vertex->SetMark(1);
      }
    }
  }

  // Create shape (triangle array)
  R3TriangleArray *shape = new R3TriangleArray(tri_verts, tris);
        
  // Insert shape
  element->InsertShape(shape);
        
  // Insert element
  node->InsertElement(element);

  // Return success
  return 1;
}



static int
ReadObjMtlFile(RNArray<R3Material *>& materials, const char *dirname, const char *mtlname)
{
  // Open file
  char filename[1024];
  sprintf(filename, "%s/%s", dirname, mtlname);
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Parse file
  char buffer[1024];
  int line_count = 0;
  R3Brdf *brdf = NULL;
  R2Texture *texture = NULL;
  R3Material *material = NULL;
  while (fgets(buffer, 1023, fp)) {
    // Increment line counter
    line_count++;

    // Skip white space
    char *bufferp = buffer;
    while (isspace(*bufferp)) bufferp++;

    // Skip blank lines and comments
    if (*bufferp == '#') continue;
    if (*bufferp == '\0') continue;

    // Get keyword
    char keyword[80];
    if (sscanf(bufferp, "%s", keyword) != 1) {
      fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
      return 0;
    }

    // Check keyword
    if (!strcmp(keyword, "newmtl")) {
      // Parse line
      char name[1024];
      if (sscanf(bufferp, "%s%s", keyword, name) != (unsigned int) 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Create new material
      texture = NULL;
      brdf = new R3Brdf();
      material = new R3Material(brdf, texture, name);
      materials.Insert(material);
    }
    else if (!strcmp(keyword, "Ka")) {
      // Parse line
      double r, g, b;
      if (sscanf(bufferp, "%s%lf%lf%lf", keyword, &r, &g, &b) != (unsigned int) 4) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Set ambient reflectance 
      if (material && brdf) {
        brdf->SetAmbient(RNRgb(r, g, b));
        material->Update();
      }
    }
    else if (!strcmp(keyword, "Kd")) {
      // Parse line
      double r, g, b;
      if (sscanf(bufferp, "%s%lf%lf%lf", keyword, &r, &g, &b) != (unsigned int) 4) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Set diffuse reflectance 
      if (material && brdf) {
        brdf->SetDiffuse(RNRgb(r, g, b));
        material->Update();
      }
    }
    else if (!strcmp(keyword, "Ks")) {
      // Parse line
      double r, g, b;
      if (sscanf(bufferp, "%s%lf%lf%lf", keyword, &r, &g, &b) != (unsigned int) 4) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Set specular reflectance 
      if (material && brdf) {
        brdf->SetSpecular(RNRgb(r, g, b));
        material->Update();
      }
    }
    else if (!strcmp(keyword, "Ns")) {
      // Parse line
      double ns;
      if (sscanf(bufferp, "%s%lf", keyword, &ns) != (unsigned int) 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Set shininess
      if (material && brdf) {
        brdf->SetShininess(ns);
        material->Update();
      }
    }
    else if (!strcmp(keyword, "Ni")) {
      // Parse line
      double index_of_refraction;
      if (sscanf(bufferp, "%s%lf", keyword, &index_of_refraction) != (unsigned int) 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Set index of refraction
      if (material && brdf) {
        brdf->SetIndexOfRefraction(index_of_refraction);
        material->Update();
      }
    }
    else if (!strcmp(keyword, "d")) {
      // Parse line
      double transparency;
      if (sscanf(bufferp, "%s%lf", keyword, &transparency) != (unsigned int) 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Set opacity
      if (material && brdf) {
        brdf->SetOpacity(1 - transparency);
        material->Update();
      }
    }
    else if (!strcmp(keyword, "map_Kd")) {
      // Parse line
      char texture_name[1024];
      if (sscanf(bufferp, "%s%s", keyword, texture_name) != (unsigned int) 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Set texture
      if (material) {
        char texture_filename[1024];
        sprintf(texture_filename, "%s/%s", dirname, texture_name);
        R2Image *image = new R2Image();
        if (!image->Read(texture_filename)) return 0;
        R2Texture *texture = new R2Texture(image);
        material->SetTexture(texture);
        material->Update();
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



static int
ReadObj(R3Scene *scene, R3SceneNode *node, const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Determine directory name (for texture image files)
  char dirname[1024];
  strncpy(dirname, filename, 1024);
  char *endp = strrchr(dirname, '/');
  if (!endp) endp = strrchr(dirname, '\\');
  if (!endp) strcpy(dirname, ".");
  else *endp = '\0';

  // Read body
  char buffer[1024];
  int line_count = 0;
  R3Material *material = &R3default_material;
  RNArray<R3Material *> materials;
  RNArray<R2Point *> texture_coords;
  RNArray<R3TriangleVertex *> verts;
  RNArray<R3Triangle *> tris;
  while (fgets(buffer, 1023, fp)) {
    // Increment line counter
    line_count++;

    // Skip white space
    char *bufferp = buffer;
    while (isspace(*bufferp)) bufferp++;

    // Skip blank lines and comments
    if (*bufferp == '#') continue;
    if (*bufferp == '\0') continue;

    // Get keyword
    char keyword[80];
    if (sscanf(bufferp, "%s", keyword) != 1) {
      fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
      return 0;
    }

    // Check keyword
    if (!strcmp(keyword, "v")) {
      // Read vertex coordinates
      double x, y, z;
      if (sscanf(bufferp, "%s%lf%lf%lf", keyword, &x, &y, &z) != 4) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Create vertex
      R3TriangleVertex *vertex = new R3TriangleVertex(R3Point(x, y, z));
      verts.Insert(vertex);
    }
    else if (!strcmp(keyword, "vt")) {
      // Read texture coordinates
      double u, v;
      if (sscanf(bufferp, "%s%lf%lf", keyword, &u, &v) != 3) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Create texture coordinates
      R2Point *vt = new R2Point(u, v);
      texture_coords.Insert(vt);
    }
    else if (!strcmp(keyword, "f")) {
      // Read vertex indices
      int quad = 1;
      char s1[128], s2[128], s3[128], s4[128] = { '\0' };
      if (sscanf(bufferp, "%s%s%s%s%s", keyword, s1, s2, s3, s4) != 5) {
        quad = 0;;
        if (sscanf(bufferp, "%s%s%s%s", keyword, s1, s2, s3) != 4) {
          fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
          return 0;
        }
      }

      // Parse vertex indices
      int vi1 = -1, vi2 = -1, vi3 = -1, vi4 = -1;
      int ti1 = -1, ti2 = -1, ti3 = -1, ti4 = -1;
      char *p1 = strchr(s1, '/'); 
      if (p1) { *p1 = 0; vi1 = atoi(s1); p1++; if (*p1) ti1 = atoi(p1); }
      else { vi1 = atoi(s1); ti1 = vi1; }
      char *p2 = strchr(s2, '/'); 
      if (p2) { *p2 = 0; vi2 = atoi(s2); p2++; if (*p2) ti2 = atoi(p2); }
      else { vi2 = atoi(s2); ti2 = vi2; }
      char *p3 = strchr(s3, '/'); 
      if (p3) { *p3 = 0; vi3 = atoi(s3); p3++; if (*p3) ti3 = atoi(p3); }
      else { vi3 = atoi(s3); ti3 = vi3; }
      if (quad) {
        char *p4 = strchr(s4, '/'); 
        if (p4) { *p4 = 0; vi4 = atoi(s4); p4++; if (*p4) ti4 = atoi(p4); }
        else { vi4 = atoi(s4); ti4 = vi4; }
      }

      // Get vertices
      R3TriangleVertex *v1 = verts.Kth(vi1-1);
      R3TriangleVertex *v2 = verts.Kth(vi2-1);
      R3TriangleVertex *v3 = verts.Kth(vi3-1);
      R3TriangleVertex *v4 = (quad) ? verts.Kth(vi4-1) : NULL;
      
      // Assign texture coordinates
      if ((ti1 >= 0) && (ti1 < texture_coords.NEntries())) v1->SetTextureCoords(*(texture_coords.Kth(ti1-1)));
      if ((ti2 >= 0) && (ti2 < texture_coords.NEntries())) v2->SetTextureCoords(*(texture_coords.Kth(ti2-1)));
      if ((ti3 >= 0) && (ti3 < texture_coords.NEntries())) v3->SetTextureCoords(*(texture_coords.Kth(ti3-1)));
      if (quad) {
        if ((ti4 >= 0) && (ti4 < texture_coords.NEntries())) v4->SetTextureCoords(*(texture_coords.Kth(ti4-1)));
      }

      // Check vertices
      if ((v1 == v2) || (v2 == v3) || (v1 == v3)) continue;
      if ((quad) && ((v4 == v1) || (v4 == v2) || (v4 == v3))) quad = 0;

      // Create first triangle
      R3Triangle *triangle = new R3Triangle(v1, v2, v3);
      tris.Insert(triangle);

      // Create second triangle
      if (quad) {
        R3Triangle *triangle = new R3Triangle(v1, v3, v4);
        tris.Insert(triangle);
      }
    }
    else if (!strcmp(keyword, "mtllib")) {
      // Read fields
      char mtlname[1024];
      if (sscanf(bufferp, "%s%s", keyword, mtlname) != 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Empty previous set of materials
      materials.Empty();

      // Read materials
      if (!ReadObjMtlFile(materials, dirname, mtlname)) return 0;
    }
    else if (!strcmp(keyword, "usemtl")) {
      // Read fields
      char mtlname[1024];
      if (sscanf(bufferp, "%s%s", keyword, mtlname) != 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Process triangles from previous material
      if (material && (verts.NEntries() > 0) && (tris.NEntries() > 0)) {
        InsertSceneElement(scene, node, material, verts, tris);
      }

      // Find material
      material = FindMaterial(materials, mtlname);
      if (!material) {
        fprintf(stderr, "Unable to find material %s at on line %d in file %s", mtlname, line_count, filename);
        return 0;
      }

      // Empty tris
      tris.Empty();
    }
  }

  // Process triangles from previous material
  if (material && (verts.NEntries() > 0) && (tris.NEntries() > 0)) {
    InsertSceneElement(scene, node, material, verts, tris);
  }

  // Delete texture coordinates
  for (int i = 0; i < texture_coords.NEntries(); i++) {
    R2Point *vt = texture_coords.Kth(i);
    delete vt;
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3Scene::
ReadObjFile(const char *filename)
{
  // Read obj file, and put contents in root node
  return ReadObj(this, root, filename);
}



int R3Scene::
WriteObjFile(const char *filename) const
{
  // Not implemented yet
  RNAbort("Not implemented");

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// MESH FILE I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

static R3TriangleArray *
ReadMesh(const char *filename)
{
  R3Mesh mesh;
  if (!mesh.ReadFile(filename)) {
    fprintf(stderr, "Unable to read mesh %s\n", filename);
    return NULL;
  }
  
  // Create array of vertices
  RNArray<R3TriangleVertex *> vertices;
  for (int i = 0; i < mesh.NVertices(); i++) {
    R3MeshVertex *mesh_vertex = mesh.Vertex(i);
    const R3Point& position = mesh.VertexPosition(mesh_vertex);
    const R3Vector& normal = mesh.VertexNormal(mesh_vertex);
    R3TriangleVertex *triangle_vertex = new R3TriangleVertex(position, normal);
    vertices.Insert(triangle_vertex);
  }

  // Create array of triangles
  RNArray<R3Triangle *> triangles;
  for (int i = 0; i < mesh.NFaces(); i++) {
    R3MeshFace *mesh_face = mesh.Face(i);
    int i0 = mesh.VertexID(mesh.VertexOnFace(mesh_face, 0));
    int i1 = mesh.VertexID(mesh.VertexOnFace(mesh_face, 1));
    int i2 = mesh.VertexID(mesh.VertexOnFace(mesh_face, 2));
    R3TriangleVertex *v0 = vertices.Kth(i0);
    R3TriangleVertex *v1 = vertices.Kth(i1);
    R3TriangleVertex *v2 = vertices.Kth(i2);
    R3Triangle *triangle = new R3Triangle(v0, v1, v2);
    triangles.Insert(triangle);
  }

  // Return triangle array
  return new R3TriangleArray(vertices, triangles);
}

 

int R3Scene::
ReadMeshFile(const char *filename)
{
  // Load triangles into scene
  R3TriangleArray *shape = ReadMesh(filename);
  if (!shape) return 0;
  R3SceneElement *element = new R3SceneElement();
  R3SceneNode *node = new R3SceneNode(this);
  element->InsertShape(shape);
  node->InsertElement(element);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// PRINCETON SCENE FILE I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

static int
ReadPrinceton(R3Scene *scene, R3SceneNode *node, const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Create array of materials
  RNArray<R3Material *> materials;

  // Create stack of group information
  const int max_depth = 1024;
  R3SceneNode *group_nodes[max_depth] = { NULL };
  R3Material *group_materials[max_depth] = { NULL };
  group_nodes[0] = (node) ? node : scene->Root();
  group_materials[0] = new R3Material(&R3default_brdf);
  int depth = 0;

  // Read body
  char cmd[128];
  int command_number = 1;
  while (fscanf(fp, "%s", cmd) == 1) {
    if (cmd[0] == '#') {
      // Comment -- read everything until end of line
      do { cmd[0] = fgetc(fp); } while ((cmd[0] >= 0) && (cmd[0] != '\n'));
    }
    else if (!strcmp(cmd, "tri")) {
      // Read data
      int m;
      R3Point p1, p2, p3;
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf", &m, 
        &p1[0], &p1[1], &p1[2], &p2[0], &p2[1], &p2[2], &p3[0], &p3[1], &p3[2]) != 10) {
        fprintf(stderr, "Unable to read triangle at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.NEntries()) {
          material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at tri command %d in file %s\n", command_number, filename);
          return 0;
        }
      }

      // Create triangle
      R3TriangleVertex *v1 = new R3TriangleVertex(p1);
      R3TriangleVertex *v2 = new R3TriangleVertex(p2);
      R3TriangleVertex *v3 = new R3TriangleVertex(p3);
      R3Triangle *triangle = new R3Triangle(v1, v2, v3);

      // Insert triangle into scene
      R3SceneElement *element = new R3SceneElement(material);
      element->InsertShape(triangle);
      R3SceneNode *node = new R3SceneNode(scene);
      node->InsertElement(element);
      group_nodes[depth]->InsertChild(node);
    }
    else if (!strcmp(cmd, "box")) {
      // Read data
      int m;
      R3Point p1, p2;
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf%lf", &m, &p1[0], &p1[1], &p1[2], &p2[0], &p2[1], &p2[2]) != 7) {
        fprintf(stderr, "Unable to read box at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.NEntries()) {
          material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at box command %d in file %s\n", command_number, filename);
          return 0;
        }
      }

      // Sort coordinates
      if (p1[0] > p2[0]) { RNCoord swap = p1[0]; p1[0] = p2[0]; p2[0] = swap; }
      if (p1[1] > p2[1]) { RNCoord swap = p1[1]; p1[1] = p2[1]; p2[1] = swap; }
      if (p1[2] > p2[2]) { RNCoord swap = p1[2]; p1[2] = p2[2]; p2[2] = swap; }

      // Create box
      R3Box *box = new R3Box(p1, p2);

      // Insert box into scene
      R3SceneElement *element = new R3SceneElement(material);
      element->InsertShape(box);
      R3SceneNode *node = new R3SceneNode(scene);
      node->InsertElement(element);
      group_nodes[depth]->InsertChild(node);
    }
    else if (!strcmp(cmd, "sphere")) {
      // Read data
      int m;
      R3Point c;
      double r;
      if (fscanf(fp, "%d%lf%lf%lf%lf", &m, &c[0], &c[1], &c[2], &r) != 5) {
        fprintf(stderr, "Unable to read sphere at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.NEntries()) {
          material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at sphere command %d in file %s\n", command_number, filename);
          return 0;
        }
      }

      // Create sphere
      R3Sphere *sphere = new R3Sphere(c, r);

      // Insert sphere into scene
      R3SceneElement *element = new R3SceneElement(material);
      element->InsertShape(sphere);
      R3SceneNode *node = new R3SceneNode(scene);
      node->InsertElement(element);
      group_nodes[depth]->InsertChild(node);
    }
    else if (!strcmp(cmd, "cylinder")) {
      // Read data
      int m;
      R3Point c;
      double r, h;
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf", &m, &c[0], &c[1], &c[2], &r, &h) != 6) {
        fprintf(stderr, "Unable to read cylinder at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.NEntries()) {
          material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at cyl command %d in file %s\n", command_number, filename);
          return 0;
        }
      }

      // Create cylinder
      R3Point p1 = c - 0.5 * h * R3posy_vector;
      R3Point p2 = c + 0.5 * h * R3posy_vector;
      R3Cylinder *cylinder = new R3Cylinder(p1, p2, r);

      // Insert cylinder into scene
      R3SceneElement *element = new R3SceneElement(material);
      element->InsertShape(cylinder);
      R3SceneNode *node = new R3SceneNode(scene);
      node->InsertElement(element);
      group_nodes[depth]->InsertChild(node);
    }
    else if (!strcmp(cmd, "cone")) {
      // Read data
      int m;
      R3Point c;
      double r, h;
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf", &m, &c[0], &c[1], &c[2], &r, &h) != 6) {
        fprintf(stderr, "Unable to read cone at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.NEntries()) {
          material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at cone command %d in file %s\n", command_number, filename);
          return 0;
        }
      }

      // Create cone
      R3Point p1 = c - 0.5 * h * R3posy_vector;
      R3Point p2 = c + 0.5 * h * R3posy_vector;
      R3Cone *cone = new R3Cone(p1, p2, r);

      // Insert cone into scene
      R3SceneElement *element = new R3SceneElement(material);
      element->InsertShape(cone);
      R3SceneNode *node = new R3SceneNode(scene);
      node->InsertElement(element);
      group_nodes[depth]->InsertChild(node);
    }
    else if (!strcmp(cmd, "mesh")) {
      // Read data
      int m;
      char meshname[256];
      if (fscanf(fp, "%d%s", &m, meshname) != 2) {
        fprintf(stderr, "Unable to parse mesh command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.NEntries()) {
          material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at cone command %d in file %s\n", command_number, filename);
          return 0;
        }
      }

      // Get mesh filename
      char buffer[2048];
      strcpy(buffer, filename);
      char *bufferp = strrchr(buffer, '/');
      if (bufferp) *(bufferp+1) = '\0';
      else buffer[0] = '\0';
      strcat(buffer, meshname);

      // Read mesh
      R3TriangleArray *mesh = ReadMesh(buffer);
      if (!mesh) return 0;

      // Insert mesh into scene
      R3SceneElement *element = new R3SceneElement(material);
      element->InsertShape(mesh);
      R3SceneNode *node = new R3SceneNode(scene);
      node->InsertElement(element);
      group_nodes[depth]->InsertChild(node);
    }
    else if (!strcmp(cmd, "line")) {
      // Read data
      int m;
      double x1, y1, z1, x2, y2, z2;
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf%lf", &m, &x1, &y1, &z1, &x2, &y2, &z2) != 7) {
        fprintf(stderr, "Unable to read line at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.NEntries()) {
          material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at cyl command %d in file %s\n", command_number, filename);
          return 0;
        }
      }

      // Create cylinder representing line
      R3Point p1(x1, y1, z1);
      R3Point p2(x2, y2, z2);
      R3Cylinder *cylinder = new R3Cylinder(p1, p2, RN_BIG_EPSILON);

      // Insert cylinder into scene
      R3SceneElement *element = new R3SceneElement(material);
      element->InsertShape(cylinder);
      R3SceneNode *node = new R3SceneNode(scene);
      node->InsertElement(element);
      group_nodes[depth]->InsertChild(node);
    }
    else if (!strcmp(cmd, "begin")) {
      // Read data
      int m;
      double matrix[16];
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &m, 
        &matrix[0], &matrix[1], &matrix[2], &matrix[3], 
        &matrix[4], &matrix[5], &matrix[6], &matrix[7], 
        &matrix[8], &matrix[9], &matrix[10], &matrix[11], 
        &matrix[12], &matrix[13], &matrix[14], &matrix[15]) != 17) {
        fprintf(stderr, "Unable to read begin at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Get material
      R3Material *material = group_materials[depth];
      if (m >= 0) {
        if (m < (int) materials.NEntries()) {
          material = materials[m];
        }
        else {
          fprintf(stderr, "Invalid material id at cone command %d in file %s\n", command_number, filename);
          return 0;
        }
      }

      // Create new node
      R3SceneNode *node = new R3SceneNode(scene);
      node->SetTransformation(R3Affine(R4Matrix(matrix)));
      group_nodes[depth]->InsertChild(node);

      // Push node onto stack
      depth++;
      group_nodes[depth] = node;
      group_materials[depth] = material;
    }
    else if (!strcmp(cmd, "end")) {
      // Check stack depth
      if (depth <= 0) {
        fprintf(stderr, "Extra end statement at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Pop node from stack
      depth--;
    }
    else if (!strcmp(cmd, "material")) {
      // Read data
      RNRgb ka, kd, ks, kt, e;
      double n, ir;
      char texture_name[256];
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%s", 
          &ka[0], &ka[1], &ka[2], &kd[0], &kd[1], &kd[2], &ks[0], &ks[1], &ks[2], &kt[0], &kt[1], &kt[2], 
          &e[0], &e[1], &e[2], &n, &ir, texture_name) != 18) {
        fprintf(stderr, "Unable to read material at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Create brdf
      R3Brdf *brdf = new R3Brdf(ka, kd, ks, kt, e, n, ir);

      // Create texture
      R2Texture *texture = NULL;
      if (strcmp(texture_name, "0")) {
        // Get texture filename
        char buffer[2048];
        strcpy(buffer, filename);
        char *bufferp = strrchr(buffer, '/');
        if (bufferp) *(bufferp+1) = '\0';
        else buffer[0] = '\0';
        strcat(buffer, texture_name);

        // Read texture file
        R2Image *image = new R2Image();
        if (!image->Read(buffer)) {
          fprintf(stderr, "Unable to read texture from %s at command %d in file %s\n", buffer, command_number, filename);
          return 0;
        }
        
        // Create texture
        texture = new R2Texture(image);
      }

      // Create material
      R3Material *material = new R3Material(brdf, texture);
      materials.Insert(material);
    }
    else if (!strcmp(cmd, "dir_light")) {
      // Read data
      RNRgb c;
      R3Vector d;
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf", 
        &c[0], &c[1], &c[2], &d[0], &d[1], &d[2]) != 6) {
        fprintf(stderr, "Unable to read directional light at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Normalize direction
      d.Normalize();

      // Create directional light
      R3DirectionalLight *light = new R3DirectionalLight(d, c);
      scene->InsertLight(light);
    }
    else if (!strcmp(cmd, "point_light")) {
      // Read data
      RNRgb c;
      R3Point p;
      double ca, la, qa;
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf", &c[0], &c[1], &c[2], &p[0], &p[1], &p[2], &ca, &la, &qa) != 9) {
        fprintf(stderr, "Unable to read point light at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Create point light
      R3PointLight *light = new R3PointLight(p, c, 1, TRUE, ca, la, qa);
      scene->InsertLight(light);
    }
    else if (!strcmp(cmd, "spot_light")) {
      // Read data
      RNRgb c;
      R3Point p;
      R3Vector d;
      double ca, la, qa, sc, sd;
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
        &c[0], &c[1], &c[2], &p[0], &p[1], &p[2], &d[0], &d[1], &d[2], &ca, &la, &qa, &sc, &sd) != 14) {
        fprintf(stderr, "Unable to read point light at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Normalize direction
      d.Normalize();

      // Create spot light
      R3SpotLight *light = new R3SpotLight(p, d, c, sd, sc, 1, TRUE, ca, la, qa);
      scene->InsertLight(light);
    }
    else if (!strcmp(cmd, "area_light")) {
      // Read data
      RNRgb c;
      R3Point p;
      R3Vector d;
      double radius, ca, la, qa;
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
        &c[0], &c[1], &c[2], &p[0], &p[1], &p[2], &d[0], &d[1], &d[2], &radius, &ca, &la, &qa) != 13) {
        fprintf(stderr, "Unable to read area light at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Normalize direction
      d.Normalize();

      // Create spot light
      R3AreaLight *light = new R3AreaLight(p, radius, d, c, 1, TRUE, ca, la, qa);
      scene->InsertLight(light);
    }
    else if (!strcmp(cmd, "camera")) {
      // Read data
      R3Point e;
      R3Vector t, u;
      RNScalar xfov, neardist, fardist;
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &e[0], &e[1], &e[2], &t[0], &t[1], &t[2], &u[0], &u[1], &u[2], &xfov, &neardist, &fardist) != 12) {
        fprintf(stderr, "Unable to read camera at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Assign camera
      R3Camera camera(e, t, u, xfov, xfov, neardist, fardist);
      scene->SetCamera(camera);
    }
    else if (!strcmp(cmd, "include")) {
      // Read data
      char scenename[256];
      if (fscanf(fp, "%s", scenename) != 1) {
        fprintf(stderr, "Unable to read include command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Get scene filename
      char buffer[2048];
      strcpy(buffer, filename);
      char *bufferp = strrchr(buffer, '/');
      if (bufferp) *(bufferp+1) = '\0';
      else buffer[0] = '\0';
      strcat(buffer, scenename);

      // Read scene from included file
      if (!ReadPrinceton(scene, group_nodes[depth], buffer)) {
        fprintf(stderr, "Unable to read included scene: %s\n", buffer);
        return 0;
      }
    }
    else if (!strcmp(cmd, "background")) {
      // Read data
      double r, g, b;
      if (fscanf(fp, "%lf%lf%lf", &r, &g, &b) != 3) {
        fprintf(stderr, "Unable to read background at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Assign background color
      scene->SetBackground(RNRgb(r, g, b));
    }
    else if (!strcmp(cmd, "ambient")) {
      // Read data
      double r, g, b;
      if (fscanf(fp, "%lf%lf%lf", &r, &g, &b) != 3) {
        fprintf(stderr, "Unable to read ambient at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Assign ambient color
      scene->SetAmbient(RNRgb(r, g, b));
    }
    else {
      fprintf(stderr, "Unrecognized command %d in file %s: %s\n", command_number, filename, cmd);
      return 0;
    }
	
    // Increment command number
    command_number++;
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3Scene::
ReadPrincetonFile(const char *filename)
{
  // Read princeton file and insert contents into root node
  return ReadPrinceton(this, root, filename);
}



int R3Scene::
WritePrincetonFile(const char *filename) const
{
  // Not implemented yet
  RNAbort("Not implemented");

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// SUPPORT HIERARCHY FILE I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3Scene::
ReadSupportHierarchyFile(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Read scene
  char buffer[1024];
  int line_count = 0;
  RNArray<R3SceneNode *> nodes;
  while (fgets(buffer, 1023, fp)) {
    // Increment line counter
    line_count++;

    // Skip white space
    char *bufferp = buffer;
    while (isspace(*bufferp)) bufferp++;

    // Skip blank lines and comments
    if (*bufferp == '#') continue;
    if (*bufferp == '\0') continue;

    // Get keyword
    char keyword[256];
    if (sscanf(bufferp, "%s", keyword) != 1) {
      fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
      return 0;
    }

    // Check keyword
    if (!strcmp(keyword, "newModel")) {
      // Read fields
      int model_index;
      char model_name[1024];
      if (sscanf(bufferp, "%s%d%s", keyword, &model_index, model_name) != (unsigned int) 3) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Create node
      assert(nodes.NEntries() == model_index);
      R3SceneNode *node = (model_index == 0) ? Root() : new R3SceneNode(this);
      node->SetName(model_name);
      nodes.Insert(node);

      // Read obj file
      char model_filename[1024];
      sprintf(model_filename, "models/%s.obj", model_name);
      if (!ReadObj(this, node, model_filename)) return 0;
    }
    else if (!strcmp(keyword, "parentIndex")) {
      // Read fields
      int parent_index;
      if (sscanf(bufferp, "%s%d", keyword, &parent_index) != (unsigned int) 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Check parent index
      if (parent_index < 0) {
        if (nodes.NEntries() != 1) {
          fprintf(stderr, "Root node was not first in file %s", filename);
          return 0;
        }
      }
      else {
        // Just checking
        if (parent_index >= nodes.NEntries()) {
          fprintf(stderr, "Invalid parent node index %d on line %d in file %s", parent_index, line_count, filename);
          return 0;
        }

        // Set last node's parent
        R3SceneNode *node = nodes.Tail();
        R3SceneNode *parent = nodes.Kth(parent_index);
        parent->InsertChild(node);
      }
    }
    else if (!strcmp(keyword, "transform")) {
      // Read fields
      double m[16];
      if (sscanf(bufferp, "%s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
        keyword, &m[0], &m[1], &m[2], &m[3], &m[4], &m[5], &m[6], &m[7], 
        &m[8], &m[9], &m[10], &m[11], &m[12], &m[13], &m[14], &m[15]) != (unsigned int) 17) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // NOTE: Our transformation is stored relative to parent's coordinate system
      // So, we must compute transformation from inverse of transform from parent node to world coordinates so that can 
      // convert file's absolute transform (which goes from node's coordinates to world coordinates)
      // to our transform (which goes from node's coordinates to parent node's coordinates)
      R3Affine transformation = R3identity_affine;
      R3SceneNode *node = nodes.Tail();
      R3SceneNode *ancestor = node->Parent();
      while (ancestor) {
        transformation.InverseTransform(ancestor->Transformation());
        ancestor = ancestor->Parent();
      }

      // Set last node's transformation
      // Note that file's matrix is for post-multiplication, while ours is for pre-multiplication, so need flip
      R4Matrix matrix(m); matrix.Flip();
      transformation.Transform(R3Affine(matrix, 0));
      node->SetTransformation(transformation);
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3Scene::
WriteSupportHierarchyFile(const char *filename) const
{
  // Not implemented yet
  RNAbort("Not implemented");

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// GRAMMAR HIERARCHY FILE I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3Scene::
ReadGrammarHierarchyFile(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Extract base name
  char basename[4096] = { '\0' };
  const char *startp = strrchr(filename, '/');
  if (!startp) startp = filename;
  else startp++;
  const char *endp = strrchr(filename, '.');
  if (!endp) endp = startp + strlen(startp);
  int basename_length = endp - startp;
  if (basename_length > 4095) basename_length = 4095;
  strncpy(basename, startp, basename_length);

  // Read scene
  char buffer[1024];
  int line_count = 0;
  RNBoolean leaf = TRUE;
  RNArray<R3SceneNode *> parsed_nodes;
  R3SceneNode *current_node = NULL;
  while (fgets(buffer, 1023, fp)) {
    // Increment line counter
    line_count++;

    // Skip white space
    char *bufferp = buffer;
    while (isspace(*bufferp)) bufferp++;

    // Skip blank lines and comments
    if (*bufferp == '#') continue;
    if (*bufferp == '\0') continue;

    // Get keyword
    char keyword[256];
    if (sscanf(bufferp, "%s", keyword) != 1) {
      fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
      return 0;
    }

    // Check keyword
    if (!strcmp(keyword, "newModel")) {
      // Read fields
      int model_index;
      if (sscanf(bufferp, "%s%d", keyword, &model_index) != (unsigned int) 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Create node
      while (parsed_nodes.NEntries() <= model_index) {
        char node_name[4096];
        sprintf(node_name, "%d", parsed_nodes.NEntries());
        R3SceneNode *node = new R3SceneNode(this);
        node->SetName(node_name);
        parsed_nodes.Insert(node);
      }

      // Remember current node
      current_node = parsed_nodes.Kth(model_index);
    }
    else if (!strcmp(keyword, "root")) {
      // Read fields
      int root_index;
      if (sscanf(bufferp, "%s%d", keyword, &root_index) != (unsigned int) 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Create node
      while (parsed_nodes.NEntries() <= root_index) {
        char node_name[4096];
        sprintf(node_name, "%d", parsed_nodes.NEntries());
        R3SceneNode *node = new R3SceneNode(this);
        node->SetName(node_name);
        parsed_nodes.Insert(node);
      }

      // Insert root of this parse as child of root
      R3SceneNode *node = parsed_nodes.Kth(root_index);
      root->InsertChild(node);
    }
    else if (!strcmp(keyword, "parent")) {
      // Read fields
      int parent_index;
      if (sscanf(bufferp, "%s%d", keyword, &parent_index) != (unsigned int) 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Check parent index
      if (parent_index >= 0) {
        // Create parent node
        while (parsed_nodes.NEntries() <= parent_index) {
          char node_name[4096];
          sprintf(node_name, "I%d", parsed_nodes.NEntries());
          R3SceneNode *node = new R3SceneNode(this);
          node->SetName(node_name);
          parsed_nodes.Insert(node);
        }
      
        // Set last node's parent
        R3SceneNode *parent = parsed_nodes.Kth(parent_index);
        parent->InsertChild(current_node);
      }
    }
    else if (!strcmp(keyword, "children")) {
      const char *token = strtok(bufferp, " \n\t");
      assert(token && (!strcmp(token, "children")));
      token = strtok(NULL, " \n\t");
      leaf = (token) ? FALSE : TRUE;
    }
    else if (!strcmp(keyword, "leaf_group")) {
      // Check current node
      if (!current_node) {
        fprintf(stderr, "leaf_group before first newModel at line %d in %s\n", line_count, filename);
        return 0;
      }

      // Read models
      if (leaf) {
        const char *token = strtok(bufferp, " \n\t");
        assert(token && (!strcmp(token, "leaf_group")));
        while (TRUE) {
          token = strtok(NULL, " \n\t");
          if (!token) break;
          int model_index = atoi(token);
          // char node_name[4096];
          // sprintf(node_name, "L%d", model_index);
          // R3SceneNode *node = new R3SceneNode(this);
          // node->SetName(node_name);
          // current_node->InsertChild(node);
          char model_filename[1024];
          sprintf(model_filename, "models/%s/%d.off", basename, model_index);
          R3TriangleArray *shape = ReadMesh(model_filename);
          if (shape) {
            R3SceneElement *element = new R3SceneElement();
            element->SetMaterial(&R3default_material);
            element->InsertShape(shape);
            current_node->InsertElement(element);
          }
        }    
      }    
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3Scene::
WriteGrammarHierarchyFile(const char *filename) const
{
  // Not implemented yet
  RNAbort("Not implemented");

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// PARSE FILE I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3Scene::
ReadParseFile(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Read header
  char buffer[4096];
  if (!fgets(buffer, 4096, fp)) {
    fprintf(stderr, "Unable to read object parse file %s\n", filename);
    return 0;
  }

  // Check header
  if (strncmp(buffer, "OBJECT PARSE 1.0", 16)) {
    fprintf(stderr, "Error in header of oject parse file %s\n", filename);
    return 0;
  }

  // Read file
  int line_number = 0;
  int assignment_index = 0;
  RNArray<R3Shape *> shapes;
  char mesh_directory[4096] = { '.', '\0' };
  while (fgets(buffer, 4096, fp)) {
    // Check line
    line_number++;
    char *bufferp = buffer;
    while (*bufferp && isspace(*bufferp)) bufferp++;
    if (!bufferp) continue;
    if (*bufferp == '#') continue;

    // Parse line
    char keyword[4096];
    if (sscanf(bufferp, "%s", keyword) == (unsigned int) 1) {
      if (!strcmp(keyword, "A")) {
        // Parse assignment
        double score, m[16];
        int segmentation_index, model_index, dummy;
        if (sscanf(bufferp, "%s%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%d%d%d%d", keyword, 
          &segmentation_index, &model_index, 
          &m[0], &m[1], &m[2], &m[3], &m[4], &m[5], &m[6], &m[7], 
          &m[8], &m[9], &m[10], &m[11], &m[12], &m[13], &m[14], &m[15], 
          &score, &dummy, &dummy, &dummy, &dummy) != (unsigned int) 24) {
          fprintf(stderr, "Error parsing assignment at line %d of %s\n", line_number, filename);
          return 0;
        }

        // Create node
        R3SceneNode *node = new R3SceneNode(this);
        root->InsertChild(node);

        // Create shape element
        R3Shape *shape = shapes.Kth(model_index);
        R3Brdf *brdf = new R3Brdf(RNRgb(0, 0.25 + score, 0));
        R3Material *material = new R3Material(brdf);
        R3SceneElement *element = new R3SceneElement();
        element->SetMaterial(material);
        element->InsertShape(shape);
        node->InsertElement(element);

        // Set node name
        char node_name[1024];
        sprintf(node_name, "A%d_M%d_S%03d", assignment_index++, model_index, (int) (1000*score));
        node->SetName(node_name);

        // Set node transformation
        R4Matrix matrix(m); 
        R3Affine affine(matrix, 0);
        node->SetTransformation(affine);
      }
      else if (!strcmp(keyword, "M")) {
        // Parse model
        int dummy;
        char mesh_name[4096];
        double cx, cy, cz, r, h;
        if (sscanf(bufferp, "%s%d%lf%lf%lf%lf%lf%d%s%d%d%d%d%d", keyword, 
          &dummy, &cx, &cy, &cz, &r, &h, &dummy, mesh_name, &dummy, &dummy, &dummy, &dummy, &dummy) != (unsigned int) 14) {
          fprintf(stderr, "Error parsing model at line %d of %s\n", line_number, filename);
          return 0;
        }

        // Read mesh
        char mesh_filename[4096];
        sprintf(mesh_filename, "%s/%s", mesh_directory, mesh_name);
        R3Shape *shape = ReadMesh(mesh_filename);
        if (!shape) return 0;

        // Insert shape into list
        shapes.Insert(shape);
      }
      else if (!strcmp(keyword, "MD")) {
        // Parse model directory
        if (sscanf(bufferp, "%s%s", keyword, mesh_directory) != (unsigned int) 2) {
          fprintf(stderr, "Error parsing model directory at line %d of %s\n", line_number, filename);
          return 0;
        }
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3Scene::
WriteParseFile(const char *filename) const
{
  // Not implemented yet
  RNAbort("Not implemented");

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// RECTANGLE FILE I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3Scene::
ReadRectangleFile(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Read file
  char buffer[4096];
  int line_number = 0;
  int assignment_index = 0;
  while (fgets(buffer, 4096, fp)) {
    // Check line
    line_number++;
    char *bufferp = buffer;
    while (*bufferp && isspace(*bufferp)) bufferp++;
    if (!bufferp) continue;
    if (*bufferp == '#') continue;

    // Parse line
    char name[4096];
    double x1, y1, x2, y2, x3, y3, x4, y4, zmin, zmax, score;
    if (sscanf(bufferp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%s", 
       &x1, &y1, &x2, &y2, &x3, &y3, &x4, &y4, &zmin, &zmax, &score, name) != (unsigned int) 12) {
      fprintf(stderr, "Error parsing line %d of %s\n", line_number, filename);
      return 0;
    }

    // Create node
    R3SceneNode *node = new R3SceneNode(this);
    root->InsertChild(node);

    // Create shape element
    RNArray<R3Triangle *> triangles;
    RNArray<R3TriangleVertex *> vertices;
    R3TriangleVertex *v1 = new R3TriangleVertex(R3Point(x1, y1, zmin)); vertices.Insert(v1);
    R3TriangleVertex *v2 = new R3TriangleVertex(R3Point(x2, y2, zmin)); vertices.Insert(v2);
    R3TriangleVertex *v3 = new R3TriangleVertex(R3Point(x3, y3, zmin)); vertices.Insert(v3);
    R3TriangleVertex *v4 = new R3TriangleVertex(R3Point(x4, y4, zmin)); vertices.Insert(v4);
    R3Triangle *t1 = new R3Triangle(v1, v2, v3); triangles.Insert(t1);
    R3Triangle *t2 = new R3Triangle(v1, v3, v4); triangles.Insert(t2);
    R3TriangleArray *base = new R3TriangleArray(vertices, triangles);
    R3Point base_centroid = base->BBox().Centroid();
    R3Point top_centroid = base_centroid + (zmax - zmin) * R3posz_vector;
    R3Cylinder *marker = new R3Cylinder(base_centroid, top_centroid, 0.01 * base->BBox().DiagonalRadius());
    R3Brdf *brdf = new R3Brdf(RNRgb(0, 0.25 + score, 0));
    R3Material *material = new R3Material(brdf);
    R3SceneElement *element = new R3SceneElement();
    element->SetMaterial(material);
    element->InsertShape(base);
    element->InsertShape(marker);
    node->InsertElement(element);
    
    // Set node name
    char node_name[1024];
    sprintf(node_name, "A%d_S%03d", assignment_index++, (int) (1000*score));
    node->SetName(node_name);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}





