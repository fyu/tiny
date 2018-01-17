// Source file for the scene information program



// Include files 

#include "R3Graphics/R3Graphics.h"



// Program variables

static char *input_scene_name = NULL;
static int print_nodes = 0;
static int print_elements = 0;
static int print_lights = 0;



static R3Scene *
ReadScene(char *filename)
{
  // Allocate scene
  R3Scene *scene = new R3Scene();
  if (!scene) {
    fprintf(stderr, "Unable to allocate scene for %s\n", filename);
    return NULL;
  }

  // Read scene from file
  if (!scene->ReadFile(filename)) {
    delete scene;
    return NULL;
  }

  // Return scene
  return scene;
}



static int
PrintNodes(R3Scene *scene, R3SceneNode *node, int level)
{
  // Construct indent
  char indent[4096] = { '\0' };
  for (int i = 0; i < level; i++) {
    strcat(indent, "  ");
  }

  // Print node
  const char *name = node->Name();
  R3Box bbox = node->BBox();
  R3Point centroid = node->Centroid();
  R4Matrix m = node->Transformation().Matrix();
  printf("%sNode %s\n", indent, (name) ? name : "Null");
  printf("%s  Parent  = %s\n", indent, (node->Parent()) ? node->Parent()->Name() : "None");
  printf("%s  Parent index = %d\n", indent, node->ParentIndex());
  printf("%s  Scene index = %d\n", indent, node->SceneIndex());
  printf("%s  # Elements = %d\n", indent, node->NElements());
  printf("%s  BBox = ( %g %g %g ) ( %g %g %g )\n", indent, bbox[0][0], bbox[0][1], bbox[0][2], bbox[1][0], bbox[1][1], bbox[1][2]);
  printf("%s  Centroid = %g %g %g\n", indent, centroid[0], centroid[1], centroid[2]);
  printf("%s  Transformation matrix = \n", indent);
  printf("%s    %12.3g %12.3g %12.3g %12.3g\n", indent, m[0][0], m[0][1], m[0][2], m[0][3]);
  printf("%s    %12.3g %12.3g %12.3g %12.3g\n", indent, m[1][0], m[1][1], m[1][2], m[1][3]);
  printf("%s    %12.3g %12.3g %12.3g %12.3g\n", indent, m[2][0], m[2][1], m[2][2], m[2][3]);
  printf("%s    %12.3g %12.3g %12.3g %12.3g\n", indent, m[3][0], m[3][1], m[3][2], m[3][3]);
  printf("\n");

  // Print elements
  if (print_elements) {
  }

  // Print chldren
  for (int i = 0; i < node->NChildren(); i++) {
    R3SceneNode *child = node->Child(i);
    if (!PrintNodes(scene, child, level+1)) return 0;
  }

  // Return success
  return 1;
}



static int
PrintLights(R3Scene *scene)
{
  // Check lights
  if (scene->NLights() == 0) return 1;

  printf("Lights ...\n");
  for (int i = 0; i < scene->NLights(); i++) {
    R3Light *light = scene->Light(i);
    const RNRgb& color = light->Color();
    printf("  Light %d\n", i);
    printf("    Intensity = %g\n", light->Intensity());
    printf("    Color = %g %g %g\n", color.R(), color.G(), color.B());
    printf("\n");
  }

  // Return success
  return 1;
}



static int
PrintScene(R3Scene *scene)
{
  // Print scene stuff
  const R3Camera& camera = scene->Camera();
  const R3Box& bbox = scene->BBox();
  R3Point centroid = scene->Centroid();
  printf("Scene ...\n");
  printf("  # Nodes = %d\n", scene->NNodes());
  // printf("  # Elements = %d\n", scene->NElements());
  printf("  # Lights = %d\n", scene->NLights());
  printf("  Camera = %g %g %g   %g %g %g   %g %g %g  %g\n", 
    camera.Origin().X(), camera.Origin().Y(), camera.Origin().Z(),
    camera.Towards().X(), camera.Towards().Y(), camera.Towards().Z(),
    camera.Up().X(), camera.Up().Y(), camera.Up().Z(),
    camera.YFOV());
  printf("  BBox = ( %g %g %g ) ( %g %g %g )\n", 
    bbox[0][0], bbox[0][1], bbox[0][2], 
    bbox[1][0], bbox[1][1], bbox[1][2]);
  printf("  Centroid = %g %g %g\n", 
    centroid[0], centroid[1], centroid[2]);
  printf("\n");

  // Print nodes and elements
  if (print_nodes || print_elements) {
    PrintNodes(scene, scene->Root(), 0);
  }

  // Print lights
  if (print_lights) {
    PrintLights(scene);
  }

  // Return success
  return 1;
}



static int 
ParseArgs(int argc, char **argv)
{
  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-nodes")) print_nodes = 1; 
      else if (!strcmp(*argv, "-nodes")) print_nodes = 1; 
      else if (!strcmp(*argv, "-elements")) print_elements = 1; 
      else if (!strcmp(*argv, "-lights")) print_lights = 1; 
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
    else {
      if (!input_scene_name) input_scene_name = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check scene filename
  if (!input_scene_name) {
    fprintf(stderr, "Usage: scninfo inputscenefile [options]\n");
    return 0;
  }

  // Return OK status 
  return 1;
}



int main(int argc, char **argv)
{
  // Parse program arguments
  if (!ParseArgs(argc, argv)) exit(-1);

  // Read scene
  R3Scene *scene = ReadScene(input_scene_name);
  if (!scene) exit(-1);

  // Print scene
  if (!PrintScene(scene)) exit(-1);

  // Return success 
  return 0;
}

















