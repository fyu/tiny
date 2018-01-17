// Source file for the google info program



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "GSV/GSV.h"



////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////

// Program arguments

static const char *input_scene_name = NULL;
static int print_runs = 0;
static int print_cameras = 0;
static int print_lasers = 0;
static int print_segments = 0;
static int print_scans = 0;
static int print_scanlines = 0;
static int print_points = 0;
static int print_tapestries = 0;
static int print_panoramas = 0;
static int print_images = 0;
static int print_verbose = 0;
static int print_debug = 0;



////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

static GSVScene *
ReadScene(const char *filename)
{
  // Allocate scene
  GSVScene *scene = new GSVScene();
  if (!scene) {
    fprintf(stderr, "Unable to allocate scene\n");
    return NULL;
  }

  // Open google scene files
  if (!scene->ReadFile(filename, FALSE)) {
    delete scene;
    return NULL;
  }

  // Return scene
  return scene;
}



////////////////////////////////////////////////////////////////////////
// Print functions
////////////////////////////////////////////////////////////////////////

static void 
PrintDebug(GSVScene *scene)
{
#if 0
  // Write tapestry images
  for (int r = 0; r < scene->NRuns(); r++) {
    GSVRun *run = scene->Run(r);
    for (int s = 0; s < run->NSegments(); s++) {
      GSVSegment *segment = run->Segment(s);
      for (int it = 2; it <= 6; it += 4) {
        GSVTapestry *tapestry = segment->Tapestry(tt);
        GSVCamera *camera = tapestry->Camera();

        // Create output image
        int width = camera->ImageWidth();
        int height = camera->ImageHeight();
        const int max_columns = 16 * 1024;
        int total_columns = width * tapestry->NImages();
        if (total_columns > max_columns) total_columns = max_columns;
        R2Image output_image(total_columns, height, 3);

        // Fill output image
        int column_index = 0;
        for (int ii = 0; ii < tapestry->NImages(); ii++) {
          GSVImage *image = tapestry->Image(ii);
          GSVPanorama *panorama = image->Panorama();

          // Construct image name
          char input_image_name[1024];
          sprintf(input_image_name, "%s/segment_%02d/unstitched_%06d_%02d.jpg", 
                  run->Name(), segment->RunIndex(), panorama->SegmentIndex(), image->PanoramaIndex());
          
          // Read input image
          R2Image input_image;
          if (!input_image.Read(input_image_name)) {
            fprintf(stderr, "Unable to read image file %s\n", input_image_name);
            return;
          }
          
          // Copy pixels from input image to output image
          for (int x = width-1; x >= 0; x--) {
            for (int y = 0; y < height; y++) {
              RNRgb rgb = input_image.PixelRGB(x, y);
              output_image.SetPixelRGB(column_index, y, rgb);
            }
            column_index++;
            if (column_index >= max_columns) break;
          }
        }

        // Write output image
        char output_image_name[1024];
        sprintf(output_image_name, "t%d.jpg", t);      
        output_image.Write(output_image_name);
      }
    }
  }
#endif
}



static void 
PrintScene(GSVScene *scene)
{
  // Print scene stuff
  R3Box scene_bbox = scene->BBox();
  const char *scene_filename = scene->Filename();
  if (!scene_filename) scene_filename = "None";
  printf("Scene ...\n");
  printf("  # Runs = %d\n", scene->NRuns());
  printf("  # Cameras = %d\n", scene->NCameras());
  printf("  # Lasers = %d\n", scene->NLasers());
  printf("  # Segments = %d\n", scene->NSegments());
  printf("  # Scans = %d\n", scene->NScans());
  printf("  # Tapestries = %d\n", scene->NTapestries());
  printf("  # Panoramas = %d\n", scene->NPanoramas());
  printf("  # Images = %d\n", scene->NImages());
  printf("  # Scanlines = %d\n", scene->NScanlines());
  printf("  Filename = %s\n", scene_filename);
  printf("  BBox = ( %g %g %g ) ( %g %g %g )\n", 
    scene_bbox[0][0], scene_bbox[0][1], scene_bbox[0][2], 
    scene_bbox[1][0], scene_bbox[1][1], scene_bbox[1][2]);
  printf("\n");

  // Print run stuff
  if (print_runs || print_cameras || print_lasers || print_segments ||
      print_scans || print_scanlines || print_points || 
      print_tapestries || print_panoramas || print_images) {
    for (int r = 0; r < scene->NRuns(); r++) {
      GSVRun *run = scene->Run(r);
      R3Box run_bbox = run->BBox();
      printf("  Run %d ...\n", r);
      printf("    Name = %s\n", run->Name());
      printf("    Scene index = %d\n", run->SceneIndex());
      printf("    # Caneras = %d\n", run->NCameras());
      printf("    # Lasers = %d\n", run->NLasers());
      printf("    # Segments = %d\n", run->NSegments());
      printf("    BBox = ( %g %g %g ) ( %g %g %g )\n", 
        run_bbox[0][0], run_bbox[0][1], run_bbox[0][2], 
        run_bbox[1][0], run_bbox[1][1], run_bbox[1][2]);
      printf("\n");

      // Print cameras
      if (print_cameras) {
        for (int m = 0; m < run->NCameras(); m++) {
          GSVCamera *camera = run->Camera(m);
          R3Box camera_bbox = camera->BBox();
          printf("    Camera %d ...\n", m);
          printf("      Scene index = %d\n", camera->SceneIndex());
          printf("      Run index = %d\n", camera->RunIndex());
          printf("      # Tapestries = %d\n", camera->NTapestries());
          printf("      DistortionType = %d\n", camera->DistortionType());
          printf("      XFocal = %g\n", camera->XFocal());
          printf("      YFocal = %g\n", camera->YFocal());
          printf("      XCenter = %g\n", camera->XCenter());
          printf("      YCenter = %g\n", camera->YCenter());
          printf("      K1 = %g\n", camera->K1());
          printf("      K2 = %g\n", camera->K2());
          printf("      K3 = %g\n", camera->K3());
          printf("      P1 = %g\n", camera->P1());
          printf("      P2 = %g\n", camera->P2());
          printf("      MaxFov = %g\n", camera->MaxFov());
          printf("      BBox = ( %g %g %g ) ( %g %g %g )\n", 
            camera_bbox[0][0], camera_bbox[0][1], camera_bbox[0][2], 
            camera_bbox[1][0], camera_bbox[1][1], camera_bbox[1][2]);
          printf("\n");
        }
      }

      // Print lasers
      if (print_lasers) {
        for (int e = 0; e < run->NLasers(); e++) {
          GSVLaser *laser = run->Laser(e);
          R3Box laser_bbox = laser->BBox();
          printf("    Laser %d ...\n", e);
          printf("      Scene index = %d\n", laser->SceneIndex());
          printf("      Run index = %d\n", laser->RunIndex());
          printf("      # Scans = %d\n", laser->NScans());
          printf("      BBox = ( %g %g %g ) ( %g %g %g )\n", 
            laser_bbox[0][0], laser_bbox[0][1], laser_bbox[0][2], 
            laser_bbox[1][0], laser_bbox[1][1], laser_bbox[1][2]);
          printf("\n");
        }
      }

      // Print segments
      if (print_segments || print_scans || print_scanlines || print_points ||
          print_tapestries || print_panoramas || print_images) {
        for (int s = 0; s < run->NSegments(); s++) {
          GSVSegment *segment = run->Segment(s);
          R3Box segment_bbox = segment->BBox();
          printf("    Segment %d ...\n", s);
          printf("      Scene index = %d\n", segment->SceneIndex());
          printf("      Run index = %d\n", segment->RunIndex());
          printf("      # Scans = %d\n", segment->NScans());
          printf("      # Tapestries = %d\n", segment->NTapestries());
          printf("      # Panorama = %d\n", segment->NPanoramas());
          printf("      BBox = ( %g %g %g ) ( %g %g %g )\n", 
            segment_bbox[0][0], segment_bbox[0][1], segment_bbox[0][2], 
            segment_bbox[1][0], segment_bbox[1][1], segment_bbox[1][2]);
          printf("\n");

          // Print scans
          if (print_scans || print_scanlines || print_points) {
            for (int a = 0; a < segment->NScans(); a++) {
              GSVScan *scan = segment->Scan(a);
              int scan_laser_index = (scan->Laser()) ? scan->Laser()->RunIndex() : -1;
              R3Box scan_bbox = scan->BBox();
              printf("      Scan %d ...\n", a);
              printf("        Scene index = %d\n", scan->SceneIndex());
              printf("        Run index = %d\n", scan->RunIndex());
              printf("        Segment index = %d\n", scan->SegmentIndex());
              printf("        Laser = %d\n", scan_laser_index);
              printf("        # Scanlines = %d\n", scan->NScanlines());
              printf("        BBox = ( %g %g %g ) ( %g %g %g )\n", 
                scan_bbox[0][0], scan_bbox[0][1], scan_bbox[0][2], 
                scan_bbox[1][0], scan_bbox[1][1], scan_bbox[1][2]);
              printf("\n");

              // Print scanlines
              if (print_scanlines || print_points) {
                // Read points
                if (print_points) scan->ReadPoints();

                // Print scanlines
                for (int n = 0; n < scan->NScanlines(); n++) {
                  GSVScanline *scanline = scan->Scanline(n);
                  R3Point scanline_viewpoint = scanline->Pose().Viewpoint();
                  R3Quaternion scanline_orientation = scanline->Pose().Orientation();
                  R3Box scanline_bbox = scanline->BBox();
                  printf("        Scanline %d ...\n", n);
                  printf("          Scene index = %d\n", scanline->SceneIndex());
                  printf("          Run index = %d\n", scanline->RunIndex());
                  printf("          Segment index = %d\n", scanline->SegmentIndex());
                  printf("          Scan index = %d\n", scanline->ScanIndex());
                  printf("          # Points = %d\n", scanline->NPoints());
                  printf("          Timestamp = %g\n", scanline->Timestamp());
                  printf("          Viewpoint = %g %g %g\n", 
                    scanline_viewpoint[0], scanline_viewpoint[1], scanline_viewpoint[2]);
                  printf("          Orientation = %g %g %g %g\n", 
                    scanline_orientation[0], scanline_orientation[1], scanline_orientation[2], scanline_orientation[3]);
                  printf("          BBox = ( %g %g %g ) ( %g %g %g )\n", 
                    scanline_bbox[0][0], scanline_bbox[0][1], scanline_bbox[0][2], 
                    scanline_bbox[1][0], scanline_bbox[1][1], scanline_bbox[1][2]);
                  if (print_points && (scanline->NPoints() > 0)) {
                    printf("          Points:\n");
                    for (int p = 0; p < scanline->NPoints(); p++) {
                      R3Point position = scanline->PointPosition(p);
                      printf("            %g %g %g\n", position[0], position[1], position[2]);
                    }
                  }
                  printf("\n");
                }

                // Release points
                if (print_points) scan->ReleasePoints();
              }
            }
          }

          // Print tapestries
          if (print_tapestries) {
            for (int t = 0; t < segment->NTapestries(); t++) {
              GSVTapestry *tapestry = segment->Tapestry(t);
              int tapestry_camera_index = (tapestry->Camera()) ? tapestry->Camera()->RunIndex() : -1;
              R3Box tapestry_bbox = segment->BBox();
              printf("      Tapestry %d ...\n", t);
              printf("        Scene index = %d\n", tapestry->SceneIndex());
              printf("        Run index = %d\n", tapestry->RunIndex());
              printf("        Segment index = %d\n", tapestry->SegmentIndex());
              printf("        Camera = %d\n", tapestry_camera_index);
              printf("        # Images = %d\n", tapestry->NImages());
              printf("        BBox = ( %g %g %g ) ( %g %g %g )\n", 
                tapestry_bbox[0][0], tapestry_bbox[0][1], tapestry_bbox[0][2], 
                tapestry_bbox[1][0], tapestry_bbox[1][1], tapestry_bbox[1][2]);
              printf("\n");
            }
          }

          // Print panoramas
          if (print_panoramas || print_images) {
            for (int p = 0; p < segment->NPanoramas(); p++) {
              GSVPanorama *panorama = segment->Panorama(p);
              R3Point panorama_viewpoint = panorama->Viewpoint();
              printf("      Panorama %d ...\n", p);
              printf("        Scene index = %d\n", panorama->SceneIndex());
              printf("        Run index = %d\n", panorama->RunIndex());
              printf("        Segment index = %d\n", panorama->SegmentIndex());
              printf("        # Images = %d\n", panorama->NImages());
              printf("        Timestamp = %g\n", panorama->Timestamp());
              printf("        Viewpoint = %g %g %g\n", 
                panorama_viewpoint[0], panorama_viewpoint[1], panorama_viewpoint[2]);
              printf("\n");

              // Print images
              if (print_images) {
                for (int g = 0; g < panorama->NImages(); g++) {
                  GSVImage *image = panorama->Image(g);
                  int image_tapestry_index = (image->Tapestry()) ? image->Tapestry()->SegmentIndex() : -1;
                  R3Point image_viewpoint = image->Pose().Viewpoint();
                  R3Quaternion image_orientation = image->Pose().Orientation();
                  printf("        Image %d ...\n", g);
                  printf("          Scene index = %d\n", image->SceneIndex());
                  printf("          Run index = %d\n", image->RunIndex());
                  printf("          Segment index = %d\n", image->SegmentIndex());
                  printf("          Panorama index = %d\n", image->PanoramaIndex());
                  printf("          Tapestry = %d\n", image_tapestry_index);
                  printf("          Width = %d\n", image->Width());
                  printf("          Height = %d\n", image->Height());
                  printf("          Timestamp = %g\n", image->Timestamp());
                  printf("          Viewpoint = %g %g %g\n", 
                    image_viewpoint[0], image_viewpoint[1], image_viewpoint[2]);
                  printf("          Orientation = %g %g %g %g\n", 
                    image_orientation[0], image_orientation[1], image_orientation[2], image_orientation[3]);
                  printf("\n");
                }
              }
            }
          }
        }
      }
    }
  }
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
      else if (!strcmp(*argv, "-runs")) print_runs = 1;
      else if (!strcmp(*argv, "-cameras")) print_cameras = 1;
      else if (!strcmp(*argv, "-lasers")) print_lasers = 1;
      else if (!strcmp(*argv, "-segments")) print_segments = 1;
      else if (!strcmp(*argv, "-scans")) print_scans = 1;
      else if (!strcmp(*argv, "-scanlines")) print_scanlines = 1;
      else if (!strcmp(*argv, "-points")) print_points = 1;
      else if (!strcmp(*argv, "-tapestries")) print_tapestries = 1;
      else if (!strcmp(*argv, "-panoramas")) print_panoramas = 1;
      else if (!strcmp(*argv, "-images")) print_images = 1;
      else if (!strcmp(*argv, "-all")) { 
        print_runs = 1;
        print_cameras = 1;
        print_lasers = 1;
        print_segments = 1;
        print_scans = 1;
        print_scanlines = 1;
        print_points = 1;
        print_tapestries = 1;
        print_panoramas = 1;
        print_images = 1;
      }
      else { 
        fprintf(stderr, "Invalid program argument: %s", *argv); 
        exit(1); 
      }
      argv++; argc--;
    }
    else {
      if (!input_scene_name) input_scene_name = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check scene name
  if (!input_scene_name) {
    fprintf(stderr, "Usage: gsvinfo scenefile [options]\n");
    return FALSE;
  }

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

  // Read scene
  GSVScene *scene = ReadScene(input_scene_name);
  if (!scene) exit(-1);

  // Print scene information
  PrintScene(scene);

  // Print debug information
  if (print_debug) PrintDebug(scene);

  // Return success 
  return 0;
}

















