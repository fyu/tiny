// Include file for GSV classes

#ifndef __GSV__H__
#define __GSV__H__


// Dependency includes

#include "R3Shapes/R3Shapes.h"



// Constants

#define GSV_DRAW_POINTS_WITHOUT_COLOR                 0x00000000
#define GSV_DRAW_POINTS_WITH_LASER_INDEX_COLOR        0x00000001
#define GSV_DRAW_POINTS_WITH_VIEWPOINT_DISTANCE_COLOR 0x00000002
#define GSV_DRAW_POINTS_WITH_HEIGHT_COLOR             0x00000004
#define GSV_DRAW_SCANLINES_WITH_POINT_INDEX_COLOR     0x00000100
#define GSV_DRAW_MESHES_WITHOUT_COLOR                 0x00010000
#define GSV_DEFAULT_DRAW_FLAGS                        GSV_DRAW_POINTS_WITH_VIEWPOINT_DISTANCE_COLOR



// Class declarations

class GSVPose;
class GSVMesh;
class GSVScanline;
class GSVScan;
class GSVImage;
class GSVPanorama;
class GSVTapestry;
class GSVLaser;
class GSVCamera;
class GSVSegment;
class GSVRun;
class GSVScene;



// GSV includes

#include "GSVPose.h"
#include "GSVMesh.h"
#include "GSVScanline.h"
#include "GSVScan.h"
#include "GSVTapestry.h"
#include "GSVImage.h"
#include "GSVPanorama.h"
#include "GSVLaser.h"
#include "GSVCamera.h"
#include "GSVSegment.h"
#include "GSVRun.h"
#include "GSVScene.h"



#endif
