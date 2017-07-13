link_directories(/usr/local/lib/ ${TBB_LIBRARY_DIR})

# find_package(Qt4 4.8 COMPONENTS QtCore QtGui Qt3Support QtSvg QtScript QtTest QtUiTools
 #                QtHelp QtWebKit QtXmlPatterns   QtNetwork
  #                QtOpenGL QtSql QtXml QtDesigner  REQUIRED)
# find_package(Qt4 4.8 COMPONENTS QtCore QtGui QtTest QtUiTools QtXml QtOpenGL QtHelp QtDBus QtDesigner QtScript QtScriptTools QtNetwork QtXmlPatterns REQUIRED)

# Looking for generated protobuf headers
include_directories(${CMAKE_CURRENT_BINARY_DIR})

include_directories(${PROJECT_SOURCE_DIR}/3rdparty/vcglib/)

qt4_add_resources(WINDOW_QT_RCCS window_QT.qrc)
qt4_wrap_cpp(WINDOW_QT_MOC window_QT.h)

protobuf_generate_cpp(BA_SETTINGS_PROTO_SRC BA_SETTINGS_PROTO_HDR ba_settings.proto)
protobuf_generate_cpp(MRF_SETTINGS_PROTO_SRC MRF_SETTINGS_PROTO_HDR mrf_settings.proto)

add_library(model model.cc)
target_link_libraries(model furry_common ${OPENCV_LIBRARIES} glog gflags)
add_library(mvs_util cost_volume.cc depth_map.cc plane_sweeping.cc)
#add_library(mvs_ui ${WINDOW_QT_MOC} ${WINDOW_QT_RCCS} window_QT.cpp util_ui.cc)
add_library(mvs_ui util_ui.cc)
target_link_libraries(mvs_ui furry_gui_cvwindow)
add_library(tiny_bundle track.cc point_relation.cc numeric.cc tiny_bundle.cc tiny_bundle_one.cc tiny_bundle_util.cc ${BA_SETTINGS_PROTO_SRC} ${BA_SETTINGS_PROTO_HDR})

add_library(mvs_opt mvs_opt_disc.cc mvs_opt_cont.cc densecrf.cc ${MRF_SETTINGS_PROTO_SRC} ${MRF_SETTINGS_PROTO_HDR})
#target_link_libraries(mvs_opt furry_3rdparty_densecrf)

qt4_wrap_cpp(MVS_UI_MOC mvs_ui.h)
add_executable(mvs mvs.cc mvs_impl.cc ${MVS_UI_MOC} mvs_ui.cc parallel.cc)
target_link_libraries(mvs model mvs_ui mvs_util furry_common  ${GRAPHICS_LIBRARIES} ${OPENCV_LIBRARIES} glog ${GFLAGS_LIBRARIES} ${TBB_LIBRARIES} ${QT_LIBRARIES})

add_executable(tiny_solve_depth solve_depth.cc smooth.cc)
target_link_libraries(tiny_solve_depth mvs_opt mvs_util furry_3rdparty_ply furry_3rdparty_mrf furry_3rdparty_densecrf furry_common protobuf ${OPENCV_LIBRARIES} ${OPTIMIZATION_LIBRARIES})

add_executable(convert_ims_cost_volume convert_ims_cost_volume.cc)
target_link_libraries(convert_ims_cost_volume mvs_util furry_common ${OPENCV_LIBRARIES} glog gflags)

add_executable(tiny_bundle_adjust bundle_adjust.cc)
target_link_libraries(tiny_bundle_adjust model mvs_ui tiny_bundle furry_common  ${GRAPHICS_LIBRARIES} ${OPENCV_LIBRARIES} glog ${GFLAGS_LIBRARIES} ${TBB_LIBRARIES} ${QT_LIBRARIES} ${OPTIMIZATION_LIBRARIES} protobuf)

add_executable(depth_uncertainty depth_uncertainty.cc)
target_link_libraries(depth_uncertainty model tiny_bundle furry_common  ${GRAPHICS_LIBRARIES} ${OPENCV_LIBRARIES} glog ${GFLAGS_LIBRARIES} ${TBB_LIBRARIES} ${QT_LIBRARIES} ${OPTIMIZATION_LIBRARIES} protobuf)

add_executable(tiny_detect_tracks detect_tracks.cc)
target_link_libraries(tiny_detect_tracks model mvs_ui furry_common ${OPENCV_LIBRARIES} glog ${GFLAGS_LIBRARIES} protobuf ${QT_LIBRARIES})

add_executable(detect_all_tracks detect_all_tracks.cc)
target_link_libraries(detect_all_tracks model mvs_ui furry_common ${OPENCV_LIBRARIES} glog ${GFLAGS_LIBRARIES} protobuf ${QT_LIBRARIES})

add_executable(stereo stereo_cli.cc)
target_link_libraries(stereo model mvs_ui furry_common ${QT_LIBRARIES} ${GRAPHICS_LIBRARIES} ${OPENCV_LIBRARIES} glog ${GFLAGS_LIBRARIES} ${TBB_LIBRARIES})

add_executable(flow flow.cc)
target_link_libraries(flow model mvs_ui furry_common ${QT_LIBRARIES} ${GRAPHICS_LIBRARIES} ${OPENCV_LIBRARIES} glog ${GFLAGS_LIBRARIES} ${TBB_LIBRARIES})

add_executable(tiny_klt klt.cc)
target_link_libraries(tiny_klt model mvs_ui furry_common ${QT_LIBRARIES} ${GRAPHICS_LIBRARIES} ${OPENCV_LIBRARIES} glog ${GFLAGS_LIBRARIES} ${TBB_LIBRARIES})

add_executable(model2ply model2ply.cc)
target_link_libraries(model2ply model furry_common ${OPENCV_LIBRARIES} glog ${GFLAGS_LIBRARIES})

add_executable(model2pmvs model2pmvs.cc)
target_link_libraries(model2pmvs model furry_common ${OPENCV_LIBRARIES} glog ${GFLAGS_LIBRARIES})

add_executable(tiny_sweep_planes sweep_planes.cc)
target_link_libraries(tiny_sweep_planes model mvs_util furry_common ${OPENCV_LIBRARIES} glog ${GFLAGS_LIBRARIES})

add_executable(tiny_test_segmentation test_segmentation.cc)
target_link_libraries(tiny_test_segmentation furry_common furry_3rdparty_slic ${OPENCV_LIBRARIES} glog gflags)

add_executable(tiny_calc_parallax calc_parallax.cc)
target_link_libraries(tiny_calc_parallax model furry_common ${OPENCV_LIBRARIES} glog gflags)

qt4_wrap_cpp(MODEL_VIEWER_MOC model_viewer.h)
add_executable(tiny_view_model view_model.cc model_viewer.cc ${MODEL_VIEWER_MOC})
target_link_libraries(tiny_view_model model furry_3rdparty_vcg ${QT_LIBRARIES} ${GRAPHICS_LIBRARIES})

# add_executable(syn_scene syn_scene.cc)
# target_link_libraries(syn_scene mvs_util furry_common ${QT_LIBRARIES} ${GRAPHICS_LIBRARIES} glog ${GFLAGS_LIBRARIES} ${TBB_LIBRARIES})

# install(TARGETS
# 		tiny_detect_tracks
# 		tiny_bundle_adjust
# 		RUNTIME DESTINATION bin/furry)
