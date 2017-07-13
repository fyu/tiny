cmake_minimum_required(VERSION 2.6)

project(tiny)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/${CMAKE_SYSTEM_NAME})

if (APPLE)
set(CMAKE_Fortran_COMPILER "gfortran")
set(CMAKE_CXX_COMPILER "clang++")

#set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wno-reorder -Wno-enum-compare -frounding-math")
set(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++ -Wall -Wno-reorder -Wno-mismatched-tags")
set(CMAKE_EXE_LINKER_FLAGS "-stdlib=libc++")
set(CMAKE_SHARED_LINKER_FLAGS "-stdlib=libc++")
set(CMAKE_MODULE_LINKER_FLAGS "-stdlib=libc++")
set(CMAKE_Fortran_FLAGS "-O3 -ffast-math")
else (APPLE)
# set(CMAKE_CXX_COMPILER "clang++")
 set(CMAKE_CXX_COMPILER "/u/fy/usr/gcc-4.8.2/bin/g++")
#set(CMAKE_CXX_COMPILER "g++")
#set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wno-reorder -Wno-overloaded-virtual -Wno-unused-private-field -Wno-uninitialized -Wno-mismatched-tags")
set(CMAKE_CXX_FLAGS "-std=c++11 -Wno-reorder -fpermissive -fopenmp")
# set(GHOME /n/fs/graphicslab/fy/usr/)
# set(GHOME /usr/local/)
set(GHOME $ENV{GHOME}/usr/)
message(STATUS "GHOME: " $ENV{GHOME})
message(STATUS "GAPS: " $ENV{GAPS_PATH})
set(CMAKE_LIBRARY_PATH ${GHOME}/lib/ ${CMAKE_LIBRARY_PATH})
endif (APPLE)

#### gaps warning ####
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sometimes-uninitialized -Wno-overloaded-virtual")

if (CMAKE_BUILD_TYPE STREQUAL "Release")
   message("*** Release ***")
   add_definitions("-DEIGEN_NO_DEBUG")
endif()

# add_definitions("-DVCG_USE_EIGEN")

if (APPLE)
else (APPLE)
set(CMAKE_INSTALL_PREFIX ${GHOME})
endif (APPLE)

#add_definitions("-DCGAL_CFG_NO_CPP0X_ARRAY=1")

find_package(PkgConfig)
set(GFLAGS_LIBRARIES gflags)

#pkg_check_modules(GSL gsl)

# message(STATUS ${OPENCV_LIBRARY_DIR})
# message(STATUS ${GFLAGS_LDFLAGS})
# message(STATUS ${GFLAGS_LIBRARY_DIRS})
# message(STATUS ${GFLAGS_LIBRARIES})
# message(STATUS ${GSL_LIBRARY_DIRS})


if (APPLE)
#set(Boost_INCLUDE_DIR /Users/fy/install/boost_1_51_0/)
set(Boost_INCLUDE_DIR /usr/local/include/)
#set(Boost_LIBRARY_DIR /Users/fy/install/boost_1_51_0/stage/lib)
set(Boost_LIBRARY_DIR /usr/local/lib)
set(Boost_LIBRARIES ${Boost_LIBRARY_DIR}/libboost_serialization.dylib
					${Boost_LIBRARY_DIR}/libboost_thread-mt.dylib
					${Boost_LIBRARY_DIR}/libboost_filesystem.dylib
					${Boost_LIBRARY_DIR}/libboost_system.dylib)


else (APPLE)
set(Boost_INCLUDE_DIR ${GHOME}/include)
set(Boost_LIBRARY_DIR ${GHOME}/lib)
#set(Boost_LIBRARY_DIR /usr/lib64)
set(Boost_LIBRARIES ${Boost_LIBRARY_DIR}/libboost_serialization.so
#${Boost_LIBRARY_DIR}/libboost_chrono.so
#${Boost_LIBRARY_DIR}/libboost_context.so
${Boost_LIBRARY_DIR}/libboost_date_time.so
${Boost_LIBRARY_DIR}/libboost_filesystem.so
${Boost_LIBRARY_DIR}/libboost_graph.so
#${Boost_LIBRARY_DIR}/libboost_iostreams.so
#${Boost_LIBRARY_DIR}/libboost_locale.so
#${Boost_LIBRARY_DIR}/libboost_math_c99f.so
#${Boost_LIBRARY_DIR}/libboost_math_c99l.so
#${Boost_LIBRARY_DIR}/libboost_math_c99.so
#${Boost_LIBRARY_DIR}/libboost_math_tr1f.so
#${Boost_LIBRARY_DIR}/libboost_math_tr1l.so
#${Boost_LIBRARY_DIR}/libboost_math_tr1.so
${Boost_LIBRARY_DIR}/libboost_prg_exec_monitor.so
${Boost_LIBRARY_DIR}/libboost_program_options.so
#${Boost_LIBRARY_DIR}/libboost_python.so
#${Boost_LIBRARY_DIR}/libboost_random.so
${Boost_LIBRARY_DIR}/libboost_regex.so
${Boost_LIBRARY_DIR}/libboost_serialization.so
${Boost_LIBRARY_DIR}/libboost_signals.so
${Boost_LIBRARY_DIR}/libboost_system.so
${Boost_LIBRARY_DIR}/libboost_thread.so
#${Boost_LIBRARY_DIR}/libboost_timer.so
${Boost_LIBRARY_DIR}/libboost_unit_test_framework.so
${Boost_LIBRARY_DIR}/libboost_wave.so)

message(STATUS "boost include: " ${Boost_INCLUDE_DIR})
# message(STATUS "boost: " ${Boost_LIBRARIES})

endif (APPLE)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

if (APPLE)
set(EIGEN3_INCLUDE_DIR /usr/local/include/eigen3/)
set(EIGEN_INCLUDE_DIRS /usr/local/include/eigen3/)
set(EIGEN_INCLUDE_DIR /usr/local/include/eigen3/)
else ()
set(EIGEN3_INCLUDE_DIR /usr/include/eigen3)
endif()

pkg_check_modules(OPENCV opencv)
set(OPENCV_INCLUDE_DIR ${OPENCV_INCLUDE_DIR})
set(OPENCV_LIBRARIES opencv_core opencv_calib3d opencv_imgproc opencv_highgui opencv_features2d opencv_flann opencv_nonfree opencv_video)
if (APPLE)
set(OPENCV_LIBRARY_DIR ${OPENCV_LIBRARY_DIR})
else ()
set(OPENCV_LIBRARY_DIR /u/fy/usr/lib/)
endif ()

# message(STATUS "opencv library dir: " ${OPENCV_LIBRARIES})

# find_package(GMP)

find_package(SuiteSparse)
message(STATUS "find suitesparse: "  ${SUITESPARSE_FOUND})
message(STATUS "suitesparse libraries: " ${SUITESPARSE_LIBRARIES})

set(MINPACK_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/3rdparty/minpack/)
set(CMINPACK_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/3rdparty/cminpack/)

if (APPLE)
find_library(OPENGL_LIBRARY opengl)
find_library(GLUT_LIBRARY GLUT)
else (APPLE)
set(OPENGL_LIBRARY GLU GL m)
set(GLUT_LIBRARY glut)
endif(APPLE)

set(GRAPHICS_LIBRARIES furry_3rdparty_glew ${OPENGL_LIBRARY} ${GLUT_LIBRARY})
# set(GRAPHICS_LIBRARIES ${OPENGL_LIBRARY} ${GLUT_LIBRARY})
set(GL_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/3rdparty/glew-1.10.0/include)

############################################################
## Ceres
############################################################

find_package(Ceres REQUIRED)

############################################################
## Protobuf
############################################################

find_package(Protobuf REQUIRED)

############################################################
## gaps
############################################################

set(GAPS_PATH $ENV{GAPS_PATH})

set(GAPS_LIBRARY_DIR ${GAPS_PATH}/lib/x86_64)
set(GAPS_INCLUDE_DIR ${GAPS_PATH}/pkgs)
message(STATUS "gaps include dir:" ${GAPS_INCLUDE_DIR})
set(GAPS_LIBRARIES fglut ${OPENGL_LIBRARY} CSparse GSV R3Graphics R3Shapes R2Shapes R3Surfels RNMath RNBasics png jpeg)

############################################################
## TBB
############################################################

set(TBB_LIBRARIES tbb)

############################################################
## Optimization Library Settings
############################################################

if (APPLE)
set(OPTIMIZATION_LIBRARIES ${CERES_LIBRARIES} ${SUITESPARSE_LIBRARIES} amd gmp mpfr glog ${GFLAGS_LIBRARIES} lapack blas protobuf)
else (APPLE)
set(OPT_LIBRARY_DIR /u/fy/usr/lib)
set(OPTIMIZATION_LIBRARIES ${CERES_LIBRARIES} ${SUITESPARSE_LIBRARIES} cholmod amd gomp
# ${GHOME}/lib/libgomp.so 
${GHOME}/lib/libmpfr.so libsuitesparseconfig.a  
glog ${GFLAGS_LIBRARIES} lapack blas protobuf)
endif (APPLE)

############################################################
## QT
############################################################

set(QT_QMAKE_EXECUTABLE /usr/local/bin/qmake)
find_package(Qt4 4.6 COMPONENTS QtCore QtGui QtTest QtUiTools QtXml QtOpenGL REQUIRED)
include(${QT_USE_FILE})
add_definitions(${QT_DEFINITIONS})


include_directories(
					"/usr/local/include/"
					#"/usr/include/suitesparse/"
					${PROJECT_SOURCE_DIR}/../
					${PROJECT_SOURCE_DIR}/../../
					${GL_INCLUDE_DIR}
					${EIGEN3_INCLUDE_DIR}
					${GFLAGS_INCLUDE_DIR}
					${OPENCV_INCLUDE_DIR}
					${Boost_INCLUDE_DIR}
					${GMP_INCLUDE_DIR}
					${GAPS_INCLUDE_DIR}
					${TBB_INCLUDE_DIR}
					${SUITESPARSE_INCLUDE_DIRS}
					${PROTOBUF_INCLUDE_DIRS})


message(STATUS ${PROJECT_SOURCE_DIR})
add_subdirectory(3rdparty)
add_subdirectory(common)
add_subdirectory(gui)
add_subdirectory(src)