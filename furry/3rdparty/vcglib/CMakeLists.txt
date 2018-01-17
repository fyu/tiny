include_directories(${PROJECT_SOURCE_DIR}/furry/3rdparty/vcglib)

add_library(furry_3rdparty_vcg
			wrap/ply/plylib.cpp
			wrap/gui/trackball.cpp
			wrap/gui/trackmode.cpp)


