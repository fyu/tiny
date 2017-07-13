include_directories("/opt/local/include/")
add_library(furry_3rdparty_flow
						flowIO.cpp
						colorcode.cpp
						Convert.cpp
						Convolve.cpp
						Image.cpp
						ImageIO.cpp
						ImageIOpng.cpp
						RefCntMem.cpp)

target_link_libraries(furry_3rdparty_flow png z)
