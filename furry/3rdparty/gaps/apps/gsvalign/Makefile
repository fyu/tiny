#
# Application name and list of source files.
#

NAME=gsvalign
CCSRCS=$(NAME).cpp poseoptimization.cpp R3CatmullRomSpline.cpp


#
# Dependency libraries
#

USER_CFLAGS=-DRN_USE_CSPARSE -DRN_NO_SPLM
PKG_LIBS=-lGSV -lR3Graphics -lR3Shapes -lR2Shapes -lRNMath -lRNBasics -lCSparse -ljpeg -lpng -lklt




#
# Ceres
#

USER_CFLAGS=-Ic:/Funk/downloads/ceres/ceres-solver-1.3.0/include -I/usr/local/include/eigen3 -DRN_USE_CERES -I../../pkgs
PKG_LIBS=-lGSV -lR3Graphics -lR3Shapes -lR2Shapes -lRNMath -lRNBasics -ljpeg -lpng -lceres -lsplm -lCSparse -lminpack -lccolamd -lcamd -lcolamd -lamd -lsuitesparseconfig -lcholmod -lccolamd -lcamd -lcolamd -lamd -llapack -lblas -lcxsparse  /usr/local/lib/libgflags.a /usr/local/lib/libglog.a -lprotobuf -lklt 



#
# R3 application makefile
#

include ../../makefiles/Makefile.apps



