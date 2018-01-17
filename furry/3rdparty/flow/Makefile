# Makefile for flow evaluation code

SRC = flowIO.cpp colorcode.cpp colortest.cpp color_flow.cpp
BIN = colortest color_flow

IMGLIB = imageLib

CC = g++
WARN = -W -Wall
OPT ?= -O3
CPPFLAGS = $(OPT) $(WARN) -I$(IMGLIB)
LDLIBS = -L$(IMGLIB) -lImg -lpng -lz
EXE = $(SRC:.cpp=.exe)

all: $(BIN)

colortest: colortest.cpp colorcode.cpp
color_flow: color_flow.cpp flowIO.cpp colorcode.cpp

clean: 
	rm -f core *.stackdump

allclean: clean
	rm -f $(BIN) $(EXE)
