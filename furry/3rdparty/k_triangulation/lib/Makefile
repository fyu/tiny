CC		= g++
CFLAGS          = -O3

OBJS		= bidiag.o \
		  eigen_okumura.o \
		  extensions.o \
		  givens.o \
		  hh_solve.o \
		  householder.o \
		  matrix.o \
		  matutil.o \
		  svd.o \
		  svdstep.o \
		  vector.o

libmatrix.a : $(OBJS)
	      ar r $@ $(OBJS)
	      - ranlib $@

clean:;	      rm -f $(OBJS)
