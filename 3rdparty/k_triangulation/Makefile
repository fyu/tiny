CC		= g++ -O3

CFLAGS          = -Wall \
		  -I./lib

LDFLAGS         = -L./lib

LIBS            = -lmatrix \
		  -lm

SRCS		= main.c \
		  twoview_triangulation.c \
		  misc.c \
		  getargs.c

HDRS		= twoview_triangulation.h \
		  misc.h \
		  getargs.h

OBJS		= $(SRCS:.c=.o)

PROGRAM		= triangulation2

all:		$(PROGRAM)

$(PROGRAM):	$(OBJS) $(HDRS) 
		$(CC) $(OBJS) $(LDFLAGS) $(LIBS) -o $(PROGRAM)

clean:;		rm -f *.o *~ $(PROGRAM)
