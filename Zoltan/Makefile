
SYSTYPE = "CRAY"
ZOLTANDIR = /apps/zoltan/3.83/INTEL/16.0
SPRNGDIR = /apps/sprng/2.0/INTEL/16.0

ifeq ($(SYSTYPE),"CRAY")
CC = cc
OPT = -O3 -qopenmp
INCLUDES = -I. -I$(ZOLTANDIR)/include -I$(SPRNGDIR)/include
LIBS = -lm -L$(ZOLTANDIR)/lib -lzoltan -L$(SPRNGDIR)/lib -lsprng -lgmp
endif

EXEC = particles

OPTIONS = $(OPT)

OBJS = main.o init.o decompose.o io.o

INCL = particles.h

CFLAGS = $(OPTIONS) $(INCLUDES) 

$(EXEC): $(OBJS)
	$(CC) -dynamic $(OPTIONS) $(OBJS) $(LIBS) -o $(EXEC)

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC) results/*


