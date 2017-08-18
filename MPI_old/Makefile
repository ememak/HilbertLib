HilbertLibVT: HilbertLib.c
	/home/administrator/vt/bin/vtcc -vt:cc cc -o HilbertLib.o HilbertLib.c 
CC=cc
OBJ=MDPoint.o MyTree.o PtrVector.o Pair.o HilbertLib.o Main.o ToGrid.o
DEPS=AxesTranspose.h HilbertLib.h MDPoint.h MyTree.h Pair.h PtrVector.h ToGrid.h

CFLAGS=-O3 -g
%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
main: $(OBJ)	
	$(CC) $(OBJ) -o Main.x $(CFLAGS)

MyTreeTestOBJ=MDPoint.o MyTree.o PtrVector.o tests/MyTreeTest.o
MyTree: $(MyTreeTestOBJ)
	$(CC) $(MyTreeTestOBJ) -o MyTreeTest.x $(CFLAGS)
MyTreeAutoTestOBJ=MDPoint.o MyTree.o PtrVector.o MyTreeAutoTest.o
MyTreeAuto: $(MyTreeAutoTestOBJ)
	$(CC) $(MyTreeAutoTestOBJ) -o MyTreeTest.x $(CFLAGS)

	
clean:
	rm -rf *.o
	rm -rf *.x
