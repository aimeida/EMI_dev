CC = g++
CFLAGS = -Wno-deprecated -O3
TARGET = emi
SRCS = emi.cpp FibonacciHeap.cpp FastGraphCluster.cpp DegreeArray.cpp DebugFunc.cpp Misc.cpp
OBJS = emi.o FibonacciHeap.o FastGraphCluster.o DegreeArray.o DebugFunc.o Misc.o

${TARGET}: ${OBJS}
	$(CC) $(CFLAGS) -o $@ $^ 

$(OBJS): $(SRCS)
	$(CC) $(CFLAGS) -c $*.cpp

.PHONY: clean

clean: 
	-rm -f *.o ${TARGET}