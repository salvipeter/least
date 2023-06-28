all: least-test

EIGEN=/usr/include/eigen3
LIBGEOM=../libgeom
INCLUDES=-I$(LIBGEOM) -I$(EIGEN)
CXXFLAGS=-std=c++20 -pedantic -Wall -O3 $(INCLUDES)
LIBS=-L$(LIBGEOM)/release -lgeom

least-test: least-test.o least.o
	$(CXX) -o $@ $^ $(LIBS)

