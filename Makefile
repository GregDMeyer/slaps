
# for NERSC
CXX = CC
CXXFLAGS =
DEBUGFLAGS = -g -O0

EXE_TARGETS = matrix-tests vector-tests

# add in the flags for UPC++
CXXFLAGS += `upcxx-meta PPFLAGS` `upcxx-meta LDFLAGS`
LDFLAGS += `upcxx-meta LIBFLAGS`

all: $(EXE_TARGETS)

matrix-tests: test-main.o test-matrix.o catch.hpp
	$(CXX) -o $@ $(LIBS) test-main.o test-matrix.o $(CXXFLAGS) $(LDFLAGS)

vector-tests: test-main.o test-vector.o catch.hpp
	$(CXX) -o $@ $(LIBS) test-main.o test-vector.o $(CXXFLAGS) $(LDFLAGS)

test-main.o: test-main.cpp catch.hpp

CXXFLAGS += $(DEBUGFLAGS)

test-vector.o: test-vector.cpp test-vector-template.cpp vector.hpp catch.hpp

test-matrix.o: test-matrix.cpp matrix.hpp vector.hpp catch.hpp

clean:
	$(RM) *.o $(EXE_TARGETS)
