
# for NERSC
CXX = CC
CXXFLAGS = -Wall
DEBUGFLAGS = -g -O0

EXE_TARGETS = matrix-tests vector-tests utils-tests

# add in the flags for UPC++
CXXFLAGS += `upcxx-meta PPFLAGS` `upcxx-meta LDFLAGS`
LDFLAGS += `upcxx-meta LIBFLAGS`

all: $(EXE_TARGETS)

test-main.o: test-main.cpp catch.hpp

matrix-tests: test-main.o matrix-tests.o catch.hpp
	$(CXX) -o $@ $(LIBS) test-main.o matrix-tests.o $(CXXFLAGS) $(LDFLAGS)

vector-tests: test-main.o vector-tests.o catch.hpp
	$(CXX) -o $@ $(LIBS) test-main.o vector-tests.o $(CXXFLAGS) $(LDFLAGS)

utils-tests: test-main.o utils-tests.o catch.hpp
	$(CXX) -o $@ $(LIBS) test-main.o utils-tests.o $(CXXFLAGS) $(LDFLAGS)

CXXFLAGS += $(DEBUGFLAGS)

vector-tests.o: vector-tests.cpp vector-tests-template.cpp vector.hpp catch.hpp utils.hpp

utils-tests.o: utils-tests.cpp utils-tests-template.cpp catch.hpp utils.hpp

matrix-tests.o: matrix-tests.cpp matrix.hpp vector.hpp catch.hpp

clean:
	$(RM) *.o $(EXE_TARGETS)
