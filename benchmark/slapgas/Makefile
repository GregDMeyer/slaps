
# for NERSC
CXX = CC
CXXFLAGS = -Wall

ifdef DEBUG
	OPTFLAGS = -O0 -DDEBUG -g
else
	OPTFLAGS = -O3
endif

INCLUDE = -I../../include

EXE_TARGETS = benchmark

# add in the flags for UPC++
CXXFLAGS += `upcxx-meta PPFLAGS` `upcxx-meta LDFLAGS` $(INCLUDE) $(OPTFLAGS)
LDFLAGS += `upcxx-meta LIBFLAGS`

all: $(EXE_TARGETS)

benchmark: benchmark.o
	$(CXX) -o $@ $(LIBS) benchmark.o $(CXXFLAGS) $(LDFLAGS)

CXXFLAGS += $(DEBUGFLAGS)

benchmark.o: benchmark.cpp ../../include/*.hpp

clean:
	$(RM) *.o $(EXE_TARGETS)
