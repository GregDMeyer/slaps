
tests: test-main.o test-vector.o test-matrix.o catch.hpp
	$(CXX) -o $@ $(LIBS) test-main.o test-vector.o test-matrix.o

test-main.o: test-main.cpp catch.hpp

test-vector.o: test-vector.cpp vector.hpp catch.hpp

test-matrix.o: test-matrix.cpp matrix.hpp vector.hpp catch.hpp
