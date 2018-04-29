/*
 *  This file is part of GASMat
 *  (C) Greg Meyer, 2018
 */

/*
 * This file gives a generic set of tests that is
 * oblivious to the data types. The data types are #define'd
 * and then this file is included in test-utils.cpp to generate
 * the actual test cases.
 */

/* macros to turn the data types into strings for the test case name */
#define _STR(x) #x
#define TO_STR(x) _STR(x)
#define TYPE_STR " \tidx_t=" TO_STR(IDX_T) " \tdata_t=" TO_STR(DATA_T)

TEST_CASE( "constructor/destructor " TYPE_STR, "" ) {

  SECTION( "default" ) {
    Mat<IDX_T, DATA_T> m;
  }

  SECTION( "default + set_dimensions" ) {
    Mat<IDX_T, DATA_T> m;
    m.set_dimensions(10, 12);

    IDX_T M, N;
    m.get_dimensions(M, N);
    REQUIRE(M == 10);
    REQUIRE(N == 12);
  }

  SECTION( "size" ) {
    Mat<IDX_T, DATA_T> m(10, 12);

    IDX_T M, N;
    m.get_dimensions(M, N);
    REQUIRE(M == 10);
    REQUIRE(N == 12);
  }

}

TEST_CASE( "set values" TYPE_STR, "" ) {

  Mat<IDX_T, DATA_T> m;
  IDX_T start, end;

  m.set_dimensions(10, 15);
  m.reserve(10, 5);
  m.get_local_range(start, end);
  for (IDX_T i = 0; i < 45; ++i) {
    IDX_T idx = i%10;
    IDX_T idy = (i^5) % 15;
    if (idx >= start && idx < end) {
      m.set_value(idx, idy, std::sin(i));
    }
  }

  m.shrink_extra();

}

TEST_CASE( "dot" TYPE_STR, "" ) {

  Mat<IDX_T, DATA_T> m;
  Vec<IDX_T, DATA_T> x, y;
  IDX_T start, end;

  m.set_dimensions(10, 15);
  m.get_local_range(start, end);
  for (IDX_T i = 0; i < 45; ++i) {
    IDX_T idx = i%10;
    IDX_T idy = (i^5) % 15;
    if (idx >= start && idx < end) {
      m.set_value(idx, idy, std::sin(i));
    }
  }

  x.allocate_elements(15);
  y.allocate_elements(10);

  IDX_T xstart, xend;
  x.get_local_range(xstart, xend);
  auto xarr = x.get_local_array();

  for (IDX_T i = xstart; i < xend; ++i) {
    xarr[i - xstart] = i + 1;
  }
  upcxx::barrier();

  m.dot(x, y);

  upcxx::barrier();

  std::vector<DATA_T> correct = { -9.90448331, -16.34697867, 3.52624249, 6.72029531,
          2.29325781, 0.48767344, -5.57498372, -4.77166363, 11.85030538, 4.73919607};

  for (IDX_T i = 0; i < 10; ++i) {
    REQUIRE(y[i] == Approx(correct[i]));
  }

}
