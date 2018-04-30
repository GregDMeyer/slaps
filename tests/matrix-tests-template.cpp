/*
 *  This file is part of SLAPGAS
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

  SECTION( "scattered" ) {
    m.set_dimensions(10, 15);
    m.reserve(10, 5);
    m.get_local_rows(start, end);
    for (IDX_T i = 0; i < 45; ++i) {
      IDX_T idx = i%10;
      IDX_T idy = (i^5) % 15;
      if (idx >= start && idx < end) {
        m.set_value(idx, idy, std::sin(i));
      }
    }
    m.shrink_extra();
  }

  SECTION( "diag 11x23" ) {
    IDX_T row_start, row_end;
    m.set_dimensions(11, 23);
    m.get_diag_cols(start, end);
    m.get_local_rows(row_start, row_end);
    for (IDX_T i = start; i < end; ++i) {
      m.set_diag_value(row_start + ((i-row_start)%(row_end-row_start)), i, std::sin(i));
    }
    m.shrink_extra();
  }

  SECTION( "diag 23x11" ) {
    IDX_T row_start, row_end;
    m.set_dimensions(23, 11);
    m.get_diag_cols(start, end);
    m.get_local_rows(row_start, row_end);
    for (IDX_T i = row_start; i < row_end; ++i) {
      m.set_diag_value(i, start + ((i-start)%(end-start)), std::sin(i));
    }
    m.shrink_extra();
  }

}

TEST_CASE( "15x15 dot" TYPE_STR, "" ) {

  Mat<IDX_T, DATA_T> m;
  Vec<IDX_T, DATA_T> x, y;
  IDX_T start, end;

  std::vector<DATA_T> correct;

  m.set_dimensions(15, 15);
  m.get_local_rows(start, end);

  x.allocate_elements(15);
  y.allocate_elements(15);

  auto xarr = x.get_local_array();

  SECTION("diagonal identity") {
    for (IDX_T i = start; i < end; ++i) {
      m.set_value(i, i, 1);
    }

    x.set_all(0);
    correct.resize(15);
    for (IDX_T i = start; i < end; ++i) {
      xarr[i-start] = i;
      correct[i] = i;
    }
    upcxx::barrier();

    m.dot(x, y);
    upcxx::barrier();

    for (IDX_T i = start; i < end; ++i) {
      CHECK(y[i] == Approx(correct[i]));
    }
  }

  SECTION("diagonal index") {
    for (IDX_T i = start; i < end; ++i) {
      m.set_value(i, i, i);
    }

    x.set_all(0);
    correct.resize(15);
    for (IDX_T i = start; i < end; ++i) {
      xarr[i-start] = i;
      correct[i] = i*i;
    }
    upcxx::barrier();

    m.dot(x, y);
    upcxx::barrier();

    for (IDX_T i = start; i < end; ++i) {
      CHECK(y[i] == Approx(correct[i]));
    }
  }

  SECTION ("off-diagonal") {
    for (IDX_T i = start; i < end; ++i) {
      m.set_value(i, (end-start + i) % 15, 1);
    }

    x.set_all(0);
    correct.resize(15);
    for (IDX_T i = start; i < end; ++i) {
      xarr[i-start] = i;
      correct[i] = (i + (end-start))%15;
    }
    upcxx::barrier();

    m.dot(x, y);
    upcxx::barrier();

    for (IDX_T i = start; i < end; ++i) {
      CHECK(y[i] == Approx(correct[i]));
    }
  }

  SECTION ("off-diagonal index") {
    for (IDX_T i = start; i < end; ++i) {
      m.set_value(i, (end-start + i) % 15, (end-start + i) % 15);
    }

    x.set_all(0);
    correct.resize(15);
    for (IDX_T i = start; i < end; ++i) {
      xarr[i-start] = i;
      correct[i] = (i + (end-start))%15;
      correct[i] *= correct[i];
    }
    upcxx::barrier();

    m.dot(x, y);
    upcxx::barrier();

    for (IDX_T i = start; i < end; ++i) {
      CHECK(y[i] == Approx(correct[i]));
    }
  }

  SECTION ("full") {
    for (IDX_T i = start; i < end; ++i) {
      m.set_value(i, i, 1);
      m.set_value(i, (end-start + i) % 15, 1);
    }

    x.set_all(0);
    correct.resize(15);
    for (IDX_T i = start; i < end; ++i) {
      xarr[i-start] = i;
      correct[i] = (i + (end-start))%15;
      correct[i] += i;
    }
    upcxx::barrier();

    m.dot(x, y);
    upcxx::barrier();

    for (IDX_T i = start; i < end; ++i) {
      CHECK(y[i] == Approx(correct[i]));
    }
  }

}

TEST_CASE( "11x13 dot" TYPE_STR, "" ) {

  IDX_T M = 11, N = 13;

  Mat<IDX_T, DATA_T> m;
  Vec<IDX_T, DATA_T> x, y;
  IDX_T start, end;

  std::vector<DATA_T> correct;

  m.set_dimensions(M, N);
  m.get_local_rows(start, end);

  x.allocate_elements(N);
  y.allocate_elements(M);

  IDX_T xstart, xend;
  x.get_local_range(xstart, xend);
  auto xarr = x.get_local_array();

  SECTION("diagonal_identity") {
    for (IDX_T i = 0; i < end-start; ++i) {
      m.set_value(start + i, xstart + i, 1);
    }

    x.set_all(0);
    for (IDX_T i = 0; i < xend - xstart; ++i) {
      xarr[i] = xstart + i;
    }
    upcxx::barrier();

    correct.resize(M);
    for (IDX_T i = 0; i < end - start; ++i) {
      correct[start+i] = xstart + i;
    }

    m.dot(x, y);
    upcxx::barrier();

    for (IDX_T i = start; i < end; ++i) {
      CHECK(y[i] == Approx(correct[i]));
    }
  }

  SECTION("diagonal_index") {
    for (IDX_T i = 0; i < end-start; ++i) {
      m.set_value(start + i, xstart + i, start + i);
    }

    x.set_all(0);
    for (IDX_T i = 0; i < xend - xstart; ++i) {
      xarr[i] = 1;
    }
    upcxx::barrier();

    correct.resize(M);
    for (IDX_T i = start; i < end; ++i) {
      correct[i] = i;
    }

    m.dot(x, y);
    upcxx::barrier();

    for (IDX_T i = start; i < end; ++i) {
      CHECK(y[i] == Approx(correct[i]));
    }
  }

  /* TODO: fix this test */
  // SECTION ("off-diagonal") {
  //   for (IDX_T i = 0; i < end - start; ++i) {
  //     m.set_value(start + i, (xend + i) % N, 1);
  //   }
  //
  //   x.set_all(0);
  //   correct.resize(M);
  //   for (IDX_T i = 0; i < end - start; ++i) {
  //     xarr[i] = xstart + i;
  //     correct[start+i] = (xend+i) % N;
  //   }
  //   upcxx::barrier();
  //
  //   m.dot(x, y);
  //   upcxx::barrier();
  //
  //   for (IDX_T i = start; i < end; ++i) {
  //     CHECK(y[i] == Approx(correct[i]));
  //   }
  // }

  SECTION ("full") {
    for (IDX_T i = 0; i < 45; ++i) {
      IDX_T idx = i%M;
      IDX_T idy = (i^5) % N;
      if (idx >= start && idx < end) {
        m.set_value(idx, idy, 1);
      }
    }

    for (IDX_T i = xstart; i < xend; ++i) {
      xarr[i - xstart] = i+1;
    }

    correct.resize(M);
    for (IDX_T i = start; i < end; ++i) {
      IDX_T ii = i;
      correct[i] = 0;
      while (ii < 45) {
        correct[i] += ((ii^5) % N)+1;
        ii += M;
      }
    }
    upcxx::barrier();

    m.dot(x, y);
    upcxx::barrier();

    for (IDX_T i = start; i < end; ++i) {
      CHECK(y[i] == Approx(correct[i]));
    }
  }

}

TEST_CASE( "13x11 dot" TYPE_STR, "" ) {

  IDX_T M = 13, N = 11;

  Mat<IDX_T, DATA_T> m;
  Vec<IDX_T, DATA_T> x, y;
  IDX_T start, end;

  std::vector<DATA_T> correct;

  m.set_dimensions(M, N);
  m.get_local_rows(start, end);

  x.allocate_elements(N);
  y.allocate_elements(M);

  IDX_T xstart, xend;
  x.get_local_range(xstart, xend);
  auto xarr = x.get_local_array();

  SECTION("diagonal_identity") {
    for (IDX_T i = 0; i < xend-xstart; ++i) {
      m.set_value(start + i, xstart + i, 1);
    }

    x.set_all(0);
    for (IDX_T i = 0; i < xend - xstart; ++i) {
      xarr[i] = xstart + i;
    }
    upcxx::barrier();

    correct.resize(M);
    for (IDX_T i = 0; i < xend - xstart; ++i) {
      correct[start+i] = xstart + i;
    }

    m.dot(x, y);
    upcxx::barrier();

    for (IDX_T i = start; i < end; ++i) {
      CHECK(y[i] == Approx(correct[i]));
    }
  }

  SECTION("diagonal_index") {
    for (IDX_T i = 0; i < xend-xstart; ++i) {
      m.set_value(start + i, xstart + i, start + i);
    }

    x.set_all(0);
    for (IDX_T i = 0; i < xend - xstart; ++i) {
      xarr[i] = 1;
    }
    upcxx::barrier();

    correct.resize(M);
    for (IDX_T i = 0; i < xend-xstart; ++i) {
      correct[start+i] = start + i;
    }

    m.dot(x, y);
    upcxx::barrier();

    for (IDX_T i = start; i < end; ++i) {
      CHECK(y[i] == Approx(correct[i]));
    }
  }

  SECTION ("off-diagonal") {
    for (IDX_T i = 0; i < xend - xstart; ++i) {
      m.set_value(start + i, (xend + i) % N, 1);
    }

    x.set_all(0);
    correct.resize(M);
    for (IDX_T i = xstart; i < xend; ++i) {
      xarr[i-xstart] = i;
    }
    for (IDX_T i = 0; i < xend - xstart; ++i) {
      correct[start+i] = (xend+i) % N;
    }

    upcxx::barrier();

    m.dot(x, y);
    upcxx::barrier();

    for (IDX_T i = start; i < end; ++i) {
      CHECK(y[i] == Approx(correct[i]));
    }
  }

}

TEST_CASE( "dot exceptions" TYPE_STR, "" ) {

  Mat<IDX_T, DATA_T> m;
  Vec<IDX_T, DATA_T> x, y;

  SECTION( "bad x dim" ) {
    m.set_dimensions(10,20);
    x.allocate_elements(12);
    y.allocate_elements(10);
    REQUIRE_THROWS_AS( m.dot(x, y), std::invalid_argument );
  }

  SECTION( "bad y dim" ) {
    m.set_dimensions(10,20);
    x.allocate_elements(20);
    y.allocate_elements(12);
    REQUIRE_THROWS_AS( m.dot(x, y), std::invalid_argument );
  }

  SECTION( "bad x and y dim" ) {
    m.set_dimensions(10,20);
    x.allocate_elements(12);
    y.allocate_elements(13);
    REQUIRE_THROWS_AS( m.dot(x, y), std::invalid_argument );
  }

}
