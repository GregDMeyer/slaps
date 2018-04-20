/*
 *  This file is part of GASMat
 *  (C) Greg Meyer, 2018
 */

/*
 * This file gives a generic set of tests that is
 * oblivious to the data types. The data types are #define'd
 * and then this file is included in test-vector.cpp to generate
 * the actual test cases.
 */

/* macros to turn the data types into strings for the test case name */
#define _STR(x) #x
#define TO_STR(x) _STR(x)
#define TYPE_STR " idx_t=" TO_STR(IDX_T) " data_t=" TO_STR(DATA_T)

TEST_CASE( "constructor and destructor" TYPE_STR, "" ) {

  SECTION( "default constructor" ) {
    std::unique_ptr<Vec<IDX_T, DATA_T>> vp(new Vec<IDX_T, DATA_T>());
  }

  SECTION( "size constructor" ) {
    std::unique_ptr<Vec<IDX_T, DATA_T>> vp(new Vec<IDX_T, DATA_T>(10));
  }

}

TEST_CASE( "size getters and setters" TYPE_STR, "" ) {
  Vec<IDX_T, DATA_T>* v;

  SECTION( "default constructor" ) {
    v = new Vec<IDX_T, DATA_T>();
    REQUIRE(v->get_size() == 0);
    v->allocate_elements(10);
  }

  SECTION( "size constructor" ) {
    v = new Vec<IDX_T, DATA_T>(10);
  }

  REQUIRE(v->get_size() == 10);

  delete v;
}

TEST_CASE( "attempt to change size" TYPE_STR, "" ) {
  Vec<IDX_T, DATA_T>* v;

  SECTION( "default constructor" ) {
    v = new Vec<IDX_T, DATA_T>();
    v->allocate_elements(10);
  }

  SECTION( "size constructor" ) {
    v = new Vec<IDX_T, DATA_T>(10);
  }

  REQUIRE_THROWS_AS(v->allocate_elements(20), std::logic_error);

  delete v;
}

TEST_CASE( "set size to bad values" TYPE_STR, "" ) {
  Vec<IDX_T, DATA_T> v;

  SECTION( "set to 0" ) {
    REQUIRE_THROWS_AS(v.allocate_elements(0), std::length_error);
  }

  /* the following tests only make sense if we are using signed types for IDX_T
   * otherwise you just end up with a bad_malloc because it's identical to asking
   * for an obscenely large vector
   */

  if (std::numeric_limits<IDX_T>::is_signed) {
    SECTION( "set to -1" ) {
      REQUIRE_THROWS_AS(v.allocate_elements(-1), std::length_error);
    }

    SECTION( "set to -10" ) {
      REQUIRE_THROWS_AS(v.allocate_elements(-10), std::length_error);
    }
  }
}

TEST_CASE( "operator[] method" TYPE_STR, "" ) {

  IDX_T len(100);

  Vec<IDX_T, DATA_T> v(len);
  IDX_T start, end;
  v.get_local_range(start, end);

  DATA_T e = 2.718;

  SECTION("start") {
    v[start] = e;
    REQUIRE(v[start] == e);
  }

  SECTION("middle") {
    v[start + (end-start)/2] = e;
    REQUIRE(v[start + (end-start)/2] == e);
  }

  SECTION("end") {
    v[end-1] = e;
    REQUIRE(v[end-1] == e);
  }

  /* exceptions */
  SECTION("access globally low") {
    REQUIRE_THROWS_AS(v[-1], std::out_of_range);
  }

  SECTION("access locally low") {
    REQUIRE_THROWS_AS(v[start-1], std::out_of_range);
  }

  SECTION("access locally high") {
    REQUIRE_THROWS_AS(v[end], std::out_of_range);
  }

  SECTION("access globally high") {
    REQUIRE_THROWS_AS(v[len+1], std::out_of_range);
  }
}

TEST_CASE( "set_all method" TYPE_STR, "" ) {

  /* TODO: can we orthogonalize this test from the operator[] method? */

  Vec<IDX_T, DATA_T> v(100);
  IDX_T start, end;
  v.get_local_range(start, end);

  /* some subset of indices to check */
  std::vector<IDX_T> test_idxs = {0, 1, 5, 21, 40, 99};
  test_idxs.push_back(v.get_local_start());
  test_idxs.push_back(v.get_local_end()-1);

  SECTION( "value = 0" ) {
    v.set_all(0);
    for (auto idx : test_idxs) {
      if (idx >= start && idx < end) {
        REQUIRE(v[idx] == 0);
      }
    }
  }

  SECTION( "value = 3.14" ) {
    v.set_all(3.14);
    for (auto idx : test_idxs) {
      if (idx >= start && idx < end) {
        /* need to cast to DATA_T to avoid rounding errors */
        REQUIRE(v[idx] == DATA_T(3.14));
      }
    }
  }

}

TEST_CASE( "norm method" TYPE_STR, "" ) {

  Vec<IDX_T, DATA_T> v(100);
  IDX_T start, end;
  v.get_local_range(start, end);

  for (IDX_T i = start; i < end; ++i) {
    /* put some values in */
    v[i] = i / 3.2;
  }

  REQUIRE(v.norm() == Approx(179.0682263482274));
  
}
