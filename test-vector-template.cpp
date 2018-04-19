
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
