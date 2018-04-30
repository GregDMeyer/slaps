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

TEST_CASE( "partition array. n=size, r=ranks. " TYPE_STR, "" ) {
  SECTION( "n=4, r=1" ) {
    auto p = _partition_array(4, 1);

    REQUIRE(p.size() == 2);

    REQUIRE(p[0] == 0);
    REQUIRE(p[1] == 4);
  }

  SECTION( "n=10, r=3" ) {
    auto p = _partition_array(10, 3);

    REQUIRE(p.size() == 4);

    REQUIRE(p[0] == 0);
    REQUIRE(p[1] == 4);
    REQUIRE(p[2] == 7);
    REQUIRE(p[3] == 10);
  }

  SECTION( "n=500, r=16" ) {
    auto p = _partition_array(500, 16);

    REQUIRE(p.size() == 17);

    REQUIRE(p[0] == 0);
    for (int i = 0; i < 4; ++i) {
      REQUIRE(p[i+1] - p[i] == 32);
    }
    for (int i = 4; i < 16; ++i) {
      REQUIRE(p[i+1] - p[i] == 31);
    }
  }

  SECTION( "n=2, r=2" ) {
    auto p = _partition_array(2, 2);

    REQUIRE(p.size() == 3);

    REQUIRE(p[0] == 0);
    REQUIRE(p[1] == 1);
    REQUIRE(p[2] == 2);
  }

  SECTION( "n=2, r=3" ) {
    auto p = _partition_array(2, 3);

    REQUIRE(p.size() == 4);

    REQUIRE(p[0] == 0);
    REQUIRE(p[1] == 1);
    REQUIRE(p[2] == 2);
    REQUIRE(p[3] == 2);
  }
}

TEST_CASE( "idx_to_proc. n=size, r=ranks. " TYPE_STR, "" ) {
  SECTION( "n=4, r=1" ) {
    for (int i = 0; i < 4; ++i) {
      REQUIRE(_idx_to_proc(i, 4, 1) == 0);
    }
  }

  SECTION( "n=500, r=16" ) {
    for (int i = 0; i < 128; ++i) {
      REQUIRE(_idx_to_proc(i, 500, 16) == i/32);
    }
    for (int i = 0; i < 500-128; ++i) {
      REQUIRE(_idx_to_proc(i+128, 500, 16) == 4 + i/31);
    }
  }

  SECTION( "n=2, r=2" ) {
    REQUIRE(_idx_to_proc(0,2,2) == 0);
    REQUIRE(_idx_to_proc(1,2,2) == 1);
  }

  SECTION( "n=2, r=3" ) {
    REQUIRE(_idx_to_proc(0,2,3) == 0);
    REQUIRE(_idx_to_proc(1,2,3) == 1);
  }
}
