/*
 *  This file is part of SLAPS
 *  (C) Greg Meyer, 2018
 */

#include <upcxx/upcxx.hpp>
#include <string>

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

int main( int argc, char* argv[] ) {

  upcxx::init();

  /* if we have multiple ranks, set up stdout to go to various files */
  std::ofstream out;
  if (upcxx::rank_n() > 1) {
    std::string fname = "test_out_" + std::to_string(upcxx::rank_me()) + ".txt";
    out.open(fname);
    std::cout.rdbuf(out.rdbuf());
  }

  /* actually run catch */
  int result = Catch::Session().run( argc, argv );

  if (upcxx::rank_n() > 1) {
    out.close();
  }

  upcxx::finalize();

  return result;

}
