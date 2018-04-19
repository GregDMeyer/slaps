/*
 *  This file is part of GASMat
 *  (C) Greg Meyer, 2018
 */
 
#include <upcxx/upcxx.hpp>
#include <string>

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

int main( int argc, char* argv[] ) {

  upcxx::init();

  /* if we have multiple ranks, set up stdout to go to various files */
  if (upcxx::rank_n() > 1) {
    std::string fname = "test_out_" + std::to_string(upcxx::rank_me()) + ".txt";
    std::ofstream out(fname);
    std::cout.rdbuf(out.rdbuf());
  }

  /* actually run catch */
  int result = Catch::Session().run( argc, argv );

  upcxx::finalize();

  return result;

}
