
#include <upcxx/upcxx.hpp>
#include <string>

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

int main( int argc, char* argv[] ) {

  upcxx::init();

  /* set up stdout from ranks to go to various files */
  std::string fname = "test_out_" + std::to_string(upcxx::rank_me()) + ".txt";

  std::ofstream out(fname);
  std::cout.rdbuf(out.rdbuf());

  /* actually run catch */
  int result = Catch::Session().run( argc, argv );

  upcxx::finalize();

  return result;
  
}
