
#include <slapgas.hpp>
#include <cstring>
#include <iostream>
#include <chrono>
#include <upcxx/upcxx.hpp>

typedef unsigned int I;
typedef double D;
typedef SingleCSRMat<I,D> mat_t;

#define LARGE_PRIME 1046527

void parse_args(int argc, char* argv[],
                I& dim,
                I& sparsity,
                I& iterations,
                bool& quiet);

int main(int argc, char* argv[])
{
  mat_t m;
  Vec<I,D> x, y;
  double timedelta;
  bool do_print;
  I row_start, row_end;

  I dim = 100, sparsity = 10, iterations = 100;
  bool quiet = false;

  upcxx::init();
  if (upcxx::rank_me() == 0) do_print = true;
  else do_print = false;

  parse_args(argc, argv, dim, sparsity, iterations, quiet);

  if (!quiet && do_print) {
    std::cout << "Timing SLAPGAS MatVec." << std::endl;
    std::cout << " dim = " << dim << std::endl;
    std::cout << " sparsity = " << sparsity << std::endl;
    std::cout << " iterations = " << iterations << std::endl;
  }

  m.set_dimensions(dim, dim);
  x.allocate_elements(dim);
  y.allocate_elements(dim);

  /* set 1's in the correct sparsity pattern in the matrix */
  I tot_nz = dim/sparsity + 1;
  m.get_local_rows(row_start, row_end);

  for (I i = row_start; i < row_end; ++i) {
    /* start from a different spot each time */
    for (I j = (LARGE_PRIME*i) % sparsity; j < dim; j += sparsity) {
      m.set_value(i, j, 1);
    }
  }

  m.setup(tot_nz/upcxx::rank_n() + 1, tot_nz - tot_nz/upcxx::rank_n() + 1);

  x.set_all(1);

  upcxx::barrier();

  auto tick = std::chrono::system_clock::now();
  for (I i = 0; i < iterations; ++i) {
    m.dot(x, y);
  }
  /* make sure we're all done */
  upcxx::barrier();
  auto tock = std::chrono::system_clock::now();

  if (do_print) {
    if (!quiet) std::cout << "Time: ";
    auto td = std::chrono::duration_cast<std::chrono::microseconds>(tock-tick).count();
    std::cout << td / 1000000. << std::endl;
  }

  upcxx::finalize();

  return 0;

}

void parse_args(int argc, char* argv[],
                I& dim,
                I& sparsity,
                I& iterations,
                bool& quiet) {

  bool recognized = true;

  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i], "-q")) {
      quiet = true;
    }
    else if (i+1 < argc) {
      if (!strcmp(argv[i], "-d")) {
        dim = atoi(argv[i+1]);
        i++;
      }
      else if (!strcmp(argv[i], "-sp")) {
        sparsity = atoi(argv[i+1]);
        i++;
      }
      else if (!strcmp(argv[i], "-it")) {
        iterations = atoi(argv[i+1]);
        i++;
      }
      else {
        recognized = false;
      }
    }
    else {
      recognized = false;
    }

    if (!recognized && upcxx::rank_me() == 0) {
      std::cout << "Unrecognized argument \"" << argv[i] << "\"" << std::endl;
    }
  }
}
