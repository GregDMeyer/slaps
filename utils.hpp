/*
 *  This file is part of GASMat
 *  (C) Greg Meyer, 2018
 */

#include <upcxx/upcxx.hpp>
#include <assert.h>

/*
Compute the number of elements to be stored locally. Simply gives an equal number
of elements to each processor, and then from the remainder R gives one element to
each of the first R processes. This is the default partitioning used by PETSc,
which defines this function in PetscSplitOwnership (see
http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscSplitOwnership.html).
*/

/* isolate this function from upc++ */
template <typename idx_t>
std::vector<idx_t> _partition_array(const idx_t size, unsigned int nranks) {
  std::vector<idx_t> p;
  idx_t idx;

  idx = 0;
  for (unsigned int rank = 0; rank < nranks; ++rank) {
    p.push_back(idx);
    idx += size/nranks + ((size % nranks) > rank);
  }
  p.push_back(size);

  return p;
}

template <typename idx_t>
std::vector<idx_t> partition_array(const idx_t size) {
  int nranks = upcxx::rank_n();

  /* otherwise, upcxx::init() was not called */
  assert(nranks > 0);

  return _partition_array(size, nranks);
}
