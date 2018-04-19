/*
 *  GASMat
 *  (C) Greg Meyer, 2018
 */

#include <upcxx/upcxx.hpp>

/*
Compute the number of elements to be stored locally. Simply gives an equal number
of elements to each processor, and then from the remainder R gives one element to
each of the first R processes. This is the default partitioning used by PETSc,
which defines this function in PetscSplitOwnership (see
http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscSplitOwnership.html).
*/
template <typename idx_t>
idx_t compute_local_size(idx_t size) {
  auto nranks = upcxx::rank_n();
  return size/nranks + ((size % nranks) > upcxx::rank_me());
}
