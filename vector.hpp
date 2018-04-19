/*
 *  GASMat
 *  (C) Greg Meyer, 2018
 */

#pragma once

#include <upcxx/upcxx.hpp>
#include <stdexcept>
#include "utils.hpp"

template <typename idx_t, typename data_t>
class Vec {
  idx_t _size = 0;
  idx_t _local_size;
  std::vector<upcxx::global_ptr<data_t>> gptrs;

public:
  Vec() {};
  Vec(idx_t size);
  ~Vec();

  idx_t get_size();
  idx_t get_local_size();
  void allocate_elements(idx_t size);
};

template <typename idx_t, typename data_t>
Vec<idx_t, data_t>::Vec(idx_t size) {
  allocate_elements(size);
}

template <typename idx_t, typename data_t>
Vec<idx_t, data_t>::~Vec() {
  if (get_size()) { /* check that we already called allocate_elements */
    upcxx::delete_array(gptrs[upcxx::rank_me()]);
  }
}

template <typename idx_t, typename data_t>
void Vec<idx_t, data_t>::allocate_elements(idx_t size) {

  /* can only set the size once */
  if (get_size()) {
    throw std::logic_error("Called allocate_elements after size has already been set.");
  }

  if (size <= 0) {
    throw std::length_error("size must be >= 0");
  }

  _size = size;
  _local_size = compute_local_size<idx_t>(get_size());

  /* allocate shared global memory and broadcast the pointers */
  gptrs.resize(upcxx::rank_n());

  /* compute the partitioning and allocate local portion */
  gptrs[upcxx::rank_me()] = upcxx::new_array<data_t>(_local_size);

  /* broadcast global pointers */
  for (int i = 0; i < upcxx::rank_n(); i++) {
    gptrs[i] = upcxx::broadcast(gptrs[i], i).wait();
  }

}

template <typename idx_t, typename data_t>
idx_t Vec<idx_t, data_t>::get_size() {
  return _size;
}

template <typename idx_t, typename data_t>
idx_t Vec<idx_t, data_t>::get_local_size() {
  return _local_size;
}
