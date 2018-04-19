/*
 *  This file is part of GASMat
 *  (C) Greg Meyer, 2018
 */

#pragma once

#include <upcxx/upcxx.hpp>
#include <stdexcept>
#include "utils.hpp"
#include <sstream>

template <typename idx_t, typename data_t>
class Vec {
  idx_t _size = 0;
  idx_t _local_size;
  std::vector<idx_t> _partitions;
  std::vector<upcxx::global_ptr<data_t>> gptrs;
  data_t* _local_data;

public:
  Vec() {};
  Vec(idx_t size);
  ~Vec();

  /* get the global dimension of the vector */
  idx_t get_size();

  /* get the size of the locally stored portion of the vector */
  idx_t get_local_size();

  /* get the start and end indices of the locally stored portion of the vector */
  idx_t get_local_start();
  idx_t get_local_end();
  void get_local_range(idx_t &start, idx_t &end);

  /*
   * allocate memory for the Vec. this only needs to be called if the vector
   * was initialized using the default constructor
   */
  void allocate_elements(idx_t size);

  /* set the whole vector the same value */
  void set_all(data_t value);

  /* get/set LOCAL data through array subscripting */
  /*
   * TODO: this is a design decision: do we want to abstract away remote data
   * setting, or do we want the programmer to be aware of the extra cost of setting
   * remote data? I think it's cool to abstract it away but it could bite anyone
   * who isn't careful in terms of performance.
   */
  data_t& operator[](idx_t index);
};

/***** implementation *****/

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
  _partitions = partition_array<idx_t>(get_size());
  _local_size = _partitions[upcxx::rank_me()+1] - _partitions[upcxx::rank_me()];

  /* allocate shared global memory and broadcast the pointers */
  gptrs.resize(upcxx::rank_n());

  /* compute the partitioning and allocate local portion */
  gptrs[upcxx::rank_me()] = upcxx::new_array<data_t>(get_local_size());
  _local_data = gptrs[upcxx::rank_me()].local();

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

template <typename idx_t, typename data_t>
idx_t Vec<idx_t, data_t>::get_local_start() {
  return _partitions[upcxx::rank_me()];
}

template <typename idx_t, typename data_t>
idx_t Vec<idx_t, data_t>::get_local_end() {
  return _partitions[upcxx::rank_me() + 1];
}

template <typename idx_t, typename data_t>
void Vec<idx_t, data_t>::get_local_range(idx_t &start, idx_t &end) {
  start = get_local_start();
  end = get_local_end();
}

template <typename idx_t, typename data_t>
void Vec<idx_t, data_t>::set_all(data_t value) {
  auto local_size = get_local_size();
  for (idx_t i = 0; i < local_size; ++i) {
    _local_data[i] = value;
  }
}

template <typename idx_t, typename data_t>
data_t& Vec<idx_t, data_t>::operator[](idx_t index) {

/* bounds check in DEBUG mode */
#ifdef DEBUG
  idx_t start, end;
  get_local_range(start, end);

  std::ostringstream out;
  if (index < start || index >= end) {
    out << "requested index " << index;
    out << " out of local range (" << start << ", " << end <<")";
    out << " for processor " << upcxx::rank_me();
    throw std::out_of_range(out.str());
  }
#endif

  idx_t local_idx = index - get_local_start();
  return _local_data[local_idx];
}
