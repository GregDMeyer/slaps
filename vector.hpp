/*
 *  This file is part of GASMat
 *  (C) Greg Meyer, 2018
 */

#pragma once

#include <upcxx/upcxx.hpp>
#include <stdexcept>
#include "utils.hpp"
#include "proxy.hpp"
#include <sstream>
#include <complex>

template <typename idx_t, typename data_t>
class Vec {
  idx_t _size = 0;
  idx_t _local_size;

  std::vector<idx_t> _partitions;
  std::vector<upcxx::global_ptr<data_t>> gptrs;
  data_t* _local_data;

  upcxx::future<> put_fut;

public:
  /*==================================*/
  /*** constructors and destructors ***/
  Vec() : put_fut(upcxx::make_future<>()) {};
  Vec(idx_t size);
  ~Vec();

  /*
   * allocate memory for the Vec. this only needs to be called if the vector
   * was initialized using the default constructor
   */
  void allocate_elements(idx_t size);

  /*===================================*/
  /*** vector dimensions and indices ***/

  /* get the global dimension of the vector */
  idx_t get_size();

  /* get the size of the locally stored portion of the vector */
  idx_t get_local_size();

  /* get the start and end indices of the locally stored portion of the vector */
  idx_t get_local_start();
  idx_t get_local_end();
  void get_local_range(idx_t &start, idx_t &end);

  /*================================*/
  /*** getting and setting values ***/

  /* set the whole vector the same value */
  void set_all(data_t value);

  /* get/set data through array subscripting */
  RData<idx_t, data_t> operator[](idx_t index);

  /* wait for all remote values set with Vec[] to complete on all processes */
  void set_wait();

  /*======================*/
  /*** vector functions ***/

  /* compute the 2-norm of the vector */
  data_t norm();

};

/*########################*/
/***** implementation *****/

/*==================================*/
/*** constructors and destructors ***/

template <typename idx_t, typename data_t>
Vec<idx_t, data_t>::Vec(idx_t size)
: put_fut(upcxx::make_future())
{
  allocate_elements(size);
}

template <typename idx_t, typename data_t>
Vec<idx_t, data_t>::~Vec() {
  /* make sure we don't get rid of the array before we're done writing to it */
  set_wait();
  upcxx::barrier();

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

/*===================================*/
/*** vector dimensions and indices ***/

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

/*================================*/
/*** getting and setting values ***/

template <typename idx_t, typename data_t>
void Vec<idx_t, data_t>::set_all(data_t value) {
  auto local_size = get_local_size();
  for (idx_t i = 0; i < local_size; ++i) {
    _local_data[i] = value;
  }
}

template <typename idx_t, typename data_t>
RData<idx_t, data_t> Vec<idx_t, data_t>::operator[](idx_t index) {

/* bounds check in DEBUG mode */
#ifdef DEBUG
  std::ostringstream out;
  if (index < 0 || index >= get_size()) {
    out << "requested index " << index << " out of range.";
    throw std::out_of_range(out.str());
  }
#endif

  auto source_proc = idx_to_proc(index, get_size());

  return RData<idx_t, data_t>(gptrs[source_proc] + (index-_partitions[source_proc]), put_fut);
}

template <typename idx_t, typename data_t>
void Vec<idx_t, data_t>::set_wait() {
  put_fut.wait();
  upcxx::barrier();
}

/*======================*/
/*** vector functions ***/

template <typename idx_t, typename data_t>
data_t Vec<idx_t, data_t>::norm() {
  /* first sum local values */
  data_t local_sum = 0;
  idx_t local_size = get_local_size();
  for (idx_t i = 0; i < local_size; ++i) {
    /* make sure we account for complex types */
    local_sum += std::norm(_local_data[i]);
  }
  /* now allreduce the local sums */
  data_t nrm2;
  nrm2 = upcxx::allreduce(local_sum, std::plus<data_t>()).wait();
  return std::sqrt(nrm2);
}
