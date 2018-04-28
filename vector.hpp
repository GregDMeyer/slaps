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

public:
  /*==================================*/
  /*** constructors and destructors ***/
  Vec() : _put_fut(upcxx::make_future<>()) {};
  Vec(idx_t size);
  ~Vec();

  /* rule of three/five: since we have an explicit destructor we also
   * need copy/move methods */
  Vec(const Vec& v);
  Vec(Vec&& v) noexcept;
  Vec& operator= (const Vec& v);
  Vec& operator= (Vec&& v) noexcept;

  /*
   * allocate memory for the Vec. this only needs to be called if the vector
   * was initialized using the default constructor
   */
  void allocate_elements(idx_t size);

  /* return true if Vec's memory has already been allocated */
  bool allocated() const;

  /*===================================*/
  /*** vector dimensions and indices ***/

  /* get the global dimension of the vector */
  idx_t get_size() const;

  /* get the size of the locally stored portion of the vector */
  idx_t get_local_size() const;

  /* get the start and end indices of the locally stored portion of the vector */
  idx_t get_local_start() const;
  idx_t get_local_end() const;
  void get_local_range(idx_t &start, idx_t &end) const;

  /* throw an exception if the vector sizes of this and v do not match */
  void validate_dims(const Vec& v) const;

  /*================================*/
  /*** getting and setting values ***/

  /* set the whole vector the same value */
  void set_all(data_t value);

  /* get/set data through array subscripting */
  RData<idx_t, data_t> operator[](idx_t index);

  /* wait for all remote values set with Vec[] to complete on all processes */
  void set_wait();

  /* return the future that is tracking setting of remote values */
  upcxx::future<> get_put_future() const;

  /* get a pointer to the local array holding local values */
  data_t* get_local_array();

  /* get a pointer to the local array holding local values (read-only) */
  const data_t* get_local_array_read() const;

  /* copy all values from this vector to another vector v */
  void copy(Vec& v) const;

  /*======================*/
  /*** vector functions ***/

  /* compute the 2-norm of the vector */
  data_t norm() const;

  /* compute the dot product of this and a vector b (this is conjugated and transposed) */
  data_t dot(const Vec& b) const;

private:
  idx_t _size = 0;
  idx_t _local_size;
  idx_t _allocated = false;

  std::vector<idx_t> _partitions;
  std::vector<upcxx::global_ptr<data_t>> _gptrs;
  upcxx::global_ptr<data_t> _local_gptr;
  data_t* _local_data;

  upcxx::future<> _put_fut;

};

/*########################*/
/***** implementation *****/

/*==================================*/
/*** constructors and destructors ***/

template <typename idx_t, typename data_t>
Vec<idx_t, data_t>::Vec(idx_t size)
: _put_fut(upcxx::make_future())
{
  allocate_elements(size);
}

/* copy constructor */
template <typename idx_t, typename data_t>
Vec<idx_t, data_t>::Vec(const Vec& v)
: _put_fut(v.get_put_future())
{
  if (v.allocated()) {
    allocate_elements(v.get_size());
    v.copy(*this);
  }
}

/* move constructor */
template <typename idx_t, typename data_t>
Vec<idx_t, data_t>::Vec(Vec&& v) noexcept
: _size(v._size)
, _local_size(v._local_size)
, _allocated(v._allocated)
, _partitions( std::move(v._partitions) )
, _gptrs( std::move(v._gptrs) )
, _local_gptr(v._local_gptr)
, _local_data(v._local_data)
, _put_fut(v._put_fut)
{
  /* make sure we don't free the memory we just gave to our new Vec */
  v._local_gptr = nullptr;
  v._allocated = false;
}

template <typename idx_t, typename data_t>
Vec<idx_t, data_t>& Vec<idx_t, data_t>::operator= (const Vec& v)
{
  Vec tmp(v); /* reuse copy constructor */
  *this = std::move(tmp);
  return *this;
}

template <typename idx_t, typename data_t>
Vec<idx_t, data_t>& Vec<idx_t, data_t>::operator= (Vec&& v) noexcept
{
  if (this == &v) return *this;

  /* if we already are allocated, need to free memory */
  if (allocated()) {
    upcxx::delete_array(_local_gptr);
  }

  if (v.allocated()) {
    _size = v._size;
    _local_size = v._local_size;
    _partitions = std::move(v._partitions);
    _gptrs = std::move(v._gptrs);

    _local_gptr = v._local_gptr;
    v._local_gptr = nullptr;
    v._allocated = false;

    _local_data = v._local_data;
    _put_fut = v._put_fut;

    _allocated = true;
  }
  else {
    _allocated = false;
  }

  return *this;
}

template <typename idx_t, typename data_t>
Vec<idx_t, data_t>::~Vec() {
  /* make sure we don't get rid of the array before we're done writing to it */
  set_wait();
  upcxx::barrier();

  if (allocated()) { /* check that we already called allocate_elements */
    upcxx::delete_array(_local_gptr);
  }
}

template <typename idx_t, typename data_t>
void Vec<idx_t, data_t>::allocate_elements(idx_t size) {

  /* can only set the size once */
  if (allocated()) {
    throw std::logic_error("Called allocate_elements after size has already been set.");
  }

  if (size <= 0) {
    throw std::length_error("size must be >= 0");
  }

  _size = size;
  _partitions = partition_array<idx_t>(get_size());
  _local_size = _partitions[upcxx::rank_me()+1] - _partitions[upcxx::rank_me()];

  /* allocate shared global memory and broadcast the pointers */
  _gptrs.resize(upcxx::rank_n());

  /* compute the partitioning and allocate local portion */
  _local_gptr = upcxx::new_array<data_t>(get_local_size());
  _allocated = true;
  _gptrs[upcxx::rank_me()] = _local_gptr;
  _local_data = _gptrs[upcxx::rank_me()].local();

  /* broadcast global pointers */
  for (int i = 0; i < upcxx::rank_n(); i++) {
    _gptrs[i] = upcxx::broadcast(_gptrs[i], i).wait();
  }

}

template <typename idx_t, typename data_t>
bool Vec<idx_t, data_t>::allocated() const
{
  return _allocated;
}

/*===================================*/
/*** vector dimensions and indices ***/

template <typename idx_t, typename data_t>
idx_t Vec<idx_t, data_t>::get_size() const {
  return _size;
}

template <typename idx_t, typename data_t>
idx_t Vec<idx_t, data_t>::get_local_size() const {
  return _local_size;
}

template <typename idx_t, typename data_t>
idx_t Vec<idx_t, data_t>::get_local_start() const {
  return _partitions[upcxx::rank_me()];
}

template <typename idx_t, typename data_t>
idx_t Vec<idx_t, data_t>::get_local_end() const {
  return _partitions[upcxx::rank_me() + 1];
}

template <typename idx_t, typename data_t>
void Vec<idx_t, data_t>::get_local_range(idx_t &start, idx_t &end) const {
  start = get_local_start();
  end = get_local_end();
}

template <typename idx_t, typename data_t>
void Vec<idx_t, data_t>::validate_dims(const Vec& v) const {
  if (get_size() != v.get_size()) {
    std::ostringstream out;
    out << "vector sizes " << get_size() << " " << v.get_size();
    out << " do not match.";
    throw std::invalid_argument(out.str());
  }
}

/*================================*/
/*** getting and setting values ***/

template <typename idx_t, typename data_t>
void Vec<idx_t, data_t>::set_all(data_t value) {
  auto local_size = get_local_size();
  auto local_array = get_local_array();
  for (idx_t i = 0; i < local_size; ++i) {
    local_array[i] = value;
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

  return RData<idx_t, data_t>(_gptrs[source_proc] + (index-_partitions[source_proc]), _put_fut);
}

template <typename idx_t, typename data_t>
void Vec<idx_t, data_t>::set_wait() {
  _put_fut.wait();
  upcxx::barrier();

  /* might as well reset the future now, in case anything is accumulating in it */
  _put_fut = upcxx::make_future();
}

template <typename idx_t, typename data_t>
upcxx::future<> Vec<idx_t, data_t>::get_put_future() const {
  return _put_fut;
}

template <typename idx_t, typename data_t>
data_t* Vec<idx_t, data_t>::get_local_array() {
  /* TODO: should we make sure that we get this returned before we write to it again? */
  return _local_data;
}

template <typename idx_t, typename data_t>
const data_t* Vec<idx_t, data_t>::get_local_array_read() const {
  return _local_data;
}

template <typename idx_t, typename data_t>
void Vec<idx_t, data_t>::copy(Vec& v) const {

  validate_dims(v);

  auto mine = get_local_array_read();
  auto other = v.get_local_array();

  idx_t local_size = get_local_size();

  /* hopefully the compiler will optimize this to use memcpy */
  for (idx_t i = 0; i < local_size; ++i) {
    other[i] = mine[i];
  }

}

/*======================*/
/*** vector functions ***/

template <typename idx_t, typename data_t>
data_t Vec<idx_t, data_t>::norm() const {

  /* first sum local values */
  data_t local_sum = 0;
  idx_t local_size = get_local_size();
  auto local_array = get_local_array_read();

  /* TODO: numerical errors accumulate! should sum pairwise! */
  for (idx_t i = 0; i < local_size; ++i) {
    /* make sure we account for complex types */
    local_sum += std::norm(local_array[i]);
  }

  /* now allreduce the local sums */
  data_t nrm2;
  nrm2 = upcxx::allreduce(local_sum, std::plus<data_t>()).wait();
  return std::sqrt(nrm2);
}

template <typename idx_t, typename data_t>
data_t Vec<idx_t, data_t>::dot(const Vec& b) const {

  validate_dims(b);

  /* first sum local values */
  data_t local_sum = 0;
  idx_t local_size = get_local_size();
  auto local_array = get_local_array_read();
  auto other_array = b.get_local_array_read();

  /* TODO: need to be careful about complex numbers here */
  for (idx_t i = 0; i < local_size; ++i) {
    /* TODO: make sure we account for complex types */
    local_sum += local_array[i] * other_array[i];
  }

  /* now allreduce the local sums */
  data_t sum;
  sum = upcxx::allreduce(local_sum, std::plus<data_t>()).wait();
  return sum;
}
