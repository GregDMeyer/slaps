/*
 *  This file is part of SLAPS
 *  (C) Greg Meyer, 2018
 */

#pragma once

#include <upcxx/upcxx.hpp>
#include <stdexcept>
#include "utils.hpp"
#include "proxy.hpp"
#include <sstream>
#include <complex>

template <typename I, typename D>
class Vec {

public:
  /*==================================*/
  /*** constructors and destructors ***/
  Vec() : _put_fut(upcxx::make_future<>()) {};
  Vec(I size);
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
  void allocate_elements(I size);

  /* return true if Vec's memory has already been allocated */
  bool allocated() const;

  /*===================================*/
  /*** vector dimensions and indices ***/

  /* get the global dimension of the vector */
  I get_size() const;

  /* get the size of the locally stored portion of the vector */
  I get_local_size() const;

  /* get the start and end indices of the locally stored portion of the vector */
  I get_local_start() const;
  I get_local_end() const;
  void get_local_range(I &start, I &end) const;

  /* throw an exception if the vector sizes of this and v do not match */
  void validate_dims(const Vec& v) const;

  /*================================*/
  /*** getting and setting values ***/

  /* set the whole vector the same value */
  void set_all(D value);

  /* get/set data through array subscripting */
  RData<I, D> operator[](I index);

  /* wait for all remote values set with Vec[] to complete on all processes */
  void set_wait();

  /* return the future that is tracking setting of remote values */
  upcxx::future<> get_put_future() const;

  /* get a pointer to the local array holding local values */
  D* get_local_array();

  /* get a pointer to the local array holding local values (read-only) */
  const D* get_local_array_read() const;

  /* copy all values from this vector to another vector v */
  void copy(Vec& v) const;

  /*
   * get a range of local or remote values by index, and write them into buf.
   * buf must be already allocated with room for end-start values!
   */
  void read_range(I start, I end, D* buf);

  /* asynchronous */
  void read_range_begin(I start, I end, D* buf);
  void read_range_complete();

  /*======================*/
  /*** vector functions ***/

  /* compute the 2-norm of the vector */
  D norm() const;

  /* compute the dot product of this and a vector b (this is conjugated and transposed) */
  D dot(const Vec& b) const;

private:
  I _size = 0;
  I _local_size;
  I _allocated = false;

  std::vector<I> _partitions;
  std::vector<upcxx::global_ptr<D>> _gptrs;
  upcxx::global_ptr<D> _local_gptr;
  D* _local_data;

  upcxx::future<> _put_fut, _range_get_fut;
  bool getting = false;

};

/*########################*/
/***** implementation *****/

/*==================================*/
/*** constructors and destructors ***/

template <typename I, typename D>
Vec<I, D>::Vec(I size)
: _put_fut(upcxx::make_future())
{
  allocate_elements(size);
}

/* copy constructor */
template <typename I, typename D>
Vec<I, D>::Vec(const Vec& v)
: _put_fut(v.get_put_future())
{
  if (v.allocated()) {
    allocate_elements(v.get_size());
    v.copy(*this);
  }
}

/* move constructor */
template <typename I, typename D>
Vec<I, D>::Vec(Vec&& v) noexcept
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

template <typename I, typename D>
Vec<I, D>& Vec<I, D>::operator= (const Vec& v)
{
  Vec tmp(v); /* reuse copy constructor */
  *this = std::move(tmp);
  return *this;
}

template <typename I, typename D>
Vec<I, D>& Vec<I, D>::operator= (Vec&& v) noexcept
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

template <typename I, typename D>
Vec<I, D>::~Vec() {
  /* make sure we don't get rid of the array before we're done writing to it */
  set_wait();
  upcxx::barrier();

  if (allocated()) { /* check that we already called allocate_elements */
    upcxx::delete_array(_local_gptr);
  }
}

template <typename I, typename D>
void Vec<I, D>::allocate_elements(I size) {

  /* can only set the size once */
  if (allocated()) {
    throw std::logic_error("Called allocate_elements after size has already been set.");
  }

  if (size <= 0) {
    throw std::length_error("size must be > 0");
  }

  _size = size;
  _partitions = partition_array<I>(get_size());
  _local_size = _partitions[upcxx::rank_me()+1] - _partitions[upcxx::rank_me()];

  /* allocate shared global memory and broadcast the pointers */
  _gptrs.resize(upcxx::rank_n());

  /* compute the partitioning and allocate local portion */
  _local_gptr = upcxx::new_array<D>(get_local_size());
  _allocated = true;
  _gptrs[upcxx::rank_me()] = _local_gptr;
  _local_data = _gptrs[upcxx::rank_me()].local();

  /* broadcast global pointers */
  for (int i = 0; i < upcxx::rank_n(); i++) {
    _gptrs[i] = upcxx::broadcast(_gptrs[i], i).wait();
  }

}

template <typename I, typename D>
bool Vec<I, D>::allocated() const
{
  return _allocated;
}

/*===================================*/
/*** vector dimensions and indices ***/

template <typename I, typename D>
I Vec<I, D>::get_size() const {
  return _size;
}

template <typename I, typename D>
I Vec<I, D>::get_local_size() const {
  return _local_size;
}

template <typename I, typename D>
I Vec<I, D>::get_local_start() const {
  return _partitions[upcxx::rank_me()];
}

template <typename I, typename D>
I Vec<I, D>::get_local_end() const {
  return _partitions[upcxx::rank_me() + 1];
}

template <typename I, typename D>
void Vec<I, D>::get_local_range(I &start, I &end) const {
  start = get_local_start();
  end = get_local_end();
}

template <typename I, typename D>
void Vec<I, D>::validate_dims(const Vec& v) const {
  if (get_size() != v.get_size()) {
    std::ostringstream out;
    out << "vector sizes " << get_size() << " " << v.get_size();
    out << " do not match.";
    throw std::invalid_argument(out.str());
  }
}

/*================================*/
/*** getting and setting values ***/

template <typename I, typename D>
void Vec<I, D>::set_all(D value) {
  auto local_size = get_local_size();
  auto local_array = get_local_array();
  for (I i = 0; i < local_size; ++i) {
    local_array[i] = value;
  }
}

template <typename I, typename D>
RData<I, D> Vec<I, D>::operator[](I index) {

/* bounds check in DEBUG mode */
#ifdef DEBUG
  std::ostringstream out;
  if (index < 0 || index >= get_size()) {
    out << "requested index " << index << " out of range.";
    throw std::out_of_range(out.str());
  }
#endif

  auto source_proc = idx_to_proc(index, get_size());

  return RData<I, D>(_gptrs[source_proc] + (index-_partitions[source_proc]), _put_fut);
}

template <typename I, typename D>
void Vec<I, D>::set_wait() {
  _put_fut.wait();
  upcxx::barrier();

  /* might as well reset the future now, in case anything is accumulating in it */
  _put_fut = upcxx::make_future();
}

template <typename I, typename D>
upcxx::future<> Vec<I, D>::get_put_future() const {
  return _put_fut;
}

template <typename I, typename D>
D* Vec<I, D>::get_local_array() {
  /* TODO: should we make sure that we get this returned before we write to it again? */
  return _local_data;
}

template <typename I, typename D>
const D* Vec<I, D>::get_local_array_read() const {
  return _local_data;
}

template <typename I, typename D>
void Vec<I, D>::copy(Vec& v) const {

  validate_dims(v);

  auto mine = get_local_array_read();
  auto other = v.get_local_array();

  I local_size = get_local_size();

  /* hopefully the compiler will optimize this to use memcpy */
  for (I i = 0; i < local_size; ++i) {
    other[i] = mine[i];
  }

}

/*
 * get a range of local or remote values by index, and write them into buf.
 * buf must be already allocated with room for end-start values!
 */
template <typename I, typename D>
void Vec<I, D>::read_range(I start, I end, D* buf)
{
  read_range_begin(start, end, buf);
  read_range_complete();
}

/* asynchronous */
template <typename I, typename D>
void Vec<I, D>::read_range_begin(I start, I end, D* buf)
{
  I vend;
  vend = get_size();

  /* check bounds */
  if (start < 0) {
    std::ostringstream out;
    out << "start index " << start << " out of range.";
    throw std::out_of_range(out.str());
  }
  if (end > vend) {
    std::ostringstream out;
    out << "end index " << end << " greater than vector size " << vend;
    throw std::out_of_range(out.str());
  }
  if (start > end) {
    std::ostringstream out;
    out << "start index " << start << " greater than end index " << end;
    throw std::out_of_range(out.str());
  }

  /* find starting process to get from */
  I proc = idx_to_proc(start, vend);
  I tmp_start = start;

  /* set up a future to conjoin to */
  _range_get_fut = upcxx::make_future();
  auto gptr = _gptrs[proc] + (start - _partitions[proc]);

  while (tmp_start < end) {

    _range_get_fut = upcxx::when_all(_range_get_fut,
      upcxx::rget(gptr, buf + (tmp_start-start),
                  std::min(_partitions[proc+1], end) - tmp_start)
    );

    proc++;
    tmp_start = _partitions[proc];
    gptr = _gptrs[proc];
  }

  getting = true;
}

template <typename I, typename D>
void Vec<I, D>::read_range_complete()
{
  _range_get_fut.wait();
  getting = false;
}

/*======================*/
/*** vector functions ***/

template <typename I, typename D>
D Vec<I, D>::norm() const {

  /* first sum local values */
  D local_sum = 0;
  I local_size = get_local_size();
  auto local_array = get_local_array_read();

  /* TODO: numerical errors accumulate! should sum pairwise! */
  for (I i = 0; i < local_size; ++i) {
    /* make sure we account for complex types */
    local_sum += std::norm(local_array[i]);
  }

  /* now allreduce the local sums */
  D nrm2;
  nrm2 = upcxx::allreduce(local_sum, std::plus<D>()).wait();
  return std::sqrt(nrm2);
}

template <typename I, typename D>
D Vec<I, D>::dot(const Vec& b) const {

  validate_dims(b);

  /* first sum local values */
  D local_sum = 0;
  I local_size = get_local_size();
  auto local_array = get_local_array_read();
  auto other_array = b.get_local_array_read();

  /* TODO: need to be careful about complex numbers here */
  for (I i = 0; i < local_size; ++i) {
    /* TODO: make sure we account for complex types */
    local_sum += local_array[i] * other_array[i];
  }

  /* now allreduce the local sums */
  D sum;
  sum = upcxx::allreduce(local_sum, std::plus<D>()).wait();
  return sum;
}
