/*
 *  GASMat
 *  (C) Greg Meyer, 2018
 */

#pragma once

#include <vector>
#include "vector.hpp"
#include "utils.hpp"

/*
 * This Mat class defines a Mat in CSR storage format.
 * In the future there may be a generic Mat class which is subclassed
 * for other storage formats.
 *
 * The Mat elements are not stored in globally addressable memory.
 * This makes sense, due to the fact that Mat elements don't currently need to be
 * transferred between processes. It also saves room in shared memory.
 * Implementing SpMSpM may require it.
 */

template <typename I, typename D>
class Mat
{

public:
  Mat() {};

  /* construct a Mat with dimensions M, N */
  Mat(I M, I N);

  /* set dimensions of Mat. clears all values. */
  void set_dimensions(I M, I N);

  /* get dimensions of the Mat */
  void get_dimensions(I& M, I& N) const;

  /* get the range of rows stored locally */
  void get_local_range(I& start, I& end) const;

  /* get the range of rows stored locally */
  I get_local_size() const;

  /* reserve space for elements. */
  void reserve(I nonzero_diag, I nonzero_offdiag);

  /* discard unused memory */
  void shrink_extra();

  /* set a Mat element */
  void set_value(I row, I col, D value);

  /* set a Mat element that we know is local */
  void set_diag_value(I row, I col, D value);

  /* set a Mat element that we know is off-diagonal */
  void set_offdiag_value(I row, I col, D value);

  /* TODO: reorganize our values to improve speed */
  void neaten();

  /* Mat-vector product y = A*x */
  void dot(Vec<I, D>& x, Vec<I, D>& y) const;

  /* Mat-vector sum product y = A*x + y */
  void plusdot(Vec<I, D>& x, Vec<I, D>& y) const;

private:
  I _M, _N;
  I _local_rows;

  /* vector storing local (block diagonal) Mat elements */
  std::vector< std::vector< std::pair<I, D>>> _local;

  /* vector storing remote (off-diagonal) Mat elements */
  std::vector< std::vector< std::pair<I, D>>> _remote;

  std::vector<I> _partitions;
};

template <typename I, typename D>
Mat<I, D>::Mat(I M, I N)
{
  set_dimensions(M, N);
}

template <typename I, typename D>
void Mat<I, D>::set_dimensions(I M, I N)
{

  if (M <= 0 || N <= 0) {
    throw std::length_error("Mat dimensions must be > 0");
  }

  _M = M;
  _N = N;

  /* compute local rows */
  _partitions = partition_array(M);
  _local_rows = _partitions[upcxx::rank_me()+1] - _partitions[upcxx::rank_me()];

  /* clear anything already in the Mat */
  _local.clear();
  _remote.clear();

  _local.resize(_local_rows);
  _remote.resize(_local_rows);

}

/* set dimensions of Mat, clearing all values, and reserve space for elements. */
template <typename I, typename D>
void Mat<I, D>::get_dimensions(I& M, I& N) const
{
  M = _M;
  N = _N;
}

template <typename I, typename D>
void Mat<I, D>::get_local_range(I& start, I& end) const
{
  start = _partitions[upcxx::rank_me()];
  end = _partitions[upcxx::rank_me() + 1];
}

template <typename I, typename D>
I Mat<I, D>::get_local_size() const
{
  I start, end;
  get_local_range(start, end);
  return end-start;
}

/* reserve space for elements. */
template <typename I, typename D>
void Mat<I, D>::reserve(I nonzero_diag, I nonzero_offdiag)
{

  for (I i = 0; i < get_local_size(); ++i) {
    _local[i].reserve(nonzero_diag);
    _remote[i].reserve(nonzero_offdiag);
  }
}

/* discard unused memory */
template <typename I, typename D>
void Mat<I, D>::shrink_extra()
{
  for (I i = 0; i < get_local_size(); ++i) {
    _local[i].shrink_to_fit();
    _remote[i].shrink_to_fit();
  }
}

template <typename I, typename D>
void Mat<I, D>::set_value(I row, I col, D value)
{
  I start, end;
  get_local_range(start, end);

  if (col >= start && col < end) {
    set_diag_value(row, col, value);
  }
  else {
    set_offdiag_value(row, col, value);
  }
}

template <typename I, typename D>
void Mat<I, D>::set_diag_value(I row, I col, D value)
{
  I start, end;
  get_local_range(start, end);

#ifdef DEBUG
  assert(row >= start && row < end);
  assert(col >= start && col < end);
#endif

  _local[row-start].push_back(std::make_pair(col-start, value));
}

template <typename I, typename D>
void Mat<I, D>::set_offdiag_value(I row, I col, D value)
{
  I start, end;
  get_local_range(start, end);

  #ifdef DEBUG
    assert(row >= start && row < end);
    assert(col < start || col >= end);
  #endif

  _remote[row-start].push_back(std::make_pair(col, value));
}

/* Mat-vector product y = A*x */
template <typename I, typename D>
void Mat<I, D>::dot(Vec<I, D>& x, Vec<I, D>& y) const
{
  y.set_all(0);
  plusdot(x, y);
}

/* Mat-vector sum product y = A*x + y */
template <typename I, typename D>
void Mat<I, D>::plusdot(Vec<I, D>& x, Vec<I, D>& y) const
{

  /* first do the local matvec */
  auto x_array = x.get_local_array_read();
  auto y_array = y.get_local_array();

  I local_size = _local.size();

  /* TODO: see if we can/should optimize this loop */
  for (I i = 0; i < local_size; ++i) {
    for (const auto& p : _local[i]) {
      y_array[i] += p.second * x_array[p.first];
    }
  }

  /* now remote part */

  /* TODO: this definitely will be optimized */
  for (I i = 0; i < local_size; ++i) {
    for (const auto& p : _remote[i]) {
      y_array[i] += p.second * x[p.first];
    }
  }

}
