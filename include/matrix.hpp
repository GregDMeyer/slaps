/*
 *  SLAPGAS
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

/* this enum defines different implementations of MatVec that we can use */
enum class DotImpl { Naive, Single, Block };

/* how many values to prefetch */
#define DOT_BLOCK_SIZE 2048

template <typename I, typename D>
class Mat
{

public:
  /*==================================*/
  /*** constructors and destructors ***/

  Mat() {};

  /* construct a Mat with dimensions M, N */
  Mat(I M, I N);

  /*============================*/
  /*** dimensions and indices ***/

  /* set dimensions of Mat. clears all values. */
  void set_dimensions(I M, I N);

  /* get dimensions of the Mat */
  void get_dimensions(I& M, I& N) const;

  /* get the range of rows stored locally */
  void get_local_rows(I& start, I& end) const;

  /* get the range of cols corresponding to local part of the vector */
  void get_diag_cols(I& start, I& end) const;

  /* get the range of rows stored locally */
  I get_local_rows_size() const;

  /*=========================================*/
  /*** value setting and memory allocation ***/

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
  // void neaten();

  /*============================*/
  /*** matrix-vector products ***/

  /* Mat-vector product y = A*x */
  void dot(Vec<I,D>& x, Vec<I,D>& y) const;

  /* Mat-vector product y = A*x using implementation `impl` */
  void dot(Vec<I,D>& x, Vec<I,D>& y, DotImpl impl) const;

  /* Mat-vector sum product y = A*x + y */
  void plusdot(Vec<I,D>& x, Vec<I,D>& y) const;

  /* Mat-vector sum product y = A*x + y using implementation `impl` */
  void plusdot(Vec<I,D>& x, Vec<I,D>& y, DotImpl impl) const;

  /* Naive implementation of plusdot with no prefetching */
  void _plusdot_naive(Vec<I,D>& x, Vec<I,D>& y) const;

  /* Optimized implementation of plusdot fetching blocks of x */
  void _plusdot_block(Vec<I,D>& x, Vec<I,D>& y) const;

  /* Optimized implementation of plusdot fetching individual values from x */
  void _plusdot_single(Vec<I,D>& x, Vec<I,D>& y) const;


private:
  I _M, _N;
  I _local_rows;

  /* vector storing local (block diagonal) Mat elements */
  std::vector< std::vector< std::pair<I, D>>> _local;

  /* vector storing remote (off-diagonal) Mat elements */
  std::vector< std::vector< std::pair<I, D>>> _remote;

  std::vector<I> _row_partitions, _col_partitions;
};

/*==================================*/
/*** constructors and destructors ***/

template <typename I, typename D>
Mat<I, D>::Mat(I M, I N)
{
  set_dimensions(M, N);
}

/*============================*/
/*** dimensions and indices ***/

template <typename I, typename D>
void Mat<I, D>::set_dimensions(I M, I N)
{

  if (M <= 0 || N <= 0) {
    throw std::length_error("Mat dimensions must be > 0");
  }

  _M = M;
  _N = N;

  /* compute local rows */
  _row_partitions = partition_array(M);
  _col_partitions = partition_array(N);
  _local_rows = _row_partitions[upcxx::rank_me()+1] - _row_partitions[upcxx::rank_me()];

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
void Mat<I, D>::get_local_rows(I& start, I& end) const
{
  start = _row_partitions[upcxx::rank_me()];
  end = _row_partitions[upcxx::rank_me() + 1];
}

template <typename I, typename D>
void Mat<I, D>::get_diag_cols(I& start, I& end) const
{
  start = _col_partitions[upcxx::rank_me()];
  end = _col_partitions[upcxx::rank_me() + 1];
}

template <typename I, typename D>
I Mat<I, D>::get_local_rows_size() const
{
  I start, end;
  get_local_rows(start, end);
  return end-start;
}

/* reserve space for elements. */
template <typename I, typename D>
void Mat<I, D>::reserve(I nonzero_diag, I nonzero_offdiag)
{

  for (I i = 0; i < get_local_rows_size(); ++i) {
    _local[i].reserve(nonzero_diag);
    _remote[i].reserve(nonzero_offdiag);
  }
}

/* discard unused memory */
template <typename I, typename D>
void Mat<I, D>::shrink_extra()
{
  for (I i = 0; i < get_local_rows_size(); ++i) {
    _local[i].shrink_to_fit();
    _remote[i].shrink_to_fit();
  }
}

/*=========================================*/
/*** value setting and memory allocation ***/

template <typename I, typename D>
void Mat<I, D>::set_value(I row, I col, D value)
{
  I start, end;
  get_diag_cols(start, end);

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
  I row_start, row_end, col_start, col_end;
  get_local_rows(row_start, row_end);
  get_diag_cols(col_start, col_end);

#ifdef DEBUG
  assert(row >= row_start && row < row_end);
  assert(col >= col_start && col < col_end);
#endif

  _local[row-row_start].push_back(std::make_pair(col-col_start, value));
}

template <typename I, typename D>
void Mat<I, D>::set_offdiag_value(I row, I col, D value)
{
  I row_start, row_end;
  get_local_rows(row_start, row_end);

  #ifdef DEBUG
    I col_start, col_end;
    get_diag_cols(col_start, col_end);
    assert(row >= row_start && row < row_end);
    assert(col < col_start || col >= col_end);
  #endif

  _remote[row-row_start].push_back(std::make_pair(col, value));
}

/*============================*/
/*** matrix-vector products ***/

/* Mat-vector product y = A*x */
template <typename I, typename D>
void Mat<I, D>::dot(Vec<I,D>& x, Vec<I,D>& y) const
{
  y.set_all(0);
  /* TODO: choose the right implementation automatically based on sparsity */
  plusdot(x, y, DotImpl::Naive);
}

/* Mat-vector product y = A*x */
template <typename I, typename D>
void Mat<I, D>::dot(Vec<I,D>& x, Vec<I,D>& y, DotImpl impl) const
{
  y.set_all(0);
  plusdot(x, y, impl);
}

/* Mat-vector sum product y = A*x + y */
template <typename I, typename D>
void Mat<I,D>::plusdot(Vec<I,D>& x, Vec<I,D>& y) const
{
  plusdot(x, y, DotImpl::Naive);
}

/* Mat-vector sum product y = A*x + y */
template <typename I, typename D>
void Mat<I,D>::plusdot(Vec<I,D>& x, Vec<I,D>& y, DotImpl impl) const
{

  I M, N;
  get_dimensions(M, N);
  if (x.get_size() != N) {
    std::ostringstream out;
    out << "vector x size " << x.get_size() << " does not match ";
    out << "matrix row length " << N;
    throw std::invalid_argument(out.str());
  }
  if (y.get_size() != M) {
    std::ostringstream out;
    out << "vector y size " << y.get_size() << " does not match ";
    out << "matrix column length " << M;
    throw std::invalid_argument(out.str());
  }

  if (impl == DotImpl::Naive) {
    _plusdot_naive(x, y);
  }
  else if (impl == DotImpl::Block) {
    _plusdot_block(x, y);
  }
  else if (impl == DotImpl::Single) {
    _plusdot_single(x, y);
  }

}

/* Mat-vector sum product y = A*x + y */
template <typename I, typename D>
void Mat<I,D>::_plusdot_naive(Vec<I,D>& x, Vec<I,D>& y) const
{

  /* first do the local matvec */
  auto x_array = x.get_local_array_read();
  auto y_array = y.get_local_array();

  I local_size = get_local_rows_size();

  for (I i = 0; i < local_size; ++i) {
    for (const auto& p : _local[i]) {
      y_array[i] += p.second * x_array[p.first];
    }
  }

  /* now remote part */
  for (I i = 0; i < local_size; ++i) {
    for (const auto& p : _remote[i]) {
      /* x[p.first] implicitly gets the remote value */
      y_array[i] += p.second * x[p.first];
    }
  }

}

/* find the next value to prefetch, for _plusdot_single */
template <typename I, typename D>
static bool _prefetch_seek_next(I& pfch_row, I& pfch_row_idx, I local_size,
                 const std::vector< std::vector< std::pair<I, D>>>& _remote)
{
  bool done = false;

  pfch_row_idx++;

  /* move on to the next row when we hit the end of one */
  while (pfch_row_idx >= I(_remote[pfch_row].size())) {
    pfch_row++;
    if (pfch_row >= local_size) {
      done = true;
      break;
    }

    pfch_row_idx = 0;
  }

  return done;
}

/* Mat-vector sum product y = A*x + y */
template <typename I, typename D>
void Mat<I,D>::_plusdot_single(Vec<I,D>& x, Vec<I,D>& y) const
{

  bool prefetch_done = false;

  auto x_array = x.get_local_array_read();
  auto y_array = y.get_local_array();

  I local_size = get_local_rows_size();

  /* start fetching the x values for the first row of mat */
  std::vector<RData<I,D>> prefetched;
  I pfch_idx, pfch_read_idx = 0, pfch_row = 0, pfch_row_idx = -1;

  prefetched.resize(DOT_BLOCK_SIZE);
  for (pfch_idx = 0; pfch_idx < DOT_BLOCK_SIZE; ++pfch_idx) {
    prefetch_done = _prefetch_seek_next(pfch_row, pfch_row_idx, local_size, _remote);
    if (prefetch_done) break;

    prefetched[pfch_idx].update(x[_remote[pfch_row][pfch_row_idx].first].get_address());
    prefetched[pfch_idx].prefetch();
  }
  pfch_idx %= DOT_BLOCK_SIZE;

  /* do the local matvec while those values are on their way */

  for (I i = 0; i < local_size; ++i) {
    for (const auto& p : _local[i]) {
      y_array[i] += p.second * x_array[p.first];
    }
  }

  /* now remote part */
  D val;
  for (I i = 0; i < local_size; ++i) {
    I max_idx = _remote[i].size();
    for (I j = 0; j < max_idx; ++j) {

      /* get the prefetched value */
      val = prefetched[pfch_read_idx].get();

      /* prefetch the next one */
      if (!prefetch_done) {
        prefetch_done = _prefetch_seek_next(pfch_row, pfch_row_idx, local_size, _remote);
        if (!prefetch_done) {
          /* update it to point at the new address */
          prefetched[pfch_idx].update(x[_remote[pfch_row][pfch_row_idx].first].get_address());
          prefetched[pfch_idx].prefetch();
        }
      }

      y_array[i] += _remote[i][j].second * val;

      ++pfch_idx;
      ++pfch_read_idx;
      pfch_idx %= DOT_BLOCK_SIZE;
      pfch_read_idx %= DOT_BLOCK_SIZE;
    }
  }

}

/* Mat-vector sum product y = A*x + y */
template <typename I, typename D>
void Mat<I,D>::_plusdot_block(Vec<I,D>& x, Vec<I,D>& y) const
{

  /* TODO */

}
