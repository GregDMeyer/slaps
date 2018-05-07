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

/* how many values to prefetch */
#define DOT_BLOCK_SIZE 2048
#define NBUFS 2

template <typename I, typename D>
class Mat
{

public:
  /*==================================*/
  /*** constructors and destructors ***/

  Mat() {};

  /* construct a Mat with dimensions M, N */
  Mat(I M, I N) { set_dimensions(M, N); };

  /*============================*/
  /*** dimensions and indices ***/

  /* set dimensions of Mat. clears all values. */
  void set_dimensions(I M, I N);

  /* get dimensions of the Mat */
  void get_dimensions(I& M, I& N) const;

  /* check that x and y are compatible with the matrix */
  void check_dimensions(const Vec<I,D>& x, const Vec<I,D>& y) const;

  /* get the range of rows stored locally */
  void get_local_rows(I& start, I& end) const;

  /* get the range of cols corresponding to local part of the vector */
  void get_diag_cols(I& start, I& end) const;

  /* get the range of rows stored locally */
  I get_local_rows_size() const;

  /*=========================================*/
  /*** value setting and memory allocation ***/

  /* set a Mat element */
  void set_value(I row, I col, D value);

protected:
  I _M, _N;
  I _local_rows;

  bool size_set = false;

  /* vector storing COO (coordinate format) Mat elements, possibly out of order */
  std::vector< std::pair< std::pair<I,I>, D> > _elements;

  std::vector<I> _row_partitions, _col_partitions;
};

/* child classes */
/*
 * Matrix types defined here:
 *  - CSRMat
 *    -> NaiveCSRMat
 *    -> SingleCSRMat
 *    -> BlockCSRMat
 *  - RCMat
 */

template <typename I, typename D>
class CSRMat : public Mat<I,D>
{

public:
  /*==================================*/
  /*** constructors and destructors ***/

  CSRMat() {};

  /* construct a CSRMat with dimensions M, N */
  CSRMat(I M, I N) { this->set_dimensions(M, N); };

  /*=========================================*/
  /*** value setting and memory allocation ***/

  /* set up CSR storage format */
  /* optional arguments reserve memory for matrix elements:
   *  - dnz : number of diagonal nonzeros per row
   *  - onz : number of off-diagonal nonzeros per row
   */
  void setup(I dnz = 0, I onz = 0);

protected:

  /* vector storing local (block diagonal) Mat elements */
  std::vector< std::vector< std::pair<I, D>>> _local;

  /* vector storing remote (off-diagonal) Mat elements */
  std::vector< std::vector< std::pair<I, D>>> _remote;

  bool is_set_up = false;

};

template <typename I, typename D>
class NaiveCSRMat : public CSRMat<I,D>
{

public:
  /*==================================*/
  /*** constructors and destructors ***/

  NaiveCSRMat() {};

  /* construct a CSRMat with dimensions M, N */
  NaiveCSRMat(I M, I N) { this->set_dimensions(M, N); };

  /*============================*/
  /*** matrix-vector products ***/

  /* Mat-vector product y = A*x */
  void dot(Vec<I,D>& x, Vec<I,D>& y) const;

  /* Mat-vector sum product y = A*x + y */
  void plusdot(Vec<I,D>& x, Vec<I,D>& y) const;

};

template <typename I, typename D>
class SingleCSRMat : public CSRMat<I,D>
{

public:
  /*==================================*/
  /*** constructors and destructors ***/

  SingleCSRMat() {};

  /* construct a CSRMat with dimensions M, N */
  SingleCSRMat(I M, I N) { this->set_dimensions(M, N); };

  /*============================*/
  /*** matrix-vector products ***/

  /* Mat-vector product y = A*x */
  void dot(Vec<I,D>& x, Vec<I,D>& y) const;

  /* Mat-vector sum product y = A*x + y */
  void plusdot(Vec<I,D>& x, Vec<I,D>& y) const;

};

template <typename I, typename D>
class BlockCSRMat : public CSRMat<I,D>
{

public:
  /*==================================*/
  /*** constructors and destructors ***/

  BlockCSRMat() {};

  /* construct a CSRMat with dimensions M, N */
  BlockCSRMat(I M, I N) { this->set_dimensions(M, N); };

  /*============================*/
  /*** matrix-vector products ***/

  /* Mat-vector product y = A*x */
  void dot(Vec<I,D>& x, Vec<I,D>& y) const;

  /* Mat-vector sum product y = A*x + y */
  void plusdot(Vec<I,D>& x, Vec<I,D>& y) const;

};

/*
 * RCMat is a "row-partition column matrix". It's like CSC format, but the
 * matrix is still partitioned across processors by row.
 */

template <typename I, typename D>
class RCMat : public Mat<I,D>
{

public:
  /*==================================*/
  /*** constructors and destructors ***/

  RCMat() {};

  /* construct an RCMat with dimensions M, N */
  RCMat(I M, I N) { this->set_dimensions(M, N); };

  /*=========================================*/
  /*** value setting and memory allocation ***/

  /* set up CSR storage format */
  /* optional arguments reserve memory for matrix elements:
   *  - nnz : expected number of nonzeros per column (will be divided by # procs)
   */
  void setup(I nnz = 0);

  /*============================*/
  /*** matrix-vector products ***/

  /* Mat-vector product y = A*x */
  void dot(Vec<I,D>& x, Vec<I,D>& y) const;

  /* Mat-vector sum product y = A*x + y */
  void plusdot(Vec<I,D>& x, Vec<I,D>& y) const;

private:

  /* a vector storing the columns! */
  std::vector< std::pair<I, std::vector< std::pair<I, D>>>> _cols;

  bool is_set_up = false;

};

/* MAT */
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

  size_set = true;

}

/* set dimensions of Mat, clearing all values, and reserve space for elements. */
template <typename I, typename D>
void Mat<I, D>::get_dimensions(I& M, I& N) const
{
  M = _M;
  N = _N;
}

template <typename I, typename D>
void Mat<I, D>::check_dimensions(const Vec<I,D>& x, const Vec<I,D>& y) const
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

/*=========================================*/
/*** value setting and memory allocation ***/

template <typename I, typename D>
void Mat<I, D>::set_value(I row, I col, D value)
{
#ifdef DEBUG
  I rstart, rend;
  get_local_rows(rstart, rend);
  assert(row >= rstart && row < rend);
  assert(col >= 0 && col < _N);
#endif

  _elements.push_back(std::make_pair( std::make_pair(row, col), value ));
}

/*=====================*/
/* CSR MATRIX          */
/*=====================*/

template <typename I, typename D>
void CSRMat<I, D>::setup(I dnz, I onz)
{
  if (!this->size_set) {
    throw std::logic_error("Must set size before calling setup()");
  }
  if (is_set_up) {
    throw std::logic_error("Matrix already set up");
  }

  I rstart, rend;
  I cstart, cend;

  this->get_local_rows(rstart, rend);
  this->get_diag_cols(cstart, cend);

  _local.resize(rend-rstart);
  _remote.resize(rend-rstart);

  for (auto& row: _local) {
    row.reserve(dnz);
  }
  for (auto& row: _remote) {
    row.reserve(onz);
  }

  /* the plan: put the elements into _local and _remote, then sort the rows */
  /* this is more efficient than sorting first */

  for (auto& e: this->_elements) {
    I row = e.first.first, col = e.first.second;
    D val = e.second;
    if (col >= cstart && col < cend) {
      _local[row-rstart].push_back( std::make_pair(col-cstart, val) );
    }
    else {
      _remote[row-rstart].push_back( std::make_pair(col, val) );
    }
  }

  /* now sort them */
  for (auto& row: _local) {
    row.shrink_to_fit();
    std::sort(row.begin(), row.end(),
              [] (std::pair<I,D> const& a, std::pair<I,D> const& b) { return a.first < b.first; });
  }
  for (auto& row: _remote) {
    row.shrink_to_fit();
    std::sort(row.begin(), row.end(),
              [] (std::pair<I,D> const& a, std::pair<I,D> const& b) { return a.first < b.first; });
  }

  is_set_up = true;
}

/*=====================*/
/* CSRNAIVE MATRIX     */
/*=====================*/

/* Mat-vector product y = A*x */
template <typename I, typename D>
void NaiveCSRMat<I, D>::dot(Vec<I,D>& x, Vec<I,D>& y) const
{
  y.set_all(0);
  plusdot(x, y);
}

/* Mat-vector sum product y = A*x + y */
template <typename I, typename D>
void NaiveCSRMat<I,D>::plusdot(Vec<I,D>& x, Vec<I,D>& y) const
{

  this->check_dimensions(x, y);
  if (!this->is_set_up) {
    throw std::logic_error("Must set up matrix with ::setup() before calling ::plusdot");
  }

  /* first do the local matvec */
  auto x_array = x.get_local_array_read();
  auto y_array = y.get_local_array();

  I local_size = this->get_local_rows_size();

  for (I i = 0; i < local_size; ++i) {
    for (const auto& p : this->_local[i]) {
      y_array[i] += p.second * x_array[p.first];
    }
  }

  /* now remote part */
  for (I i = 0; i < local_size; ++i) {
    for (const auto& p : this->_remote[i]) {
      /* x[p.first] implicitly gets the remote value */
      y_array[i] += p.second * x[p.first].get();
    }
  }

}

/*=====================*/
/* CSRBLOCK MATRIX     */
/*=====================*/

/* Mat-vector product y = A*x */
template <typename I, typename D>
void BlockCSRMat<I, D>::dot(Vec<I,D>& x, Vec<I,D>& y) const
{
  y.set_all(0);
  plusdot(x, y);
}

/* Mat-vector sum product y = A*x + y */
template <typename I, typename D>
void BlockCSRMat<I,D>::plusdot(Vec<I,D>& x, Vec<I,D>& y) const
{

  this->check_dimensions(x, y);

  auto x_array = x.get_local_array_read();
  auto y_array = y.get_local_array();

  I local_size = this->get_local_rows_size();

  std::vector< std::vector<D> > bufs(NBUFS, std::vector<D>(DOT_BLOCK_SIZE));

  /* TODO: could probably come up with a better way to do this */
  std::vector<I> row_starts(this->_M, 0);
  I buf_start_idx = 0;
  int which_buf = 0;

  /* get the first block */
  x.read_range_begin(buf_start_idx, std::min(this->_N, buf_start_idx + DOT_BLOCK_SIZE), bufs[which_buf].data());

  /* do the local matvec while those values are on their way */
  for (I i = 0; i < local_size; ++i) {
    for (const auto& p : this->_local[i]) {
      y_array[i] += p.second * x_array[p.first];
    }
  }

  /* now remote part */
  while (buf_start_idx < this->_N) {

    /* finish getting previous */
    x.read_range_complete();

    /* fetch next block */
    if (buf_start_idx+DOT_BLOCK_SIZE < this->_N) {
      x.read_range_begin(buf_start_idx+DOT_BLOCK_SIZE,
        std::min(this->_N, buf_start_idx + 2*DOT_BLOCK_SIZE), bufs[(which_buf+1)%NBUFS].data());
    }

    /* our data is in bufs[which_buf] */

    for (I i = 0; i < local_size; ++i) {
      while (row_starts[i] < I(this->_remote[i].size()) && \
             this->_remote[i][row_starts[i]].first < buf_start_idx + DOT_BLOCK_SIZE) {
        I buf_idx = this->_remote[i][row_starts[i]].first - buf_start_idx;
        y_array[i] += this->_remote[i][row_starts[i]].second * bufs[which_buf][buf_idx];
        row_starts[i]++;
      }
    }

    which_buf++;
    which_buf %= NBUFS;

    buf_start_idx += DOT_BLOCK_SIZE;

  }
}

/*=====================*/
/* CSRSINGLE MATRIX    */
/*=====================*/

/* Mat-vector product y = A*x */
template <typename I, typename D>
void SingleCSRMat<I, D>::dot(Vec<I,D>& x, Vec<I,D>& y) const
{
  y.set_all(0);
  plusdot(x, y);
}

/* find the next value to prefetch, for BlockCSRMat::dot */
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
void SingleCSRMat<I,D>::plusdot(Vec<I,D>& x, Vec<I,D>& y) const
{

  this->check_dimensions(x, y);

  bool prefetch_done = false;

  auto x_array = x.get_local_array_read();
  auto y_array = y.get_local_array();

  I local_size = this->get_local_rows_size();

  /* start fetching the x values for the first row of mat */
  std::vector<RData<I,D>> prefetched;
  I pfch_idx, pfch_read_idx = 0, pfch_row = 0, pfch_row_idx = -1;

  prefetched.resize(DOT_BLOCK_SIZE);
  for (pfch_idx = 0; pfch_idx < DOT_BLOCK_SIZE; ++pfch_idx) {
    prefetch_done = _prefetch_seek_next(pfch_row, pfch_row_idx, local_size, this->_remote);
    if (prefetch_done) break;

    prefetched[pfch_idx].update(x[this->_remote[pfch_row][pfch_row_idx].first].get_address());
    prefetched[pfch_idx].prefetch();
  }
  pfch_idx %= DOT_BLOCK_SIZE;

  /* do the local matvec while those values are on their way */

  for (I i = 0; i < local_size; ++i) {
    for (const auto& p : this->_local[i]) {
      y_array[i] += p.second * x_array[p.first];
    }
  }

  /* now remote part */
  D val;
  for (I i = 0; i < local_size; ++i) {
    I max_idx = this->_remote[i].size();
    for (I j = 0; j < max_idx; ++j) {

      /* get the prefetched value */
      val = prefetched[pfch_read_idx].get();

      /* prefetch the next one */
      if (!prefetch_done) {
        prefetch_done = _prefetch_seek_next(pfch_row, pfch_row_idx, local_size, this->_remote);
        if (!prefetch_done) {
          /* update it to point at the new address */
          prefetched[pfch_idx].update(x[this->_remote[pfch_row][pfch_row_idx].first].get_address());
          prefetched[pfch_idx].prefetch();
        }
      }

      y_array[i] += this->_remote[i][j].second * val;

      ++pfch_idx;
      ++pfch_read_idx;
      pfch_idx %= DOT_BLOCK_SIZE;
      pfch_read_idx %= DOT_BLOCK_SIZE;
    }
  }
}

/*=====================*/
/* RC MATRIX           */
/*=====================*/

template <typename I, typename D>
void RCMat<I, D>::setup(I nnz)
{
  if (!this->size_set) {
    throw std::logic_error("Must set size before calling setup()");
  }
  if (is_set_up) {
    throw std::logic_error("Matrix already set up");
  }

  I rstart, rend;
  I cstart, cend;

  this->get_local_rows(rstart, rend);
  this->get_diag_cols(cstart, cend);

  /* guess half-filling to start */
  //_cols.reserve(this->_N/2);

  /* we need to sort all the accumulated elements first */
  /* sort by column, then by row, BUT shifted by the column start for this process */
  /* that way all the processes don't spam process 0 with requests at once */
  std::sort(this->_elements.begin(), this->_elements.end(),
            [&] (std::pair<std::pair<I,I>,D> const& a, std::pair<std::pair<I,I>,D> const& b)
               { return ((a.first.second+cstart)%this->_N) < ((b.first.second+cstart)%this->_N) || \
                 (a.first.second == b.first.second && a.first.first < b.first.first); });

  /* now they're sorted, so we can fill our thing */
  I cur_col = -1;
  I col_idx = -1;
  for (auto& e: this->_elements) {

    I row = e.first.first, col = e.first.second;
    D val = e.second;

    if (col != cur_col) {
      cur_col = col;

      if (col_idx != I(-1)) {
        _cols[col_idx].second.shrink_to_fit();
      }

      _cols.emplace_back();
      col_idx++;
      _cols[col_idx].first = col;
      _cols[col_idx].second.reserve(nnz/upcxx::rank_n());
    }

    /* TODO: if the user adds two identical elements, sum them */

    _cols[col_idx].second.push_back(std::make_pair(row-rstart, val));

  }

  if (col_idx != I(-1)) {
    _cols[col_idx].second.shrink_to_fit();
  }
  _cols.shrink_to_fit();

  is_set_up = true;
}

/* Mat-vector product y = A*x */
template <typename I, typename D>
void RCMat<I, D>::dot(Vec<I,D>& x, Vec<I,D>& y) const
{
  y.set_all(0);
  plusdot(x, y);
}

/* Mat-vector sum product y = A*x + y */
template <typename I, typename D>
void RCMat<I,D>::plusdot(Vec<I,D>& x, Vec<I,D>& y) const
{

  this->check_dimensions(x, y);
  if (!this->is_set_up) {
    throw std::logic_error("Must set up matrix with ::setup() before calling ::plusdot");
  }

  auto y_array = y.get_local_array();
  I local_size = this->get_local_rows_size();

  /* keep an array of prefetched values */
  std::vector<RData<I,D>> prefetched(DOT_BLOCK_SIZE);
  I pfch_idx, get_idx = 0;
  bool prefetch_done = false;

  /* get the initial set of prefetch values */
  for (pfch_idx = 0; pfch_idx < std::min(I(_cols.size()), I(DOT_BLOCK_SIZE)); ++pfch_idx) {
    prefetched[pfch_idx].update(x[_cols[pfch_idx].first].get_address());
    prefetched[pfch_idx].prefetch();
  }

  if (_cols.size() <= DOT_BLOCK_SIZE) {
    prefetch_done = true;
  }

  for (const auto& e: _cols) {

    D val = prefetched[get_idx].get();
    get_idx++;
    get_idx %= DOT_BLOCK_SIZE;

    if (!prefetch_done) {
      /* prefetch the next one */
      if (pfch_idx == I(_cols.size())) {
        prefetch_done = true;
      }
      else {
        prefetched[pfch_idx % DOT_BLOCK_SIZE].update(x[_cols[pfch_idx].first].get_address());
        prefetched[pfch_idx % DOT_BLOCK_SIZE].prefetch();
      }
      pfch_idx++;
    }

    /* actually do the multiplication for this value */
    for (const auto& p : e.second) {
      y_array[p.first] += p.second * val;
    }
  }

}
