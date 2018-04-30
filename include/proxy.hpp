/*
 *  This file is part of SLAPGAS
 *  (C) Greg Meyer, 2018
 */

/*
 * This file defines a remote class, which acts to
 * the user as if it is a data_t. However it abstracts
 * communication operations so that they
 * become implicit.
 */

#include <upcxx/upcxx.hpp>

template <typename idx_t, typename data_t>
class RData
{
  upcxx::global_ptr<data_t> addr;

  /* store a future for the get operation */
  upcxx::future<data_t> get_fut;
  bool fetched;

  /* future from the vector class that chains all put operations together */
  upcxx::future<> &put_fut;

public:
  RData(upcxx::global_ptr<data_t> addr,
        upcxx::future<> &put_fut)
        : addr(addr)
        , put_fut(put_fut)
          {};

  /* asynchronously request the data */
  void prefetch() {
    get_fut = upcxx::rget(addr);
    fetched = true;
  }

  /* return the data when it comes in */
  data_t get() {
    if (!fetched) {
      get_fut = upcxx::rget(addr);
      fetched = true;
    }

    return get_fut.wait();
  }

  /* implicitly retrieve remote data whenever we are cast to data_t */
  /* TODO: is there a way to make this not const, so I can use get_fut.wait() ? */
  operator data_t() const {
    return upcxx::rget(addr).wait();
  }

  /* set remote data with = */
  RData& operator= (const data_t &val) {
    /* chain to put_fut, which came from our vec */
    put_fut = upcxx::when_all(put_fut, upcxx::rput(val, addr));
    return *this;
  }
};
