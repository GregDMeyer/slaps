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
  bool fetched = false;

  /* future from the vector class that chains all put operations together */
  upcxx::future<>* put_fut_p = nullptr;

public:
  /* read-only with no put_fut */
  RData() {};

  RData(upcxx::global_ptr<data_t> addr) : addr(addr) {};

  RData(upcxx::global_ptr<data_t> addr,
        upcxx::future<> &put_fut)
        : addr(addr)
        , put_fut_p(&put_fut)
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

  /* set remote data with = */
  RData& operator= (const data_t &val) {

    assert(put_fut_p);

    /* chain to put_fut, which came from our vec */
    *put_fut_p = upcxx::when_all(*put_fut_p, upcxx::rput(val, addr));
    return *this;
  }

  /* get global address this points to */
  upcxx::global_ptr<data_t> get_address() {
    return addr;
  }

  /* update to point to a new address */
  void update (upcxx::global_ptr<data_t> new_addr) {
    addr = new_addr;
    fetched = false;
  }
};
