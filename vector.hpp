/*
 *  GASMat
 *  (C) Greg Meyer, 2018
 */

#include <upcxx>
#include <stdexcept>

template <typename idx_t, typename data_t>
class Vec {
  std::vector<upcxx::global_ptr> gptrs;

public:
  Vec();
  Vec(idx_t size);
  ~Vec();

  idx_t get_size();
  void set_size(idx_t size);

}

template <typename idx_t, typename data_t>
Vec<idx_t, data_t>::Vec(idx_t size) {
  set_size(size);
}

template <typename idx_t, typename data_t>
Vec<idx_t, data_t>::~Vec() {
  if (!gptrs.empty()) /* we didn't call set_size */
    upcxx::_delete(gptrs[upcxx::rank_me()]);
  }
}

template <typename idx_t, typename data_t>
void Vec<idx_t, data_t>::set_size(idx_t size) {
  /* can only set the size once */
  if (!gptrs.empty()) {
    throw std::logic_error("Called set_size after size has already been set.");
  }
}
