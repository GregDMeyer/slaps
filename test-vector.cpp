
#include "gasmat.hpp"
#include "catch.hpp"
#include <memory>

#define IDX_T int
#define DATA_T float
#include "test-vector-template.cpp"
#undef IDX_T
#undef DATA_T

#define IDX_T unsigned long
#define DATA_T double
#include "test-vector-template.cpp"
#undef IDX_T
#undef DATA_T
