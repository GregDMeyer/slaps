/*
 *  This file is part of GASMat
 *  (C) Greg Meyer, 2018
 */

#include "gasmat.hpp"
#include "catch.hpp"
#include <memory>

#define IDX_T int
#define DATA_T float
#include "vector-tests-template.cpp"
#undef IDX_T
#undef DATA_T

#define IDX_T unsigned long
#define DATA_T double
#include "vector-tests-template.cpp"
#undef IDX_T
#undef DATA_T