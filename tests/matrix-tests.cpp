/*
 *  This file is part of SLAPS
 *  (C) Greg Meyer, 2018
 */

#include "slaps.hpp"
#include "catch.hpp"
#include <cmath>

#define IDX_T int
#define DATA_T float

#define MAT_T NaiveCSRMat
#include "matrix-tests-template.cpp"
#undef MAT_T

#define MAT_T SingleCSRMat
#include "matrix-tests-template.cpp"
#undef MAT_T

#define MAT_T BlockCSRMat
#include "matrix-tests-template.cpp"
#undef MAT_T

#define MAT_T RCMat
#include "matrix-tests-template.cpp"
#undef MAT_T

#undef IDX_T
#undef DATA_T

/******/

#define IDX_T unsigned long
#define DATA_T double

#define MAT_T NaiveCSRMat
#include "matrix-tests-template.cpp"
#undef MAT_T

#define MAT_T SingleCSRMat
#include "matrix-tests-template.cpp"
#undef MAT_T

#define MAT_T BlockCSRMat
#include "matrix-tests-template.cpp"
#undef MAT_T

#define MAT_T RCMat
#include "matrix-tests-template.cpp"
#undef MAT_T

#undef IDX_T
#undef DATA_T
