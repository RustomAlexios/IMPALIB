// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
// https://opensource.org/licenses/MIT)

#pragma once

#include <limits>
#include <string>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numeric>
#include <set>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;

// specify type of variables (double or float)
#define IMPALIB_TYPE double

#ifdef IMPALIB_TYPE
typedef IMPALIB_TYPE impalib_type;
#else
typedef float impalib_type;
#endif

const impalib_type zero_value = 0.0;
const impalib_type value_inf = std::numeric_limits<impalib_type>::infinity();
