// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include <limits>
#include <string>
//#include <bits/stdc++.h>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <typeinfo>
#include <vector>

#include <array>
#include <cmath>
#include <iomanip>
#include <iterator>
#include <set>
#include <unordered_map>
#include <unordered_set>

using namespace std;

// specify type of variables (double or float)
#define IMPALIB_TYPE double

#ifdef IMPALIB_TYPE
typedef IMPALIB_TYPE impalib_type;
#else
typedef float impalib_type;
#endif

const impalib_type zero_value = 0.0;
const impalib_type value_inf  = std::numeric_limits<impalib_type>::infinity();
