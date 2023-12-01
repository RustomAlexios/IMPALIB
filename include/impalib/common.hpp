// Copyright 2023, Alexios Rustom.
// https://github.com/RustomAlexios/IMPALIB
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include <limits>
#include <string>
//#include <bits/stdc++.h>
#include <typeinfo>
#include <chrono>
#include <iostream>
#include <vector>
#include<fstream>
#include <algorithm>
#include <numeric>
#include <functional>

#include <cmath>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <iterator>
#include <iomanip> 

using namespace std;

#define IMPALIB_TYPE double



#ifdef IMPALIB_TYPE
    typedef IMPALIB_TYPE impalib_type;
#else
    typedef float impalib_type;
#endif

impalib_type zero_value = 0.0;
impalib_type value_inf = std::numeric_limits<impalib_type>::infinity();

