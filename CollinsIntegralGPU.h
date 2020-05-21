#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include <Windows.h>
#include <vector>
#include <complex>
#include "cuda_runtime.h"
#include "BMP.h"

constexpr auto PI = 3.1415926535897932384626433832795;

using namespace std;

inline void error(const string& s)
{
	throw runtime_error(s);
}