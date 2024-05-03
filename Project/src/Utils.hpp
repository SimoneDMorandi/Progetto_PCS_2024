#pragma once

#include<iostream>
#include<string>
#include<vector>

#include "Traccia.hpp"

using namespace std;
using vec3 = vector<vector<vector<double>>>;

bool importFractures(const string& path, vec3& v);
