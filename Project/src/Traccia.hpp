#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Eigen>

// NOME DA CAMBIARE VISTO CHE INLCUDIAMO ANCHE LA STRUCT PER LE FRATTURE

using namespace std;
using namespace Eigen;

struct Fractures
{
    unsigned int N_frac; // Numero delle fratture
    vector<int> frac_id;
    vector<int> N_vert;
    vector<vector<vector<double>>> frac_vertices;  //Vector3d
    vector<int> N_traces = {};
};

// Definisco la struttura che contiene le informazioni della traccia
struct Traces
{
    vector<int> traces_id;
    vector<vector<int>> traces_gen; // generatori della traccia
    vector<vector<Vector3d>> traces_points;
    vector<pair<int,double>> pass;
    vector<pair<int,double>> not_pass;
    vector<double> traces_length;
};
