#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include <utility>

using namespace std;
using namespace Eigen;

// Definisco la struttura che contiene le informazioni delle fratture.
struct Fractures
{
    unsigned int N_frac;
    vector<unsigned int> frac_id;
    vector<unsigned int> N_vert;
    vector<vector<Vector3d>> frac_vertices;
    vector<pair<vector<unsigned int>,vector<unsigned int>>> trace_type;
        // Struttura salvavita:
        //per ogni posizione del vettore esterno, che corrisponde all'id della frattura
        //salvo un vettore con le id delle tracce e i corrispetivi tips
};

// Definisco la struttura che contiene le informazioni delle tracce
struct Traces
{
    vector<unsigned int> traces_id;
    vector<vector<unsigned int>> traces_gen; // Generatori della traccia.
    vector<vector<Vector3d>> traces_points;
    vector<double> traces_length;
};
