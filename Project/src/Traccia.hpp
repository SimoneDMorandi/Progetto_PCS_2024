#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include <utility>

// NOME DA CAMBIARE VISTO CHE INLCUDIAMO ANCHE LA STRUCT PER LE FRATTURE

using namespace std;
using namespace Eigen;

struct Fractures
{
    unsigned int N_frac; // Numero delle fratture
    vector<unsigned int> frac_id;
    vector<unsigned int> N_vert;
    vector<vector<Vector3d>> frac_vertices;
    map<unsigned int, pair<vector<unsigned int>,vector<unsigned int>>> trace_type;
        // Struttura salvavita:
       //per ogni id del poligono salvo un vettore con le id delle tracce e i corrispetivi tips
        //
};

// Definisco la struttura che contiene le informazioni della traccia
struct Traces
{
    vector<unsigned int> traces_id;
    vector<vector<unsigned int>> traces_gen; // generatori della traccia
    vector<vector<Vector3d>> traces_points;
    vector<double> traces_length;
};
