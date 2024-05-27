#pragma once

#include<iostream>
#include<string>
#include<vector>
#include "Structures.hpp"
#include <cassert> // MergeSort

using namespace std;

// Funzione che importa i dati delle fratture.
bool importFractures(const string& path, Fractures& fractures_list);

// Funzione che trova la traccia di due fratture.
void Find_Traces(Fractures& fractures_list, Traces& traces_list);

// Funzione che calcola il Bounding Box, date le coordinate di una frattura.
vector<Vector3d> Calculate_Bounding_Box(vector<Vector3d> &polygon);

// Funzione che calcola la retta passante per il lato della frattura in forma cartesiana, dati due vertici.
pair<Vector4d, Vector4d> equazioneRetta(const Vector3d& v1, const Vector3d& v2);

// Funzione che calcola il piano passante per una frattura, dati 3 vertici.
Vector4d pianoFrattura(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3);

// Funzione che verificase lè passante per una frattura, data la retta del lato e un punto.
bool check_pass(const Vector4d& pi1, const Vector4d& pi2, const Vector3d& point, const double &eps);

// Funzione che esporta le informazioni delle tracce.
bool Export_traces_Info(Traces& traces_list);

// Funzione che esporta le tipologie delle tracce.
bool Export_traces_Type(Fractures &f, Traces& traces_list);

// Funzione di ordinamento della struttura salvavita.
void Sort_Traces_Type(Fractures&f, Traces& t);

// Funzione che calcola il centroide
Vector3d calculateCentroid(const vector<Vector3d>& points);

// Funzione che calcola l'angolo polare
Vector2d calculatePolarAngles(const Vector3d& vec);

// Funzione che prolunga una traccia fino ad incontrare i lati della frattura
vector<Vector3d> extendTraceToEdges(vector<Vector3d>& frac_vertices, vector<Vector3d>& traces_points);
