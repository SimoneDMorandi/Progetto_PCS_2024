#pragma once

#include<iostream>
#include<string>
#include<vector>
#include "Structures.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Parte 1
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Funzione che importa i dati delle fratture.
bool importFractures(const string& path, Fractures& fractures_list);

// Funzione che trova la traccia di due fratture.
void Find_Traces(Fractures& fractures_list, Traces& traces_list);

// Funzione che calcola il Bounding Box, date le coordinate di una frattura.
vector<Vector3d> Calculate_Bounding_Box(vector<Vector3d>& polygon);

// Funzione che calcola la retta passante per il lato della frattura in forma cartesiana, dati due vertici.
pair<Vector4d, Vector4d> equazioneRetta(const Vector3d& v1, const Vector3d& v2);

// Funzione che calcola il piano passante per una frattura, dati 3 vertici.
Vector4d pianoFrattura(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3);

// Funzione che verifica se la traccia è passante per una frattura, data la retta del lato e un punto.
bool check_pass(const Vector4d& pi1, const Vector4d& pi2, const Vector3d& point);

// Funzione che esporta le informazioni delle tracce.
bool Export_traces_Info(Traces& traces_list);

// Funzione che esporta le tipologie delle tracce.
bool Export_traces_Type(Fractures& f, Traces& traces_list);

// Funzione di ordinamento della struttura salvavita.
void Sort_Traces_Type(Fractures& f, Traces& t);

// Funzione che ordina un vettore fino ad un certo indice e cambia l'altro.
template<typename T>
void sort_pair(vector<T>& vec1, vector<unsigned int>& vec2);

// Paraview parte 1.
void Export_Paraview(Fractures& f, Traces& t);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Parte 2.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Funzione che prolunga una traccia fino ad incontrare i lati della frattura.
vector<Vector3d> extendTraceToEdges(vector<Vector3d>& frac_vertices, vector<Vector3d>& traces_points);

pair<vector<Vector3d>,vector<Vector3d>> subPolygons(vector<Vector3d> frac_vertices,
                                                     vector<Vector3d> traces_points,
                                                     unsigned int tip);

// Funzioni che tagliano i poligoni.
bool cutPolygons(Fractures& f, Traces& t, vector<vector<Vector3d>>& found_polygons);

// Paraview parte 2.
void Export_Paraview(vector<vector<Vector3d> > &subPolygons);
