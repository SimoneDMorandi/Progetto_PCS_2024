#pragma once

#include<iostream>
#include<string>
#include<vector>

#include "Traccia.hpp"

    using namespace std;

// Funzione che importa i dati delle fratture
bool importFractures(const string& path, Fractures& fractures_list);

// Funzione che trova la traccia di due fratture
void Find_Traces(Fractures& fractures_list, Traces& traces_list);

// FUnzione che esporta le informazioni delle tracce
bool Export_traces_Info(Traces& traces_list);

// FUnzione che esporta le tipologie delle tracce
bool Export_traces_Type(Fractures &f, Traces& traces_list);

// Funzione che calcola la retta passante per il lato in forma cartesiana.
void equazioneRetta(const vector<double>& v1, const vector<double>& v2, vector<double>& pi1, vector<double>& pi2);

vector<double> pianoFrattura(const vector<double>& v1, const vector<double>& v2, const vector<double>& v3);

vector<double> crossProduct(const vector<double>& u, const vector<double>& v);

double dotProduct(const vector<double>& v1, const vector<double>& v2);

vector<double> sottrazione(const vector<double>& v1, const vector<double>& v2);

// Funzione che calcola il Bounding Box, date le coordinate di un poligono.
vector<vector<double>> Calculate_Bounding_Box(vector<vector<double>>& polygon);

// Funzione di ordinamento della struttura salvavita.
bool Sort_Traces_Type(Fractures&fractures_list);
