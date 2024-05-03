#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "Traccia.hpp"
#include "Utils.hpp"

using namespace std;
using vec3 = vector<vector<vector<double>>>;


int main()
{
    vec3 fracture;
    string path = "C:/Users/ASUS/Desktop/esercitazioni/Progetto_PCS_2024/Project/DATA/FR3_data.txt";
    if(!importFractures(path, fracture)){
        cerr << "Errore nell'import" << endl;
        return false;
    }
    else
        cout << "Import effettuato con successo" << endl;

    return 0;
}

