#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include "Traccia.hpp"
#include "Utils.hpp"

using namespace std;
using namespace Eigen;
using vec3 = vector<vector<Vector3d>>;

int main()
{
    // Definizione del contenitore delle fratture e delle tracce.
    Fractures fractures_list;
    Traces traces_list;

    // Impostazione precisione dei dati in uscita.
    cout.precision(16);
    cout << scientific << endl;

    // Lettura del file.
    string path = "./DATA/FR3_data.txt";
    if(!importFractures(path, fractures_list)){
        cerr << "Errore nell'import." << endl;
        return 1;
    }
    cout << "Fratture importate correttamente." << endl;

    // Calcolo delle tracce
    Find_Traces(fractures_list, traces_list);

    // Esportazione tracce.
    bool result_info =  Export_traces_Info(traces_list);
    if(!result_info)
    {
        return 1;
    }
    cout << "Tracce esportate correttamente." << endl;

    /* Ordinamento della struttura salvavita E DELLA LUNGHEZZA
    Sort_Traces_Type(fractures_list, traces_list);

    bool result_type = Export_traces_Type(fractures_list,traces_list);
    if(!result_type)
    {
        return 1;
    } */
    return 0;
}


