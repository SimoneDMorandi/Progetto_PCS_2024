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
using vec3 = vector<vector<vector<double>>>;

int main()
{
    // Definizione del contenitore delle fratture e lettura del file con stampa.
    Fractures fractures_list;
    Traces traces_list;
    string path = "./DATA/FR82_data.txt";
    //string path = "C:/Users/ASUS/Desktop/esercitazioni/Progetto_PCS_2024/Project/DATA/FR3_data.txt";

    cout.precision(16);
    cout << scientific << endl;
    if(!importFractures(path, fractures_list)){
        cerr << "Errore nell'import." << endl;
        return false;
    }
    else
    {
        cout << "Import effettuato con successo." << endl;
    }

    Find_Traces(fractures_list, traces_list);
    /*
    bool result_info =  Export_traces_Info(traces_list);
    if(!result_info)
    {
        return 1;
    }
    cout << "Tracce esportate correttamente." << endl;

    // Ordinamento della struttura salvavita E DELLA LUNGHEZZA
    Sort_Traces_Type(fractures_list, traces_list);

    bool result_type = Export_traces_Type(fractures_list,traces_list);
    if(!result_type)
    {
        return 1;
    }
    cout << "Tipologia di tracce esportate correttamente." << endl; */
    return 0;
}


