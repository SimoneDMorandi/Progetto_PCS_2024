#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include "Structures.hpp"
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

    //Export_Paraview(fractures_list);

    // Calcolo delle tracce
    Find_Traces(fractures_list, traces_list);

    // Esportazione tracce.
    bool result_info =  Export_traces_Info(traces_list);
    if(!result_info)
    {
        return 1;
    }
    cout << "Tracce esportate correttamente." << endl;

    Sort_Traces_Type(fractures_list,traces_list);
    bool result_type = Export_traces_Type(fractures_list,traces_list);
    if(!result_type)
    {
        return 1;
    }

    //Export_Paraview(fractures_list,traces_list);

    Fractures found_polygons;
    bool result_cut = cutPolygons(fractures_list,traces_list, found_polygons);
    if(!result_cut)
    {
        return 1;
    }

    Export_Paraview(found_polygons);

    return 0;
}


