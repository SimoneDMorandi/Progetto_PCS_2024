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

    // Lettura del file.
    string path = "./DATA/FR3_data.txt";
    if(!importFractures(path, fractures_list)){
        cerr << "Errore nell'import." << endl;
        return 1;
    }
    cout << "Fratture importate correttamente." << endl;

    // Calcolo delle tracce.
    Find_Traces(fractures_list, traces_list);

    // Esportazione tracce.
    bool result_info =  Export_traces_Info(traces_list);
    if(!result_info)
    {
        return 1;
    }
    cout << "Tracce esportate correttamente." << endl;

    // Esportazione tipologia di tracce.
    Sort_Traces_Type(fractures_list,traces_list);
    bool result_type = Export_traces_Type(fractures_list,traces_list);
    if(!result_type)
    {
        return 1;
    }
    cout << "Tipologie di tracce esportate correttamente." << endl;

    // Esportazione delle fratture e delle tracce su paraview.
    Export_Paraview(fractures_list,traces_list);

    // PARTE 2.

    // Taglio dei poligoni rispetto alla tracce ordinate.
    vector<vector<Vector3d>> found_polygons;
    bool result_cut = cutPolygons(fractures_list,traces_list, found_polygons);
    if(!result_cut)
    {
        return 1;
    }
    cout << "Sottopoligoni esportati correttamente" << endl;

    // Esportazione dei sottopoligoni triangolati.
    Export_Paraview(found_polygons);

    return 0;
}


