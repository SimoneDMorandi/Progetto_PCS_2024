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

//**************************************************************************************************************

// Funzione che calcola il piano che contiene una certa frattura, prende in input 3 punti e restituisce un vettore 4x1
// piano identificato da (a,b,c,d) => ax+by+cz+d=0
vector<double> pianoFrattura(const vector<double>& v1, const vector<double>& v2, const vector<double>& v3)
{
    vector<double> AB = sottrazione(v1, v2);
    vector<double> AC = sottrazione(v1, v3);
    vector <double> n1 = crossProduct(AB, AC); // vettore normale al piano

    double d = - dotProduct(n1,v1);
    vector <double> piano = n1;
    piano.push_back(d);

    return piano;
}

// Funzione che dati 4 vettori che rappresentano i piano, restituisce il punto di intersezione
Vector3d soluzione(const vector<double>& piano1, const vector<double>& piano2, const vector<double>& pi1, const vector<double>& pi2)
{
    Matrix<double, 4, 3> coeff;
    coeff.row(0) << piano1[0], piano1[1], piano1[2];
    coeff.row(1) << piano2[0], piano2[1], piano2[2];
    coeff.row(2) << pi1[0], pi1[1], pi1[2];
    coeff.row(3) << pi2[0], pi2[1], pi2[2];

    Vector4d termineNoto;
    termineNoto[0] = -piano1[3];
    termineNoto[1] = -piano2[3];
    termineNoto[2] = -pi1[3];
    termineNoto[3] = -pi2[3];

    HouseholderQR<Matrix<double, 4, 3>> qr(coeff);
    Vector3d sol = qr.solve(termineNoto);

    return sol;
}


//>>>>>>> Stashed changes

int main()
{
    // Definizione del contenitore delle fratture e lettura del file con stampa.
    Fractures fractures_list;
    Traces traces_list;
    string path = "C:/Users/simod/Dropbox (Politecnico Di Torino Studenti)/PC/Desktop/Progetto_PCS_2024/Project/DATA/FR3_data.txt";
    if(!importFractures(path, fractures_list)){
        cerr << "Errore nell'import." << endl;
        return false;
    }
    else
        cout << "Import effettuato con successo." << endl;

    Find_Traces(fractures_list, traces_list);
    bool result_info =  Export_traces_Info(traces_list);
    if(!result_info)
    {
        return 1;
    }
    cout << "Tracce esportate correttamente." << endl;

    bool result_type = Export_traces_Type(fractures_list,traces_list);
    if(!result_type)
    {
        return 1;
    }
    cout << "Tipologia di tracce esportate correttamente." << endl;
    return 0;
}


