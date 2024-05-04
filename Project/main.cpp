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

<<<<<<< Updated upstream
// Funzione che controlla che due fratture si intersechino e calcola i punti della traccia esistente
// Fracture è il nome dell'oggetto che contiene tutte le fratture.
void Find_Traces(Fracture& frac)
{
    // n è la dimensione letta nel file di avvio, al più tutto si interseca e ho n tracce
    Traces.traces_id.reserve(n);
    double epsilon = numeric_limits<decltype(epsilon)>::epsilon();
    // Prendo ogni singola frattura e calcolo il BoundBox
    // Boundbox è un vettore di due vettori di coordinate
    for(unsigned int i = 0; i < frac.size()-1; i++)
    {
        vector<vector<double>> Bounding_Box_1 = Calculate_Bounding_Box(frac(i));

        // Calcolo il secondo Bounding box e verifico l'intersezione
        for(unsigned int j = i+1; ;j++)
        {
            vector<vector<double>> Bounding_Box_2 = Calculate_Bounding_Box(frac(j));
            bool intersection =  CheckBbox(Bounding_Box_1,Bounding_Box_2);

            // Se c'è intersezione, procedo col calcolare i piani dei poligoni i e j e trovare i punti di intersezione.
            if (intersection)
            {
                Vector4d plane_parameters_1 = Calculate_Plane(frac(i)); //[a,b,c,-d]
                Vector4d plane_parameters_2 = Calculate_Plane(frac(j));
                // Sapendo che al più ho due punti di intersezione, so la dimensione del risultato
                vector<Vector3d> trace_points = solve_linear(plane_parameters_1,plane_parameters_2);
                // Controllo che l'intersezione non sia un punto
                if(trace_points[0] - trace_points[1] < epsilon)
                {
                    break;
                }
                vector<int> gen = {i,j};
                Traces.traces_id.append(i); // Numero della traccia
                Traces.traces_gen.append(gen);
                Traces.traces_points.append(trace_points);

                // Forse conviene memorizzare la distanza in un vettore.
            }
        }
    }
}

// Funzione che determina se la traccia è passante
void Check_Passing (Traces& traces, Fracture& fractures)
{
    for(unsigned int k = 0; k<traces.size(); k++)
    {
        // Verifica passante-non passante
        // Devo calcolare le rette passanti per due punti e vedere se un punto della traccia appartiene
        bool Tips;
        //if fracture[vector[0]] vector(i,j) = traces_gen[k] appartiene alle rette)
        // calcolo per l'altro punto
        // iffracture[vector[1]] vector(i,j) = traces_gen[k] appartiene alle rette)
        Tips = true;
        Traces.passing_traces.append(i,Tips);
        // altrimenti
        Tips = false;
        Traces.passing_traces.append(i,Tips);
    }

}


// Funzione che calcola il Bounding Box prendendo le coordinate di un poligono.
vector Calculate_Bounding_Box(vector& polygon)
{
    // Inizializzo le coordinate del punto massimo e minimo, la prima colonna contiene il punto minimo
    vector<vector<double>> Bbox = {polygon[0],polygon[0]};

    // Itero per conoscere le coordinate del bounding box
    for (const auto& coord : polygon )
    {
        (Bbox[0])[0] = min((Bbox[0])[0], coord[0]);
        (Bbox[1])[0] = max((Bbox[1])[0], coord[0]);

        (Bbox[0])[1] = min((Bbox[0])[1], coord[0]);
        (Bbox[1])[1] = max((Bbox[1])[1], coord[0]);

        (Bbox[0])[2] = min((Bbox[0])[2], coord[0]);
        (Bbox[1])[2] = max((Bbox[1])[2], coord[0]);
    }

    return Bbox;
}

// Funzione che controlla l'intersezione di due BoundingBox
bool CheckBbox( vector&  Bbox_1, vector Bbox_2, )
{
    bool overlap_x = (Bbox_1[1])[0] >= (Bbox_2[0])[0] && (Bbox_2[1])[0] >= (Bbox_1[0])[0];
    bool overlap_y= (Bbox_1[1])[1] >= (Bbox_2[0])[1] && (Bbox_2[1])[1] >= (Bbox_1[0])[1];
    bool overlap_z = (Bbox_1[1])[2] >= (Bbox_2[0])[2] && (Bbox_2[1])[2] >= (Bbox_1[0])[2];
    if (overlap_x && overlap_y && overlap_z)
    {
        return true;
    }

    return false;
}

// Funzione che calcola i parametri dell'equazione del piano contenente un poligono
Vector4d Calculate_Plane(vector& polygon)
{
    VectorXd cross = polygon[0].cross(polygon[1]);
    cross[3] = -polygon[0].lpNorm<1>();
    return cross;
}

// Funzione che risolve il sistema lineare tra due piani contenenti i poligoni.

vector solve_linear(Vector4d& coeff1, Vector4d& coeff2)
{
    direction = coeff1[03].cross(coeff2[03]); // direzione della retta di intersezione

}

struct Traces
{
    vector<int> traces_id;
    vector<vector<int>> traces_gen; // generatori della traccia
    vector<vector<Vector3d>> traces_points;
    map<int,bool> passing_traces; // id_traccia tips 1/0

};
=======
vector<double> vettore(vec3 fracture, vector <double> v, unsigned int i, unsigned int j) // i indica quale frattura considero
{
    v = { fracture[i][j][0] - fracture[i][j-1][0],
          fracture[i][j][1] - fracture[i][j-1][1],
          fracture[i][j][2] - fracture[i][j-1][2] };

    return v;
}

vector<double> crossProduct(const vector<double>& u, const vector<double>& v) {
    vector<double> w(3);
    w[0] = u[1] * v[2] - u[2] * v[1];
    w[1] = u[2] * v[0] - u[0] * v[2];
    w[2] = u[0] * v[1] - u[1] * v[0];
    return w;
}

double dotProduct(const vector<double>& v1, const vector<double>& v2)
{
    double ris = 0;
    for (int i = 0; i < 3; ++i) {
        ris += v1[i] * v2[i];
    }
    return -ris;
}

>>>>>>> Stashed changes

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

    // vediamo come calcolare i punti di intersezione tra due poligoni
    // scelgo frattura 0 e frattura 2

    vector <double> AB(3);
    vector <double> BC(3);
    vector <double> LM(3);
    vector <double> MN(3);

    AB = vettore(fracture, AB, 0, 1);
    BC = vettore(fracture, BC, 0, 2);
    LM = vettore(fracture, LM, 2, 1);
    MN = vettore(fracture, MN, 2, 2);

    vector <double> n1 = crossProduct(AB, BC);
    vector <double> n2 = crossProduct(LM, MN);



    // piano per frattura 0 identificato da (a,b,c,d)
    double d1 = dotProduct(n1,fracture[0][0]);
    vector <double> piano0 = n1;
    piano0.push_back(d1);

    // piano per frattura 2 identificato da (a,b,c,d)
    double d2 = dotProduct(n2,fracture[2][0]);
    vector <double> piano2 = n2;
    piano2.push_back(d2);

    // retta per AD, identificata dal punto A e il vettore AD
    vector <double> AD = { fracture[0][3][0] - fracture[0][0][0],
                         fracture[0][3][1] - fracture[0][0][1],
                         fracture[0][3][2] - fracture[0][0][2] };

    // converto la forma parametrica di AD in forma cartesiana, quindi trovo i due piani nella forma di prima
    // x = 0, z = 0 è la forma cartesiana di questa retta
    vector <double> pi1(4);
    vector <double> pi2(4);
    vector <double> pi3(4);


    pi1[0] = AD[1];
    pi1[1] = -AD[0];
    pi1[3] = AD[0]*fracture[0][0][1] - AD[1]*fracture[0][0][0];


    pi2[0] = AD[2];
    pi2[2] = -AD[0];
    pi2[3] = AD[0]*fracture[0][0][2] - AD[2]*fracture[0][0][0];

    pi3[1] = AD[2];
    pi3[2] = -AD[1];
    pi3[3] = AD[1]*fracture[0][0][2] - AD[2]*fracture[0][0][1];

    // ora scivo la matrice dei coefficienti fatta da piano0, piano 2, pi1, pi3 e vettore delle incognite (x,y,z,1),
    // l'intersezione sarà (x,y,z)

    Matrix<double, 4, 3> coeff;
    coeff.row(0) << piano0[0], piano0[1], piano0[2];
    coeff.row(1) << piano2[0], piano2[1], piano2[2];
    coeff.row(2) << pi1[0], pi1[1], pi1[2];
    coeff.row(3) << pi3[0], pi3[1], pi3[2];

    Vector4d termineNoto;
    termineNoto[0] = piano0[3];
    termineNoto[1] = piano2[3];
    termineNoto[2] = pi1[3];
    termineNoto[3] = pi3[3];

    HouseholderQR<Eigen::Matrix<double, 4, 3>> qr(coeff);
    Vector3d soluzione = qr.solve(termineNoto);
    for(unsigned int i = 0; i < 3; i++)
        cout << soluzione[i] << " ";
    cout << endl;












    return 0;

