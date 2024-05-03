#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "Traccia.hpp"
#include "Utils.hpp"

using namespace std;
using vec3 = vector<vector<vector<double>>>;

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

