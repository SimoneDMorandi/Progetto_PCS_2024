#include "Utils.hpp"
#include "Structures.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <iomanip>
#include <utility> // Pairs
#include <stdexcept> // Throw
#include <numeric> // iota
#include <algorithm>
#include <queue>

// Paraview
#include "UCDUtilities.hpp"
using namespace Gedim;

using namespace std;
using namespace Eigen;
using vec2 = vector<Vector3d>;
using vec3 = vector<vector<Vector3d>>;


//////////////////////////////////////////////////////////////////////////////////////
// PARTE 1
/////////////////////////////////////////////////////////////////////////////////////

// Funzione di lettura del file delle tracce.
bool importFractures(const string& path, Fractures &fractures_list) //OK.
{
    ifstream file;
    file.open(path);

    if(file.fail())
        return false;

    string line;
    getline(file, line);

    getline(file, line);
    unsigned int N = 0;

    N = stoi(line);

    unsigned int idFrac = 0;
    char delimiter = ';';
    unsigned int numVertices = 0;

    fractures_list.N_frac = N;
    fractures_list.frac_vertices.resize(N);
    fractures_list.frac_id.resize(N);
    fractures_list.N_vert.resize(N);
    fractures_list.trace_type.resize(N);
    for(unsigned int i = 0; i < N; i++)
    {
        fractures_list.trace_type[i].first.reserve(N - 1);
        fractures_list.trace_type[i].second.reserve(N - 1);
        getline(file, line);
        getline(file, line);

        istringstream converter(line);
        converter >> idFrac >> delimiter >> numVertices;

        fractures_list.frac_id[i] = idFrac;
        fractures_list.N_vert[i] = numVertices;
        fractures_list.frac_vertices[i].resize(numVertices);

        getline(file, line);


        for(unsigned int k = 0; k < 3; k++)
        {
            getline(file, line);

            istringstream converter(line);
            for(unsigned int j = 0; j < numVertices; j++)
            {
                double val;
                converter >> val >> delimiter;
                if(k == 0)
                {
                    fractures_list.frac_vertices[i][j].x() = val;
                }
                else if (k==1)
                {
                    fractures_list.frac_vertices[i][j].y() = val;
                }
                else
                {
                    fractures_list.frac_vertices[i][j].z() = val;
                }
            }
        }
    }
    file.close();
    return true;
}

///////////////////////////////////////////////////

// Funzione che calcola le tracce.
void Find_Traces(Fractures &fractures_list, Traces& traces_list)
{
    // N_frac è la dimensione letta nel file di avvio, al più tutto si interseca e ho (N_frac^2)/2 tracce.
    unsigned int max_N = (fractures_list.N_frac*fractures_list.N_frac)/2;
    traces_list.traces_id.reserve(max_N);
    traces_list.traces_gen.reserve(max_N);
    traces_list.traces_points.reserve(max_N);
    traces_list.traces_length.reserve(max_N);

    // Utility per id tracce.
    unsigned int count_traces = 0;

    // Definisco la precisione di macchina.
    double eps = numeric_limits<decltype(eps)>::epsilon(); // Precisione 1D
    double tau = 1e-12; // Precisione 2D

    // Prendo ogni singola frattura e calcolo il Bounding Box.
    for(unsigned int i = 0; i < fractures_list.N_frac - 1 ; i++)
    {
        vec2 Bbox_1 = Calculate_Bounding_Box(fractures_list.frac_vertices[i]);

        // Calcolo il secondo Bounding Box, per ogni frattura successiva a "i".
        for(unsigned int j = i+1; j < fractures_list.N_frac ;j++)
        {
            vec2 Bbox_2 = Calculate_Bounding_Box(fractures_list.frac_vertices[j]);

            // Verifico l'intersezione tra i Bounding Box.
            bool overlap_x = Bbox_1[1][0]>Bbox_2[0][0] && Bbox_2[1][0]>Bbox_1[0][0];
            bool overlap_y = Bbox_1[1][1]>Bbox_2[0][1] && Bbox_2[1][1]>Bbox_1[0][1];
            bool overlap_z = Bbox_1[1][2]>Bbox_2[0][2] && Bbox_2[1][2]>Bbox_1[0][2];

            // Se c'è intersezione, procedo col calcolare i piani delle fratture i e j e trovare la retta di intersezione.
            if(overlap_x && overlap_y && overlap_z)
            {
                // Calcolo i piani delle fratture.
                Vector4d plane_1 = pianoFrattura(fractures_list.frac_vertices[i][0], fractures_list.frac_vertices[i][1],
                                                        fractures_list.frac_vertices[i][2]);
                Vector4d plane_2 = pianoFrattura(fractures_list.frac_vertices[j][0], fractures_list.frac_vertices[j][1],
                                                        fractures_list.frac_vertices[j][2]);

                // Controllo che i piani non siano paralleli. In caso contrario, ignoro.
                if((plane_1.head<3>().cross(plane_2.head<3>())).lpNorm<1>()<tau)
                {
                    continue;
                }

                // Calcolo i punti di intersezione tra i lati delle fratture e la retta contenente la traccia.
                vector<Vector3d> points;
                unsigned int count = 0;
                points.resize(4);

                // Ricavo segmento della frattura 1: intersezione tra la retta contenente la traccia e due suoi lati.
                for (unsigned int k = 0; k < fractures_list.N_vert[i]; k++)
                {
                    // Calcolo l'equazione del singolo lato.
                    pair<Vector4d,Vector4d> planes;
                    if(k == fractures_list.N_vert[i] -1)
                    {
                        try
                        {
                            planes = equazioneRetta(fractures_list.frac_vertices[i][0] ,fractures_list.frac_vertices[i][k]);
                        }catch(...)
                        {
                                continue;
                        }
                    }
                    else{
                        try
                        {
                            planes = equazioneRetta(fractures_list.frac_vertices[i][k] ,fractures_list.frac_vertices[i][k+1]);
                        }catch(...)
                        {
                            continue;
                        }
                    }
                    Matrix<double,4,3> coeff;
                    coeff.row(0) << plane_1[0], plane_1[1], plane_1[2];
                    coeff.row(1) << plane_2[0], plane_2[1], plane_2[2];
                    coeff.row(2) << planes.first[0], planes.first[1], planes.first[2];
                    coeff.row(3) << planes.second[0], planes.second[1],planes.second[2];
                    FullPivLU<MatrixXd> lu_decomp(coeff);
                    lu_decomp.setThreshold(eps);
                    if (lu_decomp.rank() == 3)
                    {
                        // Se c'è intersezione con la retta della frattura, salvo il punto risolvendo il sistema lineare.
                        Vector4d termineNoto;
                        termineNoto[0] = -plane_1[3];
                        termineNoto[1] = -plane_2[3];
                        termineNoto[2] = -planes.first[3];
                        termineNoto[3] = -planes.second[3];
                        HouseholderQR<MatrixXd> qr(coeff);
                        // Controllo che il punto sia dentro la frattura.
                        Vector3d sol = qr.solve(termineNoto);
                        bool overLap_x = Bbox_1[0][0] - sol[0] <= eps && sol[0]-Bbox_1[1][0] <= eps;
                        bool overLap_y = Bbox_1[0][1] - sol[1] <= eps && sol[1]-Bbox_1[1][1] <= eps;
                        bool overLap_z = Bbox_1[0][2] - sol[2] <= eps && sol[2]-Bbox_1[1][2] <= eps;
                        if(overLap_x && overLap_y && overLap_z)
                        {
                            points[count] = sol;
                            count++;
                        }
                        else
                        {
                            continue;
                        }
                    }
                    else
                    {
                        continue;
                    }
                } // Se due punti uguali, considerane solo uno
                // se i due punti stanno sulla stessa retta, ignora entrambi
                // Ricavo segmento della frattura 2: intersezione tra la retta contenente la traccia e due suoi lati.
                for (unsigned int k = 0; k < fractures_list.N_vert[j]; k++)
                {
                    pair<Vector4d, Vector4d> planes;
                    if(k == fractures_list.N_vert[j] -1)
                    {
                        try
                        {
                        planes = equazioneRetta(fractures_list.frac_vertices[j][0] ,fractures_list.frac_vertices[j][k]);

                        }catch(...)
                        {
                            continue;
                        }
                    }
                    else{
                        try{
                        planes = equazioneRetta(fractures_list.frac_vertices[j][k] ,fractures_list.frac_vertices[j][k+1]);    
                        }catch(...)
                        {
                            continue;
                        }
                    }

                    // Intersezione tra i 3 piani
                    Matrix<double,4,3> coeff;
                    coeff.row(0) << plane_1[0], plane_1[1], plane_1[2];
                    coeff.row(1) << plane_2[0], plane_2[1], plane_2[2];
                    coeff.row(2) << planes.first[0],planes.first[1],planes.first[2];
                    coeff.row(3) <<  planes.second[0], planes.second[1], planes.second[2];
                    FullPivLU<MatrixXd> lu_decomp(coeff);
                    lu_decomp.setThreshold(eps);
                    if (lu_decomp.rank() == 3)
                    {
                        Vector4d termineNoto;
                        termineNoto[0] = -plane_1[3];
                        termineNoto[1] = -plane_2[3];
                        termineNoto[2] = -planes.first[3];
                        termineNoto[3] = -planes.second[3];
                        HouseholderQR<MatrixXd> qr(coeff);
                        // Controllo che il punto sia dentro la frattura.
                        Vector3d sol = qr.solve(termineNoto);
                        bool overLap_x = Bbox_2[0][0] - sol[0] <= eps && sol[0]-Bbox_2[1][0] <= eps;
                        bool overLap_y = Bbox_2[0][1] - sol[1] <= eps && sol[1]-Bbox_2[1][1] <= eps;
                        bool overLap_z = Bbox_2[0][2] - sol[2] <= eps && sol[2]-Bbox_2[1][2] <= eps;
                        if(overLap_x && overLap_y && overLap_z)
                        {
                            points[count] = sol;
                            count++;
                        }
                        else
                        {
                            continue;
                        }
                    }
                    else
                    {
                        continue;
                    }
                }

                // Caso di traccia passante per la prima e per la seconda frattura.
                if( count == 4 && ((points[1]-points[3]).lpNorm<1>() < tau && (points[0]-points[2]).lpNorm<1>() < tau)
                     ||((points[1]-points[2]).lpNorm<1>() < tau && (points[0]-points[3]).lpNorm<1>() < tau))
                {
                    // Completo la struttura TRACES.
                    traces_list.traces_id.push_back(count_traces);
                    traces_list.traces_points.push_back({points[0], points[1]});
                    traces_list.traces_length.push_back((points[1]-points[0]).lpNorm<1>());
                    traces_list.traces_gen.push_back({i,j});
                    // Completo la struttura FRACTURES salvavita per entrambe le fratture.
                    fractures_list.trace_type[i].first.push_back(traces_list.traces_id[count_traces]);
                    fractures_list.trace_type[i].second.push_back(0);
                    fractures_list.trace_type[j].first.push_back(traces_list.traces_id[count_traces]);
                    fractures_list.trace_type[j].second.push_back(0);
                    count_traces++;
                }

                // Caso di traccia non passante per almeno una frattura
                 /* Calcolo del segmento corrispondente alla traccia, avendo i 4 punti [A,C//B,D]
                 Visualmente, se ho [A,C] e [B,C] tutti sulla stessa retta, i primi appartenenti alla traccia
                 1 e gli altri alla traccia 2, se calcolo il minimo tra gli inizi trovo l'inizio della traccia
                  e calcolando il massimo tra le code trovo la fine della traccia, rappresentata da [Start,Finish]*/
                else if (count == 4)
                {
                    sort(points.begin(),points.end(), [eps](const Vector3d& a, const Vector3d& b){
                        if(abs(a.x() -b.x()) > eps)
                            return a.x() < b.x();
                        else if (abs(a.y()-b.y()) > eps)
                            return a.y() < b.y();
                        else
                            return a.z() < b.z();
                    });

                    if((points[2]-points[1]).lpNorm<1>() > eps)
                    {
                        // Completo la struttura TRACES
                        traces_list.traces_id.push_back(count_traces);
                        traces_list.traces_points.push_back({points[1], points[2]});
                        traces_list.traces_length.push_back((points[2]-points[1]).lpNorm<1>());
                        traces_list.traces_gen.push_back({i,j});

                        // Controllo se la traccia è passante per almeno una frattura.
                        // Prima frattura
                        unsigned int count_pass = 0;
                        for (unsigned int k = 0; k < fractures_list.N_vert[i]; k++)
                        {
                            pair<Vector4d, Vector4d> planes;
                            if(k == fractures_list.frac_vertices[i].size() -1)
                            {
                                try
                                {
                                    planes = equazioneRetta(fractures_list.frac_vertices[i][0] ,fractures_list.frac_vertices[i][k]);
                                }catch(...)
                                {
                                    continue;
                                }
                            }
                            else{
                                try
                                {
                                    planes = equazioneRetta(fractures_list.frac_vertices[i][k] ,fractures_list.frac_vertices[i][k+1]);
                                }catch(...)
                                {
                                    continue;
                                }
                            }

                            if(check_pass(planes.first,planes.second,points[1],eps))
                            {
                                count_pass ++;
                            }
                            if(check_pass(planes.first,planes.second,points[2],eps))
                            {
                                count_pass ++;
                            }
                        }
                        if(count_pass == 2) // Caso di traccia passante per il primo poligono e non per il secondo
                        {
                            fractures_list.trace_type[i].first.push_back(traces_list.traces_id[count_traces]);
                            fractures_list.trace_type[i].second.push_back(0);
                            fractures_list.trace_type[j].first.push_back(traces_list.traces_id[count_traces]);
                            fractures_list.trace_type[j].second.push_back(1);
                            count_traces++;
                        }
                        else  // Caso di traccia non passante per il primo ma forse per il secondo
                        {
                            // Riempio Fractures per il poligono i
                            fractures_list.trace_type[i].first.push_back(traces_list.traces_id[count_traces]);
                            fractures_list.trace_type[i].second.push_back(1);

                            // Controllo se passante per j
                            count_pass = 0;
                            for (unsigned int k = 0; k < fractures_list.N_vert[j]; k++)
                            {
                                pair<Vector4d, Vector4d> planes;
                                if(k == fractures_list.frac_vertices[j].size() -1)
                                {
                                    try
                                    {
                                        planes = equazioneRetta(fractures_list.frac_vertices[j][0] ,fractures_list.frac_vertices[j][k]);
                                    }catch(...){
                                        continue;
                                    }
                                }
                                else{
                                    try
                                    {
                                        planes = equazioneRetta(fractures_list.frac_vertices[j][k] ,fractures_list.frac_vertices[j][k+1]);
                                    }catch(...){
                                        continue;
                                    }
                                }
                                if(check_pass(planes.first,planes.second,points[1],eps))
                                {
                                    count_pass ++;
                                }
                                if(check_pass(planes.first,planes.second,points[2],eps))
                                {
                                    count_pass ++;
                                }
                            }
                            if(count_pass == 2) // Caso di traccia passante per j
                            {
                                fractures_list.trace_type[j].first.push_back(traces_list.traces_id[count_traces]);
                                count_traces++;
                                fractures_list.trace_type[j].second.push_back(0);
                            }
                            else // Caso di traccia non passante per j
                            {
                                fractures_list.trace_type[j].first.push_back(traces_list.traces_id[count_traces]);
                                count_traces++;
                                fractures_list.trace_type[j].second.push_back(1);
                            }
                        }
                    }
                }
                // Fine classificazione della traccia.
            }
            // Fine processo calcolo traccia tra i e j.
        }
        // Fine scorrimento fratture successive ad i.
    }
    // Fine scorrimento elenco fratture.

    /* Test stampa struttura salvavita
    for (unsigned int n=0; n < fractures_list.N_frac; n++)
    {
        cout << "Key: " << n << endl;
        cout << "Id traccia: ";
        for (unsigned int val : fractures_list.trace_type[n].first)
        {
            cout << val << " ";
        }
        cout << endl;

        cout << "Tips:";
        for (int val : fractures_list.trace_type[n].second)
        {
            cout << val << " ";
        }
        cout << endl;
        cout << endl;
    }*/
}

////////////////////////////////////////////////////////////////

/*
Funzione che calcola la forma parametrica di una retta e la trasforma in cartesiana.
Prende in input 2 punti e restituisce due vettori 4x1 che identificano la retta mediante il seguente calcolo
(x-x0)/a  = (y-y0)/b e (x-x0)/a  = (z-z0)/c dove (a,b,c) coincide con la direzione della retta.
I piani vanno definiti fuori dalla funzione in modo da poter essere restituiti in modo corretto.
*/
pair<Vector4d, Vector4d> equazioneRetta(const Vector3d& v1, const Vector3d& v2)
{
    if ((v1 - v2).lpNorm<1>() < 1e-9) // bisogna usare una tolleranza altrimenti non funzionano i test
    {
        throw invalid_argument("v1 and v2 must be different vectors");
    }
    Vector4d pi1, pi2, pi3;

    Vector3d n = v1-v2; // Direzione retta, retta: v1+t*n (P0+t*n, n reale)

    // Converto in forma cartesiana.
    pi1[0] = n[1];
    pi1[1] = -n[0];
    pi1[3] = n[0]*v1[1] - n[1]*v1[0];
    pi1[2] = 0.0;

    pi2[0] = n[2];
    pi2[2] = -n[0];
    pi2[3] = n[0]*v1[2] - n[2]*v1[0];
    pi2[1] = 0.0;

    pi3[0] = n[0];
    pi3[2] = -n[1];
    pi3[3] = n[1]*v1[2] - n[2]*v1[1];
    pi3[1] = n[2];

    if (pi1.lpNorm<1>() < 1e-9) {
        return make_pair(pi3, pi2);
    }
    else if (pi2.lpNorm<1>() < 1e-9) {
        return make_pair(pi1, pi3);
    }
    else {
        return make_pair(pi1, pi2);
    }
}

/////////////////////////////////

// Funzione che calcola il piano passante per un poligono.
Vector4d pianoFrattura(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3)
{
    double eps = numeric_limits<decltype(eps)>::epsilon();
    Vector3d AB = v2-v1;
    Vector3d AC = v3-v1;
    Vector3d n1 = AB.cross(AC); // Vettore normale al piano.
    // Verifica collinearità usando il prodotto vettoriale -> se l'area del triangolo è nulla i punti sono allineati
    if (n1.lpNorm<1>() < eps)
    {
        throw invalid_argument("The points do not define a valid plane (points are collinear or coincident).");
    }
    double d = - n1.dot(v1);
    Vector4d piano;
    piano << n1, d;
    return piano;
}

/////////////////////////////////

// Funzione che stampa le informazioni della traccia sul file di output.
bool Export_traces_Info(Traces& t)
{
    ofstream of("traces_info.txt", ios::out);
    if(!of.is_open())
    {
        cerr << "Errore nell'apertura del file di Output per le tracce." << endl;
        return false;
    }

    of.precision(16);
    of << scientific;
    // Verifica di coerenza delle dimensioni dei vettori
    if (t.traces_id.size() != t.traces_gen.size() || t.traces_id.size() != t.traces_points.size())
    {
        cerr << "Errore: Le dimensioni dei vettori di tracce non sono coerenti." << endl;
        return false;
    }
    of << "# Number of Traces" << "\n";
    of << t.traces_id.size() << "\n";
    //of << "# TraceId; FractureId1; Fracture Id2; X1; Y1; Z1; X2; Y2; Z2 \n";
    for(unsigned int i = 0; i < t.traces_id.size(); i++)
    {
        of << "# TraceId; FractureId1; Fracture Id2; X1; Y1; Z1; X2; Y2; Z2 \n";
        of << t.traces_id[i] << ";" << t.traces_gen[i][0] << ";" << t.traces_gen[i][1] << ";";

        if (t.traces_points[i][0].size() != 3 || t.traces_points[i][1].size() != 3)
        {
            cerr << "Errore: Le coordinate della traccia non sono nel formato corretto." << endl;
            return false;
        }
        of << t.traces_points[i][0][0] << ";" << t.traces_points[i][0][1] << ";" <<t.traces_points[i][0][2] << ";"
        << t.traces_points[i][1][0] << ";" << t.traces_points[i][1][1] << ";" << t.traces_points[i][1][2] << "\n";
    }
    of.close();
    return true;
}

// Funzione che verifica se una traccia è passante per una frattura.
bool check_pass(const Vector4d& pi1, const Vector4d& pi2, const Vector3d& point,const double& tau)
{
    // Calcola la distanza del punto dal primo piano
    double dist1 = pi1.head<3>().dot(point) + pi1[3];
    // Calcola la distanza del punto dal secondo piano
    double dist2 = pi2.head<3>().dot(point) + pi2[3];

    // Verifica se il punto si trova abbastanza vicino a uno dei piani entro una tolleranza eps
    return (abs(dist1) < tau) && (abs(dist2) < tau);
}

/////////////////////////////////

// Funzione che stampa le informazioni della traccia sul file di output
bool Export_traces_Type(Fractures& f,Traces& t)
{
    ofstream of("traces_type.txt", ios::out);
    if(!of.is_open())
    {
        cerr << "Errore nell'apertura del file di output per il tipo di tracce." << endl;
        return false;
    }
    of.precision(16);
    of << scientific;
    for(unsigned int j = 0; j < f.N_frac; j++)
    {
        if(f.trace_type[j].first.size() != 0)
        {
            of << "# FractureId; NumTraces \n";
            of << j << ";" << f.trace_type[j].first.size() << "\n";
            if(f.trace_type[j].first.size() == 0)
            {
                continue;
            }
            of << "# TraceId; Tips; Length \n";
            for (unsigned int i = 0; i<f.trace_type[j].first.size(); i++)
            {
                of << f.trace_type[j].first[i] << ";" << f.trace_type[j].second[i] << ";" << t.traces_length[f.trace_type[j].first[i]] <<"\n";
            }
        }
    }
    of.close();
    return true;
}

/////////////////////////////////

// Funzione che calcola il Bounding Box, date le coordinate di un poligono.
vector<Vector3d> Calculate_Bounding_Box(vector<Vector3d>& polygon)
{
    if (polygon.empty()) {
        throw invalid_argument("Il vettore di input è vuoto.");
    }
    // Inizializzo le coordinate del punto massimo e minimo, la prima colonna contiene il punto minimo.
    vector<Vector3d> Bbox = {polygon[0],polygon[0]};

    // Itero per conoscere le coordinate del bounding box.
    for (const auto& coord : polygon )
    {
        (Bbox[0])[0] = min((Bbox[0])[0], coord[0]);
        (Bbox[1])[0] = max((Bbox[1])[0], coord[0]);

        (Bbox[0])[1] = min((Bbox[0])[1], coord[1]);
        (Bbox[1])[1] = max((Bbox[1])[1], coord[1]);

        (Bbox[0])[2] = min((Bbox[0])[2], coord[2]);
        (Bbox[1])[2] = max((Bbox[1])[2], coord[2]);
    }

    return Bbox;
}

/////////////////////////////////

// Funzione di ordinamento della struttura salvavita.
void Sort_Traces_Type(Fractures& f, Traces &t)
{
    for(auto& fracture_pair : f.trace_type)
    {
        // Evito il calcolo se il vettore ha un solo elemento
        if(fracture_pair.first.size() == 1 || fracture_pair.first.size() == 0 )
        {
            continue;
        }
        // Evito il calcolo se il vettore ha una sola traccia passante e una non passante
        if(fracture_pair.second.size() == 2 && fracture_pair.second[0] == 0 &&
            fracture_pair.second[1] == 1)
        {
            continue;
        }
        // Vettori di prova
        // Ordino in base a passante e non passante.
        if (!is_sorted(fracture_pair.second.begin(), fracture_pair.second.end()))
        {
            sort_pair(fracture_pair.second,fracture_pair.first);
        }

        // Ricavo l'indice i di ultima traccia passante con ricerca binaria
        int left = 0;
        int right = fracture_pair.second.size() - 1;
        int lastZeroIndex = -1;
        while (left <= right)
        {
            int mid = left + (right - left) / 2;
            if (fracture_pair.second[mid] == 0)
            {
                lastZeroIndex = mid;
                left = mid + 1;
            } else
            {
                right = mid - 1;
            }
        }
        if(lastZeroIndex == -1 || lastZeroIndex == fracture_pair.first.size()) // Tutte 0 o tutte 1
        {
            // Creo il vettore delle lunghezze per poter ordinare gli ID.
            vector<double> lengths;
            lengths.reserve(fracture_pair.first.size());
            for(auto& el: fracture_pair.first)
            {
                lengths.push_back(t.traces_length[el]);
            }
            //vector<double> vettore temporaneo delle lunghezze
            if (!is_sorted(lengths.begin(), lengths.end()))
            {
                sort_pair(lengths,fracture_pair.first);
            }
        }
        else
        {
            vector<double> lengths;
            lengths.reserve(fracture_pair.first.size());
            for(auto& el: fracture_pair.first)
            {
                lengths.push_back(t.traces_length[el]);
            }
            if (!is_sorted(lengths.begin(), lengths.end()))
            {

                // vector<double> vettore temporaneo con solo metà lunghezze ordinate per tips
                // vector<double> lunghezze parziali dal lastzerIndex fino alla fine
                vector<double> first_half(lengths.begin(), lengths.begin() + lastZeroIndex + 1);
                vector<unsigned int> first_half_ids(fracture_pair.first.begin(), fracture_pair.first.begin() + lastZeroIndex + 1);

                vector<double> second_half(lengths.begin() + lastZeroIndex + 1, lengths.end());
                vector<unsigned int> second_half_ids(fracture_pair.first.begin() + lastZeroIndex + 1, fracture_pair.first.end());

                sort_pair(first_half , first_half_ids);
                sort_pair(second_half, second_half_ids);

                copy(first_half.begin(), first_half.end(), lengths.begin());
                copy(second_half.begin(), second_half.end(), lengths.begin() + lastZeroIndex + 1);

                copy(first_half_ids.begin(), first_half_ids.end(), fracture_pair.first.begin());
                copy(second_half_ids.begin(), second_half_ids.end(), fracture_pair.first.begin() + lastZeroIndex + 1);
            }
        }
    }
}

template<typename T>
// Funzione che ordina vec1 fino all'index e applica gli scambi a vec2
void sort_pair(vector<T>& vec1, vector<unsigned int>& vec2)
{
    vector<unsigned int> indices(vec1.size());
    iota(indices.begin(), indices.end(), 0);

    // Sort indices based on values in vec1
    sort(indices.begin(), indices.end(),
              [&vec1](unsigned int i, unsigned int j) { return vec1[i] < vec1[j]; });

    // Reorder vec1 and vec2 based on sorted indices
    vector<T> sorted_vec1(vec1.size());
    vector<unsigned int> sorted_vec2(vec2.size());

    for (size_t i = 0; i < vec1.size(); ++i) {
        sorted_vec1[i] = vec1[indices[i]];
        sorted_vec2[i] = vec2[indices[i]];
    }

    // Update vec1 and vec2 with sorted values
    vec1 = move(sorted_vec1);
    vec2 = move(sorted_vec2);
}




// Paraview

void Export_Paraview(Fractures &f)
{
    int N = 0;
    for (const auto& vec : f.frac_vertices) {
        N += vec.size();
    }

    MatrixXd points(3,N); // OK
    int col = 0;
    for (const auto& vec : f.frac_vertices) {
        for (const auto& point : vec) {
            points(0, col) = point.x();
            points(1, col) = point.y();
            points(2, col) = point.z();
            ++col;
        }
    }

    vector<vector<unsigned int >> polygon_vertices; // OK
    int start_index = 0;
    for (unsigned int n : f.N_vert) {
        vector<unsigned int> polygon_ids;
        for (unsigned int i = 0; i < n; ++i) {
            polygon_ids.push_back(start_index + i);
        }
        polygon_vertices.push_back(polygon_ids);
        start_index += n;
    }

    vector<UCDProperty<double>> points_properties;
    vector<UCDProperty<double>> polygons_properties;
    VectorXi material;
    Gedim::UCDUtilities UCD;
    UCD.ExportPoints("points_paraview.inp",points,points_properties,material);
    UCD.ExportPolygons("polygons_paraview.inp", points, polygon_vertices, points_properties, polygons_properties, material);
}

/*void Export_Paraview(vector<vector<Vector3d>>& found_polygons)
{
    int N = 0;
    for (const auto& vec : found_polygons) {
        N += vec.size();
    }

    MatrixXd points(3,N); // OK
    int col = 0;
    for (const auto& vec : found_polygons) {
        for (const auto& point : vec) {
            points(0, col) = point.x();
            points(1, col) = point.y();
            points(2, col) = point.z();
            ++col;
        }
    }

    vector<vector<unsigned int >> polygon_vertices; // OK
    int start_index = 0;
    for (unsigned int n = 0; n < 28; n++) {
        vector<unsigned int> polygon_ids;
        for (unsigned int i = 0; i < n; ++i) {
            polygon_ids.push_back(start_index + i);
        }
        polygon_vertices.push_back(polygon_ids);
        start_index += n;
    }



    vector<UCDProperty<double>> points_properties;
    vector<UCDProperty<double>> polygons_properties;
    VectorXi material;
    Gedim::UCDUtilities UCD;
    UCD.ExportPolygons("polygons_paraview.inp", points, polygon_vertices, points_properties, polygons_properties, material);
}*/

////////////////////////////////////////////////////////////////////////////////////
// PARTE 2
////////////////////////////////////////////////////////////////////////////////////

// Funzione che prolunga una traccia fino ad incontrare i lati della frattura
// Calcolo la retta passante per la traccia, e per ogni lato della frattura, trovo la loro intersezione
vector<Vector3d> extendTraceToEdges(vector<Vector3d>& frac_vertices,
                                      vector<Vector3d>& traces_points, unsigned int & tip)
{
    pair<Vector4d, Vector4d> points = equazioneRetta(traces_points[0], traces_points[1]);
    Vector3d sol;
    vector<Vector3d> result;
    result.reserve(6);
    result.push_back(traces_points[0]);
    result.push_back(traces_points[1]);
    double eps = 1e-9;

    unsigned int count = 0;

    for(unsigned int i = 0; i < frac_vertices.size(); i++)
    {
        pair<Vector4d, Vector4d> ver;
        Vector3d coor1, coor2;
        vector<Vector3d> bBox = Calculate_Bounding_Box(frac_vertices);
        if(i == 3)
        {
            if ((frac_vertices[0] - frac_vertices[3]).lpNorm<1>() >= eps)
            {
                ver = equazioneRetta(frac_vertices[0], frac_vertices[3]);
                coor1 = frac_vertices[0];
                coor2 = frac_vertices[3];
            }
            else
                break;
        }
        else
        {
            if ((frac_vertices[i] - frac_vertices[i+1]).lpNorm<1>() >= eps)
            {
                ver = equazioneRetta(frac_vertices[i], frac_vertices[i+1]);
                coor1 = frac_vertices[i];
                coor2 = frac_vertices[i+1];
            }
            else
                break;
        }
        Matrix<double,4,3> coeff;
        coeff.row(0) << points.first[0], points.first[1], points.first[2];
        coeff.row(1) << points.second[0], points.second[1], points.second[2];
        coeff.row(2) << ver.first[0],ver.first[1],ver.first[2];
        coeff.row(3) <<  ver.second[0], ver.second[1], ver.second[2];
        JacobiSVD<MatrixXd> svd(coeff);
        double cond =svd.singularValues().minCoeff();
        if (cond > 1e-9) // evita il caso di lati paralleli
        {
            if(tip == 0)
            {
                if(find(result.begin(), result.end(), coor1) == result.end()) // controllo se ho già trovato gli estremi del segmento
                    result.push_back(coor1);
                if(find(result.begin(), result.end(), coor2) == result.end())
                    result.push_back(coor2);
            }
            else
            {
                Vector4d termineNoto;
                termineNoto[0] = -points.first[3];
                termineNoto[1] = -points.second[3];
                termineNoto[2] = -ver.first[3];
                termineNoto[3] = -ver.second[3];
                HouseholderQR<MatrixXd> qr(coeff);
                sol = qr.solve(termineNoto);

                bool overlap_x = (sol[0] >= bBox[0][0] - eps) && (sol[0] <= bBox[1][0] + eps);
                bool overlap_y = (sol[1] >= bBox[0][1] - eps) && (sol[1] <= bBox[1][1] + eps);
                bool overlap_z = (sol[2] >= bBox[0][2] - eps) && (sol[2] <= bBox[1][2] + eps);

                if(overlap_x && overlap_y && overlap_z)
                {
                    result[count] = sol;
                    result.push_back(coor1);
                    result.push_back(coor2);
                    count++ ;
                }
            }
        }
        else
        {
            continue;
        }
    }
    /*cout << "vettore result" << endl;
    for(auto& el : result) {
        cout << el.transpose() << endl;
    }*/
    return result;
}

bool cutPolygons(Fractures& f, Traces& t, Fractures &final_pol)
{
    vector<PolygonalMesh> result;
    vector<vector<Vector3d>> found_polygons;
    double eps = 1e-9;

    for (unsigned int i = 0; i < f.N_frac; i++)
    {
        PolygonalMesh mesh;
        // Inizializza la coda con il poligono iniziale
        queue<vector<Vector3d>> polygon_queue;
        polygon_queue.push(f.frac_vertices[i]);

        for (unsigned int k = 0; k < f.trace_type[i].first.size(); k++)
        {
            int trace_id = f.trace_type[i].first[k];
            unsigned int tip = f.trace_type[i].second[k];
            unsigned int size_queue = polygon_queue.size();

            for (unsigned int j = 0; j < size_queue; j++)
            {
                vector<Vector3d> current_polygon = polygon_queue.front();
                polygon_queue.pop();

                unsigned int flag = 0;

                // Il tip è 1 se considero le tracce successive alla prima
                if(k != 0)
                    tip = 1;

                // Controllo che la traccia non coincida con un lato del poligono
                for(auto& el : current_polygon)
                {
                    unsigned int count = 0;
                    for(unsigned int n = 0; n < 3; n++)
                        if (abs(el[n] - t.traces_points[trace_id][0][n]) <= eps || abs(el[n] - t.traces_points[trace_id][1][n]) <= eps)
                        {
                            count++;
                        }
                    if(count == 3)
                        flag++;
                }
                // Controllo che la traccia sia interna al poligono
                vector<Vector3d> bBox1 = Calculate_Bounding_Box(current_polygon);

                for (auto& el : t.traces_points[trace_id])
                {
                    bool overlap_x = (el[0] >= bBox1[0][0] - eps) && (el[0] <= bBox1[1][0] + eps);
                    bool overlap_y = (el[1] >= bBox1[0][1] - eps) && (el[1] <= bBox1[1][1] + eps);
                    bool overlap_z = (el[2] >= bBox1[0][2] - eps) && (el[2] <= bBox1[1][2] + eps);
                    if(!overlap_x || !overlap_y || !overlap_z)
                    {
                        flag = 1; // flag non considerare poligono
                        break;
                    }
                }
                if (flag >= 1)
                {
                    found_polygons.push_back(current_polygon);
                    continue;
                }
                else
                {
                    pair<vector<Vector3d>, vector<Vector3d>> polygons;

                    if(t.traces_points[trace_id][0] != t.traces_points[trace_id][1])
                    {
                        polygons = subPolygons(current_polygon, t.traces_points[trace_id], tip);
                        if (!polygons.first.empty())
                        {
                            polygon_queue.push(polygons.first);
                        }
                        if (!polygons.second.empty())
                        {
                            polygon_queue.push(polygons.second);
                        }
                    }


                }
            } // chiudo ciclo su queue
        } // chiudo ciclo sulle tracce

        while(!polygon_queue.empty())
        {
            found_polygons.push_back(polygon_queue.front());
            polygon_queue.pop();
        }

        // Riempio la struttura dati richiesta
        /*for(unsigned int i = result.end()->NumberOfCell2Ds; i < found_polygons.size(); i++)
        {
            auto polygons = found_polygons[i];
            vector<unsigned int> vec(polygons.size());
            iota(vec.begin(), vec.end(), mesh.NumberOfCell0Ds);

            // Celle 0D
            mesh.NumberOfCell0Ds += polygons.size();
            mesh.IdCell0Ds.insert(mesh.IdCell0Ds.end(), vec.begin(), vec.end());
            mesh.CoordinatesCell0Ds.insert(mesh.CoordinatesCell0Ds.end(), polygons.begin(), polygons.end());

            // Cell1D
            mesh.NumberOfCell1Ds += polygons.size();
            mesh.IdCell1Ds.insert(mesh.IdCell1Ds.end(), vec.begin(), vec.end());
            mesh.VerticesCell1Ds.push_back(vec);

            mesh.NumberOfVertices.push_back(polygons.size());
            mesh.NumberOfEdges.push_back(polygons.size());

            // Cell2D
            mesh.NumberOfCell2Ds += 1;
            //mesh.IdCell2Ds.insert(mesh.IdCell2Ds.end(), {0, 1});
            mesh.VerticesCell2Ds.push_back(vec);
            mesh.EdgesCell2Ds.push_back(vec);
        }
        result.push_back(mesh);*/


    } // chiudo ciclo sui poligoni

    final_pol.N_frac = found_polygons.size();
    vector<unsigned int> vec(found_polygons.size());
    iota(vec.begin(), vec.end(), found_polygons.size());

    final_pol.frac_id.insert(final_pol.frac_id.end(), vec.begin(), vec.end());
    for (auto& pol : found_polygons)
    {
        final_pol.N_vert.push_back(pol.size());
    }
    final_pol.frac_vertices = found_polygons;


    // Stampa dei poligoni trovati
    for (const auto& polygons : found_polygons)
    {
        cout << "Poligono risultante:" << endl;
        for (const auto& el : polygons) {
            cout << el.transpose() << endl;
        }
    }

    printSubPolygons(result);

    return true;
}



pair<vector<Vector3d>,vector<Vector3d>> subPolygons(vector<Vector3d> frac_vertices,
                                                     vector<Vector3d> traces_points,
                                                     unsigned int tip)
{
    vector<Vector3d> result = extendTraceToEdges(frac_vertices,
                                                 traces_points,
                                                 tip);

    unsigned int index1 = 0;
    unsigned int index2 = 0;
    if (result.size() == 6)
    {
        auto it = find(frac_vertices.begin(), frac_vertices.end(), result[2]);
        index1 = distance(frac_vertices.begin(), it);
        it = find(frac_vertices.begin(), frac_vertices.end(), result[4]);
        index2 = distance(frac_vertices.begin(), it);
    }
    //else
    //{
    //    return false;
    //}
    // al posto di questo lanciare eccezione
    vector<Vector3d> pol1, pol2;

    bool newPol = false;
    for (unsigned int j = 0; j < frac_vertices.size(); j++)
    {
        if(newPol == false)
        {
            pol1.push_back(frac_vertices[j]);
        }
        else
        {
            pol2.push_back(frac_vertices[j]);
        }


        if (j == index1)
        {
            newPol = true;
            pol1.push_back(result[0]);
            pol1.push_back(result[1]);
        }

        if (j == index2)
        {
            newPol = false;
            pol2.push_back(result[1]);
            pol2.push_back(result[0]);
        }
    }

    return make_pair(pol1, pol2);
}

void printSubPolygons(vector<PolygonalMesh>& sub_polygons)
{
    for (size_t i = 0; i < sub_polygons.size(); ++i)
    {
        cout << "PolygonalMesh " << i << ":\n";
        const auto& mesh = sub_polygons[i];

        cout << "NumberOfCell0Ds: " << mesh.NumberOfCell0Ds << "\n";
        cout << "CoordinatesCell0Ds:\n";
        for (auto coord : mesh.CoordinatesCell0Ds) {
            cout << coord.transpose() << "\n";
        }

        cout << "NumberOfCell1Ds: " << mesh.NumberOfCell1Ds << "\n";
        cout << "VerticesCell1Ds:\n";
        for (auto& edge : mesh.VerticesCell1Ds) {
            for (auto& e : edge)
                cout << e << " ";
        }

        cout << "NumberOfCell2Ds: " << mesh.NumberOfCell2Ds << "\n";
        cout << "VerticesCell2Ds:\n";
        for (auto& vertices : mesh.VerticesCell2Ds) {
            for (auto& v : vertices) {
                cout << v << " ";
            }
            cout << "\n";
        }

        cout << "EdgesCell2Ds:\n";
        for (auto& edges : mesh.EdgesCell2Ds) {
            for (auto& e : edges) {
                cout << e << " ";
            }
            cout << "\n";
        }

        cout << "\n";
    }
}
