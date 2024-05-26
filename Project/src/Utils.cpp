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
#include <algorithm> // stable_sort

using namespace std;
using namespace Eigen;
using vec2 = vector<Vector3d>;
using vec3 = vector<vector<Vector3d>>;

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
    for(unsigned int i = 0; i < N; i++) // per ogni frattura salvo le coordinate dei vertici
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
    double eps = 1e-9;//numeric_limits<decltype(eps)>::epsilon();

    // Prendo ogni singola frattura e calcolo il Bounding Box.
    for(unsigned int i = 0; i < fractures_list.N_frac - 1 ; i++)
    {
        vec2 Bbox_1 = Calculate_Bounding_Box(fractures_list.frac_vertices[i]);

        // Calcolo il secondo Bounding Box, per ogni frattura successiva a "i".
        for(unsigned int j = i+1; j < fractures_list.N_frac ;j++)
        {
            vec2 Bbox_2 = Calculate_Bounding_Box(fractures_list.frac_vertices[j]);

            // Verifico l'intersezione tra i Bounding Box.
            bool overlap_x = (Bbox_1[1])[0]-(Bbox_2[0])[0]>eps && (Bbox_2[1])[0]-(Bbox_1[0])[0]>eps;
            bool overlap_y = (Bbox_1[1])[1]-(Bbox_2[0])[1]>eps && (Bbox_2[1])[1]-(Bbox_1[0])[1]>eps;
            bool overlap_z = (Bbox_1[1])[2]-(Bbox_2[0])[2]>eps && (Bbox_2[1])[2]-(Bbox_1[0])[2]>eps;

            // Se c'è intersezione, procedo col calcolare i piani delle fratture i e j e trovare la retta di intersezione.
            if(overlap_x && overlap_y && overlap_z)
            {
                // Calcolo i piani delle fratture.
                Vector4d plane_1 = pianoFrattura(fractures_list.frac_vertices[i][0], fractures_list.frac_vertices[i][1],
                                                        fractures_list.frac_vertices[i][2]);
                Vector4d plane_2 = pianoFrattura(fractures_list.frac_vertices[j][0], fractures_list.frac_vertices[j][1],
                                                        fractures_list.frac_vertices[j][2]);

                // Controllo che i piani non siano paralleli. In caso contrario, ignoro.
                if((plane_1.head<3>().cross(plane_2.head<3>())).lpNorm<1>()<eps)
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
                        try{
                            planes = equazioneRetta(fractures_list.frac_vertices[i][0] ,fractures_list.frac_vertices[i][k]);
                        }catch(...)
                        {
                                continue;
                            }
                    }
                    else{
                        try{
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
                    JacobiSVD<MatrixXd> svd(coeff);
                    double cond = svd.singularValues().minCoeff();
                    if (cond  >= eps)
                    {
                        // Se c'è intersezione con la retta della frattura, salvo il punto risolvendo il sistema lineare.
                        Vector4d termineNoto;
                        termineNoto[0] = -plane_1[3];
                        termineNoto[1] = -plane_2[3];
                        termineNoto[2] = -planes.first[3];
                        termineNoto[3] = -planes.second[3];
                        HouseholderQR<MatrixXd> qr(coeff);
                        points[count] = qr.solve(termineNoto);
                        count++;
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
                        try{
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
                    JacobiSVD<MatrixXd> svd(coeff);
                    double cond =svd.singularValues().minCoeff();
                    if (cond>= eps)
                    {
                        Vector4d termineNoto;
                        termineNoto[0] = -plane_1[3];
                        termineNoto[1] = -plane_2[3];
                        termineNoto[2] = -planes.first[3];
                        termineNoto[3] = -planes.second[3];
                        HouseholderQR<MatrixXd> qr(coeff);
                        points[count] = qr.solve(termineNoto);
                        count++ ;
                    }
                    else
                    {
                        continue;
                    }
                }

                // Caso di traccia passante per la prima e per la seconda frattura.
                if(count == 2 && (points[1]-points[0]).lpNorm<1>() > eps)
                {
                    // Completo la struttura TRACES.
                    traces_list.traces_id.push_back(count_traces);
                    count_traces++;
                    traces_list.traces_points.push_back({points[0], points[1]});
                    traces_list.traces_length.push_back((points[1]-points[0]).lpNorm<1>());
                    traces_list.traces_gen.push_back({i,j});
                    // Completo la struttura FRACTURES salvavita per entrambe le fratture.
                    fractures_list.trace_type[i].first.push_back(traces_list.traces_id[i]);
                    fractures_list.trace_type[i].second.push_back(0);
                    fractures_list.trace_type[j].first.push_back(traces_list.traces_id[i]);
                    fractures_list.trace_type[j].second.push_back(0);
                }

                // Caso in cui la traccia è solo un punto, non lo considero.
                else if (count == 1)
                {
                    continue;
                }

                // Caso di traccia non passante per almeno una frattura
                 /* Calcolo del segmento corrispondente alla traccia, avendo i 4 punti [A,C//B,D]
                 Visualmente, se ho [A,C] e [B,C] tutti sulla stessa retta, i primi appartenenti alla traccia
                 1 e gli altri alla traccia 2, se calcolo il minimo tra gli inizi trovo l'inizio della traccia
                  e calcolando il massimo tra le code trovo la fine della traccia, rappresentata da [Start,Finish]*/
                else if (count > 1 && (points[1]-points[0]).lpNorm<1>() > eps)
                {
                    Vector3d start;
                    Vector3d finish;
                    start = points[0].array().min(points[2].array());
                    finish = points[1].array().max(points[3].array());
                    if((start-finish).lpNorm<1>() > eps)
                    {
                        // Completo la struttura TRACES
                        traces_list.traces_id.push_back(count_traces);
                        traces_list.traces_points.push_back({start, finish});
                        traces_list.traces_length.push_back((finish-start).lpNorm<1>());
                        traces_list.traces_gen.push_back({i,j});

                        // Controllo se la traccia è passante per almeno una frattura.
                        // Prima frattura
                        unsigned int count_pass = 0;
                        for (unsigned int k = 0; k < fractures_list.N_vert[i]; k++)
                        {
                            pair<Vector4d, Vector4d> planes;
                            if(k == fractures_list.frac_vertices[i].size() -1)
                            {
                                try{
                                planes = equazioneRetta(fractures_list.frac_vertices[i][0] ,fractures_list.frac_vertices[i][k]);
                                }catch(...)
                                {
                                    continue;
                                }
                            }
                            else{
                                try{
                                planes = equazioneRetta(fractures_list.frac_vertices[i][k] ,fractures_list.frac_vertices[i][k+1]);
                                }catch(...)
                                {
                                    continue;
                                }
                            }

                            if(check_pass(planes.first,planes.second,start,eps))
                            {
                                count_pass ++;
                            }
                            if(check_pass(planes.first,planes.second,finish,eps))
                            {
                                count_pass ++;
                            }
                        }
                        if(count_pass == 2) // Caso di traccia passante per il primo poligono e non per il secondo
                        {
                            fractures_list.trace_type[i].first.push_back(traces_list.traces_id[count_traces]);
                            count_traces++;
                            fractures_list.trace_type[i].second.push_back(0);
                            fractures_list.trace_type[j].first.push_back(traces_list.traces_id[i]);
                            fractures_list.trace_type[j].second.push_back(1);
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
                                    try{
                                    planes = equazioneRetta(fractures_list.frac_vertices[j][0] ,fractures_list.frac_vertices[j][k]);
                                    }catch(...){
                                        continue;
                                    }
                                }
                                else{
                                    try{
                                    planes = equazioneRetta(fractures_list.frac_vertices[j][k] ,fractures_list.frac_vertices[j][k+1]);
                                    }catch(...){
                                        continue;
                                    }
                                }
                                if(check_pass(planes.first,planes.second,start,eps))
                                {
                                    count_pass ++;
                                }
                                if(check_pass(planes.first,planes.second,finish,eps))
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

    // Test stampa struttura salvavita
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
    }
    /* Test stampa lunghezza OK
    cout << "Lunghezze" << endl;
    for ( const auto& length : traces_list.traces_length)
    {
        cout << length << endl;
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
    Vector4d pi1, pi2;

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

    return make_pair(pi1, pi2);
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
    of << scientific << endl;
    // Verifica di coerenza delle dimensioni dei vettori
    if (t.traces_id.size() != t.traces_gen.size() || t.traces_id.size() != t.traces_points.size())
    {
        cerr << "Errore: Le dimensioni dei vettori di tracce non sono coerenti." << endl;
        return false;
    }
    of << "# Number of Traces" << "\n";
    of << t.traces_id.size() << "\n";
    of << "# TraceId; FractureId1; Fracture Id2; X1; Y1; Z1; X2; Y2; Z2 \n";
    for(unsigned int i = 0; i < t.traces_id.size(); i++)
    {
        of << t.traces_id[i] << ";" << t.traces_gen[i][0] << ";" << t.traces_gen[i][1] << ";";

        if (t.traces_points[i][0].size() != 3 || t.traces_points[i][1].size() != 3)
        {
            cerr << "Errore: Le coordinate della traccia non sono nel formato corretto." << endl;
            return false;
        }

        for(auto& coord : t.traces_points)
        {
            of << coord[0][0] << ";" << coord[0][1] << ";" << coord[0][2]
            << coord[1][0] << ";" << coord[1][1] << ";" << coord[1][2] << "\n";
        }
    }
    of.close();
    return true;
}

// Funzione che verifica se una traccia è passante per una frattura.
bool check_pass(const Vector4d& pi1, const Vector4d& pi2, const Vector3d& point,const double& eps)
{
    // Assumo che sia non passante in principio, caso più probabile;
    Matrix<double,2,3> A;
    A << pi1.head<3>().transpose(), pi2.head<3>().transpose();
    Vector2d d;
    d << pi1[3], pi2[3];
    if((A*point - d).lpNorm<1>() < eps)
    {
        return true;
    }
    else
    {
        return false;
    }
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
    of << scientific << endl;
    for(unsigned int j = 0; j < f.N_frac; j++)
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
    cout << endl;
    cout << endl;
    cout << "Inizio stampa ordinata per tips";
    cout << endl;
    cout << endl;
    for(auto& fracture_pair : f.trace_type)
    {
        // Evito il calcolo se il vettore ha un solo elemento
        if(fracture_pair.first.size() == 1 || fracture_pair.first.size() == 0 )
        {
            continue;
        }
        // Caso in cui ho solo una traccia passante e una non passante
        if(fracture_pair.second.size() == 2 && fracture_pair.second[0] == 0 &&
            fracture_pair.second[1] == 1)
        {
            continue;
        }

        // Ordino in base a passante non passante. OK

        // Creo un vettore di indici.
        vector<size_t> indices(fracture_pair.second.size());
        iota(indices.begin(), indices.end(), 0);

        // Ordino rispetto alla seconda componente della pair.
        sort(indices.begin(),indices.end(),
             [&](size_t a, size_t b) {
                 return fracture_pair.second[a] < fracture_pair.second[b];
             });

        // Salvo gli scambi su due vettori.
        vector<unsigned int> sortedId(fracture_pair.first.size());
        vector<unsigned int> sortedTips(fracture_pair.second.size());

        for(size_t i = 0; i<indices.size(); i++)
        {
            sortedId[i] = fracture_pair.first[indices[i]];
            sortedTips[i] = fracture_pair.second[indices[i]];
        }

        // Sostituisco i vettori ordinati alla pair.
        fracture_pair.first = sortedId;
        fracture_pair.second = sortedTips;

        // Stampa di prova.
        cout << "ID: ";
        for (int num : fracture_pair.first) {
            cout << num << " ";
        }
        cout << endl;

        cout << "Tips: ";
        for (int num : fracture_pair.second) {
            cout << num << " ";
        }
        cout << endl;

    /////

        cout << endl;
        cout << endl;
        cout << "Inizio stampa ordinata lunghezza";
        cout << endl;
        cout << endl;
        unsigned int left = 0;
        unsigned int right = fracture_pair.second.size() - 1;
        unsigned int i = -1; // Posizione dell'ultimo 0.

        // Ricerco la posizione dell'ultimo 0 con ricerca binaria -> efficiente se ordinato e lo è.
        while (left <= right) {
            unsigned int mid = left + (right - left) / 2;
            if (fracture_pair.second[mid] == 0) {
                i = mid;
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }

        // Ora ho i.
        vector<double> partial_length;
        for (unsigned int j=0; j<=i; j++)
        {
            double value = t.traces_length[fracture_pair.first[j]];
            partial_length.push_back(value);
        }
        // Genero un vettore di indici
        vector<size_t> permutation(partial_length.size());
        iota(permutation.begin(), permutation.end(), 0);

        stable_sort(permutation.begin(), permutation.end(),
                         [&](size_t a, size_t b) { return partial_length[a] < partial_length[b]; });

        // Applico i cambiamenti fino a i al primo elemento della pair.
        vector<unsigned int> &first_vector = fracture_pair.first;
        vector<unsigned int> sorted_first_i(i+1);
        for (int k = 0; k <= i; ++k) {
            sorted_first_i[k] = first_vector[permutation[k]];
        }
        for (int k = 0; k <= i; ++k) {
            first_vector[k] = sorted_first_i[k];
        }

        partial_length.clear();
        for (unsigned int j = i + 1; j < fracture_pair.first.size(); j++) {
            double value = t.traces_length[fracture_pair.first[j]];
            partial_length.push_back(value);
        }

        // Ordino la seconda parte da i+1 fino alla fine.
        permutation.resize(partial_length.size());
        iota(permutation.begin(), permutation.end(), 0);

        stable_sort(permutation.begin(), permutation.end(),
                    [&](size_t a, size_t b) { return partial_length[a] < partial_length[b]; });

        vector<unsigned int> sorted_remaining(fracture_pair.first.size() - (i + 1));
        for (unsigned int k = 0; k < sorted_remaining.size(); ++k) {
            sorted_remaining[k] = fracture_pair.first[i + 1 + permutation[k]];
        }
        for (unsigned int k = 0; k < sorted_remaining.size(); ++k) {
            fracture_pair.first[i+1+k] = sorted_remaining[k];
        }

        // Stampa di prova
        cout << "Updated first vector (first " << i+1 << " positions): ";
        for (unsigned int num : first_vector) {
            cout << num << " ";
        }
        cout << endl;
    }
}
