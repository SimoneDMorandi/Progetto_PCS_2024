#include "Utils.hpp"
#include "Traccia.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;
using vec2 = vector<vector<double>>;
using vec3 = vector<vector<vector<double>>>;

bool importFractures(const string& path, Fractures &fractures_list)
{
    vec3 v;

    ifstream file;
    file.open(path);

    if(file.fail())
        return false;

    string line;
    getline(file, line);

    getline(file, line);
    size_t N = 0; // numero di fratture

    N = stoi(line);

    unsigned int idFrac = 0;
    char delimiter = ';';
    size_t numVertices = 0;

    v.resize(N);
    fractures_list.N_frac = N;
    fractures_list.frac_vertices.resize(fractures_list.N_frac);

    for(unsigned int i = 0; i < N; i++) // per ogni frattura salvo le coordinate dei vertici
    {
        getline(file, line);
        getline(file, line);

        istringstream converter(line);
        converter >> idFrac >> delimiter >> numVertices;
        fractures_list.frac_id.push_back(idFrac);
        fractures_list.N_vert.push_back(numVertices);

        v[i].resize(3); // riservo lo spazio in memoria per salvare le coordinate

        getline(file, line);

        for(unsigned int k = 0; k < 3; k++)
        {
            getline(file, line);
            //v[i][k].resize(numVertices);

            istringstream converter(line);
            for(unsigned int j = 0; j < numVertices; j++)
            {
                double val;
                converter >> val >> delimiter;
                v[i][k].push_back(val);
            }
        }
    }
    file.close();

    // Definisci una nuova matrice per i vertici trasposti

    // Effettua la trasposizione
    vec3 v_transposed(N, vector<vector<double>>(numVertices, vector<double>(3)));

    // Effettua la trasposizione
    for (unsigned int i = 0; i < N; i++)
        for (unsigned int k = 0; k < 3; k++)
            for (unsigned int j = 0; j < numVertices; j++)
                v_transposed[i][j][k] = v[i][k][j];

    fractures_list.frac_vertices = v_transposed;

    // QUI COPIO, SI PUO' OTTIMIZZARE
    //vec3 fractures_lists = v_transposed;

    // Stampa dati
    for(int i = 0; i < fractures_list.N_frac; i++)
    {
        cout << "Id frattura: \t"<< i << endl;
        for( int j = 0; j < fractures_list.N_vert[i]; j++)
        {
            for(unsigned int k = 0; k < 3; k++)
                cout << fractures_list.frac_vertices[i][j][k] << "\t";
            cout << endl;
        }
        cout << endl;
    }
    return true;
}

///////////////////////////////////////////////////

void Find_Traces(Fractures &fractures_list, Traces& traces_list)
{
    // n è la dimensione letta nel file di avvio, al più tutto si interseca e ho n/2 tracce
    traces_list.traces_id.reserve(fractures_list.N_frac);
    fractures_list.trace_type.reserve(fractures_list.N_frac);
    double epsilon = numeric_limits<decltype(epsilon)>::epsilon();
    // Prendo ogni singola frattura e calcolo il BoundBox
    // Boundbox è un vettore di due vettori di coordinate
    for(int i = 0; i < fractures_list.N_frac - 1; i++)
    {
        vec2 Bbox_1 = Calculate_Bounding_Box(fractures_list.frac_vertices[i]);

        // Calcolo il secondo Bounding box, per ogni frattura successiva ad i
        for(int j = i+1; j < fractures_list.N_frac  ;j++)
        {
            vec2 Bbox_2 = Calculate_Bounding_Box(fractures_list.frac_vertices[j]);
            // Verifico l'intersezione
            bool overlap_x = (Bbox_1[1])[0] >= (Bbox_2[0])[0] && (Bbox_2[1])[0] >= (Bbox_1[0])[0];
            bool overlap_y= (Bbox_1[1])[1] >= (Bbox_2[0])[1] && (Bbox_2[1])[1] >= (Bbox_1[0])[1];
            bool overlap_z = (Bbox_1[1])[2] >= (Bbox_2[0])[2] && (Bbox_2[1])[2] >= (Bbox_1[0])[2];

            // Se c'è intersezione, procedo col calcolare i piani dei poligoni i e j e trovare i punti di intersezione
            if (overlap_x && overlap_y && overlap_z)
            {
                // Calcolo i piani delle fratture
                vector<double> plane_1 = pianoFrattura(fractures_list.frac_vertices[i][0], fractures_list.frac_vertices[i][1],
                                                        fractures_list.frac_vertices[i][2]);

                vector<double> plane_2 = pianoFrattura(fractures_list.frac_vertices[j][0], fractures_list.frac_vertices[j][1],
                                                        fractures_list.frac_vertices[j][2]);

                vector<Vector3d> points;
                unsigned int count = 0;
                points.resize(4);


                // Ricavo segmento del poligono 1
                for (unsigned int k = 0; k < fractures_list.frac_vertices[i].size(); k++)
                {
                    vector<double> pi1(4);
                    vector<double> pi2(4);

                    // Equazione del singolo lato
                    if(k == fractures_list.frac_vertices[i].size() -1)
                    {
                        equazioneRetta(fractures_list.frac_vertices[i][0] ,fractures_list.frac_vertices[i][k], pi1, pi2);
                    }
                    else{
                        equazioneRetta(fractures_list.frac_vertices[i][k] ,fractures_list.frac_vertices[i][k+1], pi1,pi2);
                    }

                    Matrix<double,4,3> coeff;
                    coeff.row(0) << plane_1[0], plane_1[1], plane_1[2];
                    coeff.row(1) << plane_2[0], plane_2[1], plane_2[2];
                    coeff.row(2) << pi1[0], pi1[1], pi1[2];
                    coeff.row(3) << pi2[0], pi2[1], pi2[2];
                    FullPivLU<MatrixXd> lu_decomp(coeff);
                    if (lu_decomp.rank() == 3)
                    {
                        Vector4d termineNoto;
                        termineNoto[0] = -plane_1[3];
                        termineNoto[1] = -plane_2[3];
                        termineNoto[2] = -pi1[3];
                        termineNoto[3] = -pi2[3];
                        HouseholderQR<MatrixXd> qr(coeff);
                        // ordino [A,C]
                        points[count] = qr.solve(termineNoto);
                        count++;
                    }
                    else
                    {
                        continue;
                    }
                }
                // Ricavo segmento del poligono 2
                for (unsigned int k = 0; k < fractures_list.frac_vertices[j].size(); k++)
                {
                    vector<double> pi1(4);
                    vector<double> pi2(4);

                    // Equazione del singolo lato
                    if(k == fractures_list.frac_vertices[i].size() -1)
                    {
                        equazioneRetta(fractures_list.frac_vertices[j][0] ,fractures_list.frac_vertices[j][k], pi1, pi2);
                    }
                    else{

                        equazioneRetta(fractures_list.frac_vertices[j][k] ,fractures_list.frac_vertices[j][k+1], pi1,pi2);
                    }

                    // Intersezione tra i 3 piani
                    Matrix<double,4,3> coeff;
                    coeff.row(0) << plane_1[0], plane_1[1], plane_1[2];
                    coeff.row(1) << plane_2[0], plane_2[1], plane_2[2];
                    coeff.row(2) << pi1[0], pi1[1], pi1[2];
                    coeff.row(3) << pi2[0], pi2[1], pi2[2];
                    FullPivLU<MatrixXd> lu_decomp(coeff);
                    if (lu_decomp.rank() == 3)
                    {
                        Vector4d termineNoto;
                        termineNoto[0] = -plane_1[3];
                        termineNoto[1] = -plane_2[3];
                        termineNoto[2] = -pi1[3];
                        termineNoto[3] = -pi2[3];
                        HouseholderQR<MatrixXd> qr(coeff);
                        points[count] = qr.solve(termineNoto);
                        count++ ;
                    }
                    else
                    {
                        continue;
                    }
                }

                // Caso di traccia passante per un poligono e non passante per l'altra
                if(count == 2 && (points[1]-points[0]).lpNorm<1>() > epsilon)
                {
                    // La traccia è passante per il primo poligono
                    traces_list.traces_id.push_back(i);
                    traces_list.traces_points.push_back({points[0], points[1]});
                    traces_list.traces_length.push_back((points[1]-points[0]).lpNorm<1>());
                    traces_list.traces_gen.push_back({i,j});

                    cout << points[0] << endl;
                    cout << endl;
                    cout << points[1] << endl;
                    cout << endl;
                    // Da riempire la struttura salcavita con Tips = 1 per il poligono i e per j.
                }
                else if (count == 1)
                {
                    continue;
                }
                // Caso di traccia non passante per nessuno dei due poligoni
                 // Calcolo del segmento corrispondente alla traccia, avendo i 4 punti [A,C//B,D]

                 /* Visualmente, se ho [A,C] e [B,C] tutti sulla stessa retta, i primi appartenenti al poligono
                 1 e gli altri al poligono 2, se calcolo il massimo tra gli inizi trovo l'inizio della traccia
                  e calcolando il minimo tra le code trovo la fine della traccia, rappresentata da [Start,Finish]*/
                else
                {
                    Vector3d start;
                    Vector3d finish;
                    // Determino l'intersezione dei due segmenti trovati.
                    start = points[0].array().min(points[2].array());
                    finish = points[1].array().max(points[3].array());
                    if((start-finish).lpNorm<1>() > epsilon)
                    {
                        // Completo la struttura TRACES e la parte di struttura 'salvavita' già con tips = 1??
                         traces_list.traces_points.push_back({start, finish});
                         traces_list.traces_length.push_back((finish-start).lpNorm<1>());
                         traces_list.traces_gen.push_back({i,j});
                         cout << start << endl;
                         cout << endl;
                         cout << finish << endl;
                         cout << endl;
                    }
                }

                // Aggiungo l'id della frattura all'elenco di fratture per ogni poligono
                // fractures_list.trace_type[i] = make_pair(traces_list.traces_id[i], 1); // some_double_value is just a placeholder, replace it with the appropriate value
                // fractures_list.trace_type[j] = make_pair(traces_list.traces_id[i], 1);

                // // CONTROLLO TRACCIA PASSANTE, NON PASSANTE

                // for(int i = 0; i < fractures_list.trace_type.size(); i++) // Corrisponde alla posizione del poligono in Fractures
                // {
                //     unsigned int counter = 0; // Assumo in partenza che sia non passante, caso più probabile.
                //     for(auto& trace_id : fractures_list.trace_type[i].first)
                //     {
                //         // Calcolo l'equazione della retta passante per la traccia
                //         vector<double> pi1_trace(4);
                //         vector<double> pi2_trace(4);
                //         equazioneRetta(traces_list.traces_points[trace_id][0],
                //                        traces_list.traces_points[trace_id][1], pi1_trace, pi2_trace);

                //         // Calcolo l'equazioni delle rette del poligono e le confronto con la traccia
                //         for(j = 0; j < fractures_list.frac_vertices[i].size(); j++)
                //         {
                //             if(j < fractures_list.frac_vertices[i].size()-1)
                //             {
                //                 vector<double> pi1_edge(4);
                //                 vector<double> pi2_edge(4);
                //                 equazioneRetta(fracture_list.frac_vertices[i][j],
                //                                traces_list.traces_points[trace_id][j+1], pi1_edge, pi2_edge);

                //                 // Metto a sistema le due rette, non credo vada bene la funzione di prima.

                //                 //if(c'è soluzione)
                //                 //{
                //                 //counter ++
                //                 //}

                //             }
                //             else
                //             {
                //                 vector<double> pi1_edge(4);
                //                 vector<double> pi2_edge(4);
                //                 equazioneRetta(fracture_list.frac_vertices[i][j],
                //                                traces_list.traces_points[trace_id][0], pi1_edge, pi2_edge);

                //                 // Metto a sistema le due rette, non credo vada bene la funzione di prima.

                //                 //if(c'è soluzione)
                //                 //{
                //                     //counter ++
                //                 //}
                //             }

                //         }
                //         if(counter == 2)
                //         {
                //             // Salvo nella struttura salvavita
                //             fractures_list.trace_type[i].second.push_back(0); // Traccia passante
                //         }
                //         else
                //         {
                //             // Salvo nella struttura salvavita
                //             fractures_list.trace_type[i].second.push_back(1); // Traccia non passante
                //         }

                    // }
                // }
                // Struttura salvavita piena, bisogna però ordinarla, meglio farlo in un'altra funzione con mergesort.
            }
            // Fine processo calcolo traccia tra i e j
        }     
    }
    // Fine scorrimento elenco poligoni
}

// Funzione che calcola la forma parametrica di una retta e la trasforma in cartesiana
// prende in input 2 punti e restituisce due vettori 4x1 che identificano la retta mediante il seguente calcolo
// (x-x0)/a  = (y-y0)/b e (x-x0)/a  = (z-z0)/c dove (a,b,c) coincide con la direzione della retta

void equazioneRetta(const vector<double>& v1, const vector<double>& v2,
                    vector<double>& pi1, vector<double>& pi2)
{
    vector<double> n = sottrazione(v1,v2); // direzione retta, retta: v1+t*n (P0+t*n, n reale)

    // converto in forma cartesiana
    pi1[0] = n[1];
    pi1[1] = -n[0];
    pi1[3] = n[0]*v1[1] - n[1]*v1[0];

    pi2[0] = n[2];
    pi2[2] = -n[0];
    pi2[3] = n[0]*v1[2] - n[2]*v1[0];
}

// Funzione che calcola il piano passante per un poligono
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


// Funzione di prodotto vettoriale da eliminare
vector<double> crossProduct(const vector<double>& u, const vector<double>& v) {
    vector<double> w(3);
    w[0] = u[1] * v[2] - u[2] * v[1];
    w[1] = u[2] * v[0] - u[0] * v[2];
    w[2] = u[0] * v[1] - u[1] * v[0];
    return w;
}

// Funzione di prodotto scalare da eliminare
double dotProduct(const vector<double>& v1, const vector<double>& v2)
{
    double ris = 0;
    for (int i = 0; i < 3; ++i) {
        ris += v1[i] * v2[i];
    }
    return ris;
}

// Funzione di sottraziojne vettoriale da eliminare
vector<double> sottrazione(const vector<double>& v1, const vector<double>& v2)
{
    vector<double> ris(3);
    for (int i = 0; i < 3; ++i) {
        ris[i] = v2[i] - v1[i];
    }
    return ris;
}

/////////////////////////////////

// Funzione che stampa le informazioni della traccia sul file di output
bool Export_traces_Info(Traces& t)
{
    ofstream of("traces_info.txt");
    if(!of.is_open())
    {
        cerr << "Errore nell'apertura del file di Output per le tracce." << endl;
        return false;
    }
    of << "# Number of Traces" << "\n";
    of << t.traces_id.size() << "\n";
    for(unsigned int i = 0; i < t.traces_id.size(); i++)
    {
        of << "# TraceId; FractureId1; Fracture Id2; X1; Y1; Z1; X2; Y2; Z2 \n";
        of << t.traces_id[i] << ";" << t.traces_gen[i][0] << ";" << t.traces_gen[i][1] << ";";
        for(auto& coord : t.traces_points)
        {
            of << coord[0](0) << ";" << coord[0](1) << ";" << coord[0](2) << "\n";
        }
    }
    of.close();
    return true;
}

///////////////////

// Funzione che stampa le informazioni della traccia sul file di output
bool Export_traces_Type(Fractures& f,Traces& t)
{
    ofstream of("traces_type.txt");
    if(!of.is_open())
    {
        cerr << "Errore nell'apertura del file di Output per le tracce." << endl;
        return false;
    }
    for(unsigned int i = 0; i < t.traces_id.size(); i++)
    {
        of << "# FractureId1; NumTraces \n";
        of << f.frac_id[i] << ";" << f.trace_type[i].first.size() << "\n";
        cout << "# TraceId; Tips; Length" << "\n";
        for(unsigned int j = 0; j < f.trace_type[i].first.size(); j++)
        {
            cout << f.trace_type[i].first[j] << ";" << f.trace_type[i].second[j]
                 << ";" << t.traces_length[f.trace_type[i].first[j]] <<  "\n";
        }
    }
    of.close();
    return true;
}
////////////////////////////////

// Funzione che calcola il Bounding Box, date le coordinate di un poligono. OK
vector<vector<double>> Calculate_Bounding_Box(vector<vector<double>>& polygon)
{
    // Inizializzo le coordinate del punto massimo e minimo, la prima colonna contiene il punto minimo
    vector<vector<double>> Bbox = {polygon[0],polygon[0]};

    // Itero per conoscere le coordinate del bounding box
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


// Funzione di ordinamento della struttura salvavita.
bool Sort_Traces_Type(Fractures& f, Traces &t)
{
    // Ordinamento secondo tips
    for(auto pair : f.trace_type)
    {
        // Mergesort su pair.second
        // Quando modifico le posizioni, faccio la stessa cosa per pair.first[i] con pair.first[j]

        VectorXd temporary_length;
        unsigned int i = 0, change = 0;
        bool flag = false;
        // Ordinamento secondo length
        for(auto id : pair.first)
        {
            // Riempio un vettore che contiene le lunghezze corrispondenti alle tracce e salvo quando
            // le tracce iniziano ad essere passanti
            temporary_length(i) = t.traces_length[id];
            if(pair.second[i] != 0 && !flag)
            {
                change = i;
                flag = true;
            }
            i++;
        }
        // Applico Mergesort a temporary_length FINO A CHANGE scambiando gli elementi di pair.first
        // non è necessario scambiare gli elementi di pair.second perché sono tutti 0
        // faccio poi la stessa cosa con temporary_lenght DA CHANGE fino alla fine.
    }
    return true;
}
