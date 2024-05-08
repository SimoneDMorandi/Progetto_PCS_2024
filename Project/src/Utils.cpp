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
    ifstream file;
    file.open(path);

    if(file.fail())
        return false;

    string line;
    getline(file, line);

    getline(file, line);  // Aggiungo il numero delle fratture
    fractures_list.N_frac = stoi(line);

    unsigned int idFrac = 0;
    char delimiter = ';';
    size_t numVertices = 0;

    fractures_list.frac_vertices.resize(fractures_list.N_frac);

    for(unsigned int i = 0; i < fractures_list.N_frac; i++) // per ogni frattura salvo le coordinate dei vertici
    {
        getline(file, line);
        getline(file, line);

        istringstream converter(line);
        converter >> idFrac >> delimiter >> numVertices;
        fractures_list.frac_id.push_back(idFrac);
        fractures_list.N_vert.push_back(numVertices);

        getline(file, line);

        for(unsigned int k = 0; k < 3; k++)
        {
            getline(file, line);
            fractures_list.frac_vertices[i].reserve(numVertices);

            istringstream converter(line);
            for(unsigned int j = 0; j < numVertices; j++)
            {
                double val;
                converter >> val >> delimiter;
                fractures_list.frac_vertices[i][k].push_back(val);
            }
        }
    }
    file.close();

    // INIZIO STAMPA, NON RICHIESTA MA UTILE PER TEST
    vec3 v_transposed(fractures_list.N_frac, vec2(fractures_list.N_vert.size(), vector<double>(3)));

    // Effettua la trasposizione
    for (unsigned int i = 0; i < fractures_list.N_frac; i++)
        for (unsigned int k = 0; k < 3; k++)
            for (unsigned int j = 0; j < fractures_list.N_vert[i]; j++)
                v_transposed[i][j][k] = fractures_list.frac_vertices[i][k][j];

    // QUI COPIO, SI PUO' OTTIMIZZARE
    //vec3 fractures_lists = v_transposed;

    // Stampa dati
    for(unsigned int i = 0; i < fractures_list.N_frac; i++)
    {
        cout << "Id frattura: \t"<< i << endl;
        for(unsigned int j = 0; j < fractures_list.N_vert[i]; j++)
        {
            for(unsigned int k = 0; k < 3; k++)
                cout << v_transposed[i][j][k] << "\t";
            cout << endl;
        }
        cout << endl;
    }\
        return true;
}

///////////////////////////////////////////////////

void Find_Traces(Fractures &fractures_list, Traces& traces_list)
{
    // n è la dimensione letta nel file di avvio, al più tutto si interseca e ho n/2 tracce
    traces_list.traces_id.reserve(fractures_list.N_frac);
    double epsilon = numeric_limits<decltype(epsilon)>::epsilon();
    //vec<pair<int, double>> ns_pass;
    //vec<pair<int, double>> ns_notpass;
    // Prendo ogni singola frattura e calcolo il BoundBox
    // Boundbox è un vettore di due vettori di coordinate
    for(unsigned int i = 0; i < fractures_list.N_frac - 1; i++)
    {
        vec2 Bbox_1 = Calculate_Bounding_Box(fractures_list.frac_vertices[i]);

        // Calcolo il secondo Bounding box, per ogni frattura successiva ad i
        for(unsigned int j = i+1; ;j++)
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
                points.reserve(4);

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

                    Matrix<double, 4, 3> coeff;
                    coeff.row(0) << plane_1[0], plane_1[1], plane_1[2];
                    coeff.row(1) << plane_2[0], plane_2[1], plane_2[2];
                    coeff.row(2) << pi1[0], pi1[1], pi1[2];
                    coeff.row(3) << pi2[0], pi2[1], pi2[2];

                    Vector4d termineNoto;
                    termineNoto[0] = -plane_1[3];
                    termineNoto[1] = -plane_2[3];
                    termineNoto[2] = -pi1[3];
                    termineNoto[3] = -pi2[3];
                    HouseholderQR<Matrix<double, 4, 3>> qr(coeff);
                    if (qr.rank() == 3)
                    {
                        points.push_back(qr.Solve(termineNoto));
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
                    Matrix<double, 4, 3> coeff;
                    coeff.row(0) << plane_1[0], plane_1[1], plane_1[2];
                    coeff.row(1) << plane_2[0], plane_2[1], plane_2[2];
                    coeff.row(2) << pi1[0], pi1[1], pi1[2];
                    coeff.row(3) << pi2[0], pi2[1], pi2[2];

                    Vector4d termineNoto;
                    termineNoto[0] = -plane_1[3];
                    termineNoto[1] = -plane_2[3];
                    termineNoto[2] = -pi1[3];
                    termineNoto[3] = -pi2[3];
                    HouseholderQR<Matrix<double, 4, 3>> qr(coeff);
                    if (qr.rank() == 3)
                    {
                        points.push_back(qr.solve(termineNoto));
                    }
                    else
                    {
                        continue;
                    }
                }

                traces_list.traces_id.push_back(i);
                traces_list.traces_points.push_back((points[1]-points[0]) - (points[3]-points[2]));
                traces_list.traces_length.push_back(abs(traces_list.traces_points[i](0)-traces_list.traces_points[i](1)));
                traces_list.traces_gen.push_back({i,j});




                /* Dopo aver ricavato i punti di intersezione.*/
                /* MANCA IL CONTROLLO DELLA TRACCIA DOPPIA
                // CONTROLLO TRACCIA PASSANTE, NON PASSANTE
                unsigned int counter = 0;
                for (auto & coord : polygon(i))
                {
                    for (unsigned int k = 0; k < coord.size(); k++)
                    {// Calcolo retta tra coord[k] e coord[k+1]
                        // Calcolo intersezione tra la retta e la retta del piano
                        if(epsilon < point - Bbox_1[0] && point - Bbox_1[1])
                        {
                            counter++;
                        }
                    }
                    if(counter < 2)
                    {
                        ns_pass.push_back({traces_list.traces_id(i), length});
                    }
                    else
                    {
                        ns_notpass.push_back({traces_list.traces_id(i),length});
                    }
                }
                for(auto & coord : polygon(j))
                {
                    // Stessa  di riga 150.
                    if(counter < 2)
                    {
                        ns_pass.push_back({traces_list.traces_id(i), length});
                    }
                    else
                    {
                        ns_notpass.push_back({traces_list.traces_id(i),length});
                    }
                }

                // Ordino i vettori {id, length}
                traces_list.pass = sort(ns_pass.begin(), ns_pass.end(), [](const pair<int, double>& a, const pair<int, double>& b) {
                    return a.second > b.second;});
                traces_list.not_pass = sort(ns_notpass.begin(), ns_notpass.end(), [](const pair<int, double>& a, const pair<int, double>& b) {
                    return a.second > b.second;});*/
            }
        }
    }
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
    return ris;
}

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
    ostream of;
    of.open("traces_info.csv");
    if(!output_file.open)
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
    ostream of;
    of.open("traces_type.csv");
    if(!output_file.open)
    {
        cerr << "Errore nell'apertura del file di Output per le tracce." << endl;
        return false;
    }
    for(unsigned int i = 0; i < t.traces_id.size(); i++)
    {
        of << "# FractureId1; NumTraces \n";
        of << f.frac_id[i] << ";" << f.N_traces[i] % 2 << "\n";
        cout << "# TraceId; Tips; Length" << "\n";
    }
    of.close();
    return true;
}
////////////////////////////////

// Funzione che calcola il Bounding Box, date le coordinate di un poligono.
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
