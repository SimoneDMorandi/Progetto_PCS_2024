#include "Utils.hpp"
#include "Traccia.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Eigen>
#include <string>
#include <vector>


using namespace std;
using vec3 = vector<vector<vector<double>>>;


bool importFractures(const string& path, vec3& v)
{
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

    for(unsigned int i = 0; i < N; i++) // per ogni frattura salvo le coordinate dei vertici
    {
        getline(file, line);
        getline(file, line);

        istringstream converter(line);
        converter >> idFrac >> delimiter >> numVertices;

        v[i].resize(3); // riservo lo spazio in memoria per salvare le coordinate

        getline(file, line);

        for(unsigned int k = 0; k < 3; k++)
        {
            getline(file, line);
            v[i][k].reserve(numVertices);

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

    vector<vector<vector<double>>> v_transposed(N, vector<vector<double>>(numVertices, vector<double>(3)));

    // Effettua la trasposizione
    for (unsigned int i = 0; i < N; i++)
        for (unsigned int k = 0; k < 3; k++)
            for (unsigned int j = 0; j < numVertices; j++)
                v_transposed[i][j][k] = v[i][k][j];

    v = v_transposed;
/*
    // Stampa dati
    for(unsigned int i = 0; i < N; i++)
    {
        cout << "Id frattura: \t"<< i << endl;
        for(unsigned int j = 0; j < numVertices; j++)
        {
            for(unsigned int k = 0; k < 3; k++)
                cout << v[i][j][k] << "\t";
            cout << endl;
        }
        cout << endl;
    }*/

    return true;
}
