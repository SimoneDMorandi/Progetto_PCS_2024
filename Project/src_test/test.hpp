#pragma once

#include "Utils.hpp"
#include "Utils.cpp"
#include <gtest/gtest.h>
#include <math.h>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <stdexcept>


using namespace std;
using namespace Eigen;
using vec2 = vector<Vector3d>;
using vec3 = vector<vector<Vector3d>>;

double eps = numeric_limits<decltype(eps)>::epsilon();

///////////////////////////////////////////////////////////////////////////////////////////////
TEST(EquazioneRettaTest, BasicTest)
{
    Vector3d v1(1.0, 2.0, 3.0);
    Vector3d v2(4.0, 5.0, 6.0);

    pair<Vector4d, Vector4d> result = equazioneRetta(v1, v2);

    // Calcola i valori attesi
    Vector3d n = v1 - v2;
    Vector4d expected_pi1, expected_pi2;

    expected_pi1[0] = n[1];
    expected_pi1[1] = -n[0];
    expected_pi1[2] = 0.0;
    expected_pi1[3] = n[0] * v1[1] - n[1] * v1[0];

    expected_pi2[0] = n[2];
    expected_pi2[1] = 0.0;
    expected_pi2[2] = -n[0];
    expected_pi2[3] = n[0] * v1[2] - n[2] * v1[0];

    for (int i = 0; i < 4; ++i)
    {
        EXPECT_NEAR(result.first[i], expected_pi1[i], eps);
        EXPECT_NEAR(result.second[i], expected_pi2[i], eps);
    }
}

TEST(EquazioneRettaTest, IdenticalPoints)
{
    Vector3d v1(1.0, 2.0, 3.0);
    Vector3d v2(1.0, 2.0, 3.0);

    EXPECT_THROW(equazioneRetta(v1, v2), invalid_argument);
}

TEST(EquazioneRettaTest, PointsTooClose)
{
    Vector3d v1(1.0, 2.0, 3.0);
    Vector3d v2(1.0 + eps, 2.0 + eps, 3.0 + eps);

    EXPECT_THROW(equazioneRetta(v1, v2), invalid_argument);
}
// non posso calcolare la retta quando i vettori in input sono uguali o circa uguali

//////////////////////////////////////////////////////////////////////////////////////////////
TEST(PianoFratturaTest, ValidPoints)
{
    Vector3d v1(0.8, 0.0, -0.1);
    Vector3d v2(0.8, 0.0, 0.3);
    Vector3d v3(0.8, 1.0, 0.3);

    Vector4d expected_piano;
    Vector3d n1 = (v2 - v1).cross(v3 - v1);
    double d = -n1.dot(v1);
    expected_piano << n1, d;

    Vector4d result = pianoFrattura(v1, v2, v3);

    EXPECT_NEAR(result[0], expected_piano[0], eps);
    EXPECT_NEAR(result[1], expected_piano[1], eps);
    EXPECT_NEAR(result[2], expected_piano[2], eps);
    EXPECT_NEAR(result[3], expected_piano[3], eps);
}

TEST(PianoFratturaTest, IdenticalPoints)
{
    Vector3d v1(1.0, 2.0, 3.0);
    Vector3d v2(1.0, 2.0, 3.0);
    Vector3d v3(1.0, 2.0, 3.0);

    EXPECT_THROW(pianoFrattura(v1, v2, v3), invalid_argument);
}

TEST(PianoFratturaTest, CollinearPoints)
{
    Vector3d v1(1.0, 2.0, 3.0);
    Vector3d v2(2.0, 4.0, 6.0); // v2 = 2 * v1
    Vector3d v3(3.0, 6.0, 9.0); // v3 = 3 * v1

    EXPECT_THROW(pianoFrattura(v1, v2, v3), invalid_argument);
}
// non posso calcolare il piano quando i punti sono uguali oppure allineati

////////////////////////////////////////////////////////////////////////////////////////////

TEST(ImportFracturesTest, ValidImport)
{
    Fractures fractures_list;
    string path = "./DATA/FR3_data.txt"; // Path al file di input valido

    bool result = importFractures(path, fractures_list);

    // Verifica che la funzione ritorni true (input valido)
    ASSERT_TRUE(result);

    // Verifica che la lista delle fratture sia stata riempita correttamente
    ASSERT_EQ(fractures_list.N_frac, 3);
    ASSERT_EQ(fractures_list.frac_id[0], 0);
    ASSERT_EQ(fractures_list.N_vert[0], 4);
}

TEST(ImportFracturesTest, InvalidImport)
{
    Fractures fractures_list;
    string path = "invalid_input.txt"; // Path al file di input non valido

    bool result = importFractures(path, fractures_list);

    // Verifica che la funzione ritorni false (input non valido)
    ASSERT_FALSE(result);
}

///////////////////////////////////////////////////////////////////////////////////////////

TEST(Calculate_Bounding_Box, ValidPolygon)
{
    // da file FR3_data
    vector<Vector3d> polygon = {{0.0, 0.0, 0.0},{1.0, 0.0, 0.0},{1.0, 1.0, 0.0},{0.0, 1.0, 0.0}};
    vector<Vector3d> expected_Bbox = {{0.0, 0.0, 0.0},{1.0, 1.0, 0.0}};

    ASSERT_EQ(Calculate_Bounding_Box(polygon), expected_Bbox);

    // da file FR200_data
    polygon = {{0.75, 0.075, 0.11},{1.32, 1.03, 0.11},{0.95, 1.25, 0.72},{0.39, 0.29, 0.72}};
    expected_Bbox = {{0.39, 0.075, 0.11},{1.32, 1.25, 0.72}};

    ASSERT_EQ(Calculate_Bounding_Box(polygon), expected_Bbox);
}

TEST(Calculate_Bounding_Box, InvalidPolygon)
{
    vector<Vector3d> polygon = {};

    EXPECT_THROW(Calculate_Bounding_Box(polygon), invalid_argument);

}

///////////////////////////////////////////////////////////////////////////////////////////

TEST(Export_traces_Info, ValidExport)
{
    Traces t1 = {
        {1, 2},
        {{101, 102}, {103, 104}},
        {{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}}, {{7.0, 8.0, 9.0}, {10.0, 11.0, 12.0}}}
    };
    ASSERT_TRUE(Export_traces_Info(t1));

    // Test vettore vuoto
    Traces t2 = { {}, {}, {} };
    ASSERT_TRUE(Export_traces_Info(t2));
}

TEST(Export_traces_Info, InvalidExport)
{
    // Test dimensioni non corrispondenti
    Traces t3 = {
        {1, 2},
        {{101, 102}},
        {{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}}}
    };
    ASSERT_FALSE(Export_traces_Info(t3));
}
// Warning che ricevo: non sto inizializzando tutti i membri della struttura traces perché in questo
// export non è richiesta la stampa passante e non passante

//////////////////////////////////////////////////////////////////////////////////////////

