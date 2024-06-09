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

//double eps = numeric_limits<decltype(eps)>::epsilon();
//double eps = 1e-10;

// Custom function to compare two Vector3d objects with tolerance
bool Near(const Vector3d& a, const Vector3d& b, double eps) {
    return (a - b).norm() < eps;
}

// Custom function to compare vectors of Vector3d objects
bool VectorsNear(const vector<Vector3d>& vec1, const vector<Vector3d>& vec2, double eps) {
    if (vec1.size() != vec2.size()) return false;
    for (size_t i = 0; i < vec1.size(); ++i) {
        if (!Near(vec1[i], vec2[i], eps)) {
            return false;
        }
    }
    return true;
}

// Custom function to compare vectors of vectors of Vector3d objects
bool NestedVectorsNear(const vector<vector<Vector3d>>& vec1, const vector<vector<Vector3d>>& vec2, double eps) {
    if (vec1.size() != vec2.size()) return false;
    for (size_t i = 0; i < vec1.size(); ++i) {
        if (vec1[i].size() != vec2[i].size()) return false;
        for (size_t j = 0; j < vec1[i].size(); ++j) {
            if ((vec1[i][j] - vec2[i][j]).norm() >= eps) return false;
        }
    }
    return true;
}


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

TEST(Find_Traces, ValidTraces)
{
    // Dovrei testare i seguenti casi:
    // 1-> fratture che non si intersecano
        Fractures fractures_list;
        Traces traces_list0;

        fractures_list.N_frac = 2;
        fractures_list.frac_id = {0, 1};
        fractures_list.N_vert = {3, 3};
        fractures_list.frac_vertices = {
            { Vector3d(0.8, 0.0, -0.1), Vector3d(0.8, 0.0, 0.3), Vector3d(0.8, 1.0, 0.3), Vector3d(0.8, 1.0, -0.1) }, // Frattura 1
            { Vector3d(-0.38, 0.5, -0.34), Vector3d(0.32, 0.5, -0.34), Vector3d(0.32, 0.5, 0.45), Vector3d(-0.24, 0.5, 0.45) }  // Frattura 2
        };
        fractures_list.trace_type.resize(2);

        Find_Traces(fractures_list, traces_list0);

        ASSERT_TRUE(traces_list0.traces_id.empty());
        ASSERT_TRUE(traces_list0.traces_gen.empty());
        ASSERT_TRUE(traces_list0.traces_points.empty());

    // 2-> fratture con traccia passante per entrambe
        Traces traces_list1;
        fractures_list.N_frac = 2;
        fractures_list.frac_id = {0, 1};
        fractures_list.N_vert = { 3, 3 };
        fractures_list.frac_vertices = {
            { Vector3d(0.0, 0.0, 0.0), Vector3d(1.0, 0.0, 0.0), Vector3d(1.0, 1.0, 0.0), Vector3d(0.0, 1.0, 0.0) },
            { Vector3d(0.8, 0.0, -0.1), Vector3d(0.8, 0.0, 0.3), Vector3d(0.8, 1.0, 0.3), Vector3d(0.8, 1.0, -0.1) }
        };
        Find_Traces(fractures_list, traces_list1);

        vector<vector<unsigned int>> expected_traces_gen = {{0, 1}};
        ASSERT_EQ(traces_list1.traces_gen, expected_traces_gen);
        vector<vector<Vector3d>> expected_points = {{{0.8, 0.0, 0.0}, {0.8, 0.8, 0.0}}};
        for (const auto& sottovettore : (traces_list1.traces_points)) {
            for (const auto& punto : sottovettore) {
                cout << "(" << punto[0] << ", " << punto[1] << ", " << punto[2] << ") ";
            }
            cout << endl;
        }
        ASSERT_TRUE(NestedVectorsNear(traces_list1.traces_points, expected_points, eps));

    // 3-> traccia interna
        /*Traces traces_list2;
        fractures_list.frac_vertices = {
            { Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(1, 1, 0), Vector3d(0, 1, 0) }, // Frattura 1
            { Vector3d(0.1, 0.5, -0.34), Vector3d(0.5, 0.5, -0.34), Vector3d(0.5, 0.5, 0.45), Vector3d(0.15, 0.5, 0.45) }  // Frattura 2
        };
        Find_Traces(fractures_list, traces_list2);

        expected_traces_gen = { { 0, 1 } };
        ASSERT_EQ(traces_list2.traces_gen, expected_traces_gen);
        expected_points = {{{0.12, 0.5, 0}, {0.5, 0.5, 0}}};
        ASSERT_FALSE(NestedVectorsNear(traces_list2.traces_points, expected_points, eps));*/

        Traces traces_list2;
        fractures_list.frac_vertices = {
            { Vector3d(0.1, 0.1, 0.0), Vector3d(0.9, 0.1, 0.0), Vector3d(0.9, 0.9, 0.0), Vector3d(0.1, 0.9, 0.0) }, // Frattura 1
            { Vector3d(0.3, 0.3, 0.0), Vector3d(0.7, 0.3, 0.0), Vector3d(0.7, 0.7, 0.0), Vector3d(0.3, 0.7, 0.0) }  // Frattura 2
        };
        Find_Traces(fractures_list, traces_list2);

        expected_traces_gen = { { 0, 1 } };
        ASSERT_EQ(traces_list2.traces_gen, expected_traces_gen);
        expected_points = { { { 0.4, 0.4, 0.0 }, { 0.6, 0.6, 0.0 } } };
        ASSERT_FALSE(NestedVectorsNear(traces_list2.traces_points, expected_points, eps));


    // 4-> traccia a metà di una
        Traces traces_list3;
        fractures_list.frac_vertices = {
            { Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(1, 1, 0), Vector3d(0, 1, 0) }, // Frattura 1
            { Vector3d(-0.38, 0.5, -0.34), Vector3d(0.32, 0.5, -0.34), Vector3d(0.32, 0.5, 0.45), Vector3d(-0.24, 0.5, 0.45) }  // Frattura 2
        };
        Find_Traces(fractures_list, traces_list3);

        expected_traces_gen = { { 0, 1 } };
        ASSERT_EQ(traces_list3.traces_gen, expected_traces_gen);
        expected_points = {{{0.32, 0.5, 0}, {-0.32, 0.5, 0}}};
        ASSERT_FALSE(NestedVectorsNear(traces_list3.traces_points, expected_points, eps));

    // 5-> una traccia passante per un poligono e una non passante -> è il caso interno?
}

/*TEST(Find_Traces, InvalidTraces)
{

}*/

///////////////////////////////////////////////////////////////////////////////////////////

TEST(check_pass, ValidCheck)
{
    Vector4d pi1 = {0.0, 1.0, 0.0, 0.0};
    Vector4d pi2 = {0.0, 0.0, 1.0, 0.0};
    Vector3d point = {0.8, 0.0, 0.0};

    ASSERT_TRUE(check_pass(pi1, pi2, point));

    point = {1, 1, 0};
    ASSERT_FALSE(check_pass(pi1, pi2, point));

    /*Vector4d pi3 = {0, -0.5, 0, 0.25};
    point = {0, 0.5, 0};
    ASSERT_FALSE(check_pass(pi1, pi3, point, eps));*/

}

/*TEST(check_pass, InvalidCheck)
{

}*/

///////////////////////////////////////////////////////////////////////////////////////////

TEST(Calculate_Bounding_Box, ValidPolygon)
{
    // da file FR3_data
    vector<Vector3d> polygon = {{0.0, 0.0, 0.0},{1.0, 0.0, 0.0},{1.0, 1.0, 0.0},{0.0, 1.0, 0.0}};
    vector<Vector3d> expected_Bbox = {{0.0, 0.0, 0.0},{1.0, 1.0, 0.0}};

    ASSERT_EQ(Calculate_Bounding_Box(polygon), expected_Bbox);
     //ASSERT_TRUE(NestedVectorsNear(Calculate_Bounding_Box(polygon),  expected_Bbox, eps));

    // da file FR200_data
    polygon = {{0.75, 0.075, 0.11},{1.32, 1.03, 0.11},{0.95, 1.25, 0.72},{0.39, 0.29, 0.72}};
    expected_Bbox = {{0.39, 0.075, 0.11},{1.32, 1.25, 0.72}};

    ASSERT_EQ(Calculate_Bounding_Box(polygon), expected_Bbox);
    //ASSERT_TRUE(NestedVectorsNear(Calculate_Bounding_Box(polygon),  expected_Bbox, eps));
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
        {{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}}, {{7.0, 8.0, 9.0}, {10.0, 11.0, 12.0}}},
        {5.2, 5.2}
    };
    ASSERT_TRUE(Export_traces_Info(t1));

    // Test vettore vuoto
    Traces t2 = { {}, {}, {}, {} };
    ASSERT_TRUE(Export_traces_Info(t2));
}

TEST(Export_traces_Info, InvalidExport)
{
    // Test dimensioni non corrispondenti
    Traces t3 = {
        {1, 2},
        {{101, 102}},
        {{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}}},
        {5.2}
    };
    ASSERT_FALSE(Export_traces_Info(t3));
}
// Warning che ricevo: non sto inizializzando tutti i membri della struttura traces perché in questo
// export non è richiesta la stampa passante e non passante

//////////////////////////////////////////////////////////////////////////////////////////
/*
TEST(Export_traces_Type, ValidExport)
{
    Fractures fractures_list;
    fractures_list.N_frac = 2;
    fractures_list.frac_id = {0, 1};
    fractures_list.N_vert = {3, 3};
    fractures_list.frac_vertices = {
        { Vector3d(0.8, 0.0, -0.1), Vector3d(0.8, 0.0, 0.3), Vector3d(0.8, 1.0, 0.3), Vector3d(0.8, 1.0, -0.1) }, // Frattura 1
        { Vector3d(0.0, 0.0, 0.0), Vector3d(1.0, 0.0, 0.0), Vector3d(1.0, 1.0, 0.0), Vector3d(0.0, 1.0, .00) }  // Frattura 2
    };
    fractures_list.trace_type = {{{1}, {1}}};

    Traces t1 = {
        {0},
        {{0, 1}},
        {{{0.8, 0.0, 0.0}, {0.8, 1.0, 0.0}}},
        {1.0}
    };
    ASSERT_TRUE(Export_traces_Type(fractures_list, t1));
}

TEST(Export_traces_Type, InvalidExport)
{
    Fractures f2 = { {}, {}, {}, {}, {} };
    Traces t2 = { {}, {}, {}, {} };
    ASSERT_TRUE(Export_traces_Type(f2, t2));
}
*/
//////////////////////////////////////////////////////////////////////////////////////////

TEST(Sort_Traces_Type, ValidSort)
{
    Fractures fractures;
    fractures.trace_type.push_back({{1, 2, 3}, {0, 1, 0}});
    fractures.trace_type.push_back({{0, 5, 6}, {1, 1, 1}});
    fractures.trace_type.push_back({{7, 8, 9}, {0, 0, 0}});

    Traces traces;
    traces.traces_length = {
        5.0, 2.0, 4.0, 3.0, 1.0, 7.5, 3.6, 9.0, 8.0, 10.0
    };

    Sort_Traces_Type(fractures, traces);

    vector<unsigned int> expected_ids_1 = {1,3,2};
    vector<unsigned int> expected_tips_1 = {0,0,1};

    vector<unsigned int> expected_ids_2 = {6,0,5};
    vector<unsigned int> expected_tips_2 = {1,1,1};

    vector<unsigned int> expected_ids_3 = {8,7,9};
    vector<unsigned int> expected_tips_3 = {0,0,0};

    EXPECT_EQ(fractures.trace_type[0].first, expected_ids_1);
    EXPECT_EQ(fractures.trace_type[0].second, expected_tips_1);

    EXPECT_EQ(fractures.trace_type[1].first, expected_ids_2);
    EXPECT_EQ(fractures.trace_type[1].second, expected_tips_2);

    EXPECT_EQ(fractures.trace_type[2].first, expected_ids_3);
    EXPECT_EQ(fractures.trace_type[2].second, expected_tips_3);
}


//////////////////////////////////////////////////////////////////////////////////////////
