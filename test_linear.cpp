#include <iostream>
#include <chrono>
#include <random>
#include "Vec.h"
#include "Matrix.h"

using namespace std;
using namespace chrono;
using namespace Linear;

static const float PI = M_PI;

#define TEST_MATRIX_VECTOR
#define TEST_MATRIX_MATRIX
//#define CLASS_FIRST
#define MAT_DIM 50 // The extra dimension of the multiplied matrices
#define MAT_TEST_SIZE 2000
#define NUM_ROTATIONS 100000000
#if 1
#define VEC_TEST_SIZE 20000 // twenty thousand
#define HEIGHT 100
#define WIDTH 73
#else
#define VEC_TEST_SIZE 60
#define HEIGHT 10000 // ten thousand
#define WIDTH 1000
#endif

random_device rd;
auto seed = rd();
mt19937_64 gen;
uniform_real_distribution<> dis(-100.0, 100.0);

#ifdef TEST_MATRIX_MATRIX
static Matrix<double, HEIGHT, WIDTH> sum_mat;
static double sum_mat_arr[HEIGHT][WIDTH];
static double test_matrix_matrix_class() {
    auto t1 = high_resolution_clock::now();
    sum_mat.zero();
    for (int x = 0; x < MAT_TEST_SIZE; ++x) {
        static Matrix<double, HEIGHT, MAT_DIM> lhs;
        static Matrix<double, MAT_DIM, WIDTH> rhs;
        for (unsigned int i = 0; i < lhs.height(); ++i) {
            for (unsigned int j = 0; j < lhs.width(); ++j) {
                lhs[i][j] = dis(gen);
            }
        }
        for (unsigned int i = 0; i < rhs.height(); ++i) {
            for (unsigned int j = 0; j < rhs.width(); ++j) {
                rhs[i][j] = dis(gen);
            }
        }
        static Matrix<double, HEIGHT, WIDTH> result;
        result = lhs * rhs;
        sum_mat += result;
    }
    auto t2 = high_resolution_clock::now();
    return duration_cast<duration<double> >(t2 - t1).count();
}

static double test_matrix_matrix_no_class() {
    auto t1 = high_resolution_clock::now();
    for (unsigned int i = 0; i < HEIGHT; ++i) {
        for (unsigned int j = 0; j < WIDTH; ++j) {
            sum_mat_arr[i][j] = 0.0;
        }
    }
    for (int x = 0; x < MAT_TEST_SIZE; ++x) {
        static double lhs[HEIGHT][MAT_DIM];
        static double rhs[MAT_DIM][WIDTH];
        for (unsigned int i = 0; i < HEIGHT; ++i) {
            for (unsigned int j = 0; j < MAT_DIM; ++j) {
                lhs[i][j] = dis(gen);
            }
        }
        for (unsigned int i = 0; i < MAT_DIM; ++i) {
            for (unsigned int j = 0; j < WIDTH; ++j) {
                rhs[i][j] = dis(gen);
            }
        }
        static double result[HEIGHT][WIDTH];
        for (unsigned int i = 0; i < HEIGHT; ++i) {
            for (unsigned int j = 0; j < WIDTH; ++j) {
                result[i][j] = 0.0;
                for (int k = 0; k < MAT_DIM; ++k) {
                    result[i][j] += lhs[i][k] * rhs[k][j];
                }
            }
        }
        for (unsigned int i = 0; i < HEIGHT; ++i) {
            for (unsigned int j = 0; j < WIDTH; ++j) {
                sum_mat_arr[i][j] += result[i][j];
            }
        }
    }
    auto t2 = high_resolution_clock::now();
    return duration_cast<duration<double> >(t2 - t1).count();
}
#endif

#ifdef TEST_MATRIX_VECTOR
static Vec<double, HEIGHT> sum_vec;
static double sum_arr[HEIGHT];

static double test_matrix_vector_class() {
    auto t1 = high_resolution_clock::now();
    sum_vec.zero();
    for (int x = 0; x < VEC_TEST_SIZE; ++x) {
        static Matrix<double, HEIGHT, WIDTH> mat;
        static Vec<double, WIDTH> vec;
        for (unsigned int j = 0; j < vec.size(); ++j) {
            vec[j] = dis(gen);
        }
        for (unsigned int i = 0; i < mat.height(); ++i) {
            for (unsigned int j = 0; j < mat.width(); ++j) {
                mat[i][j] = dis(gen);
            }
        }
        static Vec<double, HEIGHT> result;
        result = mat * vec;
        sum_vec += result;
    }
    auto t2 = high_resolution_clock::now();
    return duration_cast<duration<double> >(t2 - t1).count();
}

static double test_matrix_vector_no_class() {
    auto t1 = high_resolution_clock::now();
    for (unsigned int i = 0; i < HEIGHT; ++i) {
        sum_arr[i] = 0.0;
    }
    for (int x = 0; x < VEC_TEST_SIZE; ++x) {
        static double mat[HEIGHT][WIDTH];
        static double vec[WIDTH];
        for (unsigned int j = 0; j < WIDTH; ++j) {
            vec[j] = dis(gen);
        }
        for (unsigned int i = 0; i < HEIGHT; ++i) {
            for (unsigned int j = 0; j < WIDTH; ++j) {
                mat[i][j] = dis(gen);
            }
        }
        static double result[HEIGHT];
        for (unsigned int i = 0; i < HEIGHT; ++i) {
            result[i] = 0.0;
            for (unsigned int j = 0; j < WIDTH; ++j) {
                result[i] += vec[j] * mat[i][j];
            }
        }
        for (unsigned int i = 0; i < HEIGHT; ++i) {
            sum_arr[i] += result[i];
        }
    }
    auto t2 = high_resolution_clock::now();
    return duration_cast<duration<double> >(t2 - t1).count();
}
#endif

int main() {
    Vec2i64 v1 = { 4, 5 };
    v1 = v1 * 8;
    v1 = 8 * v1;
    cout << v1 << endl;
    cout << "THETA: " << v1.theta() << endl;
    cout << "MAG: " << v1.magnitude() << endl;

    cout << endl << "MATRIX STUFF" << endl << endl;

    Matrix<double, 5, 4> mt;
    for (unsigned int i = 0; i < 5; ++i) {
        for (unsigned int j = 0; j < 4; ++j) {
            mt[i][j] = i * j;
        }
    }
    cout << "FIRST" << endl << endl;
    cout << mt << endl << endl;
    ;
    auto m = mt.transpose();
    cout << "TRANSPOSED" << endl << endl;
    cout << m << endl << endl;

    cout << "MULTIPLIED BY" << endl << endl;

    Vec<double, 5> rhs;
    for (unsigned int j = 0; j < rhs.size(); ++j) {
        rhs[j] = j;
    }
    cout << rhs << endl << endl;

    auto lhs = m * rhs;
    cout << lhs << endl << endl;

    cout << "MATRIX MULTIPLCATION" << endl << endl;

    Matrix<double, 3, 2> m1;
    Matrix<double, 2, 4> m2;
    for (unsigned int i = 0; i < m1.height(); ++i) {
        for (unsigned int j = 0; j < m1.width(); ++j) {
            m1[i][j] = i + j;
        }
    }
    for (unsigned int i = 0; i < m2.height(); ++i) {
        for (unsigned int j = 0; j < m2.width(); ++j) {
            m2[i][j] = i + j;
        }
    }

    cout << m1 << endl << endl;
    cout << m2 << endl << endl;
    cout << m1 * m2 << endl << endl;

    cout << "MATRIX IN-PLACE MULTIPLCATION" << endl << endl;

    Matrix<double, 3> m3;
    Matrix<double, 3> m4;
    for (unsigned int i = 0; i < m3.height(); ++i) {
        for (unsigned int j = 0; j < m3.width(); ++j) {
            m3[i][j] = i + j;
        }
    }
    for (unsigned int i = 0; i < m4.height(); ++i) {
        for (unsigned int j = 0; j < m4.width(); ++j) {
            m4[i][j] = i * j;
        }
    }

    cout << m3 << endl << endl;
    cout << m4 << endl << endl;
    m3 *= m4;
    cout << m3 << endl << endl;

    double factor = 4;
    cout << "MULTIPLED BY : " << factor << endl << endl;
    cout << factor * m3 << endl << endl;
    cout << m3 * factor << endl << endl;

    cout << "DIVIDED BY : " << factor << endl << endl;
    cout << m3 / factor << endl << endl;

    cout << "MULTIPLED BY IDENTITY" << endl << endl;
    cout << m3 * Matrix<double, 3>::identity() << endl << endl;

    cout << "MULTIPLED BY ZEROS" << endl << endl;
    cout << m3 * Matrix<double, 3>::zeros() << endl << endl;

    Vec3d test_vec;
    test_vec.x = 3.0;
    test_vec.y = 4.0;
    test_vec.z = 5.0;

    cout << "MULTIPLED BY VEC3D" << endl << endl;
    cout << test_vec << endl << endl;
    cout << m3 * test_vec << endl << endl;

    Matrix<double, 3> m_determine = { { 0.0, 0.1, 0.2 },
                                      { 4.0, -3.2, 0.0 },
                                      { 2.1, 2.1, 9.4 } };
    cout << "THE DETERMINANT OF" << endl << endl;

    cout << m_determine << endl << endl;

    cout << m_determine.determinant() << endl << endl;

    Matrix<double, 3> m_solve = { { 8.0, 2.0, -4.0 },
                                  { 3.0, 8.0, 6.0 },
                                  { 2.0, 9.0, 4.5 } };
    Vec<double, 3> m_rhs = { 4.0, 3.0, 4.0 };
    cout << "THE SOLUTION TO" << endl << endl;

    cout << m_solve << endl << endl;
    cout << m_rhs << endl << endl;

    cout << m_solve.slow_solve(m_rhs) << endl << endl;
    cout << m_solve.slow_solve(Matrix<double, 3, 1>(m_rhs)) << endl << endl;

    // Testing conversions
    {
        Vec<double, 1, false> v_row;
        Vec<double, 1> v_col = to_col(v_row);
        v_col = to_col(v_col);
        v_row = to_row(v_col);
        v_row = to_row(v_row);
        Matrix<double, 1, 1> col_mat = v_col;
        Matrix<double, 1, 1> row_mat = v_row;
        v_row = Vec<double, 1, false>(v_row);
        v_col = Vec<double, 1>(v_col);
        v_row = row_mat;
        v_col = col_mat;
        col_mat = row_mat;
        row_mat = col_mat;
        v_col = row_mat;
        v_row = col_mat;
        row_mat = v_col;
        col_mat = v_row;

        // ERRORS
        // v_col = v_row;
        // v_row = v_col;
        // Matrix<double, 1, 2> x_mat;
        // v_col = x_mat;
        // v_row = x_mat;
        // x_mat = v_col;
        // x_mat = v_row;
    }
    {
        Vec<double, 2, false> v_row;
        Vec<double, 2> v_col = to_col(v_row);
        v_col = to_col(v_col);
        v_row = to_row(v_col);
        v_row = to_row(v_row);
        Matrix<double, 2, 1> col_mat = v_col;
        Matrix<double, 1, 2> row_mat = v_row;
        v_row = Vec<double, 2, false>(v_row);
        v_col = Vec<double, 2>(v_col);
        v_row = row_mat;
        v_col = col_mat;

        // ERRORS
        // v_col = row_mat;
        // v_row = col_mat;
        // col_mat = row_mat;
        // row_mat = col_mat;
        // v_col = v_row;
        // v_row = v_col;
        // row_mat = v_col;
        // col_mat = v_row;
        // Matrix<double, 2, 5> x_mat;
        // Matrix<double, 5, 2> y_mat;
        // v_col = x_mat;
        // v_row = x_mat;
        // x_mat = v_col;
        // x_mat = v_row;
        // v_col = y_mat;
        // v_row = y_mat;
        // y_mat = v_col;
        // y_mat = v_row;
    }
    {
        Vec<double, 3, false> v_row;
        Vec<double, 3> v_col = to_col(v_row);
        v_col = to_col(v_col);
        v_row = to_row(v_col);
        v_row = to_row(v_row);
        Matrix<double, 3, 1> col_mat = v_col;
        Matrix<double, 1, 3> row_mat = v_row;
        v_row = Vec<double, 3, false>(v_row);
        v_col = Vec<double, 3>(v_col);
        v_row = row_mat;
        v_col = col_mat;

        // ERRORS
        // v_col = row_mat;
        // v_row = col_mat;
        // col_mat = row_mat;
        // row_mat = col_mat;
        // v_col = v_row;
        // v_row = v_col;
        // row_mat = v_col;
        // col_mat = v_row;
        // Matrix<double, 3, 2> x_mat;
        // Matrix<double, 2, 3> y_mat;
        // v_col = x_mat;
        // v_row = x_mat;
        // x_mat = v_col;
        // x_mat = v_row;
        // v_col = y_mat;
        // v_row = y_mat;
        // y_mat = v_col;
        // y_mat = v_row;
    }
    {
        Vec<double, 4, false> v_row;
        Vec<double, 4> v_col = to_col(v_row);
        v_col = to_col(v_col);
        v_row = to_row(v_col);
        v_row = to_row(v_row);
        Matrix<double, 4, 1> col_mat = v_col;
        Matrix<double, 1, 4> row_mat = v_row;
        v_row = Vec<double, 4, false>(v_row);
        v_col = Vec<double, 4>(v_col);
        v_row = row_mat;
        v_col = col_mat;

        // ERRORS
        // v_col = row_mat;
        // v_row = col_mat;
        // col_mat = row_mat;
        // row_mat = col_mat;
        // v_col = v_row;
        // v_row = v_col;
        // row_mat = v_col;
        // col_mat = v_row;
        // Matrix<double, 4, 2> x_mat;
        // Matrix<double, 2, 4> y_mat;
        // v_col = x_mat;
        // v_row = x_mat;
        // x_mat = v_col;
        // x_mat = v_row;
        // v_col = y_mat;
        // v_row = y_mat;
        // y_mat = v_col;
        // y_mat = v_row;
    }

    cout << "TESTING 3D VECTOR OPERATIONS" << endl << endl;
    {
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<float> dis(-PI, PI);
        Vec3f v = { -1, 0, -1 };
        cout << v << endl;
        cout << "PITCH: " << v.pitch() << " "
             << "YAW: " << v.yaw() << endl;
        v = { -1, 0, 1 };
        cout << v << endl;
        cout << "PITCH: " << v.pitch() << " "
             << "YAW: " << v.yaw() << endl;
        v = { 1, 0, 1 };
        cout << v << endl;
        cout << "PITCH: " << v.pitch() << " "
             << "YAW: " << v.yaw() << endl;
        v = { 1, 0, -1 };
        cout << v << endl;
        cout << "PITCH: " << v.pitch() << " "
             << "YAW: " << v.yaw() << endl;
        v = { 1, 0, 0 };
        cout << v << endl;
        cout << "PITCH: " << v.pitch() << " "
             << "YAW: " << v.yaw() << endl;
        v = { -1, 0, 0 };
        cout << v << endl;
        cout << "PITCH: " << v.pitch() << " "
             << "YAW: " << v.yaw() << endl;
        v = { 0, 0, 1 };
        cout << v << endl;
        cout << "PITCH: " << v.pitch() << " "
             << "YAW: " << v.yaw() << endl;
        v = { 0, 0, -1 };
        cout << v << endl;
        cout << "PITCH: " << v.pitch() << " "
             << "YAW: " << v.yaw() << endl;
        v = { 0, 1, 0 };
        cout << v << endl;
        cout << "PITCH: " << v.pitch() << " "
             << "YAW: " << v.yaw() << endl;
        v = { 0, -1, 0 };
        cout << v << endl;
        cout << "PITCH: " << v.pitch() << " "
             << "YAW: " << v.yaw() << endl;
        v = { 1, 1, 1 };
        cout << v << endl;
        cout << "PITCH: " << v.pitch() << " "
             << "YAW: " << v.yaw() << endl;
        v = { 1, 1, 0 };
        cout << v << endl;
        cout << "PITCH: " << v.pitch() << " "
             << "YAW: " << v.yaw() << endl;
        v = { 0, 1, 1 };
        cout << v << endl;
        cout << "PITCH: " << v.pitch() << " "
             << "YAW: " << v.yaw() << endl;

        cout << endl;
        v = { 0, 0, 1 };
        v.normalize();
        for (int i = 0; i < 10; ++i) {
            cout << "PITCH: " << v.pitch() << "   "
                 << "YAW: " << v.yaw() << "   ";
            cout << "MAG: " << v.magnitude() << endl;
            v.rotate_x(0.1);
            // v.rotate_yaw(PI / 10);
            cout << v << endl;
        }
        /*cout << endl;
        v.rotate_yaw(PI / 2);
        cout << "ROTATING YAW PI/2" << endl << endl;
        for (int i = 0; i < 100; ++i) {
            cout << "PITCH: " << v.pitch() / PI << "   " << "YAW: " << v.yaw() /
        PI
        << "   ";
            cout << "MAG: " << v.magnitude()  << endl;
            //v.rotate_pitch(PI / 10);
            v.rotate_yaw(PI / 10);
        }*/
        cout << endl;
        cout << "PITCH: " << v.pitch() << "   "
             << "YAW: " << v.yaw() << "   ";
        cout << "MAG: " << v.magnitude() << endl;
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        for (int i = 0; i < NUM_ROTATIONS; ++i) {
            v.rotate_y(0.1);
        }
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        cout << "TIME TAKEN TO DO " << NUM_ROTATIONS
             << " MANY VECTOR ROTATIONS: "
             << duration_cast<duration<float> >(t2 - t1).count() << "s" << endl
             << endl;
    }

#ifdef TEST_MATRIX_VECTOR
    cout << "CLASS VS NON-CLASS MATRIX VECTOR MULTIPLY TEST" << endl << endl;
#ifdef CLASS_FIRST
    gen = mt19937_64(seed);
    cout << "Matrix class: " << test_matrix_vector_class() << "s" << endl;
    gen = mt19937_64(seed);
    cout << "Non class: " << test_matrix_vector_no_class() << "s" << endl;
#else
    gen = mt19937_64(seed);
    cout << "Non class: " << test_matrix_vector_no_class() << "s" << endl;
    gen = mt19937_64(seed);
    cout << "Matrix class: " << test_matrix_vector_class() << "s" << endl;
#endif

    for (unsigned int i = 0; i < HEIGHT; ++i) {
        if (sum_arr[i] != sum_vec[i]) {
            cout << "arr[" << i << "] == " << sum_arr[i] << " : "
                 << "vec[" << i << "] == " << sum_vec[i] << endl;
            cout << "FAILURE" << endl;
            return 1;
        }
    }
    cout << "SUCCESS" << endl;
    cout << endl;
#endif

#ifdef TEST_MATRIX_MATRIX
    cout << "CLASS VS NON-CLASS MATRIX MATRIX MULTIPLY TEST" << endl << endl;
#ifdef CLASS_FIRST
    gen = mt19937_64(seed);
    cout << "Matrix class: " << test_matrix_matrix_class() << "s" << endl;
    gen = mt19937_64(seed);
    cout << "Non class: " << test_matrix_matrix_no_class() << "s" << endl;
#else
    gen = mt19937_64(seed);
    cout << "Non class: " << test_matrix_matrix_no_class() << "s" << endl;
    gen = mt19937_64(seed);
    cout << "Matrix class: " << test_matrix_matrix_class() << "s" << endl;
#endif

    for (unsigned int i = 0; i < HEIGHT; ++i) {
        for (unsigned int j = 0; j < WIDTH; ++j) {
            if (sum_mat[i][j] != sum_mat_arr[i][j]) {
                cout << "arr[" << i << "][" << j << "] == " << sum_mat_arr[i][j]
                     << " : ";
                cout << "mat[" << i << "][" << j << "] == " << sum_mat[i][j]
                     << endl;
                cout << "FAILURE" << endl;
                return 1;
            }
        }
    }
    cout << "SUCCESS" << endl;
    cout << endl;
#endif

    return 0;
}
