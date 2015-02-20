#include <iostream>
#include <chrono>
#include <random>
#include "Vec.h"
#include "Matrix.h"

using namespace std;
using namespace chrono;

//#define CLASS_FIRST
#if 1
const int TEST_SIZE = 100000; // hundred thousand
const int HEIGHT = 100;
const int WIDTH = 73;
#else
const int TEST_SIZE = 100;
const int HEIGHT = 10000; // ten thousand
const int WIDTH = 1000;
#endif
random_device rd;
auto seed = rd();
mt19937_64 gen;
uniform_real_distribution<> dis(-100.0, 100.0);

static Vec<double, HEIGHT> sum_vec;
static double sum_arr[HEIGHT];

static double test_class() {
    auto t1 = high_resolution_clock::now();
    sum_vec.zero();
    for (int x = 0; x < TEST_SIZE; ++x) {
        static Matrix<double, HEIGHT, WIDTH> mat;
        static Vec<double, WIDTH> vec;
        for (int j = 0; j < vec.size(); ++j) {
            vec[j] = dis(gen);
        }
        for (int i = 0; i < mat.height(); ++i) {
            for (int j = 0; j < mat.width(); ++j) {
                mat[i][j] = dis(gen);
            }
        }
        static Vec<double, HEIGHT> result;
        result = mat * vec;
        sum_vec += result;
    }
    auto t2 = high_resolution_clock::now();
    return duration_cast<duration<double>>(t2 - t1).count();
}

static double test_no_class() {
    auto t1 = high_resolution_clock::now();
    for (int i = 0; i < HEIGHT; ++i) {
        sum_arr[i] = 0.0;
    }
    for (int x = 0; x < TEST_SIZE; ++x) {
        static double mat[HEIGHT][WIDTH];
        static double vec[WIDTH];
        for (int j = 0; j < WIDTH; ++j) {
            vec[j] = dis(gen);
        }
        for (int i = 0; i < HEIGHT; ++i) {
            for (int j = 0; j < WIDTH; ++j) {
                mat[i][j] = dis(gen);
            }
        }
        static double result[HEIGHT];
        for (int i = 0; i < HEIGHT; ++i) {
            result[i] = 0.0;
            for (int j = 0; j < WIDTH; ++j) {
                result[i] += vec[j] * mat[i][j];
            }
        }
        for (int i = 0; i < HEIGHT; ++i) {
            sum_arr[i] += result[i];
        }
    }
    auto t2 = high_resolution_clock::now();
    return duration_cast<duration<double>>(t2 - t1).count();
}

int main() {
    Vec2i64 v1;
    cout << v1 << endl;
    v1.x() = 4;
    v1.y() = 5;
    v1 = v1 * 8;
    v1 = 8 * v1;
    cout << v1 << endl;
    cout << "THETA: " << v1.theta() << endl;
    cout << "MAG: " << v1.magnitude() << endl;

    cout << endl << "MATRIX STUFF" << endl << endl;

    Matrix<double, 5, 4> mt;
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 4; ++j) {
            mt[i][j] = i * j;
        }
    }
    cout << "FIRST" << endl << endl;
    cout << mt << endl << endl;;
    auto m = mt.transpose();
    cout << "TRANSPOSED" << endl << endl;
    cout << m << endl << endl;

    Vec<double, 5> rhs;
    for (int j = 0; j < rhs.size(); ++j) {
        rhs[j] = j;
    }
    cout << rhs << endl << endl;

    auto lhs = m * rhs;
    cout << lhs << endl;
    cout << "CLASS VS NON-CLASS MATRIX VECTOR MULTIPLY TEST" << endl << endl;

#ifdef CLASS_FIRST
    gen = mt19937_64(seed);
    cout << "Matrix class: " << test_class() << "s" << endl;
    gen = mt19937_64(seed);
    cout << "Non class: " << test_no_class() << "s" << endl;
#else
    gen = mt19937_64(seed);
    cout << "Non class: " << test_no_class() << "s" << endl;
    gen = mt19937_64(seed);
    cout << "Matrix class: " << test_class() << "s" << endl;
#endif

    for (int i = 0; i < HEIGHT; ++i) {
        if (sum_arr[i] != sum_vec[i]) {
            cout << "arr[" << i << "] == " << sum_arr[i] << " : " <<
                    "vec[" << i << "] == " << sum_vec[i] << endl;
            cout << "FAILURE" << endl;
            return 1;
        }
    }
    cout << "SUCCESS" << endl;
    return 0;
}

