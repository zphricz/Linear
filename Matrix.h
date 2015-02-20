#ifndef MATRIX_H
#define MATRIX_H

#include "Vec.h"

template <typename T, size_t M, size_t N>
class Matrix {
private:
    std::array<Vec<T, N>, M> v;

public:
    Matrix<T, M, N>() {
    }

    ~Matrix<T, M, N>() {
    }

    Vec<T, N>& operator[](size_t i) {
        return v[i];
    }

    const Vec<T, N>& operator[](size_t i) const {
        return v[i];
    }
    
    template <size_t U>
    Matrix<T, M, U> operator*(const Matrix<T, N, U>& rhs) {
        return *this;
    }

    Vec<T, M> operator*(const Vec<T, N>& rhs) {
        Vec<T, M> rval;
        for (size_t i = 0; i < M; ++i) {
            rval[i] = v[i].dot(rhs);
        }
        return rval;
    }

    Matrix<T, N, M> transpose() const {
        Matrix<T, N, M> rval;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                rval[j][i] = v[i][j];
            }
        }
        return rval;
    }

    size_t height() const {
        return M;
    }

    size_t width() const {
        return N;
    }
};

template <typename T, size_t M, size_t N>
std::ostream& operator<<(std::ostream& os, const Matrix<T, M, N>& rhs) {
    Vec<size_t, N> v;
    for (size_t j = 0; j < N; ++j) {
        v[j] = 0;
        for (size_t i = 0; i < M; ++i) {
            std::ostringstream s;
            s << rhs[i][j];
            std::string str = s.str();
            if (str.size() > v[j]) {
                v[j] = str.size();
            }
        }
    }
    for (size_t i = 0; i < M; ++i) {
        os << "[ ";
        for (size_t j = 0; j < N; ++j) {
            os << std::setw(v[j]) << rhs[i][j] << " ";
        }
        os << "]";
        if (i != M - 1) {
            os << std::endl;
        }
    }
    return os;
}

#endif
