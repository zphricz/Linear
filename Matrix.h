#ifndef MATRIX_H
#define MATRIX_H

#include "Vec.h"

template <typename T, size_t M, size_t N = M>
class Matrix {
private:
    std::array<std::array<T, N>, M> v;

public:
    Matrix<T, M, N>() {
    }

    Matrix<T, M, N>(const Matrix<T, M, N>& other) {
        v = other.v;
    }

    ~Matrix<T, M, N>() {
    }

    static Matrix<T, M, N> zeros() {
        Matrix<T, M, N> rval;
        rval.zero();
        return rval;
    }

    static Matrix<T, 1, N> from_row(const Vec<T, N>& other) {
        Matrix<T, 1, N> rval;
        for (size_t i = 0; i < N; ++i) {
            rval[0][i] = other[i];
        }
        return rval;
    }

    static Matrix<T, M, 1> from_col(const Vec<T, M>& other) {
        Matrix<T, M, 1> rval;
        for (size_t i = 0; i < M; ++i) {
            rval[i][0] = other[i];
        }
        return rval;
    }

    static Matrix<T, M, N> identity() {
        static_assert(M == N, "Identity matrix must be square");
        Matrix<T, M, N> rval;
        rval.zero();
        for (size_t i = 0; i < M; ++i) {
            rval[i][i] = T(1);
        }
        return rval;
    }

    void zero() {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                v[i][j] = T();
            }
        }
    }

    Vec<T, N> row(size_t m) const {
        Vec<T, N> rval;
        for (size_t i = 0; i < N; ++i) {
            rval[i] = v[m][i];
        }
        return rval;
    }

    Vec<T, M> col(size_t n) const {
        Vec<T, M> rval;
        for (size_t i = 0; i < M; ++i) {
            rval[i] = v[i][n];
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

    std::array<T, N>& operator[](size_t i) {
        return v[i];
    }

    const std::array<T, N>& operator[](size_t i) const {
        return v[i];
    }

    Matrix<T, M, N> operator+(const Matrix<T, M, N>& rhs) const {
        Matrix<T, M, N> rval;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                rval[i][j] = v[i][j] + rhs[i][j];
            }
        }
        return rval;
    }

    Matrix<T, M, N> operator-(const Matrix<T, M, N>& rhs) const {
        Matrix<T, M, N> rval;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                rval[i][j] = v[i][j] - rhs[i][j];
            }
        }
        return rval;
    }

    Matrix<T, M, N>& operator+=(const Matrix<T, M, N>& rhs) {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                v[i][j] += rhs[i][j];
            }
        }
        return *this;
    }

    Matrix<T, M, N>& operator-=(const Matrix<T, M, N>& rhs) {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                v[i][j] -= rhs[i][j];
            }
        }
        return *this;
    }

    Vec<T, M> operator*(const Vec<T, N>& rhs) const {
        Vec<T, M> rval;
        for (size_t i = 0; i < M; ++i) {
            rval[i] = T();
            for (size_t j = 0; j < N; ++j) {
                rval[i] += v[i][j] * rhs[j];
            }
        }
        return rval;
    }

    Matrix<T, M, N> operator*(T factor) const {
        Matrix<T, M, N> rval;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                rval[i][j] = v[i][j] * factor;
            }
        }
        return rval;
    }

    template <size_t O>
    Matrix<T, M, O> operator*(const Matrix<T, N, O>& rhs) const {
        Matrix<T, M, O> rval;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < O; ++j) {
                rval[i][j] = T();
                for (size_t k = 0; k < N; ++k) {
                    rval[i][j] += v[i][k] * rhs[k][j];
                }
            }
        }
        return rval;
    }

    Matrix<T, M, N> operator/(T factor) const {
        Matrix<T, M, N> rval;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                rval[i][j] = v[i][j] / factor;
            }
        }
        return rval;
    }

    Matrix<T, M, N>& operator*=(T factor) {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                v[i][j] *= factor;
            }
        }
        return *this;
    }
    
    Matrix<T, M, N>& operator*=(const Matrix<T, M, N>& rhs) {
        static_assert(M == N, "Matrix dimensions must agree for in place "
                              "matrix multiply");
        return *this = *this * rhs;
    }

    Matrix<T, M, N>& operator/=(T factor) {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                v[i][j] /= factor;
            }
        }
        return *this;
    }
};

/* This needed the extra template parameter to get the compiler to stop
 * producing errors on the potentially mismatched types of the Vector
 * components and factor */
template <typename T, typename U, size_t M, size_t N>
Matrix<T, M, N> operator*(U f, const Matrix<T, M, N>& rhs) {
    T factor = T(f);
    Matrix<T, M, N> new_mat;
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            new_mat[i][j] = factor * rhs[i][j];
        }
    }
    return new_mat;
}

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
