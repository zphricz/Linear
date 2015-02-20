#ifndef MATRIX_H
#define MATRIX_H

#include "Vec.h"

template <typename T, size_t M, size_t N = M>
class Matrix {
private:
    std::array<Vec<T, N>, M> v;

public:
    Matrix<T, M, N>() {
    }

    ~Matrix<T, M, N>() {
    }

    static Matrix<T, M, N> zeros() {
        Matrix<T, M, N> rval;
        rval.zero();
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

    const Vec<T, N>& row(size_t m) const {
        return v[m];
    }

    // Prefer not to use this function as it is not fast
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

    Vec<T, N>& operator[](size_t i) {
        return v[i];
    }

    const Vec<T, N>& operator[](size_t i) const {
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
            rval[i] = row(i).dot(rhs);
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

    template <size_t U>
    Matrix<T, M, U> operator*(const Matrix<T, N, U>& rhs) const {
        Matrix<T, M, U> rval;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < U; ++j) {
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

        /*
         * This implementation is more performant, but fails when &rhs == this
         */
        for (size_t i = 0; i < M; ++i) {
            auto r = row(i);
            for (size_t j = 0; j < M; ++j) {
                v[i][j] = T();
                for (size_t k = 0; k < M; ++k) {
                    v[i][j] += r[k] * rhs[k][j];
                }
            }
        }
        return *this;
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
