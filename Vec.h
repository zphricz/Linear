#ifndef VEC_H
#define VEC_H

#include <cmath>
#include <array>
#include <iomanip>
#include <stdint.h>
#include <sstream>

template <typename T, size_t M, size_t N>
class Matrix;

// Default value for C is column vector
template <typename T, size_t N, bool C = true>
class Vec {
private:
    std::array<T, N> v;

public:
    Vec<T, N, C>() {
        static_assert(N > 0, "Vector cannot have a zero dimension");
    }

    template<bool B>
    Vec<T, N, C>(const Vec<T, N, B>& other) {
        static_assert(N > 0, "Vector cannot have a zero dimension");
        for (size_t i = 0; i < N; ++i) {
            v[i] = other[i];
        }
    }

    template <size_t M>
    Vec<T, M>(const Matrix<T, M, 1>& other) {
        static_assert(M > 0, "Vector cannot have a zero dimension");
        for (size_t i = 0; i < M; ++i) {
            v[i] = other[i][0];
        }
    }

    template <size_t M>
    Vec<T, M, false>(const Matrix<T, 1, M>& other) {
        static_assert(M > 0, "Vector cannot have a zero dimension");
        for (size_t i = 0; i < M; ++i) {
            v[i] = other[0][i];
        }
    }

    Vec<T, 1>(const Matrix<T, 1, 1>& other) {
        v[0] = other[0][0];
    }

    static Vec<T, N, C> zeros() {
        static_assert(N > 0, "Vector cannot have a zero dimension");
        Vec<T, N, C> rval;
        rval.zero();
        return rval;
    }

    ~Vec<T, N, C>() {
    }

    auto magnitude() -> decltype(sqrt(T())) const {
        return sqrt(magnitude_square());
    }

    T magnitude_square() const {
        T sum = T();
        for (size_t i = 0; i < N; ++i) {
            sum += v[i] * v[i];
        }
        return sum;
    }

    void zero() {
        for (size_t i = 0; i < N; ++i) {
            v[i] = T();
        }
    }

    void normalize() {
        *this /= magnitude();
    }

    Vec<T, N, C>& operator+=(const Vec<T, N, C>& rhs) {
        for (size_t i = 0; i < N; ++i) {
            v[i] += rhs[i];
        }
        return *this;
    }

    Vec<T, N, C>& operator-=(const Vec<T, N, C>& rhs) {
        for (size_t i = 0; i < N; ++i) {
            v[i] -= rhs[i];
        }
        return *this;
    }

    Vec<T, N, C>& operator*=(T factor) {
        for (size_t i = 0; i < N; ++i) {
            v[i] *= factor;
        }
        return *this;
    }

    Vec<T, N, C>& operator/=(T factor) {
        for (size_t i = 0; i < N; ++i) {
            v[i] /= factor;
        }
        return *this;
    }

    Vec<T, N, C>& operator=(const Vec<T, N, C>& rhs) {
        v = rhs.v;
        return *this;
    }

    Vec<T, N, C> operator+(const Vec<T, N, C>& rhs) const {
        Vec<T, N, C> new_vec;
        for (size_t i = 0; i < N; ++i) {
            new_vec[i] = v[i] + rhs[i];
        }
        return new_vec;
    }

    Vec<T, N, C> operator-(const Vec<T, N, C>& rhs) const {
        Vec<T, N, C> new_vec;
        for (size_t i = 0; i < N; ++i) {
            new_vec[i] = v[i] - rhs[i];
        }
        return new_vec;
    }

    Vec<T, N, C> operator*(T factor) const {
        Vec<T, N, C> new_vec;
        for (size_t i = 0; i < N; ++i) {
            new_vec[i] = v[i] * factor;
        }
        return new_vec;
    }

    Vec<T, N, C> operator/(T factor) const {
        Vec<T, N, C> new_vec;
        for (size_t i = 0; i < N; ++i) {
            new_vec[i] = v[i] / factor;
        }
        return new_vec;
    }

    T& operator[](size_t index) {
        return v[index];
    }

    T operator[](size_t index) const {
        return v[index];
    }

    bool operator==(const Vec<T, N, C>& rhs) const {
        return v == rhs.v;
    }

    bool operator!=(const Vec<T, N, C>& rhs) const {
        return v != rhs.v;
    }

    T dot(const Vec<T, N, C>& other) const {
        T sum = T();
        for (size_t i = 0; i < N; ++i) {
            sum += v[i] * other[i];
        }
        return sum;
    }

    size_t size() const {
        return N;
    }
};

/* Prints out the vector as a column vector
 */
template <typename T, size_t N, bool C>
std::ostream& operator<<(std::ostream& os, const Vec<T, N, C>& rhs) {
    if (!C) {
        os << "[ ";
    }
    size_t max_size = 0;
    for (size_t i = 0; i < N; ++i) {
        std::ostringstream s;
        s << rhs[i];
        std::string str = s.str();
        if (str.size() > max_size) {
            max_size = str.size();
        }
    }
    for (size_t i = 0; i < N; ++i) {
        if (C) {
            os << "[ " << std::setw(max_size) << rhs[i] << " ]";
            if (i != N - 1) {
                os << std::endl;
            }
        } else {
            os << rhs[i] << " ";
        }
    }
    if (!C) {
        os << "]";
    }
    return os;
}

template <typename T, typename U, size_t N, bool C>
Vec<T, N, C> operator*(U f, const Vec<T, N, C>& rhs) {
    T factor = T(f);
    Vec<T, N, C> new_vec;
    for (size_t i = 0; i < N; ++i) {
        new_vec[i] = factor * rhs[i];
    }
    return new_vec;
}


template <typename T, bool C>
class Vec<T, 2, C> {
public:
    T x;
    T y;

    Vec<T, 2, C>() {
    }

    template<bool B>
    Vec<T, 2, C>(const Vec<T, 2, B>& other) :
        x(other.x),
        y(other.y) {
    }

    Vec<T, 2, C>(const Matrix<T, 1, 2>& other) :
        x(other[0][0]),
        y(other[0][1]) {
    }

    Vec<T, 2, C>(const Matrix<T, 2, 1>& other) :
        x(other[0][0]),
        y(other[1][0]) {
    }

    Vec<T, 2, C>(T angle) :
        x(cos(angle)),
        y(sin(angle)) {
    }

    Vec<T, 2, C>(T x, T y) :
        x(x),
        y(y) {
    }

    static Vec<T, 2, C> zeros() {
        Vec<T, 2, C> rval;
        rval.zero();
        return rval;
    }

    ~Vec<T, 2, C>() {
    }

    auto theta() -> decltype(atan2(T(), T())) const {
        return atan2(y, x);
    }

    auto magnitude() -> decltype(sqrt(T())) const {
        return sqrt(magnitude_square());
    }

    T magnitude_square() const {
        return x * x + y * y;
    }

    void zero() {
        x = T();
        y = T();
    }

    void normalize() {
        *this /= magnitude();
    }

    Vec<T, 2, C>& operator+=(const Vec<T, 2, C>& rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    Vec<T, 2, C>& operator-=(const Vec<T, 2, C>& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    Vec<T, 2, C>& operator*=(T factor) {
        x *= factor;
        y *= factor;
        return *this;
    }

    Vec<T, 2, C>& operator/=(T factor) {
        x /= factor;
        y /= factor;
        return *this;
    }

    Vec<T, 2, C>& operator=(const Vec<T, 2, C>& rhs) {
        x = rhs.x;
        y = rhs.y;
        return *this;
    }

    Vec<T, 2, C> operator+(const Vec<T, 2, C>& rhs) const {
        return Vec<T, 2, C>(x + rhs.x, y + rhs.y);
    }

    Vec<T, 2, C> operator-(const Vec<T, 2, C>& rhs) const {
        return Vec<T, 2, C>(x - rhs.x, y - rhs.y);
    }

    Vec<T, 2, C> operator*(T factor) const {
        return Vec<T, 2, C>(x * factor, y * factor);
    }

    Vec<T, 2, C> operator/(T factor) const {
        return Vec<T, 2, C>(x / factor, y / factor);
    }

    // Not safe, Prefer accessing x and y directly to this
    T& operator[](size_t index) {
        return *(&x + index);
    }

    // Not safe, Prefer accessing x and y directly to this
    T operator[](size_t index) const {
        return *(&x + index);
    }

    bool operator==(const Vec<T, 2, C>& rhs) const {
        return x == rhs.x && y == rhs.y;
    }

    bool operator!=(const Vec<T, 2, C>& rhs) const {
        return !(*this == rhs);
    }

    T dot(const Vec<T, 2, C>& other) const {
        return x * other.x + y * other.y;
    }

    size_t size() const {
        return 2;
    }
};

typedef Vec<float, 2> Vec2f;
typedef Vec<double, 2> Vec2d;
typedef Vec<int8_t, 2> Vec2i8;
typedef Vec<int32_t, 2> Vec2i32;
typedef Vec<int64_t, 2> Vec2i64;
typedef Vec<uint8_t, 2> Vec2u8;
typedef Vec<uint32_t, 2> Vec2u32;
typedef Vec<uint64_t, 2> Vec2u64;

template <typename T, bool C>
class Vec<T, 3, C> {
public:
    T x;
    T y;
    T z;

    Vec<T, 3, C>() {
    }

    template<bool B>
    Vec<T, 3, C>(const Vec<T, 3, C>& other) :
        x(other.x),
        y(other.y),
        z(other.z) {
    }

    Vec<T, 3, C>(const Matrix<T, 1, 3>& other) :
        x(other[0][0]),
        y(other[0][1]),
        z(other[0][2]) {
    }

    Vec<T, 3, C>(const Matrix<T, 3, 1>& other) :
        x(other[0][0]),
        y(other[1][0]),
        z(other[2][0]) {
    }

    Vec<T, 3, C>(T x, T y, T z) :
        x(x),
        y(y),
        z(z) {
    }

    static Vec<T, 3, C> zeros() {
        Vec<T, 3, C> rval;
        rval.zero();
        return rval;
    }

    ~Vec<T, 3, C>() {
    }

    auto magnitude() -> decltype(sqrt(T())) const {
        return sqrt(magnitude_square());
    }

    T magnitude_square() const {
        return x * x + y * y + z * z;
    }

    void zero() {
        x = T();
        y = T();
        z = T();
    }

    void normalize() {
        *this /= magnitude();
    }

    Vec<T, 3, C>& operator+=(const Vec<T, 3, C>& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    Vec<T, 3, C>& operator-=(const Vec<T, 3, C>& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    Vec<T, 3, C>& operator*=(T factor) {
        x *= factor;
        y *= factor;
        z *= factor;
        return *this;
    }

    Vec<T, 3, C>& operator/=(T factor) {
        x /= factor;
        y /= factor;
        z /= factor;
        return *this;
    }

    Vec<T, 3, C>& operator=(const Vec<T, 3, C>& rhs) {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        return *this;
    }

    Vec<T, 3, C> operator+(const Vec<T, 3, C>& rhs) const {
        return Vec<T, 3, C>(x + rhs.x, y + rhs.y, z + rhs.z);
    }

    Vec<T, 3, C> operator-(const Vec<T, 3, C>& rhs) const {
        return Vec<T, 3, C>(x - rhs.x, y - rhs.y, z - rhs.z);
    }

    Vec<T, 3, C> operator*(T factor) const {
        return Vec<T, 3, C>(x * factor, y * factor, z * factor);
    }

    Vec<T, 3, C> operator/(T factor) const {
        return Vec<T, 3, C>(x / factor, y / factor, z / factor);
    }

    // Not safe, Prefer accessing x, y, and z directly to this
    T& operator[](size_t index) {
        return *(&x + index);
    }

    // Not safe, Prefer accessing x, y, and z directly to this
    T operator[](size_t index) const {
        return *(&x + index);
    }

    bool operator==(const Vec<T, 3, C>& rhs) const {
        return x == rhs.x && y == rhs.y && z == rhs.z;
    }

    bool operator!=(const Vec<T, 3, C>& rhs) const {
        return !(*this == rhs);
    }

    T dot(const Vec<T, 3, C>& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    size_t size() const {
        return 3;
    }

    T pitch() const {
        return 0.0;
    }

    T roll() const {
        return 0.0;
    }

    T yaw() const {
        return 0.0;
    }

    Vec<T, 3, C> cross(const Vec<T, 3, C>& other) const {
        return Vec<T, 3, C>(y * other.z - z * other.y,
                         z * other.x - x * other.z,
                         x * other.y - y * other.x);
    }
};


typedef Vec<float, 3> Vec3f;
typedef Vec<double, 3> Vec3d;
typedef Vec<int8_t, 3> Vec3i8;
typedef Vec<int32_t, 3> Vec3i32;
typedef Vec<int64_t, 3> Vec3i64;
typedef Vec<uint8_t, 3> Vec3u8;
typedef Vec<uint32_t, 3> Vec3u32;
typedef Vec<uint64_t, 3> Vec3u64;

#endif

