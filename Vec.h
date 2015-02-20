#ifndef VEC_H
#define VEC_H

#include <cmath>
#include <array>
#include <iomanip>
#include <stdint.h>
#include <sstream>

template <typename T, size_t N>
class Vec {
private:
    std::array<T, N> v;

public:
    Vec<T, N>() {
    }

    Vec<T, N>(const Vec<T, N>& other) {
        v = other.v;
    }

    static Vec<T, N> zeros() {
        Vec<T, N> rval;
        rval.zero();
        return rval;
    }

    ~Vec<T, N>() {
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

    Vec<T, N>& operator+=(const Vec<T, N>& rhs) {
        for (size_t i = 0; i < N; ++i) {
            v[i] += rhs[i];
        }
        return *this;
    }

    Vec<T, N>& operator-=(const Vec<T, N>& rhs) {
        for (size_t i = 0; i < N; ++i) {
            v[i] -= rhs[i];
        }
        return *this;
    }

    Vec<T, N>& operator*=(T factor) {
        for (size_t i = 0; i < N; ++i) {
            v[i] *= factor;
        }
        return *this;
    }

    Vec<T, N>& operator/=(T factor) {
        for (size_t i = 0; i < N; ++i) {
            v[i] /= factor;
        }
        return *this;
    }

    Vec<T, N>& operator=(const Vec<T, N>& rhs) {
        v = rhs.v;
        return *this;
    }

    Vec<T, N> operator+(const Vec<T, N>& rhs) const {
        Vec<T, N> new_vec;
        for (size_t i = 0; i < N; ++i) {
            new_vec[i] = v[i] + rhs[i];
        }
        return new_vec;
    }

    Vec<T, N> operator-(const Vec<T, N>& rhs) const {
        Vec<T, N> new_vec;
        for (size_t i = 0; i < N; ++i) {
            new_vec[i] = v[i] - rhs[i];
        }
        return new_vec;
    }

    Vec<T, N> operator*(T factor) const {
        Vec<T, N> new_vec;
        for (size_t i = 0; i < N; ++i) {
            new_vec[i] = v[i] * factor;
        }
        return new_vec;
    }

    Vec<T, N> operator/(T factor) const {
        Vec<T, N> new_vec;
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

    bool operator==(const Vec<T, N>& rhs) const {
        return v == rhs.v;
    }

    bool operator!=(const Vec<T, N>& rhs) const {
        return v != rhs.v;
    }

    T dot(const Vec<T, N>& other) const {
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
template <typename T, size_t N>
std::ostream& operator<<(std::ostream& os, const Vec<T, N>& rhs) {
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
        os << "[ " << std::setw(max_size) << rhs[i] << " ]";
        if (i != N - 1) {
            os << std::endl;
        }
    }
    return os;
}

/* This needed the extra template parameter to get the compiler to stop
 * producing errors on the potentially mismatched types of the Vector
 * components and factor */
template <typename T, typename U, size_t N>
Vec<T, N> operator*(U f, const Vec<T, N>& rhs) {
    T factor = T(f);
    Vec<T, N> new_vec;
    for (size_t i = 0; i < N; ++i) {
        new_vec[i] = factor * rhs[i];
    }
    return new_vec;
}


template <typename T>
class Vec<T, 2> {
public:
    T x;
    T y;

    Vec<T, 2>() {
    }

    Vec<T, 2>(T angle) :
        x(cos(angle)),
        y(sin(angle)) {
    }

    Vec<T, 2>(const Vec<T, 2>& other) :
        x(other.x),
        y(other.y) {
    }

    Vec<T, 2>(T x, T y) :
        x(x),
        y(y) {
    }

    ~Vec<T, 2>() {
    }

    auto theta() -> decltype(atan2(T(), T())) const {
        return atan2(y, x);
    }

    static Vec<T, 2> zeros() {
        Vec<T, 2> rval;
        rval.zero();
        return rval;
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

    Vec<T, 2>& operator+=(const Vec<T, 2>& rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    Vec<T, 2>& operator-=(const Vec<T, 2>& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    Vec<T, 2>& operator*=(T factor) {
        x *= factor;
        y *= factor;
        return *this;
    }

    Vec<T, 2>& operator/=(T factor) {
        x /= factor;
        y /= factor;
        return *this;
    }

    Vec<T, 2>& operator=(const Vec<T, 2>& rhs) {
        x = rhs.x;
        y = rhs.y;
        return *this;
    }

    Vec<T, 2> operator+(const Vec<T, 2>& rhs) const {
        return Vec<T, 2>(x + rhs.x, y + rhs.y);
    }

    Vec<T, 2> operator-(const Vec<T, 2>& rhs) const {
        return Vec<T, 2>(x - rhs.x, y - rhs.y);
    }

    Vec<T, 2> operator*(T factor) const {
        return Vec<T, 2>(x * factor, y * factor);
    }

    Vec<T, 2> operator/(T factor) const {
        return Vec<T, 2>(x / factor, y / factor);
    }

    // Not safe, Prefer accessing x and y directly to this
    T& operator[](size_t index) {
        return *(&x + index);
    }

    // Not safe, Prefer accessing x and y directly to this
    T operator[](size_t index) const {
        return *(&x + index);
    }

    bool operator==(const Vec<T, 2>& rhs) const {
        return x == rhs.x && y == rhs.y;
    }

    bool operator!=(const Vec<T, 2>& rhs) const {
        return !(*this == rhs);
    }

    T dot(const Vec<T, 2>& other) const {
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

template <typename T>
class Vec<T, 3> {
public:
    T x;
    T y;
    T z;

    Vec<T, 3>() {
    }

    Vec<T, 3>(const Vec<T, 3>& other) :
        x(other.x),
        y(other.y),
        z(other.z) {
    }

    Vec<T, 3>(T x, T y, T z) :
        x(x),
        y(y),
        z(z) {
    }

    ~Vec<T, 3>() {
    }

    static Vec<T, 3> zeros() {
        Vec<T, 3> rval;
        rval.zero();
        return rval;
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

    Vec<T, 3>& operator+=(const Vec<T, 3>& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    Vec<T, 3>& operator-=(const Vec<T, 3>& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    Vec<T, 3>& operator*=(T factor) {
        x *= factor;
        y *= factor;
        z *= factor;
        return *this;
    }

    Vec<T, 3>& operator/=(T factor) {
        x /= factor;
        y /= factor;
        z /= factor;
        return *this;
    }

    Vec<T, 3>& operator=(const Vec<T, 3>& rhs) {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        return *this;
    }

    Vec<T, 3> operator+(const Vec<T, 3>& rhs) const {
        return Vec<T, 3>(x + rhs.x, y + rhs.y, z + rhs.z);
    }

    Vec<T, 3> operator-(const Vec<T, 3>& rhs) const {
        return Vec<T, 3>(x - rhs.x, y - rhs.y, z - rhs.z);
    }

    Vec<T, 3> operator*(T factor) const {
        return Vec<T, 3>(x * factor, y * factor, z * factor);
    }

    Vec<T, 3> operator/(T factor) const {
        return Vec<T, 3>(x / factor, y / factor, z / factor);
    }

    // Not safe, Prefer accessing x, y, and z directly to this
    T& operator[](size_t index) {
        return *(&x + index);
    }

    // Not safe, Prefer accessing x, y, and z directly to this
    T operator[](size_t index) const {
        return *(&x + index);
    }

    bool operator==(const Vec<T, 3>& rhs) const {
        return x == rhs.x && y == rhs.y && z == rhs.z;
    }

    bool operator!=(const Vec<T, 3>& rhs) const {
        return !(*this == rhs);
    }

    T dot(const Vec<T, 3>& other) const {
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

    Vec<T, 3> cross(const Vec<T, 3>& other) const {
        return Vec<T, 3>(y * other.z - z * other.y,
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

