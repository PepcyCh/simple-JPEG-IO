#pragma once

#include <cstdlib>
#include <algorithm>
#include <iostream>

template <typename T>
class Mat {
public:
    Mat(size_t rs = 1, size_t cs = 1, T filled = static_cast<T>(0)) : rs(rs), cs(cs) {
        data = new T[rs * cs];
        std::fill(data, data + rs * cs, filled);
    }
    ~Mat() {
        delete[] data;
    }

    Mat(const Mat<T>& rhs) {
        rs = rhs.rs;
        cs = rhs.cs;
        data = new T[rs * cs];
        std::copy(rhs[0], rhs[0] + rs * cs, data);
    }

    Mat &operator=(const Mat<T>& rhs) {
        rs = rhs.rs;
        cs = rhs.cs;
        T* pOld = data;
        data = new T[rs * cs];
        std::copy(rhs[0], rhs[0] + rs * cs, data);
        delete[] pOld;
        return *this;
    }

    void create(size_t rs, size_t cs, T filled = static_cast<T>(0)) {
        this->rs = rs;
        this->cs = cs;
        delete[] data;
        data = new T[rs * cs];
        std::fill(data, data + rs * cs, filled);
    }

    size_t rows() const { return rs; }
    size_t cols() const { return cs; }

    T* operator[](size_t i) { return data + i * cs; }
    const T* operator[](size_t i) const { return data + i * cs; }

private:
    size_t rs, cs;
    T* data = nullptr;
};
