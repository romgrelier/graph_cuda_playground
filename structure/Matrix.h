#pragma once

#include <valarray>

template <typename T>
class Matrix
{
private:
    const size_t row;
    const size_t col;
    std::valarray<T> matrix;

public:
    Matrix(const size_t row, const size_t col, const T init) : row{row}, col{col}, matrix(init, row * col) { }

    size_t getRow() const noexcept { return row; }

    size_t getCol() const noexcept { return col; }

    T &operator()(const size_t y, const size_t x) { return matrix[col * y + x]; }

    const T &operator()(const size_t y, const size_t x) const { return matrix[col * y + x]; }
};