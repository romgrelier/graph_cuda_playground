#pragma once

#include "../structure/Matrix.h"

#include <string>
#include <string_view>
#include <fstream>
#include <sstream>
#include <memory>

class DIMACS
{
private:
    std::unique_ptr<Matrix<bool>> matrix;
    size_t vertexCount;
    size_t edgeCount;

public:
    DIMACS(const std::string filename);

    size_t getVertexCount() const noexcept { return vertexCount; }

    size_t getEdgeCount() const noexcept { return edgeCount; }

    bool operator()(const size_t y, const size_t x) const { return (*matrix)(y, x); }

    Matrix<bool> &getGraph() const { return *matrix; }
};
