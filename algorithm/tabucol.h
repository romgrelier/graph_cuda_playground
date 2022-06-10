#include "../structure/Matrix.h"

namespace tabucol
{
    namespace cpu
    {
        void tabucol(const Matrix<bool> &graph, std::valarray<unsigned int> &solution, const unsigned int max_iter, const unsigned int k, const unsigned int vertex_count);
    }

    namespace cuda
    {
        void tabucol(const Matrix<bool> &graph, std::valarray<unsigned int> &solution, const unsigned int max_iter, const unsigned int k, const unsigned int vertex_count);
    }
}
