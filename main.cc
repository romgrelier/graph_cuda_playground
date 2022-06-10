#include <iostream>
#include <assert.h>
#include <valarray>
#include <tuple>
#include <algorithm>
#include <random>
#include <chrono>

#include "format/dimacs.h"
#include "algorithm/tabucol.h"

unsigned int conflict(Matrix<bool> &graph, unsigned int *solution, const unsigned int vertex_count)
{
    unsigned int counter{0};

    for (size_t y{0}; y < vertex_count; ++y)
    {
        for (size_t x{0}; x < y; ++x)
        {
            if (graph(y, x) && solution[x] == solution[y])
            {
                ++counter;
            }
        }
    }

    return counter;
}

int main(int argc, char *argv[])
{
    std::mt19937_64 generator(std::chrono::system_clock::now().time_since_epoch().count());

    DIMACS instance("dimacs_instances/DSJC1000.5.col");

    const unsigned int targetColorCount{85};

    std::uniform_int_distribution<unsigned int> distribution(0, targetColorCount - 1);

    // init a random solution
    std::valarray<unsigned int> solution(instance.getVertexCount());
    std::generate(begin(solution), end(solution), [&]
                  { return distribution(generator); });

    unsigned int conflictCount{conflict(instance.getGraph(), &solution[0], instance.getVertexCount())};

    tabucol::cpu::tabucol(instance.getGraph(), solution, 1'000, targetColorCount, instance.getVertexCount());

    std::cout << "tabucol conflict : " << conflictCount << " -> " << conflict(instance.getGraph(), &solution[0], instance.getVertexCount()) << '\n';

    for (auto e : solution)
    {
        std::cout << e << ' ';
    }
    std::cout << std::endl;

    return 0;
}