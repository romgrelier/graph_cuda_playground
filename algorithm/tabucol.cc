#include "tabucol.h"

#include <random>
#include <chrono>
#include <tuple>

namespace tabucol
{
    namespace cpu
    {
        void tabucol(const Matrix<bool> &graph, std::valarray<unsigned int> &solution, const unsigned int max_iter, const unsigned int k, const unsigned int vertex_count)
        {
            std::mt19937_64 generator(std::chrono::system_clock::now().time_since_epoch().count());
            std::uniform_real_distribution<double> distribution(0.0, 1.0);

            // The gamma_conflict matrix gives which color is connected to a vertex
            Matrix<unsigned int> gamma_conflict(vertex_count, k, 0);
            Matrix<unsigned int> tabu(vertex_count, k, 0);

            // init gamma
            int conflict_count{0};
            for (size_t y{0}; y < vertex_count; ++y)
            {
                for (size_t x{0}; x < y; ++x)
                {
                    if (graph(y, x))
                    {
                        gamma_conflict(y, solution[x]) += 1;
                        gamma_conflict(x, solution[y]) += 1;

                        if (solution[x] == solution[y])
                        {
                            ++conflict_count;
                        }
                    }
                }
            }

            // best solution
            std::valarray<unsigned int> best_solution(solution);
            int best_objective{conflict_count};

            unsigned int i{0};
            while (i < max_iter && best_objective != 0)
            {
                std::tuple<size_t, size_t> best_update(0, 0); // 0: vertex, 1: color
                int best_improvement = gamma_conflict(std::get<0>(best_update), std::get<1>(best_update)) - gamma_conflict(std::get<0>(best_update), solution[std::get<0>(best_update)]);

                // for each vertex we search for the best new color to reduce the count of conflicts
                for (size_t v{0}; v < vertex_count; ++v)
                {
                    // if it has at least one conflict
                    if (gamma_conflict(v, solution[v]) > 0)
                    {
                        // try each available color
                        for (size_t c{0}; c < k; ++c)
                        {
                            // if the color is different and not tabu with this vertex
                            if (c != solution[v] && !(tabu(v, c) > i))
                            {
                                // compute the improvement
                                int improvement = gamma_conflict(v, c) - gamma_conflict(v, solution[v]);

                                // update the best improvement found
                                if (improvement < best_improvement)
                                {
                                    best_update = std::make_tuple(v, c);
                                    best_improvement = improvement;
                                }
                                else if (improvement == best_improvement && distribution(generator) < 0.6)
                                {
                                    best_update = std::make_tuple(v, c);
                                    best_improvement = improvement;
                                }
                            }
                        }
                    }
                }

                // update the solution with the best improvement found
                conflict_count += best_improvement;
                const unsigned int former_color = solution[std::get<0>(best_update)];
                solution[std::get<0>(best_update)] = std::get<1>(best_update);

                // update the gamma conflict matrix
                for (size_t n{0}; n < vertex_count; ++n)
                {
                    if (graph(std::get<0>(best_update), n))
                    {
                        gamma_conflict(n, former_color) -= 1;
                        gamma_conflict(n, std::get<1>(best_update)) += 1;
                    }
                }

                // update the tabu list
                tabu(std::get<0>(best_update), former_color) = i + static_cast<unsigned int>(0.6 * conflict_count);

                // update the best solution found
                if (conflict_count < best_objective)
                {
                    best_objective = conflict_count;
                    for (size_t v{0}; v < vertex_count; ++v)
                    {
                        best_solution[v] = solution[v];
                    }
                }

                ++i;
            }

            // copy the best solution
            for (size_t v{0}; v < vertex_count; ++v)
            {
                solution[v] = best_solution[v];
            }
        }
    }

    namespace cuda
    {

    }
}
