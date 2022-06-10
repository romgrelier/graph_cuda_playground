#include "graph/format/dimacs.h"

#include <iostream>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda/std/tuple>
#include <random>
#include <chrono>

#define MAX_THREAD_PER_BLOCK 1024

__global__ void kernel_curand_init(curandState *state, const unsigned int seed)
{
    const unsigned int id = threadIdx.x + blockIdx.x * blockDim.x;

    curand_init(seed, id, 0, &state[id]);
}

// __global__ void kernel(const size_t size, float *z, curandState *state)
// {
//     __shared__ int local_array[32][32];

//     const int id = blockIdx.x * blockDim.x * blockDim.y + threadIdx.x + threadIdx.y * blockDim.x;

//     // printf("y : %d x : %d id : %d\n", threadIdx.y, threadIdx.x, id);

//     if (id < size)
//     {
//         local_array[threadIdx.x][threadIdx.y] = threadIdx.x + threadIdx.y * blockDim.x;

//         z[id] = local_array[threadIdx.x][threadIdx.y];

//         z[id] = curand(&(state[id]));
//     }
// }

__global__ void tabucol(
    curandState_t *state,
    bool *graph,
    unsigned int *solutions,
    const unsigned int max_iter,
    const unsigned int k,
    const unsigned int vertex_count)
{
    const unsigned int id = threadIdx.x + blockIdx.x * blockDim.x;

    unsigned int *solution = &solutions[id * vertex_count];

    // best solution
    unsigned int *best_solution = new unsigned int[vertex_count];
    // unsigned int *best_solution = (unsigned int *)malloc(vertex_count * sizeof(unsigned int));

    memcpy(best_solution, solution, vertex_count * sizeof(unsigned int));
    // for (size_t i{0}; i < vertex_count; ++i)
    // {
    //     best_solution[i] = solution[i];
    // }

    // The gamma_conflict matrix gives which color is connected to a vertex
    unsigned int *gamma_conflict = new unsigned int[vertex_count * k];
    // unsigned int *gamma_conflict = (unsigned int *)malloc(vertex_count * k * sizeof(unsigned int));
    unsigned int *tabu = new unsigned int[vertex_count * k];
    // unsigned int *tabu = (unsigned int *)malloc(vertex_count * k * sizeof(unsigned int));

    // printf("solution %u : %p\n", id, solution);
    // printf("best_solution %u : %p\n", id, best_solution);
    // printf("gamma_conflict %u : %p\n", id, gamma_conflict);
    // printf("tabu %d : %p\n", id, tabu);

    // init tabu and gamma
    memset(tabu, 0, vertex_count * k * sizeof(unsigned int));
    memset(gamma_conflict, 0, vertex_count * k * sizeof(unsigned int));
    // for (size_t i{0}; i < vertex_count * k; ++i)
    // {
    //     tabu[i] = 0;
    //     gamma_conflict[i] = 0;
    // }

    // init gamma
    int conflict_count{0};
    for (size_t y{0}; y < vertex_count; ++y)
    {
        for (size_t x{0}; x < y; ++x)
        {
            if (graph[x + y * vertex_count])
            {
                gamma_conflict[y * k + solution[x]] += 1;
                gamma_conflict[x * k + solution[y]] += 1;

                if (solution[x] == solution[y])
                {
                    ++conflict_count;
                }
            }
        }
    }

    int best_objective{conflict_count};

    // Main loop
    unsigned int i{0};
    while (i < max_iter && best_objective != 0)
    {
        size_t best_update_v{0};
        size_t best_update_c{0};
        int best_improvement = gamma_conflict[best_update_v * k + best_update_c] - gamma_conflict[best_update_v * k + solution[best_update_v]];

        // for each vertex we search for the best new color to reduce the count of conflicts
        for (size_t v{0}; v < vertex_count; ++v)
        {
            // if it has at least one conflict
            if (gamma_conflict[v * k + solution[v]] > 0)
            {
                // try each available color
                for (size_t c{0}; c < k; ++c)
                {
                    // if the color is different and not tabu with this vertex
                    if (c != solution[v] && !(tabu[v * k + c] > i))
                    {
                        // compute the improvement
                        int improvement = gamma_conflict[v * k + c] - gamma_conflict[v * k + solution[v]];
                        // update the best improvement found
                        if (improvement < best_improvement)
                        {
                            best_update_v = v;
                            best_update_c = c;
                            best_improvement = improvement;
                        }
                        else if (improvement == best_improvement /*&& curand_uniform(&(state[id])) < 0.6*/)
                        {
                            best_update_v = v;
                            best_update_c = c;
                            best_improvement = improvement;
                        }
                    }
                }
            }
        }

        // update the solution with the best improvement found
        conflict_count += best_improvement;
        const unsigned int former_color = solution[best_update_v];
        solution[best_update_v] = best_update_c;

        // update the gamma conflict matrix
        for (size_t n{0}; n < vertex_count; ++n)
        {
            if (graph[best_update_v * vertex_count + n])
            {
                gamma_conflict[n * k + former_color] -= 1;
                gamma_conflict[n * k + best_update_c] += 1;
            }
        }

        // update the tabu list
        tabu[best_update_v * k + former_color] = i + static_cast<unsigned int>(0.6 * conflict_count);

        // update the best solution found
        if (conflict_count < best_objective)
        {
            best_objective = conflict_count;
            memcpy(best_solution, solution, vertex_count * sizeof(unsigned int));
            // for (size_t v{0}; v < vertex_count; ++v)
            // {
            //     best_solution[v] = solution[v];
            // }
        }

        ++i;
    }

    // copy the best solution
    memcpy(solution, best_solution, vertex_count * sizeof(unsigned int));
    // for (size_t v{0}; v < vertex_count; ++v)
    // {
    //     solution[v] = best_solution[v];
    // }

    delete[] best_solution;
    // free(best_solution);
    delete[] gamma_conflict;
    // free(gamma_conflict);
    delete[] tabu;
    // free(tabu);
}

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
    const size_t pop_size{64};

    size_t cudaDeviceLimit{1024 * 1024 * 1024};
    cudaDeviceSetLimit(cudaLimitMallocHeapSize, cudaDeviceLimit);
    cudaDeviceGetLimit(&cudaDeviceLimit, cudaLimitMallocHeapSize);
    std::cout << "cudaDeviceGetLimit : " << cudaDeviceLimit / 1024 / 1024 << "MB\n";

    // random init
    curandState_t *d_state;
    cudaMalloc(&d_state, pop_size * sizeof(curandState));

    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    // seed = 0;
    int numBlock = pop_size / 64;
    int threadsPerBlock = pop_size / numBlock;

    kernel_curand_init<<<numBlock, threadsPerBlock, 0, 0>>>(d_state, seed);
    std::string error = cudaGetErrorString(cudaDeviceSynchronize());
    std::cout << error << '\n';

    // Graph
    DIMACS instance("graph/dimacs_instances/DSJC1000.5.col");

    bool *d_graph = nullptr;
    cudaMalloc(&d_graph, instance.getVertexCount() * instance.getVertexCount() * sizeof(bool));
    std::cout << "d_graph : " << d_graph << '\n';
    cudaMemcpy(
        d_graph,
        &instance.getGraph()(0, 0),
        instance.getVertexCount() * instance.getVertexCount() * sizeof(bool),
        cudaMemcpyHostToDevice);

    // Solutions
    const unsigned int targetColorCount{90};
    std::mt19937_64 generator(seed);
    std::uniform_int_distribution<unsigned int> distribution(0, targetColorCount - 1);

    std::valarray<unsigned int> solutions(instance.getVertexCount() * pop_size);
    std::generate(begin(solutions), end(solutions), [&]
                  { return distribution(generator); });

    std::valarray<unsigned int> solutions_objective(pop_size);
    for (size_t i{0}; i < pop_size; ++i)
    {
        solutions_objective[i] = conflict(instance.getGraph(), &solutions[i * instance.getVertexCount()], instance.getVertexCount());
    }

    unsigned int *d_solutions = nullptr;
    cudaMalloc(&d_solutions, solutions.size() * sizeof(unsigned int));
    std::cout << "d_solutions : " << d_solutions << '\n';
    cudaMemcpy(
        d_solutions,
        &solutions[0],
        solutions.size() * sizeof(unsigned int),
        cudaMemcpyHostToDevice);

    // Start Kernel
    tabucol<<<numBlock, threadsPerBlock, 0, 0>>>(
        d_state,
        d_graph,
        d_solutions,
        500,
        targetColorCount,
        instance.getVertexCount());

    error = cudaGetErrorString(cudaDeviceSynchronize());
    std::cout << error << '\n';
    error = cudaGetErrorString(cudaPeekAtLastError());
    std::cout << error << '\n';
    error = cudaGetErrorString(cudaThreadSynchronize());
    std::cout << error << '\n';

    cudaMemcpy(
        &solutions[0],
        d_solutions,
        solutions.size() * sizeof(unsigned int),
        cudaMemcpyDeviceToHost);

    for (size_t s{0}; s < pop_size; ++s)
    {
        std::cout << "conflicts : " << solutions_objective[s] << " -> " << conflict(instance.getGraph(), &solutions[s * instance.getVertexCount()], instance.getVertexCount()) << " : ";
        //for (size_t v{0}; v < instance.getVertexCount(); ++v)
        //{
        //    std::cout << solutions[s * instance.getVertexCount() + v] << ' ';
        //}
        std::cout << '\n';        
    }

    cudaFree(d_state);
    cudaFree(d_graph);
    cudaFree(d_solutions);

    std::cout << std::endl;

    return 0;
}
