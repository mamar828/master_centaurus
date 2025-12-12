#include <omp.h>

#include "tools.h"

using namespace std;

/**
 * \brief Regroups a vector of distance and value pairs into an unordered map containing each unique distance as key and
 * a vector of corresponding values.
 * \param single_dists_and_vals_1d The vector of (distance, value) pairs to regroup.
 * \param regrouped_vals The unordered map in which to store the regrouped values.
 */
void regroup_distance_thread_local(
    const vector<array<double, 2>>& single_dists_and_vals_1d,
    double_unordered_map& regrouped_vals
) {
    // Create thread-local storage for the unordered_map
    vector<double_unordered_map> local_regroup_vals(omp_get_max_threads());

    #pragma omp parallel for
    for (size_t i = 0; i < single_dists_and_vals_1d.size(); ++i) {
        int thread_id = omp_get_thread_num();
        const auto& dist_and_val = single_dists_and_vals_1d[i];
        local_regroup_vals[thread_id][dist_and_val[0]].push_back(dist_and_val[1]);
    }

    // Merge results from all threads into the final unordered_map
    for (const auto& local_map : local_regroup_vals) {
        for (const auto& pair : local_map) {
            regrouped_vals[pair.first].insert(regrouped_vals[pair.first].end(),
                                              pair.second.begin(), pair.second.end());
        }
    }
}

/**
 * \brief Regroups a vector of (dist_x, dist_y) and value pairs into an unordered map containing each unique (dist_x,
 * dist_y) as key and a vector of corresponding values.
 * \param single_dists_and_vals_2d The vector of (dist_x, dist_y, value) groups to cluster.
 * \param regrouped_vals The unordered map in which to store the regrouped values.
 */
void regroup_distance_thread_local(
    const vector<array<double, 3>>& single_dists_and_vals_2d,
    array_unordered_map& regrouped_vals
) {
    // Create thread-local storage for the unordered_map
    vector<array_unordered_map> local_regroup_vals(omp_get_max_threads());

    #pragma omp parallel for
    for (size_t i = 0; i < single_dists_and_vals_2d.size(); ++i) {
        int thread_id = omp_get_thread_num();
        const auto& dist_and_val = single_dists_and_vals_2d[i];
        local_regroup_vals[thread_id][{dist_and_val[0], dist_and_val[1]}].push_back(dist_and_val[2]);
    }

    // Merge results from all threads into the final unordered_map
    for (const auto& local_map : local_regroup_vals) {
        for (const auto& pair : local_map) {
            regrouped_vals[pair.first].insert(regrouped_vals[pair.first].end(),
                                              pair.second.begin(), pair.second.end());
        }
    }
}

/**
 * \brief Appends the contents of one vector of arrays of two values to another.
 * \param dest The destination vector to which data will be appended.
 * \param src The source vector from which data will be taken.
 */
void combine_vectors(vector<array<double,2>>& dest, const vector<array<double,2>>& src) {
    #pragma omp critical
    dest.insert(dest.end(), src.begin(), src.end());
}

/**
 * \brief Appends the contents of one vector of arrays of three values to another.
 * \param dest The destination vector to which data will be appended.
 * \param src The source vector from which data will be taken.
 */
void combine_vectors(vector<array<double,3>>& dest, const vector<array<double,3>>& src) {
    #pragma omp critical
    dest.insert(dest.end(), src.begin(), src.end());
}

/**
 * \brief Applies an operation between each values of a vector_2d and computes the corresponding distance between each
 * pair of points.
 * \param input_array The vector_2d on which the function must be applied for every pair of points.
 * \param function Callable of two doubles that returns a float. This function is applied to each pair of elements
 * in the input_array.
 * \return Vector of arrays of two elements: the distance between the two points and the result of the function.
 */
template <typename T>
vector<array<double, 2>> apply_vector_map(const vector_2d& input_array, const T& function) {
    const size_t height = input_array.size();
    const size_t width = input_array[0].size();
    vector<array<double, 2>> single_dists_and_vals;

    size_t max_possible_size = (height * width * (height * width)) / 2;

    #pragma omp parallel
    {
        vector<array<double, 2>> thread_single_dists_and_vals;
        // Reserve an approximate size to avoid multiple allocations
        thread_single_dists_and_vals.reserve(max_possible_size / omp_get_num_threads());

        #pragma omp for collapse(2) schedule(dynamic)
        for (size_t y = 0; y < height; ++y) {
            for (size_t x = 0; x < width; ++x) {
                if (isnan(input_array[y][x])) continue;
                for (size_t j = y; j < height; ++j) {
                    for (size_t i = (j == y) ? x : 0; i < width; ++i) {  // lag=0 is considered here
                        if (isnan(input_array[j][i])) continue;
                        double dist = sqrt((i - x) * (i - x) + (j - y) * (j - y));
                        double val = function(input_array[y][x], input_array[j][i]);
                        thread_single_dists_and_vals.push_back({dist, val});
                    }
                }
            }
        }

        // Combine the thread-local results into the global vector
        combine_vectors(single_dists_and_vals, thread_single_dists_and_vals);
    }

    // Optionally shrink to fit if memory usage is a concern
    single_dists_and_vals.shrink_to_fit();

    return single_dists_and_vals;
}

/**
 * \brief Computes the product between each pair of elements in the input array, along with their distances.
 */
vector<array<double, 2>> multiply_pairs(const vector_2d& input_array) {
    return apply_vector_map(input_array, [](double a, double b) {return a * b;});
}

/**
 * \brief Computes the absolute difference between each pair of elements in the input array, along with their distances.
 */
vector<array<double, 2>> subtract_pairs(const vector_2d& input_array) {
    return apply_vector_map(input_array, [](double a, double b) {return abs(a - b);});
}
