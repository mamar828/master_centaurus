#include <omp.h>

#include "vsf.h"
// #include "stats.h"
// #include "time.h"

using namespace std;

// Time current_time;

/**
 * \brief Calculates the nth order structure function of two-dimensional data.
 * \param input_array The input as a two-dimensional vector.
 * \param order The order of the structure function to compute. For example, order=1 will only output the average
 * difference between pairs of points (normalized by the variance) as a function of their distance.
 */
vector_2d structure_function(const vector_2d& input_array, const int order) {
    // Compute the differences between each pair of elements along with their distances
    vector<array<double, 2>> single_dists_and_vals_1d = subtract_pairs(input_array);

    // Regroup the values by their pair separation distances
    double_unordered_map regrouped_vals;
    regroup_distance_thread_local(single_dists_and_vals_1d, regrouped_vals);

    vector_2d output_array;
    output_array.reserve(regrouped_vals.size());
    double variance_val = variance(input_array);

    // Thread-local storage for results
    vector<vector_2d> thread_local_results(omp_get_max_threads());

    // Compute the structure function for each pair separation in parallel
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        vector_2d& local_output = thread_local_results[thread_id];

        #pragma omp for
        for (int i = 0; i < regrouped_vals.size(); ++i) {
            auto it = next(regrouped_vals.begin(), i);  // access ith element
            const auto& [dist, vals] = *it;
            if (dist == 0) continue;  // reject zero distances

            vector<double> pow_values = pow(vals, (double)order);
            int N = pow_values.size();
            if (N == 1) continue;  // skip if there is only one value

            double mean_val = mean(pow_values);
            double std_val = standard_deviation(pow_values);
            double structure = mean_val / variance_val;
            double structure_uncertainty = std_val / (variance_val * sqrt(N - 1));  // sample standard error

            // Store result in thread-local buffer
            local_output.push_back({dist, structure, structure_uncertainty});
        }
    }

    // Combine results from all threads
    for (const auto& local_result : thread_local_results) {
        output_array.insert(output_array.end(), local_result.begin(), local_result.end());
    }

    return output_array;
}
