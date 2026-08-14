#include <omp.h>

#include "vsf.h"

using namespace std;

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
            double structure = mean_val;
            double structure_uncertainty = std_val / (sqrt(N - 1));  // sample standard error

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

/**
 * \brief Calculates the nth order structure function of two-dimensional data.
 * \param input_array The input as a two-dimensional vector.
 * \param order The order of the structure function to compute. For example, order=1 will only output the average
 * difference between pairs of points (normalized by the variance) as a function of their distance.
 * \param log_bin_width The width of the logarithmic bins for regrouping distances. If set to 0 or negative, no
 * logarithmic binning is applied.
 * \param bin_start The starting point for the logarithmic bins. This is the lower bound of the first bin. If set to 0
 * or negative, no logarithmic binning is applied.
 */
vector_2d structure_function(
    const vector_2d& input_array,
    const int order,
    const double log_bin_width,
    const double bin_start
) {
    if (log_bin_width <= 0 || bin_start <= 0) {
        return structure_function(input_array, order);
    }

    // Compute the differences between each pair of elements along with their distances
    vector<array<double, 2>> single_dists_and_vals_1d = subtract_pairs(input_array);

    // Regroup the values by their pair separation distances
    double_unordered_map regrouped_vals;
    regroup_distance_thread_local(single_dists_and_vals_1d, regrouped_vals);

    vector_2d output_array;
    output_array.reserve(regrouped_vals.size());
    double variance_val = variance(input_array);

    // Create bin groups
    double_unordered_map bin_groups;
    for (const auto& [dist, vals] : regrouped_vals) {
        if (dist == 0) continue;  // reject zero distances
        if (dist < bin_start) {
            bin_groups[dist].assign(vals.begin(), vals.end());
        } else {  // Apply logarithmic binning
            double bin = bin_start
                * pow(10, (floor(log10(dist / bin_start) / log_bin_width) * log_bin_width) + (log_bin_width / 2));
            bin_groups[bin].insert(bin_groups[bin].end(), vals.begin(), vals.end());
        }
    }

    // Calculate the structure function for each bin group in parallel
    vector<vector_2d> thread_local_results(omp_get_max_threads());

    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        vector_2d& local_output = thread_local_results[thread_id];

        #pragma omp for
        for (int i = 0; i < bin_groups.size(); ++i) {
            auto it = next(bin_groups.begin(), i);  // access ith element
            const auto& [dist, vals] = *it;

            vector<double> pow_values = pow(vals, (double)order);
            int N = pow_values.size();
            if (N == 1) continue;  // skip if there is only one value

            double mean_val = mean(pow_values);
            double std_val = standard_deviation(pow_values);
            double structure = mean_val;
            double structure_uncertainty = std_val / (sqrt(N - 1));  // sample standard error

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
