#ifndef WEIGHTED_SCAN_ALGOS_H
#define WEIGHTED_SCAN_ALGOS_H

#include "clusteralgo.h"

void compute_modulairty(const vector<vector<int>> &clusters, double &m) {
    m = 0;
    for (auto clu: clusters) {
        for (int vertex_i: clu) {
            for (int vertex_j: clu) {
                if (vertex_i == vertex_j)continue;
                double real_weight = 0;
                if (graph.edge_weight[vertex_i].find(vertex_j) != graph.edge_weight[vertex_i].end()) {
                    if (graph.data_folder.find("graphs_coauthors") != graph.data_folder.npos) {
                        real_weight = exp(1 / graph.edge_weight[vertex_i].at(vertex_j) - 1);
                    } else {
                        real_weight = (1 - graph.edge_weight[vertex_i].at(vertex_j)) / 0.9;
                    }
                }
                m += real_weight - graph.weight_degree[vertex_i] * graph.weight_degree[vertex_j] / graph.whole_weight;
            }
        }
    }
    m /= graph.whole_weight;
}


double compute_qs(const vector<vector<int>> &clusters, const vector<int> &cluster_id, const double &ts) {
    double qs = 0;
    for (vector<int> clu: clusters) {
        double is = 0, ds = 0;
        for (int i = 0; i < clu.size(); ++i) {
            int s = clu[i];
            for (int nei: graph.d_neighbors[s]) {
                if (cluster_id[s] == cluster_id[nei]) {
                    is += graph.get_jac_res(s, nei);
                }
                ds += graph.get_jac_res(s, nei);
            }
        }
        qs += is / ts - (ds * ds / ts / ts / ts);
    }
    return qs;
}

double get_mis_labeled_rate(const vector<int> &ground_truth, const vector<int> &test_result) {
    int miss_counter = 0;
    for (int i = 0; i < ground_truth.size(); ++i) {
        long long res = ground_truth[i];
        res *= test_result[i];
        if (res < 0)
            miss_counter++;
    }
    return double(miss_counter) / ground_truth.size();
}

inline long long c_n_2(long long n) {
    return n * (n - 1) / 2;
}


double get_ari(const vector<int> &ground_truth, const vector<int> &test_result, int clu_num_gt, int clu_num_test) {
    vector<unordered_map<int, int>> contingency_table = vector<unordered_map<int, int>>(clu_num_gt + 1,
                                                                                        unordered_map<int, int>{});
    vector<unordered_map<int, int>> contingency_table2 = vector<unordered_map<int, int>>(clu_num_test + 1,
                                                                                         unordered_map<int, int>{});

    for (int i = 0; i < ground_truth.size(); ++i) {
        int gt_id = ground_truth[i], tr_id = test_result[i];
        if (gt_id < 0) gt_id = clu_num_gt;
        if (tr_id < 0) tr_id = clu_num_test;
        contingency_table[gt_id][tr_id]++;
        contingency_table2[tr_id][gt_id]++;
    }
    vector<int> a = vector<int>(clu_num_gt + 1), b = vector<int>(clu_num_test + 1);
    long long sum_nij = 0, sum_a = 0, sum_b = 0;
    for (int i = 0; i < clu_num_gt + 1; ++i) {
        for (auto item: contingency_table[i]) {
            a[i] += item.second;
            sum_nij += c_n_2(item.second);
        }
        sum_a += c_n_2(a[i]);
    }
    for (int i = 0; i < clu_num_test + 1; ++i) {
        for (auto item: contingency_table2[i]) {
            b[i] += item.second;
        }
        sum_b += c_n_2(b[i]);
    }
    assert(sum_a > 0);
    assert(sum_b > 0);
    assert(sum_nij > 0);
    long double eri = sum_a / double(c_n_2(ground_truth.size())) * sum_b;
    long double numerator = sum_nij - eri;
    long double denominator = (sum_a + sum_b) / 2.0 - eri;
    //INFO(eri, numerator, denominator);
    return numerator / denominator;
}

void serialize_clusters(const vector<vector<int>> &cluster_result) {
    string prefix;
    if (config.algo == BOTK_SCAN) {
        prefix = "_algo_" + config.algo + "_d_" + to_str(config.distance) + "_k_" + to_string(config.hash_k)
                 + "_e_" + to_string(config.epsilon) + "_u_" + to_string(config.mu);
    } else if (config.algo == MY_ADS) {
        prefix =
                "_algo_" + config.algo + "_d_" + to_str(config.distance) + "_algo_type_" + to_string(config.algo_type) +
                "_k_" + to_string(config.hash_k) + "_e_" + to_string(config.epsilon) + "_u_" + to_string(config.mu);
    } else {
        prefix = "_algo_" + config.algo + "_d_" + to_str(config.distance)
                 + "_e_" + to_string(config.epsilon) + "_u_" + to_string(config.mu);
    }
    string idx = config.graph_location + prefix + "_cluster_result.txt";
    std::ofstream info_ofs(idx);
    boost::archive::text_oarchive info_oa(info_ofs);
    info_oa << cluster_result << graph.clusterID << graph.is_core;
}

void deserialize_clusters(vector<vector<int>> &cluster_result, vector<int> &clusterid, vector<int> &iscore) {
    string prefix;
    if (config.algo == BOTK_SCAN) {
        prefix = "_algo_" + config.algo + "_d_" + to_str(config.distance) + "_k_" + to_string(config.hash_k)
                 + "_e_" + to_string(config.epsilon) + "_u_" + to_string(config.mu);
    } else if (config.algo == MY_ADS) {
        prefix =
                "_algo_" + config.algo + "_d_" + to_str(config.distance) + "_algo_type_" + to_string(config.algo_type) +
                "_k_" + to_string(config.hash_k) + "_e_" + to_string(config.epsilon) + "_u_" + to_string(config.mu);
    } else {
        prefix = "_algo_" + config.algo + "_d_" + to_str(config.distance)
                 + "_e_" + to_string(config.epsilon) + "_u_" + to_string(config.mu);
    }
    string idx = config.graph_location + prefix + "_cluster_result.txt";
    if (exists_test(idx)) {
        std::ifstream info_ifs(idx);
        boost::archive::text_iarchive info_ia(info_ifs);
        info_ia >> cluster_result >> clusterid >> iscore;
    }
}

#endif //WEIGHTED_SCAN_ALGOS_H
