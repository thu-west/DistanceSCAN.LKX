#ifndef WEIGHTED_SCAN_QUERY_H
#define WEIGHTED_SCAN_QUERY_H

#include "algos.h"

void save_clusters(vector<int> cluster_result) {
    string result_path =
            config.graph_location + "/" + config.algo + "d_" + to_str(config.distance) + "mu_" + to_str(config.mu) +
            "e" + to_str(config.epsilon) + ".csv";
    ofstream fout(result_path);
    fout << "node_id,cluster_id" << endl;
    for (int i = 0; i < cluster_result.size(); ++i) {
        fout << i << "," << cluster_result[i] << endl;
    }
    fout.close();
}

void set_result(int used_counter) {
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - result.startTime).count();
    result.whole_time_usage = duration / TIMES_PER_SEC;;
    result.total_mem_usage = get_proc_memory() / 1000.0;
    result.total_time_usage = Timer::used(used_counter);

}

void dataOutput() {
    //write
    Timer timer(123);
    char record[100];
    FILE *fs;
    string f;
    if (config.operation != QUERY) {
        f = config.operation + "_result.txt";
    } else {
        f = config.algo + "_result.txt";
    }
    fs = fopen(f.c_str(), "a");
    char fsch;
    if (fs == NULL) {
        printf("ERROR!");
        exit(1);
    } else {
        sprintf(record, "%s", config.graph_alias.c_str());
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);
        sprintf(record, "algo:%s", config.algo.c_str());
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);
        sprintf(record, "epsilon: %f", config.epsilon);
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);
        sprintf(record, "mu: %d", config.mu);
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);
        sprintf(record, "distance: %f", config.distance);
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);
        sprintf(record, "max_distance: %f", config.max_distance);
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);
        sprintf(record, "hash_k: %d", config.hash_k);
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);
        sprintf(record, "whole_time_usage: %.2lfs", result.whole_time_usage);
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);
        sprintf(record, "total_mem_usage: %.2lfMB", result.total_mem_usage);
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);
        sprintf(record, "total_time_usage: %.2lfs", result.total_time_usage);
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);
        if (config.operation == CLUSTER_VALIDATION) {
            sprintf(record, "type: %d", config.type);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "similarityType: %d", config.similarityType);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "m_scan: %f", result.m_scan);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "m_wscan: %f", result.m_wscan);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "m_dis: %f", result.m_dis);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "m_my_ads: %f", result.m_my_ads);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            /*
            sprintf(record, "qs_scan: %f", result.qs_scan);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "qs_dis: %f", result.qs_dis);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);*/
        } else if (config.operation == EXPONLFR) {
            sprintf(record, "miss_labelled_rate: %f", result.miss_labelled_rate);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "ari: %f", result.ari);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
        } else if (config.operation == QUALITY_VALIDATION) {
            sprintf(record, "miss_labelled_rate: %f", result.miss_labelled_rate);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "ari: %f", result.ari);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "topk_cluster_quality_1_min: %f", result.topk_cluster_quality_1_min);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "topk_cluster_quality_1_avg: %f", result.topk_cluster_quality_1_avg);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "topk_cluster_quality_5_min: %f", result.topk_cluster_quality_5_min);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "topk_cluster_quality_5_avg: %f", result.topk_cluster_quality_5_avg);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "topk_cluster_quality_10_min: %f", result.topk_cluster_quality_10_min);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "topk_cluster_quality_10_avg: %f", result.topk_cluster_quality_10_avg);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "topk_cluster_quality_50_min: %f", result.topk_cluster_quality_50_min);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "topk_cluster_quality_50_avg: %f", result.topk_cluster_quality_50_avg);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "topk_cluster_quality_100_min: %f", result.topk_cluster_quality_100_min);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "topk_cluster_quality_100_avg: %f", result.topk_cluster_quality_100_avg);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
        } else if (config.operation == QUERY) {
            sprintf(record, "algo_type: %d", config.algo_type);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "sketches_optimize: %d", config.sketches_optimize);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "algtime: %.6lf", timer.timeUsed[10] / TIMES_PER_SEC);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            if (config.algo == MY_ADS) {
                sprintf(record, "bin: %f", config.bin);
                fwrite(record, sizeof(*record), strlen(record), fs);
                fsch = putc('\t', fs);
                sprintf(record, "disscan: %.6lf", timer.timeUsed[9] / TIMES_PER_SEC);
                fwrite(record, sizeof(*record), strlen(record), fs);
                fsch = putc('\t', fs);
                sprintf(record, "sketch: %.6lf", timer.timeUsed[4] / TIMES_PER_SEC);
                fwrite(record, sizeof(*record), strlen(record), fs);
                fsch = putc('\t', fs);
                sprintf(record, "get_botk_time: %.6lf", timer.timeUsed[6] / TIMES_PER_SEC);
                fwrite(record, sizeof(*record), strlen(record), fs);
                fsch = putc('\t', fs);
                if (config.algo_type % 5 == 0) {
                    sprintf(record, "traverse_edge_time: %.6lf", timer.timeUsed[31] / TIMES_PER_SEC);
                    fwrite(record, sizeof(*record), strlen(record), fs);
                    fsch = putc('\t', fs);
                } else {
                    sprintf(record, "traverse_edge_sort_time: %.6lf", timer.timeUsed[32] / TIMES_PER_SEC);
                    fwrite(record, sizeof(*record), strlen(record), fs);
                    fsch = putc('\t', fs);
                }
                if (config.algo_type % 5 == 1 || config.algo_type % 5 == 2) {
                    sprintf(record, "traverse_node_time: %.6lf", timer.timeUsed[20] / TIMES_PER_SEC);
                    fwrite(record, sizeof(*record), strlen(record), fs);
                    fsch = putc('\t', fs);
                } else {
                    sprintf(record, "traverse_node_sort_time: %.6lf", timer.timeUsed[30] / TIMES_PER_SEC);
                    fwrite(record, sizeof(*record), strlen(record), fs);
                    fsch = putc('\t', fs);
                }
                if (config.algo_type < 6) {
                    sprintf(record, "traverse_node_time_botk: %.6lf", timer.timeUsed[22] / TIMES_PER_SEC);
                    fwrite(record, sizeof(*record), strlen(record), fs);
                    fsch = putc('\t', fs);
                } else if (config.algo_type < 11) {
                    sprintf(record, "traverse_node_time_botk_sort: %.6lf", timer.timeUsed[23] / TIMES_PER_SEC);
                    fwrite(record, sizeof(*record), strlen(record), fs);
                    fsch = putc('\t', fs);
                }
                sprintf(record, "traverse_node_time_dijkstra: %.6lf", timer.timeUsed[21] / TIMES_PER_SEC);
                fwrite(record, sizeof(*record), strlen(record), fs);
                fsch = putc('\t', fs);
            } else {
                sprintf(record, "scan: %.6lf", timer.timeUsed[1] / TIMES_PER_SEC);
                fwrite(record, sizeof(*record), strlen(record), fs);
                fsch = putc('\t', fs);
                sprintf(record, "basic: %.6lf", timer.timeUsed[2] / TIMES_PER_SEC);
                fwrite(record, sizeof(*record), strlen(record), fs);
                fsch = putc('\t', fs);
                sprintf(record, "pscan: %.6lf", timer.timeUsed[3] / TIMES_PER_SEC);
                fwrite(record, sizeof(*record), strlen(record), fs);
                fsch = putc('\t', fs);
                sprintf(record, "d_neighbor: %.6lf", timer.timeUsed[7] / TIMES_PER_SEC);
                fwrite(record, sizeof(*record), strlen(record), fs);
                fsch = putc('\t', fs);
            }
            sprintf(record, "read_graph: %.6lf", timer.timeUsed[5] / TIMES_PER_SEC);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);

        } else if (config.operation == CONSTRUCT_SKETCHES) {
            sprintf(record, "sketches: %.6lf", timer.timeUsed[4] / TIMES_PER_SEC);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "construct_sketch_time: %.6lf", timer.timeUsed[13] / TIMES_PER_SEC);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "compute_all_botk_time: %.6lf", timer.timeUsed[8] / TIMES_PER_SEC);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "save_all_botk_time: %.6lf", timer.timeUsed[11] / TIMES_PER_SEC);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "update_histogram_time: %.6lf", timer.timeUsed[12] / TIMES_PER_SEC);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
        } else if (config.operation == GRAPH_MAINTAIN) {
            sprintf(record, "update_edge_nums: %d", config.update_edge_nums);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "update_mode: %d", config.update_mode);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "edge_del_time: %.6lf", timer.timeUsed[14] / TIMES_PER_SEC);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "edge_ins_time: %.6lf", timer.timeUsed[15] / TIMES_PER_SEC);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "edge_update_time: %.6lf", timer.timeUsed[16] / TIMES_PER_SEC);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
            sprintf(record, "maintain_time: %.6lf", timer.timeUsed[17] / TIMES_PER_SEC);
            fwrite(record, sizeof(*record), strlen(record), fs);
            fsch = putc('\t', fs);
        }


        fsch = putc('\n', fs);
        fsch = putc('\n', fs);
    }
}

vector<vector<int>> query() {
    result.startTime = std::chrono::steady_clock::now();
    INFO(config.distance, config.epsilon, config.mu);
    int used_counter = 0;
    vector<vector<int>> cluster_result;
    cluster_result.reserve(graph.n);
    if (config.algo == SCAN || config.algo == BASIC || config.algo == BOTK_SCAN || config.algo == PSCAN_DIS ||
        config.algo == W_SCAN || config.algo == MY_ADS) {
        used_counter = SCAN_TIME;
        ClusterAlgos clusterAlgos = ClusterAlgos();
        cluster_result = clusterAlgos.get_clusters();
    } else if (config.algo == CHECK) {
        config.algo = BASIC;
        ClusterAlgos clusterAlgos = ClusterAlgos();
        vector<int> clu_res = graph.clusterID;
        graph.clusterID = vector<int>(graph.n, -1);
        graph.is_core = vector<int>(graph.n, -1);
        config.algo = BOTK_SCAN;
        ClusterAlgos clusterAlgos1 = ClusterAlgos();
        vector<int> clu_res1 = graph.clusterID;
        for (int i = 0; i < graph.n; ++i) {
            assert(clu_res[i] == clu_res1[i]);
        }
    }
    if (config.operation == QUERY) {
        set_result(WHOLE);
        dataOutput();
    }
    serialize_clusters(cluster_result);
    // save_clusters(cluster_result);
    return cluster_result;
}


void validate_modularity() {
    result.startTime = std::chrono::steady_clock::now();
    Timer timer(999);

    config.algo = SCAN;
    {
        Graph graph0;
        graph0.init(config.graph_location);
        swap(graph, graph0);
    }

    vector<vector<int>> clusters_scan = query();
    //vector<int> scan_cluster = graph.clusterID;


    config.algo = W_SCAN;
    vector<vector<int>> clusters_w_scan;
    //vector<int> w_scan_cluster;
    if (config.similarityType == config.cos) {
        Graph graph1;
        graph1.init(config.graph_location);
        swap(graph, graph1);
        clusters_w_scan = query();
        //w_scan_cluster = graph.clusterID;
    }


    config.algo = MY_ADS;
    vector<vector<int>> clusters_my_ads;
    if(config.similarityType == config.jac){
        Graph graph3;
        graph3.init(config.graph_location);
        swap(graph, graph3);
        clusters_my_ads = query();
    }


    vector<vector<int>> clusters_dis;
    if (//graph.data_folder.find("com-lj") == graph.data_folder.npos &&
        graph.data_folder.find("uk2002") == graph.data_folder.npos) {
        config.algo = PSCAN_DIS;
        Graph graph2;
        graph2.init(config.graph_location);
        swap(graph, graph2);

        clusters_dis = query();
    }


    INFO("computing mudularity");
    {
        Graph graph_tmp;
        graph_tmp.init(config.graph_location);
        swap(graph, graph_tmp);
    }
    graph.weight_degree = vector<double>(graph.n, 0);
    graph.whole_weight = 0;
    for (int i = 0; i < graph.n; ++i) {
        for (int nei: graph.adj_list[i]) {
            double real_weight = 0;
            if (graph.data_folder.find("graphs_coauthors") != graph.data_folder.npos) {
                real_weight = exp(1 / graph.edge_weight[i][nei] - 1);
            } else {
                real_weight = (1 - graph.edge_weight[i][nei]) / 0.9;
            }
            graph.weight_degree[i] += real_weight;
        }
        graph.whole_weight += graph.weight_degree[i];
    }
    compute_modulairty(clusters_scan, result.m_scan);
    if (config.similarityType == config.cos) {
        compute_modulairty(clusters_w_scan, result.m_wscan);
    }
    if(config.similarityType == config.jac){
        compute_modulairty(clusters_my_ads, result.m_my_ads);
    }
    if (//graph.data_folder.find("com-lj") == graph.data_folder.npos &&
        graph.data_folder.find("uk2002") == graph.data_folder.npos) {
        compute_modulairty(clusters_dis, result.m_dis);
    }

    set_result(999);
    dataOutput();

}

void construct_sketches() {
    result.startTime = std::chrono::steady_clock::now();
    SKETCHES sketches = SKETCHES();
    set_result(SKETCH_TIME);
    dataOutput();
}

void validate_quality(vector<int> &ground_truth, vector<int> &test, vector<vector<int>> &clus_gt,
                      vector<vector<int>> &clus_test) {
    result.startTime = std::chrono::steady_clock::now();
    Timer timer(999);
    // result.miss_labelled_rate = get_mis_labeled_rate(ground_truth, test);

    result.ari = get_ari(ground_truth, test, clus_gt.size(), clus_test.size());
    set_result(999);
    dataOutput();
}


#endif //WEIGHTED_SCAN_QUERY_H
