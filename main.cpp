#include <iostream>
#include "./include/query.h"

int main(int argc, char *argv[]) {
    ios::sync_with_stdio(false);
    program_start(argc, argv);
    Saver::init();
    srand(time(nullptr));
    config.init_config(argc, argv);
    INFO(config.operation);
    if (config.operation == CONSTRUCT_SKETCHES) {
        graph.init(config.graph_location);
        construct_sketches();
    } else if (config.operation == QUERY) {
        vector<string> possibleAlgo = {SCAN, BASIC, BOTK_SCAN, CHECK, PSCAN_DIS, W_SCAN, MY_ADS};
        if (find(possibleAlgo.begin(), possibleAlgo.end(), config.algo) == possibleAlgo.end()) {
            INFO("Wrong algo param: ", config.algo);
            exit(1);
        }
        INFO(config.algo);

        graph.init(config.graph_location);

        query();
    } else if (config.operation == CLUSTER_VALIDATION) {
        graph.init(config.graph_location);
        validate_modularity();
    } else if (config.operation == QUALITY_VALIDATION) {
        bool flag = false;
        if(config.algo == MY_ADS){
            flag = true;
        }
        config.algo = BASIC;
        vector<vector<int>> basic_clusters, pscan_clusters, botk_clusters, my_ads_clusters;
        vector<int> basic_clusterid, basic_iscore, pscan_clusterid, pscan_iscore, botk_clusterid, botk_iscore, my_ads_clusterid, my_ads_iscore;
        deserialize_clusters(basic_clusters, basic_clusterid, basic_iscore);
        config.algo = PSCAN_DIS;
        deserialize_clusters(pscan_clusters, pscan_clusterid, pscan_iscore);
        if(flag){
            config.algo = MY_ADS;
            deserialize_clusters(my_ads_clusters,my_ads_clusterid,my_ads_iscore);
        }else{
            config.algo = BOTK_SCAN;
            deserialize_clusters(botk_clusters, botk_clusterid, botk_iscore);
        }

        INFO(basic_clusters.size(),pscan_clusters.size(),botk_clusters.size(),my_ads_clusters.size());
        if (basic_clusters.empty()) {
            if (pscan_clusters.empty() || botk_clusters.empty()) {
                INFO("No files.");
            }
            if(flag){
                config.algo = PSCAN_DIS + "_" + MY_ADS + to_string(config.algo_type);
                validate_quality(pscan_clusterid,my_ads_clusterid,pscan_clusters,my_ads_clusters);
            }else{
                config.algo = PSCAN_DIS + "_" + BOTK_SCAN;
                validate_quality(pscan_clusterid, botk_clusterid, pscan_clusters, botk_clusters);
            }
        } else {
            config.algo = BASIC + "_" + PSCAN_DIS;
            validate_quality(basic_clusterid, pscan_clusterid, basic_clusters, pscan_clusters);
            if(flag){
                config.algo = BASIC + "_" + MY_ADS + to_string(config.algo_type);
                validate_quality(basic_clusterid, my_ads_clusterid, basic_clusters, my_ads_clusters);
            }else{
                config.algo = BASIC + "_" + BOTK_SCAN;
                validate_quality(basic_clusterid, botk_clusterid, basic_clusters, botk_clusters);
            }
        }
    } else if(config.operation == CONVERT_GRAPH){
        graph.init(config.graph_location);
        graph.get_undirected_similarity_weighted();
    }

    Timer::show();
    if (config.operation == QUERY) {
        Counter::show();
        auto args = combine_args(argc, argv);
        Saver::save_json(config, result, args);
    }
    return 0;
}
