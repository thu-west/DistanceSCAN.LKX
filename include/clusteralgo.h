
#ifndef WEIGHTED_SCAN_CLUSTERALGO_H
#define WEIGHTED_SCAN_CLUSTERALGO_H

#include "sketches.h"
//#include "graph.h"

class ClusterAlgos {


    // Graph graph;

    vector<vector<int>> clusters;
    vector<int> hub;
    vector<int> outlier;

    vector<int> lb,ub;
    SKETCHES sketches;




public:

    ClusterAlgos() : sketches() {
        Timer timer(ALG_TIME);
        if (config.algo == BOTK_SCAN||config.algo == MY_ADS) {
            pruned_scan();
        } else if (config.algo == BASIC || config.algo == SCAN || config.algo == W_SCAN) {
            scan_framework();
        } else if (config.algo == PSCAN_DIS){
            pscan_dis();
        }
        assign_noncore_hubs_outliers();
    }
    vector<vector<int>> get_clusters(){return clusters;}

private:

    bool test_edge_jaccard(int u, int v);

    bool check_core(int v);

    void scan_framework();

    int find_root(vector<int> & pa, int u);

    void my_union(vector<int> & pa,vector<int> & rank, int u, int v);

    void assign_clusterid(vector<int> & pa);

    void pscan_dis();

    bool test_edge_pscan(int u, int v);

    bool test_edge(int u, int v, double dis, vector<int> &parent, vector<int> &rank, bool save = false);

    void traverse_cancidate_edges(vector<int> &parent, vector<int> &rank);
    void traverse_cancidate_edges_sort(vector<int> &parent, vector<int> &rank);

    void traverse_nodes(vector<int> & parent, vector<int> & rank, int mode = 0);
    void traverse_nodes_sort(vector<int> & parent, vector<int> & rank, int mode = 0);

    void traverse_single_node_botk(int u, vector<int> & parent, vector<int> & rank);
    void traverse_single_node_botk_sort(int u, vector<int> & parent, vector<int> & rank,double bin_length, vector<vector<idpair>> & node_buckets);
    void traverse_single_node_dijkstra(int u, vector<idpair> &dis_source_vec, vector<int> & parent, vector<int> & rank);

    void refine_clusters(vector<int> & parent, vector<int> & rank);

    void pruned_scan();

    void assign_noncore_hubs_outliers();

};


#endif //WEIGHTED_SCAN_CLUSTERALGO_H
