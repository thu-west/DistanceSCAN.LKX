#ifndef WEIGHTED_SCAN_SKETCHES_H
#define WEIGHTED_SCAN_SKETCHES_H
#include "binary_tree.h"
#include "graph.h"
#include "fhq_treap.h"

class SKETCHES {

public:
    unsigned long long P = 4294967291;

    vector<int> key2value, value2key;

    vector<double> neis_in_dis;
    vector<int> neis_in_dis_lb, neis_in_dis_ub;

    vector<vector<int>> histogram; // idx * bin < dis < (idx + 1) * bin
    vector<Treap> ads_bot;
    vector<vector<int>> bot_k;
    vector<vector<double>> bot_k_dis;

    vector<int> intersection_lb;
    int intersection_ub;

    SKETCHES() {
        Timer timer(SKETCH_TIME);
        if(config.operation == CONSTRUCT_SKETCHES){
            init_hash_table();
            construct_sketches();
            serialize_sketches();
        }else if(config.operation == GRAPH_MAINTAIN){
            deserialize_sketches();
        } else if (config.algo == BOTK_SCAN||config.algo == MY_ADS) {
            //deserialize_sketches();
            deserialize_sketches2();
        }
    }

    vector<int> get_neis_in_dis_ub() {
        return neis_in_dis_ub;
    }

    static vector<unsigned long long> pick_random_coffs(int k, const unsigned long long &P);

    int double_hashing(const vector<unsigned long long int> &coff_a, const vector<unsigned long long int> &coff_b,
                       int key, int i);

    void init_hash_table();

    void construct_sketches();

    void deserialize_sketches();

    void deserialize_sketches2();

    void serialize_sketches();

    int update_histogram(int source_nid, double path_weight);

    int query_histogram(int source_nid, double dis, bool is_lb);

    int get_bin_id(double dis, bool is_lb);



    void get_approx_neis();

    double jaccard_with_botk(int u, int v);

    double jaccard_with_sketches(int u, int v, double dis);


};

#endif //WEIGHTED_SCAN_SKETCHES_H
