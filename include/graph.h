#ifndef SCAN_GRAPH_H
#define SCAN_GRAPH_H

#include "mylib.h"
#include "config.h"

class Graph {
private:
    vector<unordered_map<int, double>> path_weight;
    vector<unordered_map<int, bool>> similarity;
    vector<unordered_map<int, double>> jac_res;

public:
    int n;
    long long m;

    string data_folder;

    bool directed = false;
    bool weighted = true;

    vector<vector<int>> d_neighbors;
    vector<vector<int>> adj_list;
    vector<unordered_map<int, double>> edge_weight;

    vector<int> clusterID;//-1 not classified, -2 hub, -3 outlier.
    vector<int> is_core;//-1 init, 0 non-member, 1 core, 2 hub, 3 outlier.

    vector<double> weight_degree;
    double whole_weight;

    void init_nm();

    Graph() = default;

    explicit Graph(const string &graph_path);

    void init(const string &graph_path);

    void init_d_neighbors();

    double get_path_weight(int u, int v);

    void set_path_weight(int u, int v, double w);

    int get_similairty(int u, int v);

    void set_similarity(int u, int v, bool s);

    void init_similarity();

    double get_jac_res(int u, int v);

    void set_jac_res(int u, int v, double s);

    unordered_map<int, double>  compute_path_weight(int u, double dis_thre);

    void edge_del(int u, int v);

    void edge_ins(int u, int v, double w);

    void edge_update(int u, int v, double w);

    const vector<unordered_map<int, double>> &getPathWeight() const;

    void get_undirected_weighted_graph_snap_uniform();

    void get_undirected_weighted_graph_snap_exponent(double exp = -1);

    void get_undirected_similarity_weighted();

    void save_graph(const string &new_graph_name);

    void handle_LFR_graph(const string &graph_path, int nodes = -1);

    double jaccard_raw(int u, int v);
    double jaccard_raw_wscan(int u, int v);

    void reweighted(int type);

    ipair random_choose_edge();
    void convert_to_undirected_graph(string graph_path);

};

extern Graph graph;
#endif //SCAN_GRAPH_H
