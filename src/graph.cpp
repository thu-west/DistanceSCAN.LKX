
#include "../include/graph.h"

Graph graph;

void Graph::init_d_neighbors() {
    Timer timer(D_NEIGHBOR_TIME);
    path_weight = vector<unordered_map<int, double >>(n, unordered_map<int, double>{});
    d_neighbors = vector<vector<int >>(n, vector<int>());
    for (int n_id = 0; n_id < n; ++n_id) {
        d_neighbors[n_id].emplace_back(n_id);
        priority_queue<idpair, vector<idpair>, cmp_idpair> pq;
        pq.push(make_pair(n_id, 0));
        set_path_weight(n_id, n_id, 0);
        for (int nei_id: d_neighbors[n_id]) {
            pq.push(make_pair(nei_id, get_path_weight(n_id, nei_id)));
        }
        while (!pq.empty()) {
            auto cur_node = pq.top();
            pq.pop();
            if (cmp_double(cur_node.second, get_path_weight(n_id, cur_node.first)) == 1)continue;
            if (cur_node.first > n_id) {
                d_neighbors[n_id].emplace_back(cur_node.first);
                d_neighbors[cur_node.first].emplace_back(n_id);
            }

            for (int nei_id: adj_list[cur_node.first]) {
                if (nei_id < n_id) continue;
                double nei_dis = cur_node.second + edge_weight[cur_node.first][nei_id];
                if (cmp_double(nei_dis, config.distance) < 1) {
                    double old_dis = get_path_weight(n_id, nei_id);
                    if ((cmp_double(old_dis, -1) == 0 || cmp_double(old_dis, nei_dis) > 0)) {
                        set_path_weight(n_id, nei_id, nei_dis);
                        pq.push(make_pair(nei_id, nei_dis));
                    }
                }
            }
        }
    }
}


void Graph::init(const string &graph_path) {
    Timer timer(READ_GRAPH_TIME);
    INFO("Reading graph ...");
    this->data_folder = graph_path;
    init_nm();
    adj_list = vector<vector<int >>(n, vector<int>());
    edge_weight = vector<unordered_map<int, double >>(n, unordered_map<int, double>());
    string graph_file = data_folder;// + FILESEP;
    if (directed) {
        graph_file += "graph.txt";
    } else if (config.operation == CONVERT_GRAPH) {
        graph_file += "undirect_graph.txt";
        weighted = false;
    } else {
        if (data_folder.find("graphs_coauthors") != data_folder.npos) {
            graph_file += "undirect_graph2.txt";
        } else if (data_folder.find("LFR") != data_folder.npos || config.operation == EXPONLFR) {
            graph_file += "undirect_graph.txt";
        } else {
            //graph_file += "uniform_weighted_graph.txt";
            graph_file += "jac_graph.txt";
        }
    }

    FILE *fin = fopen(graph_file.c_str(), "r");
    if (weighted) {
        int t1, t2;
        double w;
        while (fscanf(fin, "%d%d%lf", &t1, &t2, &w) != EOF) {
            if (t1 == t2)continue;
            adj_list[t1].push_back(t2);
            edge_weight[t1][t2] = w;
        }
    } else {
        int t1, t2;
        while (fscanf(fin, "%d%d", &t1, &t2) != EOF) {
            if (t1 == t2)continue;
            adj_list[t1].push_back(t2);
        }
    }
    result.n = this->n;
    result.m = this->m;
    cout << "init graph graphn: " << this->n << " m: " << this->m << endl;

    clusterID = vector<int>(n, -1);
    is_core = vector<int>(n, -1);
    similarity = vector<unordered_map<int, bool >>(n, unordered_map<int, bool>{});
    if (config.algo == W_SCAN) {
        weight_degree = vector<double>(n, 0);
        whole_weight = 0;
        if (config.similarityType == Config::cos) {
            for (int i = 0; i < n; ++i) {
                for (int nei: adj_list[i]) {
                    weight_degree[i] += edge_weight[i][nei] * edge_weight[i][nei];
                }
                weight_degree[i] = sqrt(weight_degree[i] + 1);
                whole_weight += weight_degree[i] + 1;
            }
        }
    }
    if (data_folder.find("graphs_coauthors") != data_folder.npos || data_folder.find("LFR") != data_folder.npos) {
        reweighted(config.type);
    }
    if (config.operation == CLUSTER_VALIDATION) {
        jac_res = vector<unordered_map<int, double >>(n, unordered_map<int, double>{});
    }
}

Graph::Graph(const string &graph_path) {
    Timer timer(READ_GRAPH_TIME);
    INFO("Reading graph ...");
    this->data_folder = graph_path;
    init_nm();
    adj_list = vector<vector<int >>(n, vector<int>());
    edge_weight = vector<unordered_map<int, double >>(n, unordered_map<int, double>());
    string graph_file = data_folder;// + FILESEP;
    if (directed) {
        graph_file += "graph.txt";
    } else if (config.operation == CONVERT_GRAPH) {
        graph_file += "undirect_graph.txt";
    } else {
        if (data_folder.find("graphs_coauthors") != data_folder.npos) {
            graph_file += "undirect_graph2.txt";
        } else if (data_folder.find("LFR") != data_folder.npos || config.operation == EXPONLFR) {
            graph_file += "undirect_graph.txt";
        } else {
            //graph_file += "uniform_weighted_graph.txt";
            graph_file += "jac_graph.txt";
        }
    }

    FILE *fin = fopen(graph_file.c_str(), "r");
    if (weighted) {
        int t1, t2;
        double w;
        while (fscanf(fin, "%d%d%lf", &t1, &t2, &w) != EOF) {
            if (t1 == t2)continue;
            adj_list[t1].push_back(t2);
            edge_weight[t1][t2] = w;
        }
    } else {
        int t1, t2;
        while (fscanf(fin, "%d%d", &t1, &t2) != EOF) {
            if (t1 == t2)continue;
            adj_list[t1].push_back(t2);
        }
    }
    fclose(fin);
    result.n = this->n;
    result.m = this->m;
    cout << "init graph graph n: " << this->n << " m: " << this->m << endl;

    clusterID = vector<int>(n, -1);
    is_core = vector<int>(n, -1);
    similarity = vector<unordered_map<int, bool >>(n, unordered_map<int, bool>{});
    if (config.algo == W_SCAN) {
        weight_degree = vector<double>(n, 0);
        whole_weight = 0;
        if (config.similarityType == Config::cos) {
            for (int i = 0; i < n; ++i) {
                for (int nei: adj_list[i]) {
                    weight_degree[i] += edge_weight[i][nei] * edge_weight[i][nei];
                }
                weight_degree[i] = sqrt(weight_degree[i] + 1);
                whole_weight += weight_degree[i] + 1;
            }
        }
    }
    if (data_folder.find("graphs_coauthors") != data_folder.npos || data_folder.find("LFR") != data_folder.npos) {
        reweighted(config.type);
    }
    if (config.operation == CLUSTER_VALIDATION) {
        jac_res = vector<unordered_map<int, double >>(n, unordered_map<int, double>{});
    }
}

void Graph::init_nm() {
    string attribute_file = data_folder + "attribute.txt";
    assert_file_exist("attribute file", attribute_file);
    ifstream attr(attribute_file);
    string line;
    while (getline(attr, line)) {
        vector<string> output;
        split_string(line, output, "=");
        if (output.size() != 2)continue;
        if (output[0] == "n") {
            n = stoi(output[1]);
        } else if (output[0] == "m") {
            m = stoi(output[1]);
        } else if (output[0] == "weighted") {
            if (output[1] != "0") {
                weighted = true;
            } else {
                weighted = false;
            }
        }
    }
    attr.close();
}

double Graph::get_path_weight(int u, int v) {
    if (u > v) swap(u, v);
    if (path_weight[u].find(v) == path_weight[u].end()) {
        return -1;
    }
    return path_weight[u][v];
}

void Graph::set_path_weight(int u, int v, double w) {
    if (u > v) swap(u, v);
    path_weight[u][v] = w;
}


int Graph::get_similairty(int u, int v) {
    if (u > v) swap(u, v);
    if (similarity[u].find(v) == similarity[u].end()) {
        return -1;
    }
    return similarity[u][v];
}

void Graph::set_similarity(int u, int v, bool s) {
    if (u > v) swap(u, v);
    similarity[u][v] = s;
}

void Graph::init_similarity() {
    vector<unordered_map<int, bool >>(n, unordered_map<int, bool>{}).swap(similarity);
}


unordered_map<int, double> Graph::compute_path_weight(int u, double dis_thre) {
    unordered_map<int, double> nei_dis_table;
    nei_dis_table[u] = 0;
    priority_queue<idpair, vector<idpair>, cmp_idpair> pq;
    pq.push(make_pair(u, 0));
    while (!pq.empty()) {
        auto cur_node = pq.top();
        pq.pop();
        if (cmp_double(cur_node.second, nei_dis_table[cur_node.first]) == 1)continue;
#ifdef _DEBUG_
        assert(cmp_double(cur_node.second, nei_dis_table[cur_node.first]) == 0);
#endif
        for (int nei_id: adj_list[cur_node.first]) {
            double nei_dis = cur_node.second + edge_weight[cur_node.first][nei_id];
            if (cmp_double(nei_dis, dis_thre) < 1) {
                if (nei_dis_table.find(nei_id) == nei_dis_table.end() ||
                    cmp_double(nei_dis_table[nei_id], nei_dis) > 0) {
                    nei_dis_table[nei_id] = nei_dis;
                    pq.push(make_pair(nei_id, nei_dis));
                }
            }
        }
    }
    return nei_dis_table;
}

void Graph::edge_del(int u, int v) {
    m--;
    for (int k = 0; k < 2; ++k) {
        for (int i = adj_list[u].size() - 1; i >= 0; --i) {
            if (adj_list[u][i] == v) {
                swap(adj_list[u][i], adj_list[u].back());
                adj_list[u].pop_back();
                break;
            }
        }
        swap(u, v);
    }
    for (int k = 0; k < 2; ++k) {
        if (edge_weight[u].find(v) != edge_weight[u].end()) {
            edge_weight[u].erase(v);
        }
        swap(u, v);
    }
    if (similarity.empty())return;
    if (u < v) {
        if (similarity[u].find(v) != similarity[u].end()) {
            similarity[u].erase(v);
        }
    } else {
        if (similarity[v].find(u) != similarity[v].end()) {
            similarity[v].erase(u);
        }
    }
}

void Graph::edge_ins(int u, int v, double w) {
    m++;
    adj_list[u].emplace_back(v);
    adj_list[v].emplace_back(u);
    edge_weight[u][v] = w;
    edge_weight[v][u] = w;
}

void Graph::edge_update(int u, int v, double w) {
    edge_weight[u][v] = w;
    edge_weight[v][u] = w;
}

const vector<unordered_map<int, double>> &Graph::getPathWeight() const {
    return path_weight;
}

void Graph::save_graph(const string &new_graph_name) {
    string w_file = data_folder + new_graph_name;
    ofstream out(w_file, ios::out);
    if (out.is_open()) {
        for (int i = 0; i < n; ++i) {
            for (int nei: adj_list[i]) {
                out << i << "\t" << nei << "\t";
                if (edge_weight[i].find(nei) != edge_weight[i].end()) {
                    out << edge_weight[i][nei];
                }
                out << "\n";
            }
        }
    }
}

void Graph::get_undirected_similarity_weighted() {
    string w_file_jac = data_folder + "jac_graph.txt";
    string w_file_cos = data_folder + "cos_graph.txt";
    ofstream out_jac(w_file_jac, ios::out);
    ofstream out_cos(w_file_cos, ios::out);
    for (int i = 0; i < n; ++i) {
        for (int nei: adj_list[i]) {
            int inter = 2, i_index = 0, nei_index = 0;
            while (i_index < adj_list[i].size() && nei_index < adj_list[nei].size()) {
                if (adj_list[i][i_index] == adj_list[nei][nei_index]) {
                    inter++;
                    i_index++;
                    nei_index++;
                } else if (adj_list[i][i_index] < adj_list[nei][nei_index]) {
                    i_index++;
                } else {
                    nei_index++;
                }
            }
            double jac = double(inter) / double(adj_list[i].size() + adj_list[nei].size() + 2 - inter);
            double cos = double(inter) / sqrt(double((adj_list[i].size() + 1) * (adj_list[nei].size() + 1)));
            jac = 1 - (0.9 * jac);
            cos = 1 - (0.9 * cos);
            if (out_jac.is_open()) {
                out_jac << i << "\t" << nei << "\t" << jac << "\n";
            }
            if (out_cos.is_open()) {
                out_cos << i << "\t" << nei << "\t" << cos << "\n";
            }
        }
    }
}

void Graph::get_undirected_weighted_graph_snap_uniform() {
    static default_random_engine generator(time(NULL));
    static uniform_real_distribution<double> dis(0, 0.9);
    for (int i = 0; i < n; ++i) {
        for (int nei: adj_list[i]) {
            if (edge_weight[i].find(nei) == edge_weight[i].end()) {
                edge_weight[i][nei] = 1 - dis(generator);
                edge_weight[nei][i] = edge_weight[i][nei];
            }
        }
    }
    save_graph("uniform_weighted_graph.txt");
}

void Graph::get_undirected_weighted_graph_snap_exponent(double log_base) {
    std::default_random_engine generator(time(NULL));
    std::uniform_real_distribution<double> dis(0.1, 1);
    for (int i = 0; i < n; ++i) {
        for (int nei: adj_list[i]) {
            double w = -1 * log(1 - dis(generator));
            if (log_base > 1) {
                w /= log(log_base);
            }
            if (edge_weight[i].find(nei) == edge_weight[i].end()) {
                edge_weight[i][nei] = w;
                edge_weight[nei][i] = edge_weight[i][nei];
            }
        }
    }
    if (log_base > 1) {
        save_graph("exponent_" + to_string(log_base) + "_weighted_graph.txt");
    } else {
        save_graph("exponent_e_weighted_graph.txt");
    }
}

void Graph::handle_LFR_graph(const string &graph_path, int nodes) {
    INFO("Reading graph ...");
    this->data_folder = graph_path;
    INFO(nodes);
    if (nodes == -1) {
        string attribute_file = data_folder + "flags.dat";
        INFO(attribute_file);
        ifstream attr(attribute_file);
        string line;
        while (getline(attr, line)) {
            vector<string> output;
            split_string(line, output, " ");
            INFO(line);
            INFO(output);
            if (output.size() != 2)continue;
            if (output[0] == "-N") {
                n = stoi(output[1]) + 1;
                break;
            }
        }
        attr.close();
    } else {
        n = nodes;
    }
    INFO(n);
    adj_list = vector<vector<int >>(n, vector<int>());
    edge_weight = vector<unordered_map<int, double >>(n, unordered_map<int, double>());
    string graph_file = data_folder + "network.dat";

    m = 0;
    FILE *fin = fopen(graph_file.c_str(), "r");
    int t1, t2;
    double w;
    while (fscanf(fin, "%d%d%lf", &t1, &t2, &w) != EOF) {
        if (t1 == t2 || cmp_double(w, 0) == 0)continue;
        adj_list[t1].push_back(t2);
        if (edge_weight[t1].find(t2) == edge_weight[t1].end()) {
            edge_weight[t1][t2] = w;
            edge_weight[t2][t1] = w;
        } else {
            edge_weight[t1][t2] += w;
            edge_weight[t2][t1] += w;
            //edge_weight[t1][t2] /= 2;
            //edge_weight[t2][t1] /= 2;
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int nei: adj_list[i]) {
            if (cmp_double(edge_weight[i][nei], 0) < 1) {
                edge_del(i, nei);
            }
        }
    }
    m = 0;
    for (int i = 0; i < n; ++i) {
        m += adj_list[i].size();
    }
    fclose(fin);
    cout << "init graph graph n: " << this->n << " m: " << this->m << endl;
    save_graph("raw_undirect_graph.txt");
    for (int i = 0; i < n; ++i) {
        for (int nei: adj_list[i]) {
            edge_weight[i][nei] = 1 / edge_weight[i][nei];
        }
    }
    m = 0;
    for (int i = 0; i < n; ++i) {
        m += adj_list[i].size();
    }
    save_graph("undirect_graph.txt");
    string w_file = data_folder + "attribute.txt";
    ofstream out(w_file, ios::out);
    if (out.is_open()) {
        out << "n=" << n << "n" << endl;
        out << "m=" << m << "n" << endl;
    }
    out.close();
}

double Graph::get_jac_res(int u, int v) {
    if (u > v) swap(u, v);
    if (jac_res[u].find(v) == jac_res[u].end()) {
        return -1;
    }
    return jac_res[u][v];
}

void Graph::set_jac_res(int u, int v, double s) {
    if (u > v) swap(u, v);
    jac_res[u][v] = s;
}

double Graph::jaccard_raw_wscan(int u, int v) {
    double inter = 0;
    unordered_map<int, bool> union_table;
    for (int nei: adj_list[u]) {
        union_table[nei] = true;
    }
    for (int nei: adj_list[v]) {
        if (union_table.find(nei) != union_table.end()) {
            if (u == nei || v == nei) {
                inter += exp(1 / edge_weight[u][v] - 1);
            } else {
                double real_w1 = exp(1 / edge_weight[u][nei] - 1);
                double real_w2 = exp(1 / edge_weight[v][nei] - 1);
                inter += real_w1 * real_w2;
            }
        }
    }
    assert(config.similarityType == Config::cos);
    return inter / weight_degree[u] / weight_degree[v];
}

double Graph::jaccard_raw(int u, int v) {
    unordered_map<int, bool> union_table;
    for (int nei: adj_list[u]) {
        union_table[nei] = true;
    }
    for (int nei: adj_list[v]) {
        union_table[nei] = true;
    }
    int inter_size = adj_list[u].size() + adj_list[v].size() - union_table.size();
    double similarity_result;
    if(config.similarityType == config.jac){
        similarity_result = double(inter_size) / double(union_table.size());
    } else if (config.similarityType == config.cos) {
        similarity_result = inter_size / sqrt(adj_list[u].size()) / sqrt(adj_list[v].size());
    } else if(config.similarityType == config.set_containment1){
        similarity_result = double(inter_size) / double (adj_list[u].size());
    } else if(config.similarityType == config.set_containment2){
        similarity_result = double(inter_size) / double (adj_list[v].size());
    }
    if (config.operation == CLUSTER_VALIDATION) {
        set_jac_res(u, v, similarity_result);
    }
    return similarity_result;
}

void Graph::reweighted(int type) {
    int max_weight = 0, sum_weight = 0;
    for (int i = 0; i < n; ++i) {
        for (int nei: adj_list[i]) {
            if (max_weight < edge_weight[i][nei]) {
                max_weight = edge_weight[i][nei];
            }
            sum_weight += edge_weight[i][nei];
        }
    }
    sum_weight /= 2;
    for (int i = 0; i < n; ++i) {
        for (int nei: adj_list[i]) {
            int x = edge_weight[i][nei];
            switch (type) {
                case 0:
                    // 1/(log(x) + 1)
                    edge_weight[i][nei] = 1 / (log(x) + 1);
                    //INFO(edge_weight[i][nei]);
                    break;
                case 1:
                    // log(sum_weight/x)/log(sum_weight) 归一化后的 inverse document frequency
                    edge_weight[i][nei] = log(sum_weight / x) / log(sum_weight);
                    break;
                case 2:
                    // (log(sum_weight/(1+x))+1)/(log(sum_weight/2)+1) 归一化的 inverse document frequency smooth
                    edge_weight[i][nei] = (log(sum_weight / (1 + x)) + 1) / (log(sum_weight / 2) + 1);
                    break;
                case 3:
                    // log((max_weight+1)/x)/log(max_weight+1) 归一化的 inverse document frequency max
                    edge_weight[i][nei] = log((max_weight + 1) / x) / log(max_weight + 1);
                    break;
                case 4:
                    // log((sum_weight-x)/x)/log(max_weight-1) 归一化的 probabilistic inverse document frequency
                    edge_weight[i][nei] = log((sum_weight - x) / x) / log(sum_weight - 1);
                    break;

            }
        }
    }
}

ipair Graph::random_choose_edge() {
    static std::default_random_engine generator(time(NULL));
    std::uniform_int_distribution<int> dis(0, n - 1);
    int source = -1, target = -1;
    do {
        source = dis(generator);
    } while (adj_list[source].empty());
    int max_idx = adj_list[source].size() - 1;
    std::uniform_int_distribution<int> dis2(0, max_idx);
    target = dis2(generator);
    return make_pair(source, adj_list[source][target]);
}

void Graph::convert_to_undirected_graph(string graph_path) {
    data_folder = config.graph_location + graph_path;
    init_nm();
    vector<set<int>> undirect_adj_list(n);
    string graph_file = data_folder + FILESEP + "graph.txt";
    assert_file_exist("graph file", graph_file);
    FILE *fin = fopen(graph_file.c_str(), "r");
    char line[1000];
    long long t1, t2;
    int min_id = INT_MAX, max_id = -1, final_n = 0, final_m = 0, tmp = 0;
    unordered_map<long, int> str_index;
    fin = fopen(graph_file.c_str(), "r");
    while (fgets(line, 1000, fin) != NULL) {
        if (line[0] > '9' || line[0] < '0') { continue; }
        sscanf(line, "%lld\t%lld\n", &t1, &t2);
        if (t1 == t2)continue;
        min_id = min(min_id, (int) t1);
        min_id = min(min_id, (int) t2);
        max_id = max(max_id, (int) t1);
        max_id = max(max_id, (int) t2);
        undirect_adj_list[t1].insert(t2);
        undirect_adj_list[t2].insert(t1);
    }
    string undirect_graph_file = data_folder + FILESEP + "undirect_graph.txt";
    FILE *fout = fopen(undirect_graph_file.c_str(), "w");
    for (int j = 0; j < undirect_adj_list.size(); ++j) {
        if (!undirect_adj_list[j].empty()) {
            final_n++;
            for (int des: undirect_adj_list[j]) {
                final_m++;
                fprintf(fout, "%d\t%d\n", j, des);
            }
        }
    }
    fclose(fout);
    string undirect_graph_attr_file = data_folder + FILESEP + "undirect_graph_attribute.txt";
    fout = fopen(undirect_graph_attr_file.c_str(), "w");
    fprintf(fout, "n=%d\nm=%lld\n", n, m);
    fprintf(fout, "final_n=\t%d\nfinal_m=\t%d\n", final_n, final_m);
    fprintf(fout, "min_node=\t%d\nmax_node=\t%d\n", min_id, max_id);
    for (auto item: str_index) {
        fprintf(fout, "%ld\t%d\n", item.first, item.second);
    }
    fclose(fin);
    fclose(fout);
}