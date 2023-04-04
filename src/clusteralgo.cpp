#include "clusteralgo.h"

bool ClusterAlgos::test_edge_jaccard(int u, int v) {
    int is_sim = graph.get_similairty(u, v);
    if (is_sim != -1) return is_sim == 1;
    if (config.algo == W_SCAN) {
        is_sim = cmp_double(graph.jaccard_raw_wscan(u, v), config.epsilon) > -1;
    } else {
        is_sim = cmp_double(graph.jaccard_raw(u, v), config.epsilon) > -1;
    }
    graph.set_similarity(u, v, is_sim);
    return is_sim;
}

bool ClusterAlgos::check_core(int v) {
    int sim_neis_counter = 0;
    if (graph.is_core[v] == -1) {
        if (graph.adj_list[v].size() >= config.mu) {
            for (int nei: graph.adj_list[v]) {
                if (test_edge_jaccard(v, nei)) {
                    sim_neis_counter++;
                }
            }
        }
        if (sim_neis_counter >= config.mu) {
            graph.is_core[v] = 1;
            return true;
        } else {
            //if v is not a core, it is labeled as non-member label
            graph.is_core[v] = 0;
            return false;
        }
    } else if (graph.is_core[v] == 1) {
        return true;
    } else {
        return false;
    }
}

void ClusterAlgos::scan_framework() {
    Timer timer1(BASIC_TIME);
    clusters.clear();
    if (config.algo == BASIC) {
        graph.init_d_neighbors();
        swap(graph.d_neighbors, graph.adj_list);
    }
    for (int i = 0; i < graph.n; ++i) {
        if (graph.clusterID[i] != -1) continue;
        // check whether i is a core.
        if (check_core(i)) {
            //if v is a core, a new cluster is expanded
            int new_clusterID = int(clusters.size());
            vector<int> new_cluster;//{i};
            graph.clusterID[i] = new_clusterID;

            queue<int> q;
            q.push(i);
            while (!q.empty()) {
                int y = q.front();
                q.pop();
                int sim_neis_counter = 0;
                for (int nei: graph.adj_list[y]) {
                    if (test_edge_jaccard(y, nei)) {
                        sim_neis_counter++;
                    }
                }
                if (sim_neis_counter >= config.mu) {
                    new_cluster.emplace_back(y);
                    graph.is_core[y] = 1;
                    if (config.algo != BASIC) {
                        for (int nei: graph.adj_list[y]) {
                            if (graph.clusterID[nei] == -1 && test_edge_jaccard(y, nei)) {
                                graph.clusterID[nei] = new_clusterID;
                                //new_cluster.emplace_back(nei);
                                q.push(nei);
                            }
                        }
                    } else {
                        for (int nei: graph.d_neighbors[y]) {
                            if (cmp_double(graph.edge_weight[y].at(nei), config.distance) == 1) continue;
                            if (graph.clusterID[nei] == -1 && test_edge_jaccard(y, nei)) {
                                graph.clusterID[nei] = new_clusterID;
                                //new_cluster.emplace_back(nei);
                                q.push(nei);
                            }
                        }
                    }
                } else {
                    graph.is_core[y] = 0;
                }
            }
            clusters.push_back(new_cluster);
        }
    }
    if (config.algo == BASIC) {
        swap(graph.d_neighbors, graph.adj_list);
    }
}

void ClusterAlgos::assign_noncore_hubs_outliers() {
    vector<set<int>> non_core_cluster_members = vector<set<int>>(clusters.size(), set<int>());
    for (int i = 0; i < graph.n; ++i) {
        if (graph.is_core[i] == 1)continue;
        idpair nearest_core_last(-1, -1);
        bool is_hub = false, in_cluster = false;
        int nei_cls_id = -1;
        for (int nei: graph.adj_list[i]) {
            if (graph.is_core[nei] == 1 ) {
                if (config.algo == SCAN || config.algo == W_SCAN ||
                    (graph.edge_weight[i].find(nei) != graph.edge_weight[i].end() &&
                     cmp_double(graph.edge_weight[i].at(nei), config.distance) < 1)) {
                    idpair nei_dis(nei, -1);
                    if (config.algo == BOTK_SCAN || config.algo == MY_ADS) {
                        nei_dis.second = sketches.jaccard_with_botk(i, nei);
                    } else if (config.algo == SCAN || config.algo == W_SCAN) {
                        nei_dis.second = graph.jaccard_raw(i, nei);
                    } else {
                        swap(graph.adj_list[i], graph.d_neighbors[i]);
                        swap(graph.adj_list[nei], graph.d_neighbors[nei]);
                        nei_dis.second = graph.jaccard_raw(i, nei);
                        swap(graph.adj_list[i], graph.d_neighbors[i]);
                        swap(graph.adj_list[nei], graph.d_neighbors[nei]);
                    }
                    if (cmp_double(nei_dis.second, config.epsilon) > -1) {
                        in_cluster = true;
                        non_core_cluster_members[graph.clusterID[nei]].insert(i);
                    } else if (!in_cluster) {
                        if (nearest_core_last.first != -1 &&
                            graph.clusterID[nearest_core_last.first] != graph.clusterID[nei]) {
                            is_hub = true;
                        }
                    }
                    int cmp_res = cmp_double(nei_dis.second, nearest_core_last.second);
                    if (cmp_res == 1) {
                        nearest_core_last = nei_dis;
                    } else if (cmp_res == 0 && nei_dis.first < nearest_core_last.first) {
                        nearest_core_last = nei_dis;
                    }
                }
            }
            if(!in_cluster && graph.clusterID[nei] > -1){
                if(nei_cls_id == -1){
                    nei_cls_id = graph.clusterID[nei];
                } else if(nei_cls_id != graph.clusterID[nei]){
                    is_hub = true;
                }
            }
        }
        if (in_cluster) {
            graph.clusterID[i] = graph.clusterID[nearest_core_last.first];
            graph.is_core[i] = 4;
        } else if (is_hub) {
            graph.is_core[i] = 2;
            graph.clusterID[i] = -2;
            hub.emplace_back(i);
        } else {
            graph.is_core[i] = 3;
            graph.clusterID[i] = -3;
            outlier.emplace_back(i);
        }
    }
    for (int i = 0; i < non_core_cluster_members.size(); ++i) {
        for (int nid: non_core_cluster_members[i]) {
            clusters[i].emplace_back(nid);
        }
    }
}


bool ClusterAlgos::test_edge_pscan(int u, int v) {
    if (u == v)return false;
    int is_sim = graph.get_similairty(u, v);
    if (is_sim != -1) return is_sim == 1;
    is_sim = cmp_double(graph.jaccard_raw(u, v), config.epsilon) > -1;
    graph.set_similarity(u, v, is_sim == 1);
    if (is_sim) {
        if (graph.is_core[u] == -1 && ++lb[u] >= config.mu) {
            graph.is_core[u] = 1;
        }
        if (graph.is_core[v] == -1 && ++lb[v] >= config.mu) {
            graph.is_core[v] = 1;
        }
    } else {
        if (graph.is_core[u] == -1 && --ub[u] < config.mu) {
            graph.is_core[u] = 0;
        }
        if (graph.is_core[v] == -1 && --ub[v] < config.mu) {
            graph.is_core[v] = 0;
        }
    }
    return is_sim;
}

bool ClusterAlgos::test_edge(int u, int v, double dis, vector<int> &parent, vector<int> &rank, bool save) {
    if (u == v)return false;
    int is_sim = graph.get_similairty(u, v), is_edge = false;
    if (is_sim != -1) return is_sim == 1;
    is_sim = sketches.jaccard_with_sketches(u, v, dis);
    if (graph.edge_weight[u].find(v) != graph.edge_weight[u].end()) {
        graph.set_similarity(u, v, is_sim == 1);
        is_edge = true;
    }
    if (save) {
        graph.set_similarity(u, v, is_sim == 1);
    }
    if (is_sim) {
        if (graph.is_core[u] == 1 && graph.is_core[v] == 1 && is_edge &&
            cmp_double(graph.edge_weight[u].at(v), config.distance) < 1) {
            my_union(parent, rank, u, v);
        }

        if (graph.is_core[u] == -1 && ++lb[u] >= config.mu) {
            graph.is_core[u] = 1;
            for (auto item: graph.edge_weight[u]) {
                int nei = item.first;
#ifdef _DEBUG_
                assert(cmp_double(item.second, config.distance) < 1);
#endif
                if (graph.is_core[nei] == 1 && graph.get_similairty(u, nei) == 1) {
                    my_union(parent, rank, u, nei);
                }
            }
        }
        if ((is_edge||save) && graph.is_core[v] == -1 && ++lb[v] >= config.mu) {
            graph.is_core[v] = 1;
            for (auto item: graph.edge_weight[v]) {
                int nei = item.first;
                if (graph.is_core[nei] == 1 && graph.get_similairty(v, nei) == 1) {
                    my_union(parent, rank, v, nei);
                }
            }
        }
    } else {
        if (graph.is_core[u] == -1 && --ub[u] < config.mu) {
            graph.is_core[u] = 0;
        }
        if ((is_edge||save) && graph.is_core[v] == -1 && --ub[v] < config.mu) {
            graph.is_core[v] = 0;
        }
    }
    return is_sim;
}

void
ClusterAlgos::traverse_cancidate_edges(vector<int> &parent, vector<int> &rank) {
    Timer timer(TRAVERSE_EDGES_TIME);
    for (int u = 0; u < graph.n; ++u) {
        for (auto item : graph.edge_weight[u]){
            int v = item.first;
            if (graph.is_core[u] == 0 && graph.is_core[v] == 0) continue;
            if (graph.is_core[u] == 1 && graph.is_core[v] == 1 && find_root(parent, u) == find_root(parent, v))
                continue;
            test_edge(u, v, graph.edge_weight[u].at(v), parent, rank);
        }
    }
}

void
ClusterAlgos::traverse_cancidate_edges_sort(vector<int> &parent, vector<int> &rank) {
    Timer timer(TRAVERSE_EDGES_SORT_TIME);
    // int num_of_bins = ceil(config.distance / config.bin) + 1;
    int num_of_bins = ceil(config.distance / 0.01) + 1;
    vector<vector<ipair>> candidate_edges = vector<vector<ipair>>(num_of_bins, vector<ipair>());
    for (int i = 0; i < graph.n; ++i) {
        for (auto item: graph.edge_weight[i]) {
            int nei = item.first;
            double dis = item.second;
            // int bin_id = floor(dis / config.bin);
            int bin_id = floor(dis / 0.01);
            candidate_edges[bin_id].emplace_back(make_pair(i, nei));
        }
    }
    for (auto edge_set: candidate_edges) {
        for (ipair edge: edge_set) {
            int u = edge.first, v = edge.second;
            if (graph.is_core[u] == 0 && graph.is_core[v] == 0)continue;
            if (graph.is_core[u] == 1 && graph.is_core[v] == 1 && find_root(parent, u) == find_root(parent, v))
                continue;
#ifdef _DEBUG_
            assert(cmp_double(graph.edge_weight[u].at(v), config.distance) < 1);
#endif
            test_edge(u, v, graph.edge_weight[u].at(v), parent, rank);
        }
    }
}

void ClusterAlgos::traverse_nodes(vector<int> &parent, vector<int> &rank, int mode) {
    Timer timer(TRAVERSE_NODES_TIME);
    vector<idpair> dis_source_vec = vector<idpair>(graph.n, make_pair(-1, -1));
    double bin_length;
    vector<vector<idpair>> node_buckets;
    if(mode == 2){
        bin_length = (config.distance - MIN_DISTANCE) / BINS_NEIGHBOR;
        node_buckets = vector<vector<idpair>>(BINS_NEIGHBOR + 3, vector<idpair>());
    }
    for (int n_id = 0; n_id < graph.n; ++n_id) {
        if (graph.is_core[n_id] != -1)continue;
        // checkcore
        if (mode == 2) {
            traverse_single_node_botk_sort(n_id, parent, rank, bin_length, node_buckets);
            if (graph.is_core[n_id] == -1) {
                traverse_single_node_dijkstra(n_id, dis_source_vec, parent, rank);
            }
        } else if (mode == 1) {
            traverse_single_node_botk(n_id, parent, rank);
            if (graph.is_core[n_id] == -1) {
                traverse_single_node_dijkstra(n_id, dis_source_vec, parent, rank);
            }
        } else if (mode == 0) {
            traverse_single_node_dijkstra(n_id, dis_source_vec, parent, rank);
        }
        if (graph.is_core[n_id] == -1) {
            graph.is_core[n_id] = 0;
        }
    }
}

void ClusterAlgos::traverse_nodes_sort(vector<int> &parent, vector<int> &rank, int mode) {
    Timer timer(TRAVERSE_NODES_SORT_TIME);
    int max_ub = 0;
    vector<vector<int>> node_buckets_whole = vector<vector<int>>(graph.n, vector<int>());
    for (int i = 0; i < graph.n; ++i) {
        node_buckets_whole[ub[i]].emplace_back(i);
        max_ub = max(max_ub, ub[i]);
    }
    vector<idpair> dis_source_vec = vector<idpair>(graph.n, make_pair(-1, -1));

    double bin_length;
    vector<vector<idpair>> node_buckets;
    if(mode == 2){
        bin_length = (config.distance - MIN_DISTANCE) / BINS_NEIGHBOR;
        node_buckets = vector<vector<idpair>>(BINS_NEIGHBOR + 3, vector<idpair>());
    }
    for (int bucket_id = max_ub; bucket_id >= config.mu; --bucket_id) {
        for (int n_id: node_buckets_whole[bucket_id]) {
            if (bucket_id != ub[n_id]) {
                node_buckets_whole[ub[n_id]].emplace_back(n_id);
                continue;
            }
            if (graph.is_core[n_id] != -1)continue;
            // checkcore
            if (mode == 2) {
                traverse_single_node_botk_sort(n_id, parent, rank,bin_length, node_buckets);
                if (graph.is_core[n_id] == -1) {
                    traverse_single_node_dijkstra(n_id, dis_source_vec, parent, rank);
                }
            } else if (mode == 1) {
                traverse_single_node_botk(n_id, parent, rank);
                if (graph.is_core[n_id] == -1) {
                    traverse_single_node_dijkstra(n_id, dis_source_vec, parent, rank);
                }
            } else if (mode == 0) {
                traverse_single_node_dijkstra(n_id, dis_source_vec, parent, rank);
            }
            if (graph.is_core[n_id] == -1) {
                graph.is_core[n_id] = 0;
            }
        }
    }
}

void ClusterAlgos::refine_clusters(vector<int> &parent, vector<int> &rank) {
    // int num_of_bins = ceil(config.distance / config.bin) + 1;
    int num_of_bins = ceil(config.distance / 0.01) + 1;
    vector<vector<ipair>> candidate_edges = vector<vector<ipair>>(num_of_bins, vector<ipair>());

    for (int i = 0; i < graph.n; ++i) {
        for (auto item: graph.edge_weight[i]) {
            int nei = item.first;
            if (graph.is_core[i] == 1 && graph.is_core[nei] == 1 &&
                find_root(parent, i) != find_root(parent, nei) && i < nei &&
                graph.get_similairty(i, nei) == -1) {
                // int bin_id = floor(item.second / config.bin);
                int bin_id = floor(item.second / 0.01);
                candidate_edges[bin_id].emplace_back(make_pair(i, nei));
            }
        }
    }

    for (auto edge_set: candidate_edges) {
        for (ipair edge: edge_set) {
            int u = edge.first, v = edge.second;
            if (find_root(parent, u) == find_root(parent, v)) {
                continue;
            }
            test_edge(u, v, graph.edge_weight[u].at(v), parent, rank);
        }
    }
    assign_clusterid(parent);
}

void ClusterAlgos::traverse_single_node_botk(int u, vector<int> &parent, vector<int> &rank) {
    Timer timer(TRAVERSE_NODES_TIME_BOTK);
    for (int j = 0; j < sketches.bot_k[u].size(); ++j) {
        int nei = sketches.value2key[sketches.bot_k[u][j]];
        double dis = sketches.bot_k_dis[u][j];
        test_edge(u, nei, dis, parent, rank, true);
        if (graph.is_core[u] != -1)break;
    }
}

void ClusterAlgos::traverse_single_node_botk_sort(int u, vector<int> &parent, vector<int> &rank,double bin_length,  vector<vector<idpair>> & node_buckets) {
    Timer timer(TRAVERSE_NODES_TIME_BOTK_SORT);
    vector<int> bin_true_size = vector<int>(node_buckets.size(), 0);
    for (int i = 0; i < sketches.bot_k_dis[u].size(); ++i) {
        int bin_id;
        if (cmp_double(sketches.bot_k_dis[u][i], 0) == 0) {
            bin_id = 0;
        } else {
            bin_id = (sketches.bot_k_dis[u][i] - MIN_DISTANCE) / bin_length;
        }
        int nei = sketches.value2key[sketches.bot_k[u][i]];
        double dis = sketches.bot_k_dis[u][i];
        if(node_buckets[bin_id].size() > bin_true_size[bin_id]){
            node_buckets[bin_id][bin_true_size[bin_id]++] = make_pair(nei, dis);
        }else{
            node_buckets[bin_id].emplace_back(make_pair(nei, dis));
            bin_true_size[bin_id]++;
        }
    }
    for (int i = 0; i < node_buckets.size(); ++i) {
        for (int j = 0; j < bin_true_size[i]; ++j) {
            idpair item = node_buckets[i][j];
            test_edge(u, item.first, item.second, parent, rank, true);
            if (graph.is_core[u] != -1)break;
        }
    }
}

void ClusterAlgos::traverse_single_node_dijkstra(int u, vector<idpair> &dis_source_vec, vector<int> &parent,
                                                 vector<int> &rank) {
    Timer timer(TRAVERSE_NODES_TIME_DIJKSTRA);
    idpair source = make_pair(u, 0);
    dis_source_vec[u] = source;
    priority_queue<idpair, vector<idpair>, cmp_idpair> pq;
    pq.push(source);
    while (!pq.empty()) {
        idpair cur_node = pq.top();
        pq.pop();
        if (cmp_double(cur_node.second, dis_source_vec[cur_node.first].second) == 1)continue;
        test_edge(u, cur_node.first, cur_node.second, parent, rank);
        if (graph.is_core[u] != -1)break;

        for (auto item: graph.edge_weight[cur_node.first]) {
            int nei_id = item.first;
            double nei_dis = cur_node.second + graph.edge_weight[cur_node.first][nei_id];
            if (cmp_double(nei_dis, config.distance) < 1) {

                if (dis_source_vec[nei_id].first != u ||
                    cmp_double(dis_source_vec[nei_id].second, nei_dis) > 0) {
                    dis_source_vec[nei_id].first = u;
                    dis_source_vec[nei_id].second = nei_dis;
                    pq.push(make_pair(nei_id, nei_dis));
                }
            }
        }
    }
}
void ClusterAlgos::pruned_scan() {
    Timer timer(DIS_SCAN_TIME);
    sketches.intersection_ub = floor(config.hash_k * config.epsilon);
    for (int i = 0; i < config.hash_k; ++i) {
        sketches.intersection_lb.emplace_back(ceil(config.hash_k * config.epsilon - config.hash_k + i));
    }
    sketches.get_approx_neis();

    clusters.clear();
    lb = vector<int>(graph.n, 1), ub = sketches.get_neis_in_dis_ub();
    for (int i = 0; i < graph.n; ++i) {
        if (ub[i] < config.mu)graph.is_core[i] = 0;
        for (int nei: graph.adj_list[i]) {
            if (cmp_double(graph.edge_weight[i].at(nei), config.distance) == 1) {
                graph.edge_weight[i].erase(nei);
            }
        }
    }

    vector<int> parent = vector<int>(graph.n), rank = vector<int>(graph.n);
    iota(parent.begin(), parent.end(), 0);
    if (config.algo_type == 1) {
        traverse_nodes(parent, rank, 1);
    } else if (config.algo_type == 2) {
        traverse_cancidate_edges_sort(parent, rank);
        traverse_nodes(parent, rank, 1);
    } else if (config.algo_type == 3) {
        traverse_nodes_sort(parent, rank, 1);
    } else if (config.algo_type == 4) {
        traverse_cancidate_edges_sort(parent, rank);
        traverse_nodes_sort(parent, rank, 1);
    } else if (config.algo_type == 5) {
        traverse_cancidate_edges(parent, rank);
        traverse_nodes_sort(parent, rank, 1);
    } else if (config.algo_type == 6) {
        traverse_nodes(parent, rank, 2);
    } else if (config.algo_type == 7) {
        traverse_cancidate_edges_sort(parent, rank);
        traverse_nodes(parent, rank, 2);
    } else if (config.algo_type == 8) {
        traverse_nodes_sort(parent, rank, 2);
    } else if (config.algo_type == 9) {
        traverse_cancidate_edges_sort(parent, rank);
        traverse_nodes_sort(parent, rank, 2);
    } else if (config.algo_type == 10) {
        traverse_cancidate_edges(parent, rank);
        traverse_nodes_sort(parent, rank, 2);
    } else if (config.algo_type == 11) {
        traverse_nodes(parent, rank, 0);
    } else if (config.algo_type == 12) {
        traverse_cancidate_edges_sort(parent, rank);
        traverse_nodes(parent, rank, 0);
    } else if (config.algo_type == 13) {
        traverse_nodes_sort(parent, rank, 0);
    } else if (config.algo_type == 14) {
        traverse_cancidate_edges_sort(parent, rank);
        traverse_nodes_sort(parent, rank, 0);
    } else if (config.algo_type == 15) {
        traverse_cancidate_edges(parent, rank);
        traverse_nodes_sort(parent, rank, 0);
    }
    refine_clusters(parent, rank);
}

int ClusterAlgos::find_root(vector<int> &pa, int u) {
    int x = u;
    while (pa[x] != x) x = pa[x];

    while (pa[u] != x) {
        int tmp = pa[u];
        pa[u] = x;
        u = tmp;
    }

    return x;
}

void ClusterAlgos::my_union(vector<int> &pa, vector<int> &rank, int u, int v) {
    int ru = find_root(pa, u);
    int rv = find_root(pa, v);

    if (ru == rv) return;

    if (rank[ru] < rank[rv]) pa[ru] = rv;
    else if (rank[ru] > rank[rv]) pa[rv] = ru;
    else {
        pa[rv] = ru;
        ++rank[ru];
    }
}

void ClusterAlgos::assign_clusterid(vector<int> &pa) {
    unordered_map<int, vector<int>> tmp_clusters;
    for (int i = 0; i < graph.n; ++i) {
        if (graph.is_core[i] != 1)continue;
        int key = find_root(pa, i);
        tmp_clusters[key].emplace_back(i);
    }
    for (auto item: tmp_clusters) {
        for (int nid: item.second) {
            graph.clusterID[nid] = clusters.size();
        }
        clusters.emplace_back(item.second);
    }
}

void ClusterAlgos::pscan_dis() {
    Timer timer(PSCAN_TIME);
    clusters.clear();
    graph.init_d_neighbors();
    swap(graph.d_neighbors, graph.adj_list);
    int max_ub = 0;
    vector<vector<int>> node_bins = vector<vector<int>>(graph.n, vector<int>());
    lb = vector<int>(graph.n, 1);
    vector<int> parent = vector<int>(graph.n), rank = vector<int>(graph.n);
    for (int i = 0; i < graph.n; ++i) {
        parent[i] = i;
        ub.emplace_back(graph.adj_list[i].size());
        max_ub = max_ub < ub.back() ? ub.back() : max_ub;
        node_bins[ub.back()].emplace_back(i);
    }
    for (int i = 0; i < graph.n; ++i) {
        for (int nei: graph.adj_list[i]) {
            if(i==nei)continue;
            if (graph.get_similairty(i, nei) != -1)continue;
            int i_size = graph.adj_list[i].size() + 1, nei_size = graph.adj_list[nei].size() + 1;
            if (i_size < config.epsilon * nei_size || nei_size < config.epsilon * i_size) {
                if (graph.is_core[i] == -1 && --ub[i] < config.mu) {
                    graph.is_core[i] = 0;
                }
                if (graph.is_core[nei] == -1 && --ub[nei] < config.mu) {
                    graph.is_core[nei] = 0;
                }
                graph.set_similarity(i, nei, false);
            } else if (2 > (i_size + nei_size - 2) * config.epsilon) {
                if (graph.is_core[i] == -1 && ++lb[i] >= config.mu) {
                    graph.is_core[i] = 1;
                }
                if (i!= nei && graph.is_core[nei] == -1 && ++lb[nei] >= config.mu) {
                    graph.is_core[nei] = 1;
                }
                graph.set_similarity(i, nei, true);
            }
        }
    }
    for (int bucket_id = max_ub; bucket_id >= config.mu; --bucket_id) {
        for (int nid: node_bins[bucket_id]) {
            if (bucket_id != ub[nid]) {
                node_bins[ub[nid]].emplace_back(nid);
                continue;
            }
            // checkcore
            for (int nei: graph.adj_list[nid]) {
                if (graph.is_core[nid] != -1) break;
                if (graph.get_similairty(nid, nei) != -1)continue;

                test_edge_pscan(nid, nei);
            }
            if (graph.is_core[nid] == 1) {

                for (int nei: graph.d_neighbors[nid]) {
                    int is_sim = graph.get_similairty(nid, nei);
                    if (is_sim == -1 || cmp_double(graph.edge_weight[nid].at(nei), config.distance) == 1) continue;
                    if (graph.is_core[nei] == 1 && is_sim == 1) {
                        my_union(parent, rank, nid, nei);
                    }
                }
                for (int nei: graph.d_neighbors[nid]) {
                    int is_sim = graph.get_similairty(nid, nei);
                    if (is_sim != -1 || cmp_double(graph.edge_weight[nid].at(nei), config.distance) == 1) continue;
                    if (ub[nei] >= config.mu && find_root(parent, nid) != find_root(parent, nei)) {
                        test_edge_pscan(nid, nei);
                        if (graph.is_core[nei] == 1 && graph.get_similairty(nid, nei) == 1) {
                            my_union(parent, rank, nid, nei);
                        }
                    }
                }
            }
        }
    }
    assign_clusterid(parent);
    swap(graph.d_neighbors, graph.adj_list);
}