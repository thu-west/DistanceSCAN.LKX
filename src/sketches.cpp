
#include "sketches.h"

vector<unsigned long long> SKETCHES::pick_random_coffs(int k, const unsigned long long &P) {
    vector<unsigned long long> coffs;
    while (k--) {
        coffs.emplace_back(rand_ulong(0, P - 1));
    }
    return coffs;
}

int
SKETCHES::double_hashing(const vector<unsigned long long int> &coff_a, const vector<unsigned long long int> &coff_b,
                         int key, int i) {
    // 2-universal hash h(key) = (a * key + b) % P % n.
    unsigned long long int h_0 = (coff_a[0] * key + coff_b[0]) % P % (2 * graph.n);
    unsigned long long int h_1 = (coff_a[1] * key + coff_b[1]) % P % (2 * graph.n);
    //double hashing  h(key, i) = (h1(key) + i * h2(key)) mod n
    return int((h_0 + i * h_1) % (2 * graph.n));
}

void SKETCHES::init_hash_table() {
    vector<unsigned long long> coff_a, coff_b;
    coff_a = pick_random_coffs(2, P);
    coff_b = pick_random_coffs(2, P);

#ifdef _DEBUG_
    coff_a = {409230109, 2524636865};
    coff_b = {92027918, 2638940465};
#endif
    key2value = vector<int>(graph.n, -1);
    value2key = vector<int>(2 * graph.n, -1);

    //open address hash
    for (int key = 0; key < key2value.size(); ++key) {
        for (int j = 0; j < value2key.size(); ++j) {
            int value = double_hashing(coff_a, coff_b, key, j);
            if (value2key[value] == -1) {
                value2key[value] = key;
                key2value[key] = value;
                break;
            }
        }
    }
}

void SKETCHES::construct_sketches() {
    Timer timer(CONSTRUCT_SKETCH_TIME);
    // init
    INFO("constructing sketches...");

    string prefix = "dmax_" + to_str(config.max_distance) + "_k_" + to_string(config.hash_k);
    string idx;
    if (config.algo == BOTK_SCAN) {
        idx = config.graph_location + prefix + "_multi_ads_bot.idx";
    } else if (config.algo == MY_ADS) {
        idx = config.graph_location + prefix + "_my_ads.idx";
    }

    std::ofstream ofs(idx);
    boost::archive::binary_oarchive oa(ofs);

    int num_bins = get_bin_id(config.max_distance, false) + 1;
    histogram = vector<vector<int>>(graph.n, vector<int>(num_bins));
    //ads_bot = vector<Treap>(graph.n, Treap());
    vector<idpair> dis_source_vec = vector<idpair>(graph.n, make_pair(-1, -1));
    for (int source_nid = 0; source_nid < graph.n; ++source_nid) {
#ifdef _DEBUG_
        if (source_nid % 100 == 99) {
            INFO(source_nid);
        }
        int counter =0;
#endif
        Treap treap;
        ADS ads;
        if (config.algo == BOTK_SCAN) {
            treap = Treap();
        } else if (config.algo == MY_ADS) {
            ads = ADS();
        }
        idpair source = make_pair(source_nid, 0);
        dis_source_vec[source_nid] = source;
        priority_queue<idpair, vector<idpair>, cmp_idpair> pq;
        pq.push(source);
        while (!pq.empty()) {
            auto cur_node = pq.top();
            pq.pop();
            if (cmp_double(cur_node.second, dis_source_vec[cur_node.first].second) == 1)continue;

            update_histogram(source_nid, cur_node.second);
            if (config.algo == BOTK_SCAN) {
                treap.insert_botk(key2value[cur_node.first], cur_node.second);
            } else if (config.algo == MY_ADS) {
                ads.insert_botk(key2value[cur_node.first], cur_node.second);
            }
            //ads_bot[source_nid].insert_botk(key2value[cur_node.first], cur_node.second);
#ifdef _DEBUG_
            assert(dis_source_vec[cur_node.first].first == source_nid);
            assert(cmp_double(cur_node.second, dis_source_vec[cur_node.first].second) == 0);
            counter++;
#endif
            for (int nei_id: graph.adj_list[cur_node.first]) {
                double nei_dis = cur_node.second + graph.edge_weight[cur_node.first][nei_id];
                if (cmp_double(nei_dis, config.max_distance) < 1) {
                    if (dis_source_vec[nei_id].first != source_nid ||
                        cmp_double(dis_source_vec[nei_id].second, nei_dis) > 0) {
                        dis_source_vec[nei_id].first = source_nid;
                        dis_source_vec[nei_id].second = nei_dis;
                        pq.push(make_pair(nei_id, nei_dis));
                    }
                }
            }
        }
        {
            Timer timer(SAVE_BOTK_TIME);
            if (config.algo == BOTK_SCAN) {
                oa << treap;
            } else if (config.algo == MY_ADS) {
                oa << ads;
            }
        }
        {
            Timer timer(UPDATE_HISTO_TIME);
            for (int j = 1; j < histogram[source_nid].size(); ++j) {
                histogram[source_nid][j] += histogram[source_nid][j - 1];
            }
        }

    }
}

int SKETCHES::update_histogram(int source_nid, double path_weight) {
    int bin_id = get_bin_id(path_weight, false);
    histogram[source_nid][bin_id]++;
    return bin_id;
}

int SKETCHES::get_bin_id(double dis, bool is_lb) {
    int bin_id = floor(dis / config.bin);
    if (is_lb) {
        bin_id--;
    }
    return bin_id;
}

int SKETCHES::query_histogram(int source_nid, double dis, bool is_lb) {
    int bin_id = get_bin_id(dis, is_lb);
    if (bin_id < 0) {
        return 0;
    } else {
        return histogram[source_nid][bin_id];
    }
}

void SKETCHES::deserialize_sketches() {
    string prefix = "dmax_" + to_str(config.max_distance) + "_k_" + to_string(config.hash_k);
    string file_name = config.graph_location + prefix + "_multi_ads_bot.idx";
    if (!exists_test(file_name)) {
        file_name = config.graph_location + "dmax_" + to_str(config.max_distance) + "_k_65536_multi_ads_bot.idx";
    }
    INFO(file_name);
    assert_file_exist("index file", file_name);
    std::ifstream ifs(file_name);
    boost::archive::binary_iarchive ia(ifs);
    Treap tmp_treap;
    for (int i = 0; i < graph.n; ++i) {
        ia >> tmp_treap;
        ads_bot.emplace_back(tmp_treap);
    }

    file_name = config.graph_location + prefix + "_histogram.idx";
    assert_file_exist("index file", file_name);
    std::ifstream info_ifs(file_name);
    boost::archive::binary_iarchive info_ia(info_ifs);
    info_ia >> histogram >> key2value >> value2key;
}

void SKETCHES::deserialize_sketches2() {
    string prefix = "dmax_" + to_str(config.max_distance) + "_k_" + to_string(config.hash_k);
    string file_name;
    if (config.algo == BOTK_SCAN) {
        file_name = config.graph_location + prefix + "_multi_ads_bot.idx";
        if (!exists_test(file_name)) {
            file_name = config.graph_location + "dmax_" + to_str(config.max_distance) + "_k_65536_multi_ads_bot.idx";
        }
    } else if (config.algo == MY_ADS) {
        file_name = config.graph_location + prefix + "_my_ads.idx";
        if (!exists_test(file_name)) {
            file_name = config.graph_location + "dmax_" + to_str(config.max_distance) + "_k_65536_my_ads.idx";
        }
    }

    INFO(file_name);
    assert_file_exist("index file", file_name);
    std::ifstream ifs(file_name);
    boost::archive::binary_iarchive ia(ifs);
    Treap tmp_treap;
    ADS ads;
    bot_k = vector<vector<int>>(graph.n, vector<int>{});
    if (config.algo == MY_ADS) {
        bot_k_dis = vector<vector<double>>(graph.n, vector<double>{});
    }
    for (int i = 0; i < graph.n; ++i) {
        if (config.algo == BOTK_SCAN) {
            ia >> tmp_treap;
            bot_k[i] = tmp_treap.get_bot_k();
        } else if (config.algo == MY_ADS) {
            ia >> ads;
            bot_k[i] = ads.get_bot_k(bot_k_dis[i]);
        }
    }
    file_name = config.graph_location + prefix + "_hashmap.idx";
    if (!exists_test(file_name)) {
        file_name = config.graph_location + "dmax_" + to_str(config.max_distance) + "_k_65536_hashmap.idx";
    }
    INFO(file_name);
    std::ifstream info_ifs2(file_name);
    boost::archive::binary_iarchive info_ia2(info_ifs2);
    info_ia2 >> key2value >> value2key;

    prefix += "_bin_" + to_string(config.bin);
    file_name = config.graph_location + prefix + "_histogram.idx";
    if (!exists_test(file_name)) {
        file_name = config.graph_location + "dmax_" + to_str(config.max_distance) + "_k_65536_bin_" +
                    to_string(config.bin) + "_histogram.idx";
    }
    INFO(file_name);
    assert_file_exist("index file", file_name);
    std::ifstream info_ifs(file_name);
    boost::archive::binary_iarchive info_ia(info_ifs);
    info_ia >> histogram;
    INFO(histogram.size(), key2value.size(), value2key.size());
}

void SKETCHES::serialize_sketches() {
    string prefix = "dmax_" + to_str(config.max_distance) + "_k_" + to_string(config.hash_k);
    string idx = config.graph_location + prefix + "_hashmap.idx";
    std::ofstream info_ofs2(idx);
    boost::archive::binary_oarchive info_oa2(info_ofs2);
    info_oa2 << key2value << value2key;
    prefix += +"_bin_" + to_string(config.bin);
    idx = config.graph_location + prefix + "_histogram.idx";
    std::ofstream info_ofs(idx);
    boost::archive::binary_oarchive info_oa(info_ofs);
    info_oa << histogram;
}


void SKETCHES::get_approx_neis() {
    //neis_in_dis = vector<double>(graph.n);
    neis_in_dis_lb = vector<int>(graph.n);
    neis_in_dis_ub = vector<int>(graph.n);
    int bin_id = get_bin_id(config.distance, false);
    for (int i = 0; i < graph.n; ++i) {
        neis_in_dis_lb[i] = query_histogram(i, config.distance, true);
        if (bot_k[i].size() < config.hash_k) {
            neis_in_dis_ub[i] = bot_k[i].size();
        } else {
            neis_in_dis_ub[i] = query_histogram(i, config.distance, false);
        }
    }
}

double SKETCHES::jaccard_with_botk(int u, int v) {
    int u_cur = 0, v_cur = 0, union_size = 0, inter_size = 0;
    while (u_cur != bot_k[u].size() && v_cur != bot_k[v].size() && union_size < config.hash_k) {
        union_size++;
        if (bot_k[u][u_cur] < bot_k[v][v_cur]) {
            u_cur++;
        } else if (bot_k[u][u_cur] > bot_k[v][v_cur]) {
            v_cur++;
        } else {
            inter_size++;
            u_cur++;
            v_cur++;
        }
    }
    if (union_size < config.hash_k - 0.1) {
        union_size += bot_k[u].size() - u_cur + bot_k[v].size() - v_cur;
    }
    union_size = min(config.hash_k, union_size);
    return double(inter_size) / union_size;
}

double SKETCHES::jaccard_with_sketches(int u, int v, double dis) {
#ifdef  _DEBUG_
    assert(dis > -0.5);
#endif
    if (config.sketches_optimize > 1) {
        int u_ub_size = neis_in_dis_ub[u], u_lb_size = neis_in_dis_lb[u];
        int v_ub_size = neis_in_dis_ub[v], v_lb_size = neis_in_dis_lb[v];
        //int u_ub_size = query_histogram(u, config.distance, false), u_lb_size = query_histogram(u, config.distance, true);
        //int v_ub_size = query_histogram(v, config.distance, false), v_lb_size = query_histogram(v, config.distance, true);
        double sim_ub = 0;
        if (u_ub_size < v_lb_size) {
            sim_ub = max(sim_ub, double(u_ub_size) / v_lb_size);
        }
        if (v_ub_size < u_lb_size) {
            sim_ub = max(sim_ub, double(v_ub_size) / u_lb_size);
        }
        if (sim_ub != 0 && cmp_double(sim_ub, config.epsilon) == -1) {
            return false;
        }
        int u_intersection = query_histogram(u, config.distance - dis, true);
        int v_intersection = query_histogram(v, config.distance - dis, true);
        if (u_intersection < v_intersection) {
            swap(u, v);
            swap(u_intersection, v_intersection);
        }
        double sim_lb = double(u_intersection) / (u_ub_size + v_ub_size - u_intersection);
        if (cmp_double(sim_lb, config.epsilon) > -1) {
            return true;
        }
    }
    if (config.sketches_optimize > 0) {
        int u_cur = 0, v_cur = 0, union_size = 0, inter_size = 0;
        while (u_cur != bot_k[u].size() && v_cur != bot_k[v].size() && union_size < config.hash_k) {
            union_size++;
            if (bot_k[u][u_cur] < bot_k[v][v_cur]) {
                u_cur++;
            } else if (bot_k[u][u_cur] > bot_k[v][v_cur]) {
                v_cur++;
            } else {
                inter_size++;
                u_cur++;
                v_cur++;
            }
            if (inter_size > intersection_ub) {
                return true;
            } else if (inter_size < intersection_lb[union_size]) {
                return false;
            }
        }
        if (union_size < config.hash_k) {
            union_size += bot_k[u].size() - u_cur + bot_k[v].size() - v_cur;
        }
        union_size = min(config.hash_k, union_size);
        return cmp_double(double(inter_size) / union_size, config.epsilon) > -1;
    } else {
        return cmp_double(jaccard_with_botk(u, v), config.epsilon) > -1;
    }
}
