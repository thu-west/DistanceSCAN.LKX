#include "binary_tree.h"

void ADS::insert_botk(int value, double path_w) {
    Timer timer(ALLD_BOTK_TIME);
    if (ads_nodes.size() < config.hash_k) {
        ads_nodes.emplace_back(make_pair(value, path_w));
        cur_botk.emplace(value);
    } else {
        if(cur_botk.top()<value) return;
        ads_nodes.emplace_back(make_pair(value, path_w));
        cur_botk.pop();
        cur_botk.emplace(value);
    }
}
vector<int> ADS::get_bot_k(double dis) {
    priority_queue<int> pq;
    for (idpair item: ads_nodes) {
        if (cmp_double(item.second , dis)==1) break;
        int value = item.first;
        if (pq.size() < config.hash_k) {
            pq.emplace(value);
        } else if (pq.top() > value) {
            pq.pop();
            pq.emplace(value);
        }
    }
    vector<int> ret(pq.size());
    for (int i = pq.size() - 1; i > -1; --i) {
        ret[i] = pq.top();
        pq.pop();
    }
    return ret;
}

vector<int> ADS::get_bot_k(vector<double> &dis_vec, double dis) {
    Timer timer(GET_BOTK_TIME);
    priority_queue<idpair> pq;
    for (idpair item: ads_nodes) {
        if (cmp_double(item.second , dis)==1) break;
        if (pq.size() < config.hash_k) {
            pq.emplace(item);
        } else if (pq.top().first > item.first) {
            pq.pop();
            pq.emplace(item);
        }
    }
    vector<int> ret(pq.size());
    dis_vec = vector<double>(pq.size());
    for (int i = pq.size() - 1; i > -1; --i) {
        ret[i] = pq.top().first;
        dis_vec[i] = pq.top().second;
        pq.pop();
    }
    return ret;
}

ADS::~ADS() {

}

ADS::ADS() {
    ads_nodes = vector<idpair>();
    //cur_botk = BSTree<int>();
}
