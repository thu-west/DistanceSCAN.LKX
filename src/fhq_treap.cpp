#include "fhq_treap.h"

Treap::Treap() {

}

void Treap::update(int pos) {
    if (pos >= 0) {
        tnodes[pos].sz = 1;
        int left_id = tnodes[pos].left, right_id = tnodes[pos].right;
        if (left_id >= 0) {
            tnodes[pos].sz += tnodes[left_id].sz;
        }
        if (right_id >= 0) {
            tnodes[pos].sz += tnodes[right_id].sz;
        }
    }
}

pair<int, int> Treap::split(int tn, int k) {
    if (tn < 0) {
        return make_pair(-1, -1);
    }
    pair<int, int> o;
    if (k < tnodes[tn].key) {
        o = split(tnodes[tn].left, k);
        tnodes[tn].left = o.second;
        update(tn);
        return make_pair(o.first, tn);
    } else {
        o = split(tnodes[tn].right, k);
        tnodes[tn].right = o.first;
        update(tn);
        return make_pair(tn, o.second);
    }
}

int Treap::merge(int x, int y) {
    if (x < 0) {
        return y;
    }
    if (y < 0) {
        return x;
    }
    if (tnodes[x].priority > tnodes[y].priority) {
        tnodes[x].right = merge(tnodes[x].right, y);
        update(x);
        return x;
    } else {
        tnodes[y].left = merge(x, tnodes[y].left);
        update(y);
        return y;
    }
}

pair<int, int> Treap::split_persist(int tn, int k) {
    if (tn < 0) {
        return make_pair(-1, -1);
    } else {
        int ntn_id = (int) tnodes.size();
        tnodes.emplace_back(tnodes[tn]);
        pair<int, int> ret;
        if (tnodes[tn].key <= k) {
            ret = split_persist(tnodes[ntn_id].right, k);
            tnodes[ntn_id].right = ret.first;
            ret.first = ntn_id;
        } else {
            ret = split_persist(tnodes[ntn_id].left, k);
            tnodes[ntn_id].left = ret.second;
            ret.second = ntn_id;
        }
        update(ntn_id);
        return ret;
    }
}

void Treap::insert(int ins_id) {
    pair<int, int> o = split(roots.back(), tnodes[ins_id].key);
    o.first = merge(o.first, ins_id);
    roots.back() = merge(o.first, o.second);
}

void Treap::insert_persist(int ins_id) {
    pair<int, int> o = split_persist(roots.back(), tnodes[ins_id].key);
    o.first = merge(o.first, ins_id);
    roots.emplace_back(merge(o.first, o.second));
}

void Treap::update_persist(int max_key, int value) {
    pair<int, int> o = split_persist(roots.back(), max_key - 1);
    pair<int, int> p = split_persist(o.second, max_key);
#ifdef _DEBUG_
    assert(p.second == -1);
#endif
    int upd_nid = p.first;
    tnodes[upd_nid].left = -1;
    tnodes[upd_nid].right = -1;
    tnodes[upd_nid].key = value;
    tnodes[upd_nid].priority = rand();
    tnodes[upd_nid].sz = 1;
    pair<int, int> q = split_persist(o.first, value);
    o.first = merge(q.first, upd_nid);
    roots.emplace_back(merge(o.first, q.second));
}

int Treap::maximum(int root_id) {
    int max_id = roots[root_id];
    while (tnodes[max_id].right >= 0)
        max_id = tnodes[max_id].right;
    return tnodes[max_id].key;
}


void Treap::insert_botk(int value, double path_w) {
    Timer timer(ALLD_BOTK_TIME);
    int new_node_id = (int) tnodes.size();
    if (roots.empty()) {
        tnodes.emplace_back(TreapNode(value));
        roots.emplace_back(new_node_id);
        dis_vec.emplace_back(path_w);
        max_value.emplace_back(value);
    } else if (tnodes[roots.back()].sz < config.hash_k) {
        tnodes.emplace_back(TreapNode(value));
        insert_persist(new_node_id);
        dis_vec.emplace_back(path_w);
        max_value.emplace_back(max(value, max_value.back()));

    } else {
        //int max_key = maximum(int(roots.size()) - 1);
        if (max_value.back() < value) return;
        update_persist(max_value.back(), value);
        max_value.emplace_back(maximum(int(roots.size()) - 1));
        dis_vec.emplace_back(path_w);
#ifdef _DEBUG_
        assert(tnodes[roots.back()].sz == config.hash_k);
#endif
    }

#ifdef _DEBUG_
    assert(max_value.size() == dis_vec.size() && max_value.size() == roots.size());
#endif
}

void Treap::delete_botk(double dis) {
    while (!dis_vec.empty()) {
        if (dis_vec.back() >= dis) {
            dis_vec.pop_back();
            max_value.pop_back();
            roots.pop_back();
        } else {
            break;
        }
    }
    int maxid = get_maxid_in_tree(roots.back());
    while (tnodes.size() > maxid + 1) {
        tnodes.pop_back();
    }
}

int Treap::get_maxid_in_tree(int root_id) {
    int maxid = root_id;
    if (tnodes[root_id].left != -1) {
        maxid = max(maxid, get_maxid_in_tree(tnodes[root_id].left));
    }
    if (tnodes[root_id].right != -1) {
        maxid = max(maxid, get_maxid_in_tree(tnodes[root_id].right));
    }
    return maxid;
}

int Treap::get_root_id(double dis) {
    int left = 0, right = dis_vec.size(), mid;
    if (dis_vec.empty() || cmp_double(dis_vec[0], dis) == 1) {
        return -1;
    } else if (cmp_double(dis_vec.back(), dis) < 1) {
        return right - 1;
    }
    while (left < right) {
        mid = floor((left + right) / 2);

        if (cmp_double(dis_vec[mid], dis) < 1) {
#ifdef _DEBUG_
            assert(mid + 1 < dis_vec.size());
#endif
            if (cmp_double(dis_vec[mid + 1], dis) == 1) {
                break;
            } else {
                left = mid;
            }
        } else {
            right = mid;
        }
    }
    return mid;

}

int Treap::get_tau_k(double dis) {
    int root_id = get_root_id();
    if (root_id == -1) return -1;
    return maximum(root_id);
}

vector<int> Treap::get_bot_k(double dis) {
    Timer timer(GET_BOTK_TIME);
    int root_id = get_root_id();
    if (root_id == -1) return vector<int>();
    return inOrder(roots[root_id]);
}

void Treap::_inOrder(int tn_id, vector<int> & ret) {
    if (tn_id >= 0) {
        _inOrder(tnodes[tn_id].left,ret);
        if (ret.size() >= config.hash_k) return;
        ret.emplace_back(tnodes[tn_id].key);
        _inOrder(tnodes[tn_id].right,ret);
    }
    return;
}

vector<int> Treap::inOrder(int tn_id) {
    int tree_sz = tnodes[tn_id].sz;
#ifdef _DEBUG_
    vector<int> ret = _inOrder(tn_id);
    assert(ret.size() == tree_sz);
    return ret;
#endif
    vector<int> ret;
    _inOrder(tn_id,ret);
    return ret;
}

void Treap::print(int tn_id, int key, int direction) {
    if (tn_id >= 0) {
        if (direction == 0)
            cout << setw(2) << tnodes[tn_id].key << " is root" << endl;
        else
            cout << setw(2) << tnodes[tn_id].key << " is " << setw(2) << key << "'s " << setw(12)
                 << (direction == 1 ? "right child" : "left child") << endl;

        print(tnodes[tn_id].left, tnodes[tn_id].key, -1);
        print(tnodes[tn_id].right, tnodes[tn_id].key, 1);
    }
}

void Treap::print() {
    for (int i = 0; i < roots.size(); ++i) {
        print(roots[i], tnodes[roots[i]].key, 0);
    }
}