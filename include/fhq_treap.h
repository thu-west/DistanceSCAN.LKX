
#ifndef WEIGHTED_SCAN_FHQ_TREAP_H
#define WEIGHTED_SCAN_FHQ_TREAP_H
#include "mylib.h"
#include "config.h"
using namespace std;
class TreapNode{
public:
    int key = -1;
    int left = -1,right = -1;

    int priority=0;
    int sz = 1;
    TreapNode(){};
    TreapNode(int value): key(value){
        priority = rand();
    }

    TreapNode(int value, int l, int r):
            key(value),left(l),right(r) {
        priority = rand();
    }

    TreapNode(const TreapNode & tn){
        key = tn.key;
        left = tn.left;
        right = tn.right;
        priority = tn.priority;
        sz = tn.sz;
    }
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & key & left & right;
    }
};

class Treap {
private:
    vector<TreapNode> tnodes;
    vector<int> max_value;
    vector<double> dis_vec;
    vector<int> roots;
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & tnodes & dis_vec & roots;
    }
public:
    Treap();

    void insert_botk(int value, double path_w);

    int get_tau_k(double dis=config.distance);

    vector<int> get_bot_k(double dis=config.distance);

    void delete_botk(double dis);

    vector<int> inOrder(int tn_id);


    int maximum(int root_id);
    int get_maxid_in_tree(int root_id);

    void destroy();

    void print();
private:
    void insert_persist(int ins_id);

    void update_persist(int max_key, int value);

    int get_root_id(double dis=config.distance);

    void update(int pos);

    pair<int,int> split_persist(int tn, int k);

    pair<int, int> split(int tn, int k);

    int merge(int x, int y);

    void insert(int ins_id);


    void _inOrder(int tn_id, vector<int> & ret);


    void destroy(TreapNode* &tree);

    void print(int tn_id, int key, int direction) ;
};

#endif //WEIGHTED_SCAN_FHQ_TREAP_H
