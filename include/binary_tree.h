
#ifndef DISTANCE_SCAN_SIGMOD_BINARY_TREE_H
#define DISTANCE_SCAN_SIGMOD_BINARY_TREE_H
#include "mylib.h"
#include "config.h"
#include <iomanip>
#include <iostream>
using namespace std;

class ADS{
private:
    vector<idpair> ads_nodes;
    priority_queue<int> cur_botk;
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & ads_nodes;
    }

public:
    ADS();

    void insert_botk(int value, double path_w);

    vector<int> get_bot_k(double dis=config.distance);
    vector<int> get_bot_k(vector<double> &dis_vec, double dis=config.distance);

    virtual ~ADS();
};

#endif //DISTANCE_SCAN_SIGMOD_BINARY_TREE_H
