#include "../include/mylib.h"
#include "../include/config.h"

Config config;
Result result;

vector<unordered_map<int, bool>> sim_res1;
vector<unordered_map<int, bool>> sim_res2;

bool exists_test(const std::string &name) {
    ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    }
    else {
        f.close();
        return false;
    }
}

void assert_file_exist(string desc, string name) {

    if (!exists_test(name)) {
        cerr << desc << " " << name << " not find " << endl;
        exit(1);
    }
}

