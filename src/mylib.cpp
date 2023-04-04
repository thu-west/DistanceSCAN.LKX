
#include "../include/mylib.h"


template<class T>
string toStr(T t) {
    stringstream ss;
    ss << t;
    return ss.str();
}

string __n_variable(string t, int n) {
    t = t + ',';
    int i = 0;
    if (n) for (; i < SIZE(t) && n; i++) if (t[i] == ',') n--;
    n = i;
    for (; t[i] != ','; i++);
    t = t.substr((unsigned long) n, (unsigned long) (i - n));
    trim(t);
    if (t[0] == '"') return "";
    return t + "=";
}

int Counter::cnt[1000] = {0};

vector<double> Timer::timeUsed;
vector<string> Timer::timeUsedDesc;


bool cmp_pair(const pair<int, double> &p1, const pair<int, double> &p2) {

    if (p1.second < p2.second) return true;
    if (p1.second > p2.second) return false;
    return p1.first < p2.first;
    //return p1.second < p2.second;
}

unsigned long long rand_ulong(unsigned long long min, unsigned long long max) {
    static default_random_engine re(time(0));
    using Dist = uniform_int_distribution<unsigned long long>;
    static Dist uid {};
    return uid(re, Dist::param_type{min,max});
}
void split_string(const std::string& input_str, std::vector<std::string>& output, const char* delim)
{
    int pos = 0;
    int npos = 0;
    int regexlen = strlen(delim);
    while((npos = input_str.find(delim, pos)) != -1) {
        std::string tmp = input_str.substr(pos, npos - pos);
        output.push_back(tmp);
        pos = npos + regexlen;
    }
    output.push_back(input_str.substr(pos, input_str.length() - pos));
}

int cmp_double(const double &d1, const double & d2){
    if (d1-d2>1e-6){
        return 1;
    }else if (d2-d1>1e-6){
        return -1;
    }else{
        return 0;
    }
}