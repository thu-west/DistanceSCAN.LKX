#ifndef SCAN_MYLIB_H
#define SCAN_MYLIB_H

//#define _DEBUG_
#include <iostream>
#include <set>
#include <list>
#include <sstream>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdio>
#include <functional>
#include <algorithm>
#include <climits>
#include <cstring>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <map>
#include <queue>
#include <atomic>
#include <thread>
#include <future>
#include <queue>
#include <mutex>
#include <sys/resource.h>
#include <condition_variable>
/*
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
*/

#include <random>
#include <sys/stat.h>

using namespace std;
typedef unsigned int uint;
typedef unsigned char uint8;
typedef long long int64;
typedef unsigned long long uint64;
typedef pair<int, int> ipair;
typedef pair<double, double> dpair;
typedef pair<int, double> idpair;
/*
#define MP make_pair
#define F first
#define S second
*/
#ifndef TIMES_PER_SEC
#define TIMES_PER_SEC (1.0e9)
#endif

typedef char int8;
typedef unsigned char uint8;
typedef long long int64;
typedef unsigned long long uint64;

#define SIZE(t) (int)(t.size())
#define ALL(t) (t).begin(), (t).end()
#define FOR(i, n) for(int (i)=0; (i)<((int)(n)); (i)++)


const int MAXN = 1000000;

const int VectorDefaultSize = 20;
const int TOPNUM = 1;

inline static double drand() {
    return rand() * 1.0f / RAND_MAX;
}


static inline std::string &ltrim(std::string &str) {
    auto it2 = std::find_if(str.begin(), str.end(),
                            [](char ch) { return !std::isspace<char>(ch, std::locale::classic()); });
    str.erase(str.begin(), it2);
    return str;
}

static inline std::string &rtrim(std::string &str) {
    auto it1 = std::find_if(str.rbegin(), str.rend(),
                            [](char ch) { return !std::isspace<char>(ch, std::locale::classic()); });
    str.erase(it1.base(), str.end());
    return str;
}

static inline string &trim(string &s) { return ltrim(rtrim(s)); }

string __n_variable(string t, int n);

#define __expand_nv(x) __n_variable(t, x)<< t##x << " "


template<class T0>
void ___debug(string t, T0 t0, ostream &os) {
    os << __expand_nv(0);
}

template<class T0, class T1>
void ___debug(string t, T0 t0, T1 t1, ostream &os) {
    os << __expand_nv(0) << __expand_nv(1);
}

template<class T0, class T1, class T2>
void ___debug(string t, T0 t0, T1 t1, T2 t2, ostream &os) {
    os << __expand_nv(0) << __expand_nv(1) << __expand_nv(2);
}

template<class T0, class T1, class T2, class T3>
void ___debug(string t, T0 t0, T1 t1, T2 t2, T3 t3, ostream &os) {
    os << __expand_nv(0) << __expand_nv(1) << __expand_nv(2) << __expand_nv(3);
}

template<class T0, class T1, class T2, class T3, class T4>
void ___debug(string t, T0 t0, T1 t1, T2 t2, T3 t3, T4 t4, ostream &os) {
    os << __expand_nv(0) << __expand_nv(1) << __expand_nv(2) << __expand_nv(3) << __expand_nv(4);
}

template<class T0>
void ___debug(string t, deque<T0> t0, ostream &os) {
    os << __n_variable(t, 0);
    FOR(i, SIZE(t0))os << t0[i] << " ";
}

template<class T0>
void ___debug(string t, vector<T0> t0, ostream &os) {
    os << __n_variable(t, 0);
    FOR(i, SIZE(t0))os << t0[i] << " ";
}

template<class T0, class T1>
void ___debug(string t, vector<pair<T0, T1> > t0, ostream &os) {
    os << __n_variable(t, 0);
    FOR(i, SIZE(t0))os << t0[i].F << "," << t0[i].S << " ";
}

#define RUN_TIME(...) { int64 t=rdtsc();  __VA_ARGS__; t=rdtsc()-t; cout<<  #__VA_ARGS__ << " : " << t/TIMES_PER_SEC <<"s"<<endl;  }

#ifdef HEAD_TRACE
#define TRACE(...) {{ ___debug( #__VA_ARGS__,  __VA_ARGS__,cerr); cerr<<endl;  } }
#define IF_TRACE(args) args
#define TRACE_LINE(...) { ___debug( #__VA_ARGS__,  __VA_ARGS__,cerr); cerr<<"                    \033[100D";  }
#define TRACE_SKIP(a, ...) { static int c=-1; c++; if(c%a==0)TRACE( __VA_ARGS__); }
#define TRACE_LINE_SKIP(a, ...) { static int c=-1; c++; if(c%a==0) TRACE_LINE(__VA_ARGS__);  }
#define TRACE_LINE_END(...) {cerr<<endl; }
ofstream __HEAD_H__LOG("log.txt");
#define TRACE_LOG(...) { __HEAD_H__LOG.close(); ofstream cerr("log.txt", ofstream::out|ofstream::app); ___debug( #__VA_ARGS__,  __VA_ARGS__, cerr); cerr<<endl;  }
#else
#define TRACE(...) ;
#define IF_TRACE(args) ;
#define TRACE_LINE(...) ;
#define TRACE_SKIP(a, ...) ;
#define TRACE_LINE_SKIP(a, ...) ;
#define TRACE_LINE_END(...) ;
#define TRACE_LOG(...) ;
#endif //HEAD_TRACE


//void setInfoFile(string s) { __HEAD_H_FOUT.open(s.c_str()); }

#define ASSERT(v) {if (!(v)) {cerr<<"ASSERT FAIL @ "<<__FILE__<<":"<<__LINE__<<endl; exit(1);}}
#define INFO(...) do {\
     ___debug( #__VA_ARGS__,  __VA_ARGS__,cout); cout<<endl; \
    } while(0)
#define ASSERTT(v, ...) {if (!(v)) {cerr<<"ASSERT FAIL @ "<<__FILE__<<":"<<__LINE__<<endl; INFO(__VA_ARGS__); exit(1);}}

template<class T>
string toStr(T t);


class Counter {
public:
    static int cnt[1000];
    int myid = 0;

    Counter(int id = 0) {
        myid = id;
        cnt[id]++;
    }

    void add(int x) {
        cnt[myid] += x;
    }

    ~Counter() {
    }

    static void show() {
        for (int i = 0; i < 1000; i++)
            if (cnt[i] > 0)
                INFO("Counter", i, cnt[i]);
    }
};


uint64 rdtsc();


class Timer {
public:
    static vector<double> timeUsed;
    static vector<string> timeUsedDesc;
    int id;
    std::chrono::steady_clock::time_point startTime;
    bool showOnDestroy;

    Timer(int id, string desc = "", bool showOnDestroy = false) {
        this->id = id;
        while ((int) timeUsed.size() <= id) {
            timeUsed.push_back(0);
            timeUsedDesc.push_back("");
        }
        timeUsedDesc[id] = desc;
        startTime = std::chrono::steady_clock::now();
        this->showOnDestroy = showOnDestroy;
    }

    static double used(int id) {
        return timeUsed[id] / TIMES_PER_SEC;
    }

    ~Timer() {
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - startTime).count();
        if (showOnDestroy) {
            cout << "time spend on " << timeUsedDesc[id] << ":" << duration / TIMES_PER_SEC << "s" << endl;
        }
        timeUsed[id] += duration;
    }

    static void show(bool debug = false) {
        if (debug) { TRACE("### Timer");
        } else {
            INFO("### Timer");
        }
        for (int i = 0; i < (int) timeUsed.size(); i++) {
            if (timeUsed[i] > 0) {
                char str[100];
                sprintf(str, "%.6lf", timeUsed[i] / TIMES_PER_SEC);
                string s = str;
                if ((int) s.size() < 15) s = " " + s;
                char t[100];
                memset(t, 0, sizeof t);
                sprintf(t, "%4d %s %s", i, s.c_str(), timeUsedDesc[i].c_str());
                if (debug) { TRACE(t);
                } else {
                    INFO(t);
                }
            }
        }
    }

    static void reset(int id) {
        timeUsed[id] = 0;
    }

    static void clearAll() {
        timeUsed.clear();
        timeUsedDesc.clear();
    }
};

static vector<string> combine_args(int argc, char **argv) {
    vector<string> args;
    for (int i = 0; i < argc; ++i) {
        args.push_back(argv[i]);
    }
    return args;
}

static string to_str(double t) {
    stringstream ss;
    ss << t;
    return ss.str();
}

inline bool file_exists_test(const std::string &name) {
    ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }
}

static string replace(std::string str, const std::string &from, const std::string &to) {
    //original string will not be modified
    size_t start_pos = str.find(from);
    if (start_pos == std::string::npos)
        return str;
    str.replace(start_pos, from.length(), to);
    return str;
}

// unit is kilobyte (KB)
inline long long get_proc_memory() {
    struct rusage r_usage;
    getrusage(RUSAGE_SELF, &r_usage);
    return r_usage.ru_maxrss;
}

const string Green = "\033[0;32m";
const string Reset = "\033[0m";
const string Red = "\033[0;31m";


static string get_current_time_str() {
    time_t rawtime;
    struct tm *timeinfo;
    char buffer[80];

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", timeinfo);
    std::string str(buffer);

    return str;
}

static void program_start(int argc, char **argv) {

    cout << Green << "--------------start------------" << get_current_time_str() << Reset << endl;
    string combine = "";
    for (int i = 1; i < argc; i++) {
        combine += argv[i];
        combine += " ";
    }
    cout << Green << "args:" << combine << Reset << endl;
}

static void program_stop() {
    cout << Red << "--------------stop------------" << get_current_time_str() << Reset << endl;
    cout << endl;
    cout << endl;
    cout << endl;

}

#ifndef NDEBUG
#   define ASSERTMSG(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::exit(EXIT_FAILURE); \
        } \
    } while (false)
#else
#   define ASSERTMSG(condition, message) do { } while (false)
#endif

struct cmp_idpair{
    bool operator()(idpair p1,idpair p2){
        if (p1.second > p2.second) return true;
        if (p1.second < p2.second) return false;
        return p1.first < p2.first;
    }
};

bool cmp_pair(const pair<int, double> &p1, const pair<int, double> &p2);

int cmp_double(const double &d1, const double & d2);

unsigned long long rand_ulong(unsigned long long min, unsigned long long max);

void split_string(const std::string& input_str, std::vector<std::string>& output, const char* delim);

#endif //SCAN_MYLIB_H
