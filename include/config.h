#ifndef SCAN_CONFIG_H
#define SCAN_CONFIG_H


#ifdef WIN32
#define FILESEP "\\"
#else
#define FILESEP "/"
#endif

//#include <boost/progress.hpp>
#include <getopt.h>
#include <boost/math/distributions/normal.hpp>
#include <boost/timer/timer.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>
#include <unordered_map>
#include <list>


#include <boost/serialization/serialization.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/utility.hpp>
//using namespace boost;
//using namespace boost::property_tree;
//using boost::math::normal;




const string QUERY = "query";
const string BATCH_QUERY = "batch_query";
const string CLUSTER_VALIDATION = "cluster_validation";
const string QUALITY_VALIDATION = "quality_validation";
const string GRAPH_MAINTAIN = "graph_maintain";
const string CONVERT_GRAPH = "convert_graph";
const string CONSTRUCT_SKETCHES = "construct_sketches";
const string GETCLUSTERCOEFFICIENT = "get_cluster_coefficient";
const string EXPONLFR = "EXPONLFR";
const string GENRATE_EDGE_UPDATE = "generate_edge_update";

const string SCAN = "scan";
const string W_SCAN = "w_scan";
const string BASIC = "exact";
const string BOTK_SCAN = "distancescan_pst";
const string PSCAN_DIS = "pscan";
const string MY_ADS = "distancescan";

const string CHECK = "check";

const int BINS_NEIGHBOR = 10;
const double MIN_DISTANCE = 0.1;
//用作计时
const int WHOLE = 0;
const int SCAN_TIME = 1;
const int BASIC_TIME = 2;
const int PSCAN_TIME = 3;
const int SKETCH_TIME = 4;
const int READ_GRAPH_TIME = 5;
const int GET_BOTK_TIME = 6;
const int D_NEIGHBOR_TIME = 7;
const int ALLD_BOTK_TIME = 8;
const int DIS_SCAN_TIME = 9;
const int ALG_TIME = 10;
const int SAVE_BOTK_TIME = 11;
const int UPDATE_HISTO_TIME = 12;
const int CONSTRUCT_SKETCH_TIME = 13;

const int EDGE_DEL_TIME = 14;
const int EDGE_INS_TIME = 15;
const int EDGE_UPDATE_TIME = 16;
const int GRAPH_MAINTAIN_TIME = 17;


const int TRAVERSE_NODES_TIME = 20;
const int TRAVERSE_NODES_TIME_DIJKSTRA = 21;
const int TRAVERSE_NODES_TIME_BOTK = 22;
const int TRAVERSE_NODES_TIME_BOTK_SORT = 23;

const int TRAVERSE_NODES_SORT_TIME = 30;
const int TRAVERSE_EDGES_TIME = 31;
const int TRAVERSE_EDGES_SORT_TIME = 32;


//const int WIGHT_SCAN_TIME = 2;


extern vector<unordered_map < int, bool>>
sim_res1;
extern vector<unordered_map < int, bool>>
sim_res2;

#ifdef WIN32
const string parent_folder = "../../";
#else
const string parent_folder = string("./") + FILESEP;
#endif

class Config {
public:
    //bool use_cos_similarity = false;
    enum SimilarityType{jac,cos,set_containment1,set_containment2};
    SimilarityType similarityType = SimilarityType(0);

    string operation = QUERY;
    string graph_alias;
    string graph_location;

    string prefix = "../kcore_dataset/";
    string outfile = "result.txt";

    string exe_result_dir = parent_folder;

    string get_graph_folder() {
        return prefix + graph_alias + FILESEP;
    }

    double epsilon = 0.2;
    int mu = 5;
    double distance = 0.3;

    double max_distance = 1.4;
    double bin = 0.01;

    int hash_k = 65536;

    double base = -1;
    string algo;
    int type = 0;

    int update_edge_nums = 1000;
    int update_mode = 0;

    int sketches_optimize = 2;

    int algo_type = 6;

    void init_config(int argc, char *argv[]) {
        int opt;
        int digit_optind = 0;

        int option_index = 0;
        const char *optstring = "u:e:k:x:d:m:n:z:s:";

        static struct option long_options[] = {
                {"dataset",          required_argument, NULL, 'g'},
                {"prefix",           required_argument, NULL, 'p'},
                {"operation",        required_argument, NULL, 'o'},
                {"algo",             required_argument, NULL, 'a'},
                {"mu",               required_argument, NULL, 'u'},
                {"epsilon",          required_argument, NULL, 'e'},
                {"distance",         required_argument, NULL, 'd'},
                {"hash_k",           required_argument, NULL, 'k'},
                {"max_distance",     required_argument, NULL, 'm'},
                {"update_edge_nums", required_argument, NULL, 'n'},
                {"bin",required_argument,NULL,'b'},
                {"update_mode",      required_argument, NULL, 'z'},
                {"sketches_optimize",    required_argument, NULL, 'y'},
                {"algo_type",required_argument,NULL,'x'},
                {"similarity_type", required_argument, NULL,'s'},
                {0, 0, 0,                                     0}
        };

        while ((opt = getopt_long(argc, argv, optstring, long_options, &option_index)) != -1) {
            switch (opt) {
                case 'g':
                    graph_alias = string(optarg);
                    break;
                case 'p':
                    prefix = string(optarg);
                    break;
                case 'o':
                    operation = string(optarg);
                    break;
                case 'a':
                    algo = string(optarg);
                    break;
                case 'u':
                    mu = atoi(optarg);
                    break;
                case 'e':
                    epsilon = atof(optarg);
                    break;
                case 'd':
                    distance = atof(optarg);
                    break;
                case 'k':
                    hash_k = pow(2, atoi(optarg));
                    assert(hash_k < 70000);
                    break;
                case 'm':
                    max_distance = atof(optarg);
                    break;
                case 'n':
                    update_edge_nums = atoi(optarg);
                    break;
                case 'b':
                    bin = atof(optarg);
                    break;
                case 'z':
                    update_mode = atoi(optarg);
                    break;
                case 'y':
                    sketches_optimize = atoi(optarg);
                    break;
                case 'x':
                    algo_type = atoi(optarg);
                    break;
                case 's':
                    similarityType = SimilarityType(atoi(optarg));
            }
        }
        graph_location = get_graph_folder();
    }

    boost::property_tree::ptree get_data() {
        boost::property_tree::ptree data;
        data.put("graph_alias", graph_alias);
        data.put("operation", operation);
        data.put("algo", algo);
        data.put("mu", mu);
        data.put("epsilon", epsilon);
        data.put("distance", distance);
        data.put("hash_k", hash_k);
        data.put("result-dir", exe_result_dir);
        return data;
    }
};

class Result {
public:
    int n;
    long long m;

    double sim = 1;

    std::chrono::steady_clock::time_point startTime;
    double whole_time_usage;

    double total_mem_usage;
    double total_time_usage;

    vector<int> cluster_id;

    double m_scan = 0, m_dis = 0, m_wscan = 0, m_my_ads = 0, qs_scan = 0, qs_dis = 0;

    double miss_labelled_rate = 0, ari = 0;
    double topk_cluster_quality_1_min = 0, topk_cluster_quality_1_avg = 0;
    double topk_cluster_quality_5_min = 0, topk_cluster_quality_5_avg = 0;
    double topk_cluster_quality_10_min = 0, topk_cluster_quality_10_avg = 0;
    double topk_cluster_quality_20_min = 0, topk_cluster_quality_20_avg = 0;
    double topk_cluster_quality_50_min = 0, topk_cluster_quality_50_avg = 0;
    double topk_cluster_quality_100_min = 0, topk_cluster_quality_100_avg = 0;

    boost::property_tree::ptree get_data() {
        boost::property_tree::ptree data;
        data.put("n", n);
        data.put("m", m);

        //data.put("cluster_id", cluster_id);

        data.put("similarity", sim);

        data.put("whole time usage", whole_time_usage);
        data.put("total memory usage(MB)", total_mem_usage);

        data.put("total time usage(s)", total_time_usage);


        return data;
    }

};

extern Config config;
extern Result result;

bool exists_test(const std::string &name);

void assert_file_exist(string desc, string name);

namespace Saver {
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

    static string get_time_path() {
        // using namespace boost::posix_time;
        // auto tm = second_clock::local_time();
        if (!boost::algorithm::ends_with(config.exe_result_dir, FILESEP))
            config.exe_result_dir += FILESEP;
        config.exe_result_dir += "execution/";
        if (!boost::filesystem::exists(config.exe_result_dir)) {
            boost::filesystem::path dir(config.exe_result_dir);
            boost::filesystem::create_directories(dir);
        }

        string filename = config.graph_alias + "." + config.algo;

        filename += "u-" + to_string(config.mu) + ".";
        filename += "e-" + to_string(config.epsilon) + ".";


        return config.exe_result_dir + filename;
        // return config.exe_result_dir + to_iso_string(tm);
    }

    static boost::property_tree::ptree combine;

    static void init() {
        combine.put("start_time", get_current_time_str());
    }


    static void save_json(Config &config, Result &result, vector<string> args) {
        ofstream fout(get_time_path() + ".json");
        string command_line = "";
        for (int i = 1; i < args.size(); i++) {
            command_line += " " + args[i];
        }
        combine.put("end_time", get_current_time_str());
        combine.put("command_line", command_line);
        combine.put_child("config", config.get_data());
        combine.put_child("result", result.get_data());
        boost::property_tree::ptree timer;
        for (int i = 0; i < (int) Timer::timeUsed.size(); i++) {
            if (Timer::timeUsed[i] > 0) {
                timer.put(to_str(i), Timer::timeUsed[i] / TIMES_PER_SEC);
            }
        }
        combine.put_child("timer", timer);

        write_json(fout, combine, true);
    }
};


#endif //SCAN_CONFIG_H
