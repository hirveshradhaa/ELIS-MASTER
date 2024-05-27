
#include <cstdio>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <stdlib.h>
#include <ctime>
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <fstream>

using namespace std;
#define D_NUM 5
#define DEFAULT_PAA_DIMENSION 10
#define PAA_BUCKET 10
#define INF 200000
#define MAX_CHAR_PER_LINE 200000
#define MAX_CLASS 100

int PAA_DIMENSION[D_NUM] = {2, 5, 10, 20, 50};
int FUZZY = 1;

// train/test data
std::vector<std::vector<double>> data_x; 
double max_data = -INF, min_data = INF, bucket_w;
std::vector<int> label;
int datacnt[MAX_CLASS];
int shapelets_n[MAX_CLASS];

int ts_len; 
int data_n; 
int num_c; 

struct point {
    std::vector<int> f;
    int d;

    bool operator<(const point& other) const {
        return std::lexicographical_compare(f.begin(), f.end(), other.f.begin(), other.f.end());
    }
};

map<point, set<int>> tree[D_NUM][MAX_CLASS];

struct paaword {
    double score;
    set<int> allcov;
    vector<int> vec;
    int window;
    set<int> cov;
    

    bool operator<(const paaword& other) const {
        while (fabs(score - other.score) > 1e-8) {
            return score > other.score;
        }

        while (cov.size() != other.cov.size()) {
            return cov.size() > other.cov.size();
        }

        while (window * vec.size() != other.window * other.vec.size()) {
            return window * vec.size() < other.window * other.vec.size();
        }

        return lexicographical_compare(vec.begin(), vec.end(), other.vec.begin(), other.vec.end());
    }

    bool operator==(const paaword& other) const {
        if (window * vec.size() != other.window * other.vec.size()) {
            return false;
        }

        return equal(vec.begin(), vec.end(), other.vec.begin(), other.vec.end());
    }
};

vector<paaword> candidates[MAX_CLASS];
vector<paaword> results[MAX_CLASS];


void readfile(const char *data_file) {
    std::ifstream file(data_file);
    if (!file.is_open()) {
        // Handle file opening error
        std::cerr << "Error: Unable to open file " << data_file << std::endl;
        return;
    }

    std::string line;
    const char* delimiter = ", \r\n";

    data_x.reserve(data_n);
    std::set<int> unique_labels;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;

        // Extract label
        std::getline(iss, token, ',');
        double now_label = std::stod(token);
        label.push_back(now_label);
        unique_labels.insert(now_label);

        // Extract time series data
        data_x.emplace_back(std::istream_iterator<double>(iss), std::istream_iterator<double>());

        // Update min_data and max_data
        for (double now_data : data_x.back()) {
            max_data = std::max(max_data, now_data);
            min_data = std::min(min_data, now_data);
        }
    }

    num_c = unique_labels.size();

    std::vector<int> all_label(unique_labels.begin(), unique_labels.end());
    for (int i = 0; i < data_n; ++i) {
        label[i] = std::lower_bound(all_label.begin(), all_label.end(), label[i]) - all_label.begin();
        ++datacnt[label[i]];
    }

    bucket_w = (max_data - min_data) / PAA_BUCKET;
    max_data = max_data + (bucket_w / 100);
    bucket_w = (max_data - min_data) / PAA_BUCKET;

    file.close();
}


double paa2ts(int x) {
    return min_data + bucket_w * (0.5 + x);
}


void insertintohashmap(const vector<int>& qu, int offset, int d, int tsid, int len, int start) {
    // PAA window
    int window = len / PAA_DIMENSION[d];
    int s_pos = start + offset * window;

    // Calculate squared distances and accumulate
    double squared_distance = 0;
    for (int i = 0; i < PAA_DIMENSION[d]; ++i) {
        for (int j = 0; j < window; ++j) {
            double tmp_dis = paa2ts(qu[offset + i]) - data_x[tsid][s_pos + i * window + j];
            squared_distance = squared_distance + tmp_dis * tmp_dis;
        }
    }

    // Normalize the distance
    double normalized_distance = squared_distance / len;

    // Check if the distance satisfies the FUZZY condition
    if (normalized_distance > bucket_w * bucket_w * FUZZY) {
        return;
    }

    // Create the point and insert into the tree
    point f;
    f.d = PAA_DIMENSION[d];
    f.f.assign(qu.begin() + offset, qu.begin() + offset + PAA_DIMENSION[d]);
    tree[d][label[tsid]][f].insert(tsid);
}


set<int> queryfromhashmap(const vector<int>& qu, int offset, int d, int treeid) {
   
     point f;
      set<int> cnt;
    f.d = PAA_DIMENSION[d];
   
    f.f.assign(qu.begin() + offset, qu.begin() + offset + PAA_DIMENSION[d]);

    // Helper function to insert elements from the tree into the set
    auto insertElements = [&](const point& point) {
        if (tree[d][treeid].count(point)) {
            cnt.insert(tree[d][treeid][point].begin(), tree[d][treeid][point].end());
        }
    };

    // Insert elements for the original point
    insertElements(f);

    // Insert elements for neighbors (f.f[i] -= 1, f.f[i] += 2)
    for (int i = 0; i < PAA_DIMENSION[d]; ++i) {
        f.f[i] =   f.f[i] - 1;
        insertElements(f);
        f.f[i] = f.f[i] + 2;
        insertElements(f);
        f.f[i] =  f.f[i] - 1;
    }

    return cnt;
}
void createwords() {
    int shapelet_minlen = max(DEFAULT_PAA_DIMENSION, (ts_len / 20 + DEFAULT_PAA_DIMENSION - 1)
                                                     / DEFAULT_PAA_DIMENSION * DEFAULT_PAA_DIMENSION);

int len = shapelet_minlen;
while (len <= ts_len) {
    // Clear hash map
    int i = 0;
    while (i < D_NUM) {
        int j = 0;
        while (j <= num_c) {
            tree[i][j].clear();
            ++j;
        }
        ++i;
    }

    int d = 0;
    while (d < D_NUM) {
        if (len % PAA_DIMENSION[d] != 0) {
            ++d;
            continue;
        }

        // PAA window
        int window = len / PAA_DIMENSION[d];

        // Insert words into hashmap
        int i = 0;
        while (i < data_n) {
            vector<int> qu;

            // PAA
            int offset = 0;
            while (offset + window <= ts_len) {
                double avg = 0;

                // Calculate average within the window
                int j = offset;
                while (j < offset + window) {
                    avg += data_x[i][j];
                    ++j;
                }

                // Calculate bucket and push it to the vector
                int bucket = static_cast<int>((avg / window - min_data) / bucket_w + 0.5);
                qu.push_back(bucket);
                offset += window;
            }

            // Insert into hash map
            int j = 0;
            while (j + PAA_DIMENSION[d] <= qu.size()) {
                insertintohashmap(qu, j, d, i, len, j * window);
                ++j;
            }

            // Insert into hash map
            auto it = qu.begin();
            while (it + PAA_DIMENSION[d] <= qu.end()) {
                insertintohashmap(qu, it - qu.begin(), d, i, len, (it - qu.begin()) * window);
                ++it;
            }

            ++i;
        }

        // Query from hashmap to calculate tfidf score
         i = 0;
        while (i < data_n) {
            vector<int> qu;

            // PAA
            int offset = 0;
            while (offset + window <= ts_len) {
                double avg = 0;
                int j = offset;
                while (j < offset + window) {
                    avg += data_x[i][j];
                    ++j;
                }
                int bucket = static_cast<int>(floor((avg / window - min_data) / bucket_w));
                qu.push_back(bucket);
                offset += window;
            }

            // Query from hash map
            auto it = qu.begin();
            while (it + PAA_DIMENSION[d] <= qu.end()) {
                int inclasses = 0;
                vector<set<int>> select(num_c);
                set<int> allselect;

                int k = 0;
                while (k < num_c) {
                    select[k] = queryfromhashmap(qu, it - qu.begin(), d, k);
                    if (!select[k].empty()) {
                        ++inclasses;
                        allselect.insert(select[k].begin(), select[k].end());
                    }
                    ++k;
                }

                while (inclasses != 1) {
                    continue;
                }

                double idf = log2(2.0);

                k = 0;
                while (k < num_c) {
                    int inthisclass = static_cast<int>(select[k].size());

                    if (inthisclass == 0) {
                        ++k;
                        continue;
                    }

                    double tf = static_cast<double>(inthisclass) / datacnt[k];

                    paaword now;
                    now.score = tf * idf;
                    now.window = window;
                    now.allcov = allselect;
                    now.cov = select[k];

                    int cp = 0;
                    while (cp < PAA_DIMENSION[d]) {
                        now.vec.push_back(*(it + cp));
                        ++cp;
                    }

                    candidates[k].push_back(now);
                    ++k;
                }
                ++it;
            }
            ++i;
        }
        ++d;
    }
    len += shapelet_minlen;
}

}

double calculatedis(paaword paa, int tsid) {

double min_dis = std::numeric_limits<double>::infinity();
    int len = paa.window * paa.vec.size();
    
int s_pos = 0;
while (s_pos + len <= ts_len) {
    double distance = 0;

    int i = 0;
    while (i < paa.vec.size()) {
        double paa_value = paa2ts(paa.vec[i]);

        int j = 0;
        while (j < paa.window) {
            double tmp_dis = paa_value - data_x[tsid][s_pos + i * paa.window + j];
            distance += tmp_dis * tmp_dis;
            ++j;
        }
        ++i;
    }

    distance = distance / len;
    min_dis = std::min(min_dis, distance);

    ++s_pos;
}

    return min_dis;
}

int check(int classid, const vector<paaword>& tmp_result, int ENOUGH_COVER) {
    vector<int> cnt_cover(data_n, 0);

auto it_paa = tmp_result.begin();
while (it_paa != tmp_result.end()) {
    const paaword& paa = *it_paa;
    double max_dis = -1;

    auto it_cov = paa.cov.begin();
    while (it_cov != paa.cov.end()) {
        max_dis = std::max(max_dis, calculatedis(paa, *it_cov));
        ++cnt_cover[*it_cov];
        ++it_cov;
    }

    max_dis += 1e-8;

    int i = 0;
    while (i < data_n) {
        if (paa.cov.find(i) == paa.cov.end() && calculatedis(paa, i) < max_dis) {
            ++cnt_cover[i];
        }
        ++i;
    }

    ++it_paa;
}

    int cnt_wrong = count_if(cnt_cover.begin(), cnt_cover.end(),
                              [classid, ENOUGH_COVER](int cnt) { return cnt >= ENOUGH_COVER && label[cnt] != classid; });

    return cnt_wrong;
}




void getresults() {
   int i = 0;
while (i < num_c) {
    sort(candidates[i].begin(), candidates[i].end());

    int best_wrong = data_n + 1;
    int best_cover = 1;
    results[i].clear();

    int ENOUGH_COVER = 1;
    while (ENOUGH_COVER <= 5) {
        vector<paaword> tmp_result;
        map<int, int> cover;
        int total_cover = 0;

        auto it_now = candidates[i].begin();
        while (it_now != candidates[i].end()) {
            const paaword& now = *it_now;
            bool repeat = false;

            auto it_cov = now.cov.begin();
            while (it_cov != now.cov.end()) {
                if (cover[*it_cov] < ENOUGH_COVER) {
                    repeat = true;
                    break;
                }
                ++it_cov;
            }

            if (repeat) {
                ++it_now;
                continue;
            }

            bool chosen = false;

            auto it_tmp = tmp_result.begin();
            while (it_tmp != tmp_result.end()) {
                paaword& tmp = *it_tmp;
                if (now == tmp) {
                    auto it_cov_now = now.cov.begin();
                    while (it_cov_now != now.cov.end()) {
                        if (tmp.cov.find(*it_cov_now) == tmp.cov.end()) {
                            tmp.cov.insert(*it_cov_now);
                            if (cover[*it_cov_now] < ENOUGH_COVER) {
                                ++cover[*it_cov_now];
                                ++total_cover;
                            }
                        }
                        ++it_cov_now;
                    }
                    chosen = true;
                    break;
                }
                ++it_tmp;
            }

            if (!chosen) {
                tmp_result.push_back(now);
            }

            ++it_now;
        }

        // Further processing for the current ENOUGH_COVER value can be added here.

        ++ENOUGH_COVER;
    }

    ++i;
}

}
void findbestshape() {
    createwords();
    getresults();
}

void printshape(const char* out_file) {
    std::ofstream outFile(out_file);

    for (int i = 0; i < num_c; ++i) {
        outFile << shapelets_n[i] << '\n';

        for (const paaword& now : results[i]) {
            outFile << -now.score << ',' << now.window * now.vec.size();

            for (int k = 0; k < now.window * now.vec.size(); ++k) {
                outFile << ',' << paa2ts(now.vec[k / now.window]);
            }

            outFile << '\n';
        }
    }
}

int main(int argc, const char *argv[]) {
    // Check if the required command line arguments are provided
    if (argc < 4) {
        fprintf(stderr, "Usage: %s ts_len data_file data_n [PAA fuzzy=1] [cover=1]\n", argv[0]);
        return 1;
    }

    // Parse command line arguments
    ts_len = atoi(argv[1]);
    const char *data_file = argv[2];
    data_n = atoi(argv[3]);

    // Set default values for optional arguments
    FUZZY = 1;

    // Parse optional arguments if provided
    if (argc >= 5) {
        FUZZY = atoi(argv[4]);
        FUZZY *= FUZZY;
    }

    // Read data from the file
    readfile(data_file);
    printf("Reading data finished\n");

    // Measure the start time for shapelet discovery
    double start_time = clock();

    // Perform shapelet discovery
    findbestshape();

    // Measure the end time for shapelet discovery
    double end_time = clock();
    printf("Shapelet discovery finished\n");

    // Print the discovered shapelets to a file
    printshape("init.txt");

    // Print the total running time
    printf("Total Running Time = %.3f sec\n", (end_time - start_time) / CLOCKS_PER_SEC);

    return 0;
}