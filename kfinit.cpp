
#include <cstdio>
#include <iostream>
#include <cstring>
#include <cmath>
#include <stdlib.h>
#include <ctime>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <chrono>
#include <algorithm>
using namespace std;

const int RTREE_DIMENSION = 10;
const int PaaBucket = 10;
const int inf = 200000;


#define CoverSame 1
#define PAA_DIFF 2
#define PAA_SAME 0.1
#define PAA_LIKE 0.2
#define MAX_CHAR_PER_LINE 200000
#define MAX_CLASS 100

#define x2_0250 5.02
#define x2_0100 6.63
#define x2_0050 7.88
#define x2_0025 9.14
#define x2_0200 5.41
#define x2_0500 3.84
#define x2_0010 10.83
#define x2_0005 12.12

// train/test data
int datacnt[MAX_CLASS];
int shapelets_n[MAX_CLASS];
vector<vector<double> > data_x; // data_n * ts_len
double max_data = -inf, min_data = inf, bucket_w;
int num_c, min_class_id = inf; // number of class , min class id (init inf)
int ts_len; // time series length
int data_n; // number of data
vector<int> label;

struct point {
    int f[RTREE_DIMENSION];

    bool operator<(const point &other) const {
      int i = 0;
        while (i < RTREE_DIMENSION) {
            if (f[i] < other.f[i]) {
                return true;
            } else if (f[i] > other.f[i]) {
                return false;
            }
            ++i;
        }
        return false;
    }
};


map<point, set<int> > tree[MAX_CLASS];


struct candidatesax {
    double score;
    int sign;
    int window;
    std::vector<int> vec;
    std::set<int> cov;
    int a, b, c, d;
    double ta, tb, tc, td;

    bool operator<(const candidatesax& other) const {
        if (window != other.window) {
            return window < other.window;
        }

        int i = 0;
        while (i < RTREE_DIMENSION) {
            if (vec[i] < other.vec[i]) {
                return true;
            } else if (vec[i] > other.vec[i]) {
                return false;
            }
            ++i;
        }

        return false;
    }

    bool operator==(const candidatesax& f) const {
        if (f.window * 2 <= window) {
            return false;
        }

        int mincnt = f.window * RTREE_DIMENSION;

        auto checkPaaDiff = [&](int i, int j, int di, int dj) {
            return vec[di] < f.vec[dj] - PAA_DIFF || vec[di] > f.vec[dj] + PAA_DIFF;
        };

        auto checkPaaSame = [&](int& cnt, int num, int& i, int& j, int& di, int& dj) {
            if (vec[di] < f.vec[dj] - 1 || vec[di] > f.vec[dj] + 1) {
                if (checkPaaDiff(i, j, di, dj)) {
                    cnt = f.window * RTREE_DIMENSION;
                    return true;
                }
                cnt += num;
                return cnt > PAA_LIKE * RTREE_DIMENSION * f.window;
            }
            return false;
        };

        int offset = 0;
        while (offset <= window * RTREE_DIMENSION - f.window * RTREE_DIMENSION) {
            int cnt = 0;

            int i = window - offset % window;
            int di = offset / window;
            int j = f.window;
            int dj = 0;

            while (dj < RTREE_DIMENSION) {
                int num = std::min(i, j);

                if (checkPaaSame(cnt, num, i, j, di, dj)) {
                    break;
                }

                i = i - num;
                j =  j - num;
j = (j == 0) ? f.window : j;
dj += (j == f.window) ? 1 : 0;

i = (i == 0) ? window : i;
di += (i == window) ? 1 : 0;

               
            }

            if (cnt <= PAA_SAME * RTREE_DIMENSION * f.window) {
                return true;
            }

            mincnt = std::min(mincnt, cnt);

            ++offset;
        }

        if (mincnt > PAA_LIKE * RTREE_DIMENSION * f.window) {
            return false;
        }

        int cnt = 0;
        auto t1 = cov.begin();
        auto t2 = f.cov.begin();

        while (t1 != cov.end() && t2 != f.cov.end()) {
    
  cnt += ((*t1) == (*t2)) ? (++t1, ++t2, 1) : ((*t1) < (*t2) ? (++t1, 0) : (++t2, 0));
        
}

//
        return cnt >= cov.size() * CoverSame && cnt >= f.cov.size() * CoverSame;

       if (cnt >= std::max(f.cov.size(), cov.size()) * CoverSame) {
            printf("----------------PAA COVER LIKE----------------\n");
            printf("%f,%d", score * sign, window * RTREE_DIMENSION);

            int k = 0;
            while (k < RTREE_DIMENSION) {
                printf(",%d", vec[k]);
                ++k;
            }

            printf("\n");
            printf("-----------------------------------------------\n");

            printf("\ncover");

            auto printCov = [](const set<int>& cov) {
                auto x = cov.begin();
                while (x != cov.end()) {
                    printf(" %d", *x);
                    ++x;
                }
            };

            printCov(cov);
            printf("\n");

            printf("%f,%d", f.score * f.sign, f.window * RTREE_DIMENSION);

            k = 0;
            while (k < RTREE_DIMENSION) {
                printf(",%d", f.vec[k]);
                ++k;
            }

            printf("\ncover");
            printCov(f.cov);
            printf("\n");
            printf("-----------------------------------------------\n");
        
            return true;
        }
        return false;
    }
};
set<candidatesax>::iterator it;
set<candidatesax> que[MAX_CLASS];


vector<candidatesax> buffer;

void readfile(const char *data_file) {
    FILE *f = fopen(data_file, "r");

    char buff[MAX_CHAR_PER_LINE];

    data_x.reserve(data_n);
    label.reserve(data_n);

    num_c = 0;

   int i = 0;
while (i < data_n) {
    fgets(buff, MAX_CHAR_PER_LINE, f);
    char *tmp = strtok(buff, ", \r\n");

    label.push_back(atoi(tmp));

    min_class_id = std::min(min_class_id, label[i]);

    tmp = strtok(NULL, ", \r\n");
    data_x.emplace_back(); // Add an empty vector for each data point

    int j = 0;
    while (j < ts_len) {
        double now_data = atof(tmp);
        max_data = std::max(max_data, now_data);
        min_data = std::min(min_data, now_data);

        data_x[i].push_back(now_data);
        tmp = strtok(NULL, ", \r\n");
        ++j;
    }

    ++i;
}
   i = 0;
while (i < data_n) {
    label[i] -= min_class_id;
    num_c = std::max(num_c, label[i]);
    ++i;
}
    ++num_c;
   i = 0;
while (i < data_n) {
    ++datacnt[label[i]];
    ++i;
}


    bucket_w = (max_data - min_data) / PaaBucket;
    max_data = max_data + (bucket_w / 2);
    bucket_w = (max_data - min_data) / PaaBucket;

    fclose(f);
}

void insertintortree(vector<int> &qu, int offset, int tsid) {
    point f;
   int i = 0;
while (i < RTREE_DIMENSION) {
    f.f[i] = qu[offset + i];
    ++i;
}

    tree[label[tsid]][f].insert(tsid);
    tree[num_c][f].insert(tsid);
}


double x2test(double a, double b, double c, double d) {
   double ta = (a + c) * (a + b);
double td = (c + d) * (b + d);
double tc = (c + d) * (a + c);
double tb = (b + d) * (a + b);

if (ta < 6 || tb < 6 || tc < 6 || td < 6) {
    double tmp = (tmp < 0) ? -tmp - data_n / 2 : fabs(tmp) - data_n / 2;
    return tmp * tmp * data_n / ((a + b) * (a + c) * (c + d) * (b + d));
}

double tmp = d * a - c * b;
return tmp * fabs(tmp) * data_n / ((a + c) * (a + b) * (b + d) * (c + d));

}
//

set<int> queryfromrtree(vector<int>& qu, int offset, int treeid) {
   point f;
    set<int> cnt;
    

    auto insertIntoSet = [&](const point& f) {
        auto it = tree[treeid].find(f);
if (it != tree[treeid].end()) {
    cnt.insert(it->second.begin(), it->second.end());
}
    };

   int i = 0;
while (i < RTREE_DIMENSION) {
    f.f[i] = qu[offset + i];
    ++i;
}
    insertIntoSet(f);

  i = 0;
while (i < RTREE_DIMENSION) {
    f.f[i] -= 1;
    insertIntoSet(f);
    f.f[i] += 2;
    insertIntoSet(f);
    f.f[i] -= 1;
    ++i;
}
    return cnt;
}
 double x2table(double v) {
   std::map<double, double> thresholds = {
        {0.0005, 0.9995},
        {0.0010, 0.9990},
        {0.0025, 0.9975},
        {0.0050, 0.9950},
        {0.0100, 0.9900},
        {0.0200, 0.9800},
        {0.0250, 0.9750},
        {0.0500, 0.9500}
    };

    for (const auto& entry : thresholds) {
        if (entry.first < v) {
            return entry.second;
        }
    }

    return 0.0;
}
                  double calculateX2Score(int a, int b, int c, int d) {
    return x2test(a, b, c, d);
}
void findbestshape() {

    int sumvec = 0;

    int shapelet_minlen = max(RTREE_DIMENSION, (ts_len / 20 + RTREE_DIMENSION - 1) / RTREE_DIMENSION * RTREE_DIMENSION);

  int len = shapelet_minlen;
while (len <= ts_len) {
    printf("search length %d:\n", len);

    // clear rtree
  int i = 0;
while (i <= num_c) {
    tree[i].clear();
    ++i;
}


    printf("Hash map clear\n");

    // PAA window
    int window = len / RTREE_DIMENSION;

     i = 0;
    while (i < data_n) {
        int start = 0;
        while (start < window) {
            vector<int> qu;

            // PAA
            int offset = start;
            while (offset + window <= ts_len) {
                double avg = 0;
                int j = offset;
                while (j < offset + window) {
                    avg += data_x[i][j];
                    ++j;
                }
                int bucket = (int) floor((avg / window - min_data) / bucket_w);
                qu.push_back(bucket);
                offset += window;
            }

            // insert into hash map
            int j = 0;
            while (j + RTREE_DIMENSION <= qu.size()) {
                insertintortree(qu, j, i);
                ++j;
            }
            ++start;
        }
        ++i;
    }

    printf("insert finished\n");

     i = 0;
    while (i < data_n) {
        int start = 0;
        while (start < window) {
            vector<int> qu;

            // PAA
            int offset = start;
            while (offset + window <= ts_len) {
                double avg = 0;
                int j = offset;
                while (j < offset + window) {
                    avg += data_x[i][j];
                    ++j;
                }
                int bucket = (int) floor((avg / window - min_data) / bucket_w);
                qu.push_back(bucket);
                offset += window;
            }

            // query from hash map
            int j = 0;
            while (j + RTREE_DIMENSION <= qu.size()) {
                int seleceall = queryfromrtree(qu, j, num_c).size();
                int k = 0;
                while (k < num_c) {
                   auto covts = queryfromrtree(qu, j, k);
int selectc = covts.size();

double score = calculateX2Score(
    selectc, seleceall - selectc, datacnt[k] - selectc,
    data_n - seleceall - datacnt[k] + selectc
);


                    
                   if (fabs(score) > x2_0010) {
    candidatesax now;
    now.score = x2table(fabs(score));
    now.sign = (score < 0) ? -1 : 1;
    now.window = window;

    int cp = 0;
    while (cp < RTREE_DIMENSION) {
        now.vec.push_back(qu[j + cp]);
        ++cp;
    }

    now.cov = covts;
    now.a = selectc;
    now.c = datacnt[k] - selectc;
    now.b = seleceall - selectc;
    now.d = data_n - seleceall + selectc - datacnt[k];

    now.ta = (now.a + now.c) / data_n * (double)(now.a + now.b);
    now.td = (now.c + now.d) / data_n * (double)(now.b + now.d);
    now.tb = (now.b + now.d) / data_n * (double)(now.a + now.b);
    now.tc = (now.c + now.d) / data_n * (double)(now.a + now.c);

    que[k].insert(now);
} else {
    // Handle the else case if needed
}

                    ++k;
                }
                ++j;
            }
            ++start;
        }
        ++i;
    }

    printf("query finished\n");

    sumvec += tree[num_c].size();
    len += shapelet_minlen;
}

    printf("total distinct subsequence %d\n", sumvec);
}
//
void printshape(const char *out_file) {
    FILE *f = fopen(out_file, "w");
    
int i = 0;
while (i < num_c) {
    shapelets_n[i] = 0;
    buffer.clear();

  auto it = que[i].begin();

// Using std::copy_if and lambda function
vector<candidatesax> curlen;
std::copy_if(it, que[i].end(), std::back_inserter(curlen),
    [&it](const candidatesax& now) {
        bool condition = (it->window == now.window);
        ++it;
        return condition;
    });


        int j = 0;
        while (j < curlen.size()) {
            bool insertcandidate = true;
            auto found = std::find(buffer.begin(), buffer.end(), curlen[j]);
            if (found != buffer.end()) {
                insertcandidate = false;
            }
            if (insertcandidate) {
                ++shapelets_n[i];
                buffer.push_back(curlen[j]);
            }
            ++j;
        }
    }

    fprintf(f, "%d\n", shapelets_n[i]);
    int j = 0;
    while (j < shapelets_n[i]) {
        candidatesax &now = buffer[j];
        fprintf(f, "%f,%d", now.score * now.sign, now.window * RTREE_DIMENSION);
        int k = 0;
        while (k < RTREE_DIMENSION) {
            fprintf(f, ",%d", now.vec[k]);
            ++k;
        }
        fprintf(f, "\ncover");
     auto it = now.cov.begin();
while (it != now.cov.end()) {
    fprintf(f, " %d", *it);
    ++it;
}
        fprintf(f, "\n");
        fprintf(f, "%d(%f)\t\t%d(%f)\n", now.a, now.ta, now.b, now.tb);
        fprintf(f, "%d(%f)\t\t%d(%f)\n", now.c, now.tc, now.d, now.td);

        ++j;
    }
    ++i;
}


void printUsage(const char* programName) {
    std::cerr << "Usage: " << programName << " <ts_len> <data_file> <data_n>\n";
}

void runAlgorithm(const char* data_file) {
    readfile(data_file);

    std::cout << "Reading data finished\n";

    auto start_time = std::chrono::high_resolution_clock::now();
    findbestshape();
    auto end_time = std::chrono::high_resolution_clock::now();

    std::cout << "Algorithm finished, start outputting\n";

    printshape("init.txt");

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "Total Running Time = " << duration / 1000.0 << " sec\n";
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        printUsage(argv[0]);
        return 1;
    }

    ts_len = std::stoi(argv[1]);
    const char* data_file = argv[2];
    data_n = std::stoi(argv[3]);

    runAlgorithm(data_file);

    return 0;
}

