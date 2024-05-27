#include <cstdio>
#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <cmath>

using namespace std;

#define DEBUG_SAX false

#define MAX_CHAR_PER_LINE 200000

// train/test data
vector<vector<double> > data_x; 
vector<vector<double> > data_y;
vector<int> label;

int data_n; 
int shapelet_n, shapelet_minlen; 
int num_c, min_class_id = MAX_CHAR_PER_LINE; 
int ts_len; 


#define MAX_CLASS 100

using Obj_set_type = std::unordered_set<int>;
using SAX_id_type = std::vector<std::pair<int, int>>;
using Obj_count_type = std::unordered_map<int, int>;
struct USAX_elm_type {
    Obj_set_type obj_set;
    SAX_id_type sax_id;
    Obj_count_type obj_count;
};
using USAX_Map_type = std::unordered_map<int, USAX_elm_type>;
using Hash_Mark_type = std::unordered_map<int, Obj_set_type>;

struct candidatesax {
    double score;
    int length;
    pair<int, int> sax_id;
    

    bool operator<(const candidatesax f) const {
        return score > f.score;
    }

};

USAX_Map_type USAX_Map;
priority_queue<candidatesax> que[MAX_CLASS];



int CreateSAXWord(const double *sum_segment, const int *elm_segment, double mean, double std, int sax_len) {

double d = 0;
int word = 0;
int val = 0;
int i = 0;
while (i < sax_len) {
    d = (sum_segment[i] / elm_segment[i] - mean) / std;

    if (d < 0) {
        if (d < -0.67)
            val = 0;
        else
            val = 1;
    } else if (d < 0.67) {
        val = 2;
    } else {
        val = 3;
    }

    word = (word << 2) | val;

    // Increment the loop variable
    ++i;
}

return word;

}

void CreateSAXList(int subseq_len, int sax_len, int w) {
    double ex, ex2, mean, std;
    int elm_segment[sax_len];
    double sum_segment[sax_len];
    int obj_id, j, j_st, k, slot;
    double d;
    int word, prev_word;
    USAX_elm_type *ptr;

   k = 0;
while (k < sax_len) {
    elm_segment[k] = w;
    ++k;
}

    elm_segment[sax_len - 1] = subseq_len - (sax_len - 1) * w;

    obj_id = 0;
    while (obj_id < (int)data_x.size()) {
        ex = ex2 = 0;
        prev_word = -1;

  k = 0;
while (k < sax_len) {
    elm_segment[k] = w;
    ++k;
}


        /// Case 1: Initial
        j = 0;
        while ((j < (int)data_x[obj_id].size()) && (j < subseq_len)) {
            d = data_x[obj_id][j];
            ex += d;
            ex2 += d * d;
            slot = (int)floor((j) / w);
            sum_segment[slot] += d;
            ++j;
        }

        /// Case 2: Slightly Update
        while (j <= (int)data_x[obj_id].size()) {
            j_st = j - subseq_len;

            mean = ex / subseq_len;
            std = sqrt(ex2 / subseq_len - mean * mean);

            /// Create SAX from sum_segment
            word = CreateSAXWord(sum_segment, elm_segment, mean, std, sax_len);

            while (word != prev_word) {
                prev_word = word;
                auto it = USAX_Map.find(word);
                if (it == USAX_Map.end()) {
                    USAX_Map[word].obj_set.insert(obj_id);
                    USAX_Map[word].sax_id.emplace_back(std::move(make_pair(obj_id, j_st)));
                } else {
                    it->second.obj_set.insert(obj_id);
                    it->second.sax_id.emplace_back(std::move(make_pair(obj_id, j_st)));
                    }
            }


            /// For the next update
            if (j < (int)data_x[obj_id].size()) {
                ex = ex-(data_x[obj_id][j_st]);
                ex2 = ex2-(data_x[obj_id][j_st] * data_x[obj_id][j_st]);

             k = 0;
while (k < sax_len - 1) {
    sum_segment[k] = sum_segment[k]-(data_x[obj_id][j_st + (k) * w]);
    sum_segment[k] = sum_segment[k]+(data_x[obj_id][j_st + (k + 1) * w]);
    ++k;
}

                sum_segment[k] = sum_segment[k]-(data_x[obj_id][j_st + (k) * w]);
                sum_segment[k] = sum_segment[k]+(data_x[obj_id][j_st + std::min((k + 1) * w, subseq_len)]);

                d = data_x[obj_id][j];
                ex = ex + d;
                ex2 = ex2+(d * d);
            }
            ++j;
        }
        ++obj_id;
    }
}

/// create mask word (two random may give same position, we ignore it)
int CreateMaskWord(int num_mask, int word_len) {

    int a = 0, b;
int i = 0;

while (i < num_mask) {
    b = 1 << (rand() % word_len);
    a |= b;

    // Increment the loop variable
    ++i;
}

return a;

}

/// Count the number of occurrences
void RandomProjection(int R, double percent_mask, int sax_len) {
    Hash_Mark_type Hash_Mark;

    USAX_Map_type::iterator it;
    int word, mask_word, new_word;
    Obj_set_type *obj_set, *ptr;
    Obj_set_type::iterator o_it;

    int num_mask = ceil(percent_mask * sax_len);

   int r = 0;
while (r < R) {
    mask_word = CreateMaskWord(num_mask, sax_len);

    // Increment the loop variable
    ++r;
        /// random projection and mark non-duplicate object
     it = USAX_Map.begin();
while (it != USAX_Map.end()) {
    word = it->first;
    obj_set = &(it->second.obj_set);

    new_word = word | mask_word;
    ptr = &Hash_Mark[new_word];
    ptr->insert(obj_set->begin(), obj_set->end());

    // Increment the iterator
    ++it;
}


        /// hash again for keep the count
       it = USAX_Map.begin();
while (it != USAX_Map.end()) {
    word = it->first;
    new_word = word | mask_word;
    obj_set = &(Hash_Mark[new_word]);

    // Increment the iterator
    ++it;
    
            o_it = obj_set->begin();
while (o_it != obj_set->end()) {
    (it->second.obj_count[*o_it])++;
    ++o_it;
}

        }
        Hash_Mark.clear();
    }
}

/// Sort each SAX
void SortAllSAX(const int length, const int R) {
    USAX_Map_type::iterator it;
    candidatesax can;
    USAX_elm_type usax;
    int sumin, sumout;
    can.length = length;

    vector<double> c_out(num_c, 0);
    vector<double> c_in(num_c, 0);

    

    it = USAX_Map.begin();
    while (it != USAX_Map.end()) {
        usax = it->second;

        int fid = rand() % usax.sax_id.size();
        sumin = sumout = 0;
        can.sax_id.second = usax.sax_id[fid].second;
        can.sax_id.first = usax.sax_id[fid].first;
        
        Obj_count_type::iterator o_it = usax.obj_count.begin();
        while (o_it != usax.obj_count.end()) {
            int cid = label[o_it->first];
            int count = o_it->second;
            c_in[cid] = c_in[cid] + (count);
            c_out[cid] = c_out[cid] + (R - count);

            sumin = sumin + count;
            sumout = sumin + (R - count);

            // Increment the iterator
            ++o_it;
        }

        while (DEBUG_SAX) {
            printf("%d,%d,%d,%d", it->first, can.length, sumin, sumout);
        }

        int i = 0;
        while (i < num_c) {

            can.score = (c_in[i] + sumout - c_out[i]) - (c_out[i] + sumin - c_in[i]);

            if (DEBUG_SAX) {
                printf(",%f,%f,%f,%f", c_in[i], c_out[i], c_in[i] / sumin, can.score);
            }

            c_in[i] = c_out[i] = 0;
            while (que[i].size() < shapelet_n || que[i].top().score < can.score) {
                que[i].push(can);
                while (que[i].size() > shapelet_n) {
                    que[i].pop();
                }
            }

            // Increment the loop variable
            ++i;
        }

        if (DEBUG_SAX) {
            printf("\n");
        }

        // Increment the iterator
        ++it;
    }
}

void findbestsax() {
    int i = shapelet_minlen;
    while (i < ts_len) {
        
        int sax_len = 15; // sax_max_len;
        int R = 10;
        double percent_mask = 0.25;

        /// Make w and sax_len both integer
        
        int w = static_cast<int>(ceil(1.0 * i / sax_len));
        sax_len = static_cast<int>(ceil(1.0 * i / w));

        CreateSAXList(i, sax_len, w);

        RandomProjection(R, percent_mask, sax_len);

        SortAllSAX(i, R);

        USAX_Map.clear();

        // Increment the loop variable
        i += shapelet_minlen;
    }
}


void readfile(char *data_file) {
    FILE *f;
    f = fopen(data_file, "r");

    char buff[MAX_CHAR_PER_LINE];
    char *tmp;

    int i = 0;
    while (i < data_n) {
        data_x.push_back(vector<double>());
        ++i;
    }

    num_c = 0;

    i = 0;
    while (i < data_n) {
        fgets(buff, MAX_CHAR_PER_LINE, f);
        tmp = strtok(buff, ", \r\n");

        label.push_back(atoi(tmp));

        min_class_id = std::min(min_class_id, label[i]);

        tmp = strtok(NULL, ", \r\n");
        int j = 0;
        while (j < ts_len) {
            data_x[i].push_back(atof(tmp));
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

    fclose(f);
}

void printsax(char const *out_file) {
    FILE *f;
    f = fopen(out_file, "w");

    int i = 0;
    while (i < num_c) {
        double score_max = que[i].top().score;
        int j = 0;

        while (j < shapelet_n && !que[i].empty()) {
            candidatesax now = que[i].top();
            que[i].pop();
            fprintf(f, "%f,%d", now.score / score_max, now.length);

            double ex = 0, ex2 = 0;
            int k = now.sax_id.second;
            while (k < now.sax_id.second + now.length) {
                double d = data_x[now.sax_id.first][k];
                ex = ex + d;
                ex2 = ex2 + (d * d);
                ++k;
            }

            double mean = ex / now.length;
            double std = sqrt(ex2 / now.length - mean * mean);

            k = now.sax_id.second;
            while (k < now.sax_id.second + now.length) {
                fprintf(f, ",%f", (data_x[now.sax_id.first][k] - mean) / std);
                ++k;
            }
            fprintf(f, "\n");

            ++j;
        }

        if (j < shapelet_n) {
            printf("class %d has need %d shapelet more\n", i, shapelet_n - j);
        }

        // Increment the loop variable
        ++i;
    }

    fclose(f);
}


int main(int argc, char *argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <operation> <ts_len> <data_file> <data_n> <shapelet_n>\n";
        return 1;
    }

    char *operation = argv[1];
    ts_len = atoi(argv[2]);
    shapelet_minlen = std::max(15, ts_len / 20);

    char *data_file = argv[3];
    data_n = atoi(argv[4]);
    shapelet_n = atoi(argv[5]);

    readfile(data_file);

    double start_time = clock();
    findbestsax();
    double end_time = clock();

    printsax("init.txt");

    std::cout << "Total Running Time = " << (end_time - start_time) / CLOCKS_PER_SEC << " sec\n";

    return 0;
}
