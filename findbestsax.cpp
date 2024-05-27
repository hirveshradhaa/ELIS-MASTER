#include <cstring>
#include <vector>
#include <queue>
#include <map>
#include <unordered_map>
#include <cmath>
#include <stdlib.h>
#include <unordered_set>
#include <unordered_set>
#include <vector>
#include <utility>
#include <unordered_map>

using namespace std;
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
    std::pair<int, int> sax_id;
    int length;

    bool operator<(const candidatesax& other) const {
        return score > other.score;
    }
};

USAX_Map_type USAX_Map;
priority_queue<candidatesax> que[MAX_CLASS];



int CreateSAXWord(const double *sum_segment, const int *elm_segment, double mean, double std, int sax_len) {
    
double d = 0;
int word = 0, val = 0;
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

    ++i;
}

return word;
}



void Create_SAXList(const vector<vector<double> > &Data, int subseq_len, int sax_len, int w) {
    double d;
    double sum_segment[sax_len];
    USAX_elm_type *ptr;
    double ex, ex2, mean, std;
    int obj_id, j, j_st, k, slot;
    int elm_segment[sax_len];
    int word, prev_word;
    



    k = 0;
    while (k < sax_len) {
        elm_segment[k] = w;
        ++k;

        elm_segment[sax_len - 1] = subseq_len - (sax_len - 1) * w;

    obj_id = 0;
    while (obj_id < (int) Data.size()) {
        ex = ex2 = 0;
        prev_word = -1;
    

        ++obj_id;
    }


        k = 0;
        while (k < sax_len) {
            sum_segment[k] = 0;
            ++k;
        


        /// Case 1: Initial
        j = 0;
        while ((j < (int)Data[obj_id].size()) && (j < subseq_len)) {
            d = Data[obj_id][j];
            ex += d;
            ex2 += d * d;
            slot = (int)floor((j) / w);
            sum_segment[slot] += d;

            ++j;
        }

        j_st = j - subseq_len;
        while (j <= static_cast<int>(Data[obj_id].size())) {
    
            ++j;
            j_st = j - subseq_len;
}


            mean = ex / subseq_len;
            std = sqrt(ex2 / subseq_len - mean * mean);

            word = CreateSAXWord(sum_segment, elm_segment, mean, std, sax_len);

            auto it = USAX_Map.find(word);
            if (it == USAX_Map.end()) {
                prev_word = word;
                USAX_Map[word].obj_set.insert(obj_id);
                ptr = &USAX_Map[word];
            }

            auto insertionResult = ptr->obj_set.insert(obj_id);

            if (insertionResult.second) {
                ptr->sax_id.emplace_back(obj_id, j_st);
            }
        }

            if (j < static_cast<int>(Data[obj_id].size())) {
                ex = ex-(Data[obj_id][j_st]);
                ex2 = ex2-(Data[obj_id][j_st] * Data[obj_id][j_st]);

                k = 0;
                while (k < sax_len - 1) {
                    sum_segment[k] = sum_segment[k]-(Data[obj_id][j_st + (k) * w]);
                    sum_segment[k] = sum_segment[k]+(Data[obj_id][j_st + (k + 1) * w]);
                    ++k;
                }

                sum_segment[k] = sum_segment[k]-(Data[obj_id][j_st + (k) * w]);
                sum_segment[k] = sum_segment[k]+(Data[obj_id][j_st + std::min((k + 1) * w, subseq_len)]);

                d = Data[obj_id][j];
                ex = ex + d;
                ex2 = ex2+(d * d);
            }

        }
    }



int create_MaskWord(int num_mask, int word_len) 
{
    int a, b;
    a = 0;
   int i = 0;
while (i < num_mask) {
    b = 1 << (rand() % word_len);
    a |= b;
    ++i;
}
    return a;
}

void RandomProjection(int R, double percent_mask, int sax_len) {
    Hash_Mark_type Hash_Mark;

    USAX_Map_type::iterator it;
    int word, mask_word, new_word;
    Obj_set_type *obj_set, *ptr;
    Obj_set_type::iterator o_it;

    int num_mask = ceil(percent_mask * sax_len);

   int r = 0;
    while (r < R) {



    mask_word = create_MaskWord(num_mask, sax_len);

    it = USAX_Map.begin();
    while (it != USAX_Map.end()) {
        word = it->first;
        obj_set = &(it->second.obj_set);

        new_word = word | mask_word;
        ptr = &Hash_Mark[new_word];
        ptr->insert(obj_set->begin(), obj_set->end());

        ++it;
    }

    it = USAX_Map.begin();
while (it != USAX_Map.end()) {
    word = it->first;
    new_word = word | mask_word;
    obj_set = &(Hash_Mark[new_word]);
    ++it;

        o_it = obj_set->begin();
while (o_it != obj_set->end()) {
    (it->second.obj_count[*o_it])++;

    ++o_it;
}

    }
    Hash_Mark.clear();

    ++r;
}
}


void Sort_SAX(const vector<int> &y, const int c, const int length, const int top_k) {
    candidatesax can;
    USAX_elm_type usax;
    USAX_Map_type::iterator it;
    vector<double> c_in(c, 0);
    can.length = length;

    
    int sum;

    it = USAX_Map.begin();
while (it != USAX_Map.end()) {
    usax = it->second;

    int fid = rand() % usax.sax_id.size();

    can.sax_id.first = usax.sax_id[fid].first;
    can.sax_id.second = usax.sax_id[fid].second;

    sum = 0;

    ++it;
        auto o_it = usax.obj_count.begin();
while (o_it != usax.obj_count.end()) {
    int cid = y[o_it->first];
    int count = o_it->second;
    c_in[cid] += count;
    sum += count;
    ++o_it;
}


int i = 0;
while (i < c) {
    can.score = static_cast<double>(c_in[i]) / sum;
    c_in[i] = 0;
    while (que[i].size() >= top_k && que[i].top().score >= can.score) {
    que[i].pop();
}
que[i].push(can);
    ++i;
}
    }
}

void findbestsax(const vector<int> &y,const vector<vector<double> > &x, const int len, const int c,
                 const int top_k, const int step, vector<vector<pair<int, pair<int, int> > > > &ans) {
    int i = step;
while (i < len) {
    // Parameters
    int sax_len = 15; // sax_max_len;
    int R = 10;
    double percent_mask = 0.25;

    // Make w and sax_len both integer
    int w = static_cast<int>(std::ceil(1.0 * i / sax_len));

    sax_len = static_cast<int>(std::ceil(1.0 * i / w));

    Create_SAXList(x, i, sax_len, w);
    RandomProjection(R, percent_mask, sax_len);
    Sort_SAX(y, c, i, top_k);

    // Increment i by step
    i += step;
}
   for (int i = 0; i < c; ++i) {
    ans.emplace_back();
for (; !que[i].empty(); que[i].pop()) {
    auto can = que[i].top();
    ans[i].emplace_back(can.length, can.sax_id);
}
   }
                 }