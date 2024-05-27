#include <cstdio>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <ctime>
#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#define MAXDOUBLE 1e100
#define MAXCHARPERLINE 200000
#define MINACCURACY 1e-100


#define PRINTLOG true

vector<vector<double> > datay; 
vector<int> rightlog;
vector<vector<double> > datax; 


double alpha = -25.0;
int tslen; 
int datan; 
int numc; 
double learn_rate = 0.01;
int max_iter = 10000;
double regularization = 0.01;
int right_cnt = 0;
int test_n; 

vector<vector<vector<double> > > e_dis;
vector<vector<double> > w; 
vector<vector<int> > shapelet_length; 
vector<vector<vector<double> > > shapelet; 
vector<double> w0;  
vector<vector<double> > m; 
vector<int> shapelet_n;
vector<vector<vector<double> > > dis;  
vector<vector<double> > test_y; 
vector<vector<double> > test_x;
vector<vector<double> > sum_dis;


void readfile(char *data_file) {
    FILE *f;
    f = fopen(data_file, "r");

    char buff[MAXCHARPERLINE];
    char *tmp;
int i = 0;
while (i < datan) {
    datax.push_back(vector<double>());
    datay.push_back(vector<double>());
    ++i;
}

    int *label = new int[datan];
    int min_c = MAXCHARPERLINE;

  i = 0;
while (i < datan) {
    fgets(buff, MAXCHARPERLINE, f);
    tmp = strtok(buff, ", \r\n");

    label[i] = atoi(tmp);
    min_c = min(min_c, label[i]);

    tmp = strtok(NULL, ", \r\n");
    int j = 0;
    while (j < tslen) {
        datax[i].push_back(atof(tmp));
        tmp = strtok(NULL, ", \r\n");
        ++j;
    }

    ++i;
}

numc = 0;
   i = 0;
while (i < datan) {
    label[i] -= min_c;
    numc = max(numc, label[i]);
    ++i;
}

    ++numc;

   i = 0;
while (i < datan) {
    int j = 0;
    while (j < numc) {
        if (j == label[i]) {
            datay[i].push_back(1.0);
        } else {
            datay[i].push_back(0.0);
        }
        ++j;
    }
    ++i;
}

delete[] label;


    fclose(f);
}


void init(int op) {
    std::ifstream file;
    std::string line;

    if (op == 0) {
        file.open("init.txt");
        w0.assign(numc, 1.0);
    } else {
        file.open("learned.txt");
        std::getline(file, line);
        std::istringstream iss(line);
       int i = 0;
while (i < numc && iss >> w0[i]) {
    ++i;
}

    }

int i = 0;
while (i < numc && std::getline(file, line)) {
    std::istringstream iss(line);
    int shapelet_count;
    iss >> shapelet_count;
    shapelet_n.push_back(shapelet_count);

    shapelet_length.push_back(std::vector<int>());
    shapelet.push_back(std::vector<std::vector<double>>());
    w.push_back(std::vector<double>());

    int j = 0;
    while (j < shapelet_count && std::getline(file, line)) {
        std::istringstream iss(line);

        iss >> w[i][j];
        iss >> shapelet_length[i][j];
        shapelet[i].push_back(std::vector<double>());

        int k = 0;
        while (k < shapelet_length[i][j]) {
            double val;
            if (iss >> val) {
                shapelet[i][j].push_back(val);
            }
            ++k;
        }

        ++j;
    }

    ++i;
}

dis.resize(numc, std::vector<std::vector<double>>());
e_dis.resize(numc, std::vector<std::vector<double>>());
sum_dis.resize(numc, std::vector<double>());
m.resize(numc, std::vector<double>());

i = 0;
while (i < numc) {
    dis[i].resize(shapelet_n[i], std::vector<double>());
    e_dis[i].resize(shapelet_n[i], std::vector<double>());
    sum_dis[i].resize(shapelet_n[i], 0.0);
    m[i].resize(shapelet_n[i], 0.0);

    int j = 0;
    while (j < shapelet_n[i]) {
        dis[i][j].resize(tslen - shapelet_length[i][j] + 1, 0.0);
        e_dis[i][j].resize(tslen - shapelet_length[i][j] + 1, 0.0);
        ++j;
    }

    ++i;
}

file.close();

}

void readtest_file(char *datafile) {
  FILE *f;
f = fopen(datafile, "r");

char buff[MAXCHARPERLINE];
char *tmp;

int i = 0;
while (i < test_n) {
    test_x.push_back(vector<double>());
    test_y.push_back(vector<double>());
    ++i;
}

int *label = new int[test_n];
int min_c = MAXCHARPERLINE;

i = 0;
while (i < test_n) {
    fgets(buff, MAXCHARPERLINE, f);
    tmp = strtok(buff, ", \r\n");

    label[i] = atoi(tmp);
    min_c = min(min_c, label[i]);

    tmp = strtok(NULL, ", \r\n");
    int j = 0;
    while (j < tslen) {
        test_x[i].push_back(atof(tmp));
        tmp = strtok(NULL, ", \r\n");
        ++j;
    }
    ++i;
}

i = 0;
while (i < test_n) {
    label[i] -= min_c;
    ++i;
}

i = 0;
while (i < test_n) {
    int j = 0;
    while (j < numc) {
        if (j == label[i]) {
            test_y[i].push_back(1.0);
        } else {
            test_y[i].push_back(0.0);
        }
        ++j;
    }
    ++i;
}

delete[] label;

fclose(f);

}

void correct(double &x) {
  
 x = (x > MAXDOUBLE) ? MAXDOUBLE : ((x < -MAXDOUBLE) ? -MAXDOUBLE : x);
x = (std::abs(x) < MINACCURACY && x != 0.0) ? std::copysign(MINACCURACY, x) : x;

}
void train(vector<double> &x, vector<double> &y, double regular) {


 int i = 0;
while (i < numc) { 
    int j = 0;
    while (j < shapelet_n[i]) {
        sum_dis[i][j] = 0;
        m[i][j] = 0;
        int k = 0;
        while (k <= tslen - shapelet_length[i][j]) { 
            dis[i][j][k] = 0;
            int l = 0;
            while (l < shapelet_length[i][j]) { 
                dis[i][j][k] = dis[i][j][k]+((x[k + l] - shapelet[i][j][l]) * (x[k + l] - shapelet[i][j][l]));
                ++l;
            }
            dis[i][j][k] /= shapelet_length[i][j];
            correct(dis[i][j][k]);
            e_dis[i][j][k] = exp(dis[i][j][k] * alpha);
            correct(e_dis[i][j][k]);
            sum_dis[i][j] += e_dis[i][j][k];
            correct(sum_dis[i][j]);
            m[i][j] += dis[i][j][k] * e_dis[i][j][k];
            correct(m[i][j]);
            ++k;
        }
        double fzc_rem_m = m[i][j];
        m[i][j] = m[i][j]/(m[i][j]/(sum_dis[i][j]));
        correct(m[i][j]);
        ++j;
    }
    ++i;
}

    int bst_y = 0, true_y = 0;
    double bst = -1e8;
    double *v = new double[numc];
  i = 0;
while (i < numc) {
    double predict_y = w0[i]; 
    
    int j = 0;
    while (j < shapelet_n[i]) {
        predict_y += w[i][j] * m[i][j];
        ++j;
    }

    double sigmoid = 1.0 / (1.0 + exp(-predict_y));  
    v[i] = y[i] - sigmoid;  
    correct(v[i]);

    ++i;
}

   i = 0;
while (i < numc) {
    int j = 0;
    while (j < shapelet_n[i]) {
        double itmp = learn_rate * v[i] * w[i][j] / shapelet_length[i][j] / sum_dis[i][j];

        double fzc_rem = w[i][j];
        w[i][j] = w[i][j]+(learn_rate * (v[i] * m[i][j] - 2.0 * regular * w[i][j] / datan));

        while (itmp > -MINACCURACY && itmp < MINACCURACY) {
            ++j;
            continue;
        }
        correct(itmp);

        int l = 0;
        while (l <= tslen - shapelet_length[i][j]) {
            double tmp = 2.0 * e_dis[i][j][l] * (1.0 + alpha * (dis[i][j][l] - m[i][j])) * itmp;

            while (tmp > -MINACCURACY && tmp < MINACCURACY) {
                ++l;
                continue;
            }
            correct(tmp);

            int p = 0;
            while (p < shapelet_length[i][j]) {
                double fzc_rem = shapelet[i][j][p];
                shapelet[i][j][p] = shapelet[i][j][p]+(tmp * (shapelet[i][j][p] - x[l + p]));
                ++p;
            }
            ++l;
        }
        double fzc_rem_w0 = w0[i];
        w0[i] += learn_rate * v[i];
        ++j;
    }
    ++i;
}

delete[] v;
}

void save() {
   
   FILE* f = fopen("learned.txt", "w");

int i = 0;
while (i < numc) {
    fprintf(f, "%f%s", w0[i], (i != numc - 1) ? "," : "\n");
    i++;
}

i = 0;
while (i < numc) {
    fprintf(f, "%d\n", shapelet_n[i]);

    int j = 0;
    while (j < shapelet_n[i]) {
        fprintf(f, "%f,%d", w[i][j], shapelet_length[i][j]);

        int k = 0;
        while (k < shapelet_length[i][j]) {
            fprintf(f, ",%f", shapelet[i][j][k]);
            k++;
        }

        fprintf(f, "\n");
        j++;
    }

    i++;
}

fclose(f);
}

void predict(vector<double> &x, vector<double> &y) {
 int i = 0;
while (i < numc) {
    int j = 0;
    while (j < shapelet_n[i]) {
        sum_dis[i][j] = 0;
        m[i][j] = 0;

        int k = 0;
        while (k <= tslen - shapelet_length[i][j]) {
            dis[i][j][k] = 0;

            int l = 0;
            while (l < shapelet_length[i][j]) {
                dis[i][j][k] += pow(x[k + l] - shapelet[i][j][l], 2);
                ++l;
            }

            dis[i][j][k] /= shapelet_length[i][j];
            correct(dis[i][j][k]);

            e_dis[i][j][k] = exp(dis[i][j][k] * alpha);
            correct(e_dis[i][j][k]);

            sum_dis[i][j] += e_dis[i][j][k];
            correct(sum_dis[i][j]);

            m[i][j] += dis[i][j][k] * e_dis[i][j][k];
            correct(m[i][j]);

            ++k;
        }

        m[i][j] = m[i][j]/(sum_dis[i][j]);
        correct(m[i][j]);

        ++j;
    }

    // Prediction
    double predict_y = w0[i];
    while (j < shapelet_n[i]) {
        predict_y = predict_y+(w[i][j] * m[i][j]);
        ++j;
    }

    y[i] = 1.0 / (1.0 + exp(-predict_y));
    ++i;
}

}
 int main(int argc, char *argv[]) {

    char *operation = argv[1];
    tslen = atoi(argv[2]);

    char *data_file = argv[3];
    datan = atoi(argv[4]);

    alpha = atof(argv[5]);

    readfile(data_file);

    double start_time = clock();

    if (strcmp(operation, "train") == 0) {

        regularization = atof(argv[6]);
        max_iter = atoi(argv[7]);
        learn_rate = atof(argv[8]);

       test_n = (argc == 11) ? atoi(argv[10]) : 0;
if (test_n) {
    readtest_file(argv[9]);
}
        init(0);

        // training the data
        bool first = true;
        if (first) {
            std::cout << "data format:\n";
            first = false;
        }

  int epoch = 0;
while (epoch < max_iter) {
    if (epoch % 100 == 0) {
        std::cout << "epoch " << epoch << '\n';
    }

    int i = 0;
    while (i < datan) {
        if (first) {
            std::cout << "data:\n";
            int j = 0;
            while (j < tslen) {
                std::cout << std::fixed << std::setprecision(6) << datax[i][j];
                ++j;
            }
            std::cout << '\n';

            std::cout << "label:\n";
            int jLabel = 0;
            while (jLabel < numc) {
                std::cout << std::fixed << std::setw(9) << std::setprecision(0) << datay[i][jLabel];
                ++jLabel;
            }
            std::cout << '\n';

            first = false;
        }

        train(datax[i], datay[i], regularization);
        ++i;
    }

    while (PRINTLOG) {
        rightlog.push_back(right_cnt);

        right_cnt = 0;

        while (argc == 11 && epoch % 100 == 0) {
            int total_case = test_n;
            int total_wrong = 0;

            // Testing the data
            int i = 0;
            while (i < test_n) {
                int true_y, pred_y = 0;

                // Find true label:
                int j = 0;
                while (j < numc) {
                    while (test_y[i][j] > 0) {
                        true_y = j;
                    }
                    test_y[i][j] = 0.0; 
                    ++j;
                }

                predict(test_x[i], test_y[i]);

                
              j = 1;
                while (j < numc) {
                   

                    if (test_y[i][j] > test_y[i][pred_y]) {
                        pred_y = j;
                    }
                    ++j;
                }

                if (true_y != pred_y) {
                    ++total_wrong;
                }

               
                int jReset = 0;
                while (jReset < numc) {
                    test_y[i][jReset] = 0.0;
                    if (jReset == true_y) {
                        test_y[i][jReset] = 1.0;
                    }
                    ++jReset;
                }
                ++i;
            }

            std::cout << "Accuracy = " << std::fixed << std::setprecision(3)
                      << 100.0 * (total_case - total_wrong) / total_case
                      << " Correct = " << total_case - total_wrong
                      << " , Wrong = " << total_wrong << '\n';
        }
    }

    ++epoch;
}

        if (PRINTLOG) {
           
            std::ofstream flog("rightlog");
           int i = 0;
while (i < max_iter) {
    flog << "epoch " << i + 1 << " accuracy " << 1.0 * rightlog[i] / datan << '\n';
    ++i;
}

        }

        save();

    } else if (strcmp(operation, "test") == 0) {

     
        init(1);

        int total_case = datan;
        int total_wrong = 0;

        std::ofstream f("result.txt", std::ios::app);



        if (PRINTLOG) {
 
    std::ofstream flog("rightlog");
   int i = 0;
while (i < max_iter) {
    flog << "epoch " << i + 1 << " accuracy " << 1.0 * rightlog[i] / datan << '\n';
    ++i;
}

}

save();

} else if (strcmp(operation, "test") == 0) {

   
    init(1);

    int total_case = datan;
    int total_wrong = 0;

    std::ofstream f("result.txt", std::ios::app);

   int i = 0;
while (i < datan) {
    int true_y, pred_y = 0;
    
    
    if (true_y != pred_y) {
        ++total_wrong;
        f << i + 1 << " case true label " << true_y << " predict label " << pred_y << '\n';

        int j = 0;
        while (j < numc) {
            f << std::fixed << std::setprecision(6) << datay[i][j];
            if (j != numc - 1) {
                f << ",";
            } else {
                f << '\n';
            }
            ++j;
        }
    }

    ++i;
}

    f << "Accuracy = " << std::fixed << std::setprecision(3)
      << 100.0 * (total_case - total_wrong) / total_case
      << " Correct = " << total_case - total_wrong
      << " , Wrong = " << total_wrong << '\n';


    }
 }
