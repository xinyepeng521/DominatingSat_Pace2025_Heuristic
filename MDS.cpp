#include<vector>
#include<algorithm>
#include<ctime>
#include<set>
#include<stack>
#include<map>
#include<random>
#include<queue>
#include<cstdint>
#include <random>
#include<fstream>
#include<sstream>
#include <csignal>       
#include <atomic>        
#include <unordered_set>
#include<iostream>
using namespace std;
struct HashPair {
    size_t operator()(const pair<uint64_t, uint64_t>& p) const {
        return hash<uint64_t>()(p.first) ^ (hash<uint64_t>()(p.second));
    }
};
clock_t log3_time;
int cutoff_time = 1000;
int seed = 0;
int alpha = 90;
mt19937 gen(seed);
long long nbIter = 90000000000;
int current_X;

volatile atomic<bool> signal_received(false);

struct DegreeCounters {
    vector<int> L;
    vector<int> L2;
    vector<int> L3;

    DegreeCounters(size_t size) : L(size, 0), L2(size, 0), L3(size, 0) {}
    DegreeCounters() = default;
    void reset() {
        fill(L.begin(), L.end(), 0);
        fill(L2.begin(), L2.end(), 0);
        fill(L3.begin(), L3.end(), 0);
    }
};

class DiscretizeVector {
public:
    vector<int>All;
    vector<int>element;
    DiscretizeVector() {}
    DiscretizeVector(int n) {
        All.assign(n, -1);
    }
    auto erase(int ele) -> void {
        if (All[ele] == -1) return;
        int last_index = element.size() - 1;
        int last_ele = element[last_index];
        int del_index = All[ele];

        element[del_index] = last_ele;
        All[last_ele] = del_index;

        All[ele] = -1;
        element.pop_back();
    }
    auto push(int ele) -> void {
        All[ele] = element.size();
        element.push_back(ele);
    }
    auto pop() -> void {
        All[element.back()] = -1;
        element.pop_back();
    }
    auto count(int ele) -> bool {
        return All[ele] != -1;
    }
    auto insert(int ele) -> void {
        if (All[ele] == -1)
            push(ele);
    }

    auto clear() -> void {
        while (element.size() != 0)
            pop();
    }
};
void handle_signal(int sig) {
    signal_received.store(true);
}
class MDSP
{
public:
    vector<vector<int>>HL;
    vector<char>X;
    vector<char>X_minus;
    DiscretizeVector redundant_vertex;
    int best_current_X_M;
    DiscretizeVector X2;
    DiscretizeVector X_M;
    DegreeCounters degrees;
    DiscretizeVector first_model_remove;
    bool stop = false;
    vector<vector<int>>vertex;
    int now_time = 1;
    clock_t log2_time = clock();
    clock_t log_time = clock();
    double best_time = 0;
    clock_t curr_time;
    clock_t begin_time;
    int ini_num;
    int x_minus;
    int x;
    vector<int>age;
    // ﾹ￾ￏﾣￏ￠ﾹ￘ﾳ￉ￔﾱ
    vector<uint64_t> base1;
    vector<uint64_t> base2;
    uint64_t hash1;
    uint64_t hash2;
    unordered_set<pair<uint64_t, uint64_t>, HashPair> tabu_set;
    deque<pair<uint64_t, uint64_t>> tabu_list;
    int tabu_tenure = 7;

    vector<int> best_result;
    void output_result() {
        cout << best_result.size() << endl;
        for (auto i : best_result)
            cout << i + 1 << endl;
    }
    void ReductionRule(int i) {
        if (degrees.L[i] <= 1 && X_minus[i])
        {
            int num_x = i; bool flg = true; int times = 0;
            for (auto j : HL[i])
            {
                int score = X_minus[j] + degrees.L[j];
                if (score > 1)
                {
                    times++;
                    num_x = j;
                    if (times == 2)
                    {
                        flg = false;
                        break;
                    }
                }

            }
            if (flg)
            {
                remove_to_X(num_x);
                for (auto j : HL[num_x])
                {
                    if (!X[j])
                    {
                        for (auto k : HL[j])
                            ReductionRule(k);
                    }
                }
            }
        }
    }
    bool  SelectReductionInsert(int i, int& insert) {
        if (degrees.L[i] <= 1 && X_minus[i])
        {
            int num_x = i; bool flg = true; int times = 0;
            for (auto j : HL[i])
            {
                int score = X_minus[j] + degrees.L[j];
                if (score > 1)
                {
                    times++;
                    num_x = j;
                    if (times == 2)
                    {
                        flg = false;
                        break;
                    }
                }
            }
            if (flg)
            {
                insert = num_x;
                return true;
            }
            return false;
        }
        return false;
    }
    int vary_length = 25;
    bool model = false;
    void calc_L();
    void remove_to_X(int u);
    void ini_greed();
    void calc_L2();
    void calc_L3();
    auto improve(long long iter) -> void;
    void remove_vertex(int num, long long iter);
    int select_remove(int& a, int& b, long long& iter);
    void apply_move(int& a, int& b, long long iter);
    void DemDS();
    int select_insert(int& ins);
    void insert_vertex(int num);
    int  select_insert2(int& ins);
    void make_rand()
    {
        for (int i = ini_num; i < min(ini_num + vary_length, int(X2.element.size())); i++)
        {
            int rd = gen() % (X2.element.size() - ini_num) + ini_num;
            swap(X2.element[i], X2.element[rd]);
            swap(X2.All[X2.element[i]], X2.All[X2.element[rd]]);
        }
    }
    void f1_1(int p)
    {
        for (auto i : HL[p])
        {
            if (degrees.L3[i] == 0 && X[i] && degrees.L2[i] != 0)
                redundant_vertex.erase(i);
            degrees.L3[i]++;
        }
    }
    void f1_2(int p)
    {
        X_minus[p] = 0;
        X_M.erase(p);
        for (auto i : HL[p])
        {
            degrees.L[i]--;
        }
    }
    void f2_1(int p)
    {
        for (auto i : HL[p])
        {
            degrees.L3[i]--;
            if (degrees.L3[i] == 0 && X[i] && degrees.L2[i] != 0)
                redundant_vertex.insert(i);
        }
    }
    void f2_2(int p)
    {
        X_minus[p] = 1;
        X_M.push(p);
        for (auto i : HL[p])
        {
            degrees.L[i]++;
        }
    }
    void f1(int p)
    {
        f1_1(p);
        f1_2(p);
    }
    void f2(int p)
    {
        f2_1(p);
        f2_2(p);
    }
    void f3(int p)
    {
        X[p] = 1;
        x++;
        f1_2(p);
        for (auto i : HL[p])
        {
            if (degrees.L3[i] == 0 && X[i] && degrees.L2[i] == 0)
                redundant_vertex.insert(i);
            degrees.L2[i]++;

            if (X_minus[i])
            {
                f1(i);
            }
            if (degrees.L2[i] == 2 && !X[i])
            {
                for (auto j : HL[i])
                    if (X2.All[j] >= ini_num && j != p)
                        first_model_remove.insert(j);
                f2_1(i);
            }
        }


    }
    void f4(int p)
    {
        X[p] = 1;
        x++;
        for (auto i : HL[p])
        {
            if (degrees.L3[i] == 0 && X[i] && degrees.L2[i] == 0)
                redundant_vertex.insert(i);
            degrees.L2[i]++;
        }
        if (degrees.L2[p] == 1)
        {
            f2_1(p);
        }
        for (auto i : HL[p])
        {
            if (X_minus[i])
            {
                f1(i);
            }
            if (degrees.L2[i] == 2 && !X[i])
            {
                for (auto j : HL[i])
                    if (X2.All[j] >= ini_num && j != p)
                        first_model_remove.insert(j);
                f2_1(i);
            }
        }

    }
    void f5(int p) {
        X[p] = 0;
        x--;
        for (auto i : HL[p])
        {
            degrees.L2[i]--;
            if (degrees.L3[i] == 0 && X[i] && degrees.L2[i] == 0)
                redundant_vertex.erase(i);
        }
        if (degrees.L2[p] == 1)
        {
            f1_1(p);
        }
        for (auto i : HL[p])
        {
            if (degrees.L2[i] == 0 && !X[i] && !X_minus[i])
            {
                f2(i);
            }
            if (degrees.L2[i] == 1 && !X[i])
            {
                f1_1(i);
            }

        }

    }
    void f6(int p) {

        X[p] = 0;
        x--;
        for (auto i : HL[p])
        {
            degrees.L2[i]--;
            if (degrees.L3[i] == 0 && X[i] && degrees.L2[i] == 0)
                redundant_vertex.erase(i);
        }
        f2_2(p);
        for (auto i : HL[p])
        {
            if (degrees.L2[i] == 0 && !X[i] && !X_minus[i])
            {
                f2(i);
            }
            if (degrees.L2[i] == 1 && !X[i])
            {
                f1_1(i);
            }
            if (degrees.L2[i] == 0 && X[i])
            {
                X2.erase(i);
            }
        }
    }
    friend istream& operator >>(istream& istr, MDSP& g)
    {
        int V, value, E, value2;
        string temp;
        while (istr.peek() != 'p') {
            getline(istr, temp);
        }
        istr >> temp >> temp;
        istr >> V >> E;
        g.X_minus.assign(V, 1);
        g.X.assign(V, 0);
        g.HL.assign(V, vector<int>());
        g.degrees = DegreeCounters(V);
        g.X2 = DiscretizeVector(V);
        g.X_M = DiscretizeVector(V);
        g.x_minus = V;
        g.x = 0;
        g.age.assign(V, 0);
        g.first_model_remove = DiscretizeVector(V);
        g.redundant_vertex = DiscretizeVector(V);

        g.base1.resize(V);
        g.base2.resize(V);
        mt19937_64 rng(seed);
        for (int i = 0; i < V; i++) {
            g.base1[i] = rng();
            g.base2[i] = rng();
        }
        g.hash1 = 0;
        g.hash2 = 0;

        
        for (int i = 0; i < E; i++)
        {
            istr >> value >> value2;

            g.HL[value - 1].push_back(value2 - 1);
            g.HL[value2 - 1].push_back(value - 1);
        }
        return istr;
    }
};
void MDSP::calc_L()
{
    degrees.L.assign(degrees.L.size(), 0);
    for (int i = 0; i < degrees.L.size(); i++)
    {
        for (auto& j : HL[i])
        {
            if (X_minus[j])
                degrees.L[i]++;
        }
    }
}
void MDSP::remove_to_X(int u)
{
    X[u] = 1;
    x++;
    if (X_minus[u] == 1)
    {
        X_minus[u] = 0;
        x_minus--;
        for (auto& i : HL[u])
            degrees.L[i]--;
    }

    for (auto& i : HL[u])
    {
        if (X_minus[i] == 1)
        {
            X_minus[i] = 0;
            x_minus--;
            for (auto& j : HL[i])
                degrees.L[j]--;
        }
    }
}
void MDSP::ini_greed()
{
    calc_L();
    for (int i = 0; i < HL.size(); i++)
    {
        if (HL[i].size() == 0)
        {
            X_minus[i] = 0;
            x++;
            x_minus--;
            X[i] = 1;
        }
        if (HL[i].size() == 1 && X_minus[i] == 1)
        {
            remove_to_X(HL[i][0]);
        }
        if (HL[i].size() == 2 && X_minus[i] == 1)
        {
            int sec = HL[i][0], thd = HL[i][1];
            if (HL[sec].size() == 2 && (HL[sec][0] == thd || HL[sec][1] == thd))
                remove_to_X(thd);
            else if (HL[thd].size() == 2 && (HL[thd][0] == sec || HL[thd][1] == sec))
                remove_to_X(sec);
        }
    }
    for (int i = 0; i < HL.size(); i++)
    {
        ReductionRule(i);
    }
    ini_num = x;
    for (int i = 0; i < HL.size(); i++)
    {
        if (X[i])
        {
            X2.push(i);
        }
    }
    int sc_max = -1;
    int ins = -1;
    for (int i = 0; i < degrees.L.size(); i++)
    {
        int sc_tmp = degrees.L[i] + X_minus[i];
        if (sc_tmp > sc_max)
        {
            sc_max = sc_tmp; ins = i;
        }
    }
    if (sc_max > 0)
        remove_to_X(ins);
    while (x_minus != 0)
    {
        int i = ins + 1;
        for (; i < degrees.L.size(); i++)
        {
            int sc_tmp = degrees.L[i] + X_minus[i];
            if (sc_tmp == sc_max)
            {
                ins = i;
                break;
            }
        }
        if (i != degrees.L.size())
        {
            remove_to_X(i);

        }
        else
        {
            /* if (sc_max < 4)
             {
                 for (int i = 0; i < HL.size(); i++)
                 {
                     ReductionRule(i);
                 }
             }*/
            sc_max = -1;
            for (int i = 0; i < degrees.L.size(); i++)
            {
                int sc_tmp = degrees.L[i] + X_minus[i];
                if (sc_tmp > sc_max)
                {
                    sc_max = sc_tmp; ins = i;
                }
            }
            if (sc_max > 0)
            {
                remove_to_X(ins);
            }
        }
    }

    for (int i = 0; i < HL.size(); i++)
    {
        if (X[i])
        {
            if (X2.All[i] == -1)
            {
                X2.push(i);
            }
        }
    }
    hash1 = 0; hash2 = 0;
    for (int i = 0; i < X.size(); i++) {
        if (X[i]) {
            hash1 += base1[i];
            hash2 += base2[i];
        }
    }
}
void MDSP::calc_L2()
{
    degrees.L2.assign(degrees.L2.size(), 0);
    for (int i = 0; i < degrees.L2.size(); i++)
    {
        for (auto& j : HL[i])
            if (X[j])
                degrees.L2[i]++;
    }
}
void MDSP::calc_L3()
{
    degrees.L3.assign(degrees.L3.size(), 0);
    for (int i = 0; i < degrees.L3.size(); i++)
    {
        if (degrees.L2[i] == 1 && !X[i])
        {
            for (auto& j : HL[i])
            {
                degrees.L3[j]++;
            }
        }
    }
}
void MDSP::insert_vertex(int i)
{

    X2.push(i);

    if (X_minus[i])
        f3(i);
    else
        f4(i);
    if (degrees.L3[i] == 0 && degrees.L2[i] != 0)
        redundant_vertex.insert(i);
    hash1 += base1[i];
    hash2 += base2[i];
}
void MDSP::remove_vertex(int i, long long iter)
{
    if (degrees.L3[i] == 0 && degrees.L2[i] != 0)
        redundant_vertex.erase(i);
    X2.erase(i);

    if (degrees.L2[i] == 0)
        f6(i);
    else
        f5(i);
    age[i] = iter;
  

}
int MDSP::select_insert(int& ins) {
    first_model_remove.clear();
    int best_candidate = -1;          // 非禁忌最佳候选
    int best_score = -1;
    int best_age = INT_MAX;
    int best_tabu_candidate = -1;     // 禁忌最佳候选
    int best_tabu_score = -1;
    int best_tabu_age = INT_MAX;

    for (auto i : X_M.element) {
        for (auto j : HL[i]) {
            int score_now = degrees.L[j] + (X_minus[j] == 1 ? 1 : 0);
            uint64_t new_hash1 = hash1 + base1[j];
            uint64_t new_hash2 = hash2 + base2[j];
            bool is_tabu = (tabu_set.find({ new_hash1, new_hash2 }) != tabu_set.end());
            if (is_tabu) {
                // 禁忌候选：记录分数最高且年龄最小的
                if (score_now > best_tabu_score ||
                    (score_now == best_tabu_score && age[j] < best_tabu_age)) {
                    best_tabu_candidate = j;
                    best_tabu_score = score_now;
                    best_tabu_age = age[j];
                }
            }
            else {
                // 非禁忌候选：记录分数最高且年龄最小的
                if (score_now > best_score ||
                    (score_now == best_score && age[j] < best_age)) {
                    best_candidate = j;
                    best_score = score_now;
                    best_age = age[j];
                }
            }
        }
        // 处理当前顶点i（自身）
        int score_now_i = degrees.L[i] + (X_minus[i] == 1 ? 1 : 0);
        uint64_t new_hash1_i = hash1 + base1[i];
        uint64_t new_hash2_i = hash2 + base2[i];
        bool is_tabu_i = (tabu_set.find({ new_hash1_i, new_hash2_i }) != tabu_set.end());
        if (is_tabu_i) {
            if (score_now_i > best_tabu_score ||
                (score_now_i == best_tabu_score && age[i] < best_tabu_age)) {
                best_tabu_candidate = i;
                best_tabu_score = score_now_i;
                best_tabu_age = age[i];
            }
        }
        else {
            if (score_now_i > best_score ||
                (score_now_i == best_score && age[i] < best_age)) {
                best_candidate = i;
                best_score = score_now_i;
                best_age = age[i];
            }
        }
    }

    // 优先选择非禁忌候选，若无则选择禁忌候选
    if (best_candidate != -1) {
        ins = best_candidate;
        return best_score;
    }
    else if (best_tabu_candidate != -1) {
        ins = best_tabu_candidate;
        return best_tabu_score;
    }
    else {
        // 保底选择（理论上不会触发）
        if (X_M.element.empty()) {
            ins = -1;
            return -1;
        }
        int i = X_M.element[0];
        ins = i;
        return degrees.L[i] + (X_minus[i] == 1 ? 1 : 0);
    }
}
int MDSP::select_insert2(int& ins) {
    first_model_remove.clear();
    int Xm = X_M.element[gen() % X_M.element.size()];
    ins = Xm; int ins_score = degrees.L[Xm] + (X_minus[Xm] == 1 ? 1 : 0);
    int rd_ins = 2;
    for (auto i : HL[Xm])
    {
        if (X[i])
            continue;
        int score_now = degrees.L[i] + (X_minus[i] == 1 ? 1 : 0);
        //uint64_t new_hash1 = hash1 + base1[i];
        //uint64_t new_hash2 = hash2 + base2[i];
        //if (tabu_set.find({ new_hash1, new_hash2 }) != tabu_set.end()) {
        //    continue; // ￌ￸ﾹ�ﾽ￻ﾼ￉ﾽ￢
        //}
        if (score_now > ins_score)
        {
            ins = i;
            ins_score = score_now;
            rd_ins = 2;
        }
        else if (score_now == ins_score)
        {
            if (age[i] < age[ins])
            {
                ins = i;
                ins_score = score_now;
                rd_ins = 2;
            }
            else if (age[i] == age[ins] && !gen() % (rd_ins++))
            {
                ins = i;
                ins_score = score_now;
            }
        }
    }
    return ins_score;
}
int MDSP::select_remove(int& a, int& b, long long& iter) {
    int rd_1 = 2;
    int rem = -1;
    int rem_score = 99999999;
    int num = 0;
    if (X_M.element.size() != 0)
    {
        if (b != -1 && degrees.L2[b] != 0)
        {
            for (auto i : HL[b])
            {
                if (X2.All[i] < ini_num)
                    continue;
                if (X[i])
                {
                    int score_rem = degrees.L3[i];
                    if (score_rem + X_M.element.size() < best_current_X_M)
                    {
                        a = i;
                        return 1;
                    }
                    num++;
                    uint64_t new_hash1 = hash1 - base1[i];
                    uint64_t new_hash2 = hash2 - base2[i];
                    if (tabu_set.find({ new_hash1, new_hash2 }) != tabu_set.end()) {
                        continue; // ￌ￸ﾹ�ﾽ￻ﾼ￉ﾽ￢
                    }

                    if (score_rem < rem_score)
                    {
                        rem_score = degrees.L3[i];
                        rem = i;
                        rd_1 = 2;
                    }
                    else if (score_rem == rem_score)
                    {
                        if (age[i] < age[rem])
                        {
                            rem = i;
                            rd_1 = 2;
                        }
                        else if (age[i] == age[rem] && !gen() % (rd_1++))
                        {
                            rem = i;
                        }
                    }
                }
            }
        }
        for (auto i : first_model_remove.element)
        {
            if (X2.All[i] < ini_num)
                continue;
            if (X[i])
            {
                int score_rem = degrees.L3[i] + (degrees.L2[i] == 0 ? 1 : 0);
                if (score_rem + X_M.element.size() < best_current_X_M)
                {
                    a = i;
                    return 1;
                }
                num++;
                uint64_t new_hash1 = hash1 - base1[i];
                uint64_t new_hash2 = hash2 - base2[i];
                if (tabu_set.find({ new_hash1, new_hash2 }) != tabu_set.end()) {
                    continue; // ￌ￸ﾹ�ﾽ￻ﾼ￉ﾽ￢
                }

                if (score_rem < rem_score)
                {
                    rem_score = degrees.L3[i];
                    rem = i;
                    rd_1 = 2;
                }
                else if (score_rem == rem_score)
                {
                    if (age[i] < age[rem])
                    {
                        rem = i;
                        rd_1 = 2;
                    }
                    else if (age[i] == age[rem] && !(gen() % (rd_1++)))
                    {
                        rem = i;
                    }
                }
            }
        }
    }
    else
        num = 15;
    num = max(3, num);
    vary_length = num;
    for (int j = ini_num; j < min(ini_num + num, int(X2.element.size() - 1)); j++)
    {
        int i = X2.element[j];

        int score_rem = degrees.L3[i] + (degrees.L2[i] == 0 ? 1 : 0);

        if (score_rem < rem_score)
        {
            rem_score = degrees.L3[i];
            rem = i;
            rd_1 = 2;
        }
        else if (score_rem == rem_score)
        {
            if (age[i] < age[rem])
            {
                rem = i;
                rd_1 = 2;
            }
            else if (age[i] == age[rem] && !gen() % (rd_1++))
            {
                rem = i;
            }
        }
    }
    a = rem;
    return 1;
}
void MDSP::apply_move(int& a, int& b, long long  iter)
{
    remove_vertex(a, iter);
    age[b] = iter;
    age[a] = iter;

}
void MDSP::DemDS()
{
    ini_greed();
    calc_L2();
    calc_L3();
    best_result = X2.element;
    long long iter = 0;
    curr_time = clock();
    bool record = X2.element.size() < 50000;
    int flag = 1;
    int iter2 = 0;
    int flag3 = 1;
    int a, b = -1, c, d;
    tabu_tenure = 10/*sqrt(X2.element.size()) / 5*/;
    for (int i = 0; i < X2.element.size(); i++)
    {
        if (degrees.L3[X2.element[i]] == 0 && degrees.L2[X2.element[i]] != 0)
            remove_vertex(X2.element[i], iter);
    }

    while (iter < nbIter) {
        if (X2.element.size() == ini_num) {
            output_result();
            return;
        }
        if (X_M.element.size() == 0)
        {
            b = -1;
            flag = 1;
            if (redundant_vertex.element.size() != 0)
                improve(iter);
            if (X2.element.size() < best_result.size())
            {
                /* if (iter2 > X2.element.size() / 3)
                     improve(0);*/

                best_result = X2.element;
                best_current_X_M = X2.element.size();
            }
            if (flag3)
                iter2 = 0;
        }

        best_current_X_M = min(best_current_X_M, int(X_M.element.size()));
        select_remove(a, b, iter);
        remove_vertex(a, iter);
        hash1 -= base1[a];
        hash2 -= base2[a];
        pair<uint64_t, uint64_t> current_hash = { hash1, hash2 };
        if (tabu_set.find(current_hash) == tabu_set.end()) {
            tabu_set.insert(current_hash);
            tabu_list.push_back(current_hash);
            if (tabu_list.size() > tabu_tenure) {
                auto old_hash = tabu_list.front();
                tabu_list.pop_front();
                tabu_set.erase(old_hash);
            }
        }
        make_rand();
        if (X_M.element.size() != 0 && X2.element.size() < best_result.size() - 1)
        {
            if (X_M.element.size() > 20)
            {
                select_insert(b);
                insert_vertex(b);
                age[b] = iter;
                b = -1;
            }
            else
            {
                select_insert2(b);
                insert_vertex(b);
                age[b] = iter;
                //Tabu[b] = iter + 15 + gen() % 10;
            }
        }
        if (X_M.element.size() != 0)
        {
            if (X2.element.size() < best_result.size() - 1 || redundant_vertex.element.size() != 0/*||(X2.element.size() == best_result.size() - 1&&redundant_vertex!=0)*/)
            {
                int insert_v;
                if (select_insert(insert_v) == X_M.element.size())
                {
                    insert_vertex(insert_v);
                    age[insert_v] = iter;
                    b = insert_v;
                    /* cout << "* ";*/
                     /*Tabu[b] = iter+10+gen()%10;*/
                }

            }
        }
       
        if (redundant_vertex.element.size() > 1 && X_M.element.size() < 20)
            improve(iter);
        while (X_M.element.size() != 0 && X2.element.size() < best_result.size() - 1)
        {
            {
                select_insert(b);
                insert_vertex(b);
                age[b] = iter;
                b = -1;
                if (redundant_vertex.element.size() != 0)
                    improve(iter);
            }
        }
        if (iter % 10 == 0)
        {
               /*curr_time = clock();
               if ((curr_time - log2_time) / CLOCKS_PER_SEC > 1)
               {
                   cout << "#" << iter << "  #" << X_M.element.size() << "  " << best_result.size()<<"  #"<<redundant_vertex.element.size() << "   " << X2.element.size() << endl;
                   log2_time = clock();
               }*/

            if (signal_received.load()) break;
        }
        iter++;
        iter2++;
    }


    output_result();
}
int main()
{
    log3_time = clock();
    signal(SIGTERM, handle_signal);
    /*signal(SIGINT, handle_signal);*/

    ios::sync_with_stdio(false);
    cin.tie(0); cout.tie(0);
    seed = 1;
    alpha = 100;
    gen.seed(seed);

    MDSP a;
   /* ifstream infile("heuristic_080.gr");
    infile >> a;*/
    cin >> a;
    a.DemDS();
}
auto MDSP::improve(long long iter)->void {
    while (redundant_vertex.element.size() != 0)
    {
        int i = redundant_vertex.element[redundant_vertex.element.size() - 1];
        remove_vertex(i, iter);
    }
}
