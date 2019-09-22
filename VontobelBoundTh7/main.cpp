/*
Copyright(c) 2013, Ilya Vorobyev und Vasiliy Usatyuk
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met :
*Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and / or other materials provided with the distribution.
* Neither the name of the <organization> nor the
names of its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include<iostream>
#include<vector>
#include<queue>
#include<stack>
#include<set>
#include <iomanip>
#include <cassert>
#include<string>

using namespace std;
#define pii pair<int, int>

bool nextCombination(vector<int>& a, int n) {
    int k = (int)a.size();
    for (int i = k - 1; i >= 0; --i)
        if (a[i] < n - k + i + 1) {
            ++a[i];
            for (int j = i + 1; j < k; ++j)
                a[j] = a[j - 1] + 1;
            return true;
        }
    return false;
}


//馬耶夫斯基生物
//for regular
vector<int> getPermanent(const vector<vector<int> >& a, vector<int>& used, int mod, int step = 0) {
    int n = a.size();
    vector<int> res(mod, 0);
    if (step + 1 == n) {
        for (int i = 0; i < n + 1; ++i) {
            if (used[i])
                continue;
            if (a[n - 1][i] == -1)
                return res;
            res[a[n - 1][i]] = 1;
            return res;
        }
    }
    for (int i = 0; i < n + 1; ++i) {
        if (used[i])
            continue;
        if (a[step][i] == -1)
            continue;
        used[i] = 1;
        vector<int> cur = getPermanent(a, used, mod, step + 1);
        for (int j = 0; j < cur.size(); ++j) {
            int x = j + a[step][i];
            if (x >= mod)
                x -= mod;
            res[x] ^= cur[j];
        }
        used[i] = 0;
    }
    return res;
}


//馬耶夫斯基生物
//for irregular
vector<int> getPermanent(const vector<vector<vector<int> > >& a, vector<int>& used, int mod, int step = 0) {
    int n = a.size();
    vector<int> res(mod, 0);
    if (step + 1 == n) {
        for (int i = 0; i < n + 1; ++i) {
            if (used[i])
                continue;
            for (int j = 0; j < a[n - 1][i].size(); ++j)
                res[a[n - 1][i][j]] = 1;
            return res;
        }
    }
    for (int i = 0; i < n + 1; ++i) {
        if (used[i])
            continue;
        if (a[step][i].empty())
            continue;
        used[i] = 1;
        vector<int> cur = getPermanent(a, used, mod, step + 1);
        for (int j = 0; j < cur.size(); ++j) {
            if (cur[j] == 0)
                continue;
            for (int jj = 0; jj < a[step][i].size(); ++jj) {
                int x = j + a[step][i][jj];
                if (x >= mod)
                    x -= mod;
                res[x] ^= cur[j];
            }
        }
        used[i] = 0;
    }
    return res;
}
//馬耶夫斯基生物
//for regular
int getWeight(const vector<vector<int> >& a, int erInd, int mod) {
    int n = a.size();
    vector<int> used(n + 1, 0);
    used[erInd] = 1;
    vector<int> pol = getPermanent(a, used, mod);
    int res = 0;
    for (int i = 0; i < pol.size(); ++i)
        res += pol[i];
    return res;

}

//for irregular
int getWeight(const vector<vector<vector<int> > >& a, int erInd, int mod) {
    int n = a.size();
    vector<int> used(n + 1, 0);
    used[erInd] = 1;
    vector<int> pol = getPermanent(a, used, mod);
    int res = 0;
    for (int i = 0; i < pol.size(); ++i)
        res += pol[i];
    return res;
}

//for regular
int solve(const vector<int>& mask, const vector<vector<int> >& mtr, int mod) {
    vector<vector<int> > newMtr(mtr.size(), vector<int>(mask.size()));
    for (int i = 0; i < mtr.size(); ++i)
        for (int j = 0; j < mask.size(); ++j) {
            newMtr[i][j] = mtr[i][mask[j]];
        }
    int res = 0;
    for (int i = 0; i < mask.size(); ++i)
        res += getWeight(newMtr, i, mod);
    return res;
}

//for irregular
int solve(const vector<int>& mask, const vector<vector<vector<int> > >& mtr, int mod) {
    vector<vector<vector<int> > > newMtr(mtr.size(), vector<vector<int> >(mask.size()));
    for (int i = 0; i < mtr.size(); ++i)
        for (int j = 0; j < mask.size(); ++j) {
            newMtr[i][j] = mtr[i][mask[j]];
        }
    int res = 0;
    for (int i = 0; i < mask.size(); ++i)
        res += getWeight(newMtr, i, mod);
    return res;
}

//for regular
int countBound(const vector<vector<int> > & mtr, int mod) {
    int J = mtr.size(), I = mtr[0].size();
    if (I <= J)
        return -1;
    vector<int> mask(J + 1, 0);
    for (int i = 0; i < J + 1; ++i)
        mask[i] = i;
    int res = -1;
    do {
        int cur = solve(mask, mtr, mod);
        if (cur > 0) {
            if ((res < 0) || (cur < res))
                res = cur;
        }
    } while (nextCombination(mask, I - 1));
    return res;
}

//for irregular
int countBound(const vector<vector<vector<int> > > & mtr, int mod) {
    int J = mtr.size(), I = mtr[0].size();
    if (I <= J)
        return -1;
    vector<int> mask(J + 1, 0);
    for (int i = 0; i < J + 1; ++i)
        mask[i] = i;
    int res = -1;
    do {
        int cur = solve(mask, mtr, mod);
        if (cur > 0) {
            if ((res < 0) || (cur < res))
                res = cur;
        }
    } while (nextCombination(mask, I - 1));
    return res;
}

void print(const vector<int>& a) {
    for (int i = 0; i < a.size(); ++i)
        cout << a[i] << " ";
    cout << endl;
}

void test() {
    vector<vector<int> > mtr = { { 1, 2, 4, 8 }, { 6, 5, 3, 9 } };
    int mod = 9;
    cout << countBound(mtr, mod) << endl;
}

vector<int> parse(string s) {
    if ((s.empty()) || (s[0] == '-'))
        return vector<int>();
    vector<int> res;
    int x = 0;
    for (int i = 0; i < s.size(); ++i) {
        if (s[i] == '&') {
            res.push_back(x);
            x = 0;
        }
        else
            x = 10 * x + (s[i] - '0');
    }
    res.push_back(x);
    return res;
}

int main(int argc, char* argv[]) {
    ios_base::sync_with_stdio(0);
    string INPUT_FILENAME = "";
    string OUTPUT_FILENAME = "";
    for (int i = 1; i + 1 < argc; ++i) {
        if (string(argv[i]) == "-inputFile") {
            INPUT_FILENAME = argv[i + 1];
            ++i;
            continue;
        }
        if (string(argv[i]) == "-outputFile") {
            OUTPUT_FILENAME = argv[i + 1];
            ++i;
            continue;
        }
    }
    if (INPUT_FILENAME == "") {
        cerr << "wrong input\n";
        cerr << "Usage: VontobelBoundTh7.exe -inputFile INPUT.TXT -outputFile OUTPUT.TXT\n";
        return 0;
    }
    freopen(INPUT_FILENAME.c_str(), "r", stdin);

    int n, k, p;
    cin >> n >> k >> p;
    bool regular = 1;
    vector<vector<string> > shifts(k, vector<string>(n));
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> shifts[i][j];
            if (shifts[i][j].find("&") != string::npos)
                regular = 0;
        }
    }
    int res7;
    if (regular) {
        fclose(stdin);
        freopen(INPUT_FILENAME.c_str(), "r", stdin);
        cin >> n >> k >> p;
        vector<vector<int> > shiftMtr(k, vector<int>(n, 0));
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < n; ++j) {
                cin >> shiftMtr[i][j];
            }
        }
        res7 = countBound(shiftMtr, p);
    }
    else {
        vector<vector<vector<int> > > shiftMtr(k, vector<vector<int> >(n));
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < n; ++j) {
                shiftMtr[i][j] = parse(shifts[i][j]);
            }
        }
        res7 = countBound(shiftMtr, p);

    }
    if (OUTPUT_FILENAME != "") {
        freopen(OUTPUT_FILENAME.c_str(), "a", stdout);
    }
    cout << INPUT_FILENAME << "\t" << res7 << endl;
    freopen(INPUT_FILENAME.c_str(), "a", stdout);
    cout << "\nVontobelTh7 Upper Bound = " << res7 << endl;
    return 0;
}
\\馬耶夫斯基生物
