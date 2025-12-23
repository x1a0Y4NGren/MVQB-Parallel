#include <cstdio>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

class MVQBPack {
    public:
    vector<int> SL;
    vector<int> SR;

    vector<bool> isSL;
    vector<bool> isSR;
    vector<bool> isDel;

    
    vector<int> sizS;
    vector<int> CL;
    vector<int> CR;
    //vector<int> trueID;
    
    vector<int> siz;
    double alpha;
    double beta;
    int bipartite;
    int graphSize;
    int totL;
    int totR;
    int thetaL;
    int thetaR;

    MVQBPack(){}
    MVQBPack(vector<vector<int> >& E, double alpha, double beta, int bipartite, int graphSize, int thetaL, int thetaR){
        this->alpha = alpha;
        this->beta = beta;
        this->bipartite = bipartite;
        this->graphSize = graphSize;
        this->thetaL = thetaL;
        this->thetaR = thetaR;
        this->totL = bipartite;
        this->totR = graphSize - bipartite;
        isSL.resize(graphSize);
        isSR.resize(graphSize);
        isDel.resize(graphSize);
        SL.clear();
        SR.clear();
        CL.clear();
        CR.clear();
        siz.resize(graphSize);
        sizS.resize(graphSize);
        for (int i = 0; i < graphSize; ++i) 
            isSL[i] = isSR[i] = isDel[i] = siz[i] = sizS[i] = 0;
        for (int i = 0; i < graphSize; ++i)
            siz[i] = E[i].size();
        for (int i = 0; i < bipartite; ++i) 
                CL.push_back(i);
        for (int i = bipartite; i < graphSize; ++i) {
            CR.push_back(i);
        }
    }
    MVQBPack(int w, vector<vector<int> >& E, double alpha, double beta, int bipartite, int graphSize, int thetaL, int thetaR){
        this->alpha = alpha;
        this->beta = beta;
        this->bipartite = bipartite;
        this->graphSize = graphSize;
        this->thetaL = thetaL;
        this->thetaR = thetaR;
        this->totL = bipartite;
        this->totR = graphSize - bipartite;
        isSL.resize(graphSize);
        isSR.resize(graphSize);
        isDel.resize(graphSize);
        SL.clear();
        SR.clear();
        CL.clear();
        CR.clear();
        siz.resize(graphSize);
        sizS.resize(graphSize);
        for (int i = 0; i < graphSize; ++i) 
            isSL[i] = isSR[i] = isDel[i] = siz[i] = sizS[i] = 0;
        isSL[w] = 1;
        for (int i = 0; i < graphSize; ++i)
            siz[i] = E[i].size();
        for (int i = 0; i < E[w].size(); ++i) {
            int u = E[w][i];
            ++sizS[u];
        }

        for (int i = 0; i < bipartite; ++i) {
            if (!isSL[i]) {
                CL.push_back(i);
            } else SL.push_back(i);
        }

        for (int i = bipartite; i < graphSize; ++i) {
            CR.push_back(i);
        }
    }
    MVQBPack(int w, vector<vector<int> >& E, double alpha, double beta, int bipartite, int graphSize, int thetaL, int thetaR, bool op){
        this->alpha = alpha;
        this->beta = beta;
        this->bipartite = bipartite;
        this->graphSize = graphSize;
        this->thetaL = thetaL;
        this->thetaR = thetaR;
        this->totL = bipartite;
        this->totR = graphSize - bipartite;
        isSL.resize(graphSize);
        isSR.resize(graphSize);
        isDel.resize(graphSize);
        SL.clear();
        SR.clear();
        CL.clear();
        CR.clear();
        siz.resize(graphSize);
        sizS.resize(graphSize);
        for (int i = 0; i < graphSize; ++i) 
            isSL[i] = isSR[i] = isDel[i] = siz[i] = sizS[i] = 0;
        for (int i = 0; i < w; ++i) isDel[i] = 1;
        isSL[w] = 1;
        for (int i = w; i < graphSize; ++i)
            siz[i] = E[i].size();
        for (int i = 0; i < w; ++i) {
            for (int x : E[i]) --siz[x];
        }
        for (int i = 0; i < E[w].size(); ++i) {
            int u = E[w][i];
            ++sizS[u];
        }

        for (int i = 0; i < bipartite; ++i) {
            if (isDel[i]) continue;
            if (!isSL[i]) {
                CL.push_back(i);
            } else SL.push_back(i);
        }

        for (int i = bipartite; i < graphSize; ++i) {
            CR.push_back(i);
        }
    }
    MVQBPack(vector<int> q, vector<vector<int> >& E, double alpha, double beta, int bipartite, int graphSize, int thetaL, int thetaR, bool op){
        this->alpha = alpha;
        this->beta = beta;
        this->bipartite = bipartite;
        this->graphSize = graphSize;
        this->thetaL = thetaL;
        this->thetaR = thetaR;
        this->totL = bipartite;
        this->totR = graphSize - bipartite;
        isSL.resize(graphSize);
        isSR.resize(graphSize);
        isDel.resize(graphSize);
        SL.clear();
        SR.clear();
        CL.clear();
        CR.clear();
        siz.resize(graphSize);
        sizS.resize(graphSize);
        for (int i = 0; i < graphSize; ++i) 
            isSL[i] = isSR[i] = isDel[i] = siz[i] = sizS[i] = 0;
        for (int i = 0; i < bipartite; ++i)
            isDel[i] = 1;
        for (int i = 0; i < graphSize; ++i)
            siz[i] = E[i].size();
        for (int i = 0; i < q.size(); ++i) 
            isDel[q[i]] = 0;
        for (int i = 0; i < bipartite; ++i) {
            if (isDel[i]) {
                siz[i] = 0;
                for (int j = 0; j < E[i].size(); ++j) 
                    --siz[E[i][j]];
            }
        }
        int w = q.back(); isSL[w] = 1;
        for (int i = 0; i < E[w].size(); ++i) {
            int u = E[w][i];
            ++sizS[u];
        }

        for (int i = 0; i < bipartite; ++i) {
            if (isDel[i]) continue;
            if (!isSL[i]) {
                CL.push_back(i);
            } else SL.push_back(i);
        }

        for (int i = bipartite; i < graphSize; ++i) {
            CR.push_back(i);
        }
    }
    bool JudgeBiclique();
    int validSize() {return SL.size() + SR.size() + CL.size() + CR.size();}
    void CalcLimit(int& UB_L, int& UB_R, int& LB_L, int &LB_R, int &minSL, int &minSR);
    void FindPivot(int& wSL, int &wSR, int& wCL, int &wCR);
    void FindForce(int &wL, int &wR);
    MVQBPack Prune7_1(int UB_L, int UB_R, vector<vector<int> >& E);
    MVQBPack Prune7_2(int LB_L, int LB_R, vector<vector<int> >& E);
    MVQBPack Prune7_2_Low(int LB_L, int LB_R, vector<vector<int> >& E);
    MVQBPack Prune_Force(int LB_L, int LB_R, vector<vector<int> >& E);
    void PreBranch1L(vector<int>& vec, int &p, int &q, int w, int UBR, vector<vector<int> >& E);
    void PreBranch1R(vector<int>& vec, int &p, int &q, int w, int UBL, vector<vector<int> >& E);
    MVQBPack Branch1(vector<int>& vec, int w, int p, int q, int pos, bool isL, vector<vector<int> >& E);
    MVQBPack Branch2Add(int w, bool isL, vector<vector<int> >& E);
    MVQBPack Branch2Del(int w, bool isL, vector<vector<int> >& E);
    bool JudgeS(int UBL, int UBR);
    MVQBPack JudgeC(int LBL, int LBR,  vector<vector<int> >& E);
    void Debug();
    void Debug(CorePack T);
};

void MVQBPack::FindForce(int &wL, int &wR) {
    wL = wR = -1;
    for (int i = 0; i < CL.size(); ++i) {
        int u = CL[i];
        if (wL == -1 || siz[wL] > siz[u])
            wL = u;
    }

    for (int i = 0; i < CR.size(); ++i) {
        int u = CR[i];
        if (wR == -1 || siz[wR] > siz[u])
            wR = u;
    }
}

bool MVQBPack::JudgeBiclique() {
    bool valid = 1;
    int limitL = ceil(alpha * totR);
    int limitR = ceil(beta * totL);
    for (int i = 0; i < bipartite; ++i) if (!isDel[i])
        valid &= siz[i] >= limitL;
    for (int i = bipartite; i < graphSize; ++i) if (!isDel[i])
        valid &= siz[i] >= limitR;
    return valid;
}

bool MVQBPack::JudgeS(int UB_L, int UB_R) {
    int limitR = UB_R - (int)(ceil(alpha * UB_R * 100) / 100);
    int limitL = UB_L - (int)(ceil(beta * UB_L * 100) / 100);
    for (int i = 0; i < SL.size(); ++i) 
        if (SR.size() - sizS[SL[i]] > limitR) return 0;
    for (int i = 0; i < SR.size(); ++i)
        if (SL.size() - sizS[SR[i]] > limitL) return 0;
    return 1;
}

void MVQBPack::CalcLimit(int& UB_L, int& UB_R, int& LB_L, int &LB_R, int &minSL, int &minSR) {
    UB_L = UB_R = graphSize;
    if (SL.empty()) UB_R = SR.size() + CR.size();
    else {
        UB_R = SR.size() + CR.size();
        int w = UB_R; 
        for (int i = 0; i < SL.size(); ++i) {
            //UB_R = min(UB_R, (int)floor(E[SL[i]].size() / alpha));
            w = min(w, siz[SL[i]]);
        }
        UB_R = min(UB_R, (int)floor((100 * w) / (100 * alpha)));
    }

    if (SR.empty()) UB_L = SL.size() + CL.size();
    else {
        UB_L = SL.size() + CL.size();
        int w = UB_L;
        for (int i = 0; i < SR.size(); ++i) {
            //UB_L = min(UB_L, (int)floor(E[SR[i]].size() / beta));
            w = min(w, siz[SR[i]]);
        }
        UB_L = min(UB_L, (int)floor((100 * w) / (100 * beta)));
    }

    int nablaSL = 0, nablaSR = 0;

    for (int i = 0; i < SL.size(); ++i) {
        int count = 0;
        // for (int j = 0; j < siz[SL[i]]; ++j) {
        //     if (isSR[E[SL[i]][j]]) ++count;
        // }
        count = sizS[SL[i]];
        nablaSL = max(nablaSL, (int)SR.size() - count);
    }

    for (int i = 0; i < SR.size(); ++i) {
        int count = 0;
        // for (int j = 0; j < E[SR[i]].size(); ++j) {
        //     if (isSL[E[SR[i]][j]]) ++count;
        // }
        count = sizS[SR[i]];
        nablaSR = max(nablaSR, (int)SL.size() - count);
    }

    LB_R = ceil((100 * nablaSL) / (100 - 100 * alpha));
    LB_R = max(max(thetaR, LB_R), (int) SR.size());
    LB_L = ceil((100 * nablaSR) / (100 - 100 * beta));
    LB_L = max(max(thetaL, LB_L), (int) SL.size());

    minSL = SR.size() + CR.size();
    for (int i = 0; i < SL.size(); ++i) {
        minSL = min(minSL, siz[SL[i]]);
    }

    minSR = SL.size() + CL.size();
    for (int i = 0; i < SR.size(); ++i) {
        minSR = min(minSR, siz[SR[i]]);
    }
}

void MVQBPack::FindPivot(int& wSL, int &wSR, int& wCL, int &wCR) {
    wSL = wSR = wCL = wCR = -1;

    //printf("limit : %d\n", (int)ceil(alpha * totR));
    
    for (int i = 0; i < SL.size(); ++i) {
        int u = SL[i];
        if (siz[u] < (int)ceil(alpha * totR)) {
            if (wSL == -1 || siz[wSL] > siz[u])
                wSL = u;
        }
    }

    for (int i = 0; i < SR.size(); ++i) {
        int u = SR[i];
        if (siz[u] < (int)ceil(beta * totL)) {
            if (wSR == -1 || siz[wSR] > siz[u])
                wSR = u;
        }
    }

    for (int i = 0; i < CL.size(); ++i) {
        int u = CL[i];
        if (siz[u] < (int)ceil(alpha * totR)) {
            if (wCL == -1 || siz[wCL] > siz[u])
                wCL = u;
        }
    }

    for (int i = 0; i < CR.size(); ++i) {
        int u = CR[i];
        if (siz[u] < (int)ceil(beta * totL)) {
            if (wCR == -1 || siz[wCR] > siz[u])
                wCR = u;
        }
    }

    // for (int i = 0; i < bipartite; ++i) {
    //     if (!isSL[i]) {
    //         if (E[i].size() < (int)ceil(alpha * (graphSize - bipartite))) {
    //             if (wCL == -1 || E[wCL].size() > E[i].size())
    //                 wCL = i;
    //         }
    //     }
    // }

    // for (int i = bipartite; i < graphSize; ++i) {
    //     if (!isSR[i]) {
    //         if (E[i].size() < (int)ceil(beta * bipartite)) {
    //             if (wCR == -1 || E[wCR].size() > E[i].size())
    //                 wCR = i;
    //         }
    //     }
    // }
}

MVQBPack MVQBPack::Prune7_1(int UB_L, int UB_R, vector<vector<int> >& E) {
    int limitR = UB_R - (int)(ceil(alpha * UB_R * 100) / 100);
    int limitL = UB_L - (int)(ceil(beta * UB_L * 100) / 100);
    MVQBPack nex = *this;
    for (int i = 0; i < SL.size(); ++i) {
        int u = SL[i];
        if (SR.size() - sizS[u] == limitR) {
            vector<bool> isNe(graphSize, 0);

            for (int j = 0; j < E[u].size(); ++j) {
                int v = E[u][j];
                isNe[v] = 1;
            }

            for (int j = bipartite; j < graphSize; ++j) {
                int v = j;
                if (nex.isSR[v] || nex.isDel[v] || isNe[v]) continue;
                nex.isDel[v] = 1;
                for (int k = 0; k < E[v].size(); ++k) {
                    --nex.siz[E[v][k]];
                }
            }
        }
    }

    for (int i = 0; i < SR.size(); ++i) {
        int u = SR[i];
        if (SL.size() - sizS[u] == limitL) {
            vector<bool> isNe(graphSize, 0);
            for (int j = 0; j < E[u].size(); ++j) {
                int v = E[u][j];
                isNe[v] = 1;
            }

            for (int j = 0; j < bipartite; ++j) {
                int v = j;
                if (nex.isSL[v] || nex.isDel[v] || isNe[v]) continue;
                nex.isDel[v] = 1;
                for (int k = 0; k < E[v].size(); ++k)
                    --nex.siz[E[v][k]];
            }
        }
    }

    nex.CL.clear();
    nex.CR.clear();
    for (int i = 0; i < CL.size(); ++i)
        if (!nex.isDel[CL[i]]) nex.CL.push_back(CL[i]);
    for (int i = 0; i < CR.size(); ++i)
        if (!nex.isDel[CR[i]]) nex.CR.push_back(CR[i]);
    nex.totL = nex.CL.size() + nex.SL.size();
    nex.totR = nex.CR.size() + nex.SR.size();
    return nex;
}

MVQBPack MVQBPack::Prune7_2_Low(int LB_L, int LB_R, vector<vector<int> >& E) {
     MVQBPack nex = *this;

    int limitR = 2 * (int)(ceil(alpha * LB_R * 100) / 100) - LB_R - 1;
    for (int i = 0; i < CL.size(); ++i) {
        int u = CL[i]; 
        bool needDel = 0;
        if (!needDel) {
            vector<int> num(graphSize, 0);
            for (int j = 0; j < E[u].size(); ++j) {
                if (!isDel[E[u][j]]) ++num[E[u][j]];
            }
            for (int j = 0; j < SL.size(); ++j) {
                int v = SL[j]; 
                int count = 0;
                for(int k = 0; k < E[v].size(); ++k) {
                    if (isDel[E[v][k]]) continue;
                    count += num[E[v][k]];
                }
                if (count < limitR) {
                    needDel = 1;
                    break;
                }
            }
        }
        if (needDel) {
            nex.isDel[u] = 1;
            for (int j = 0; j < E[u].size(); ++j) {
                --nex.siz[E[u][j]];
            } 
        }
    }
    
    //puts("T1");
    int limitL = 2 * (int)(ceil(beta * LB_L * 100) / 100) - LB_L - 1;
    for (int i = 0; i < CR.size(); ++i) {
        int u = CR[i]; 
        bool needDel = 0;
        if (!needDel) {
            vector<int> num(graphSize, 0);
            for (int j = 0; j < E[u].size(); ++j) {
                if (!isDel[E[u][j]]) ++num[E[u][j]];
            }
            for (int j = 0; j < SR.size(); ++j) {
                int v = SR[j]; 
                int count = 0;
                for(int k = 0; k < E[v].size(); ++k) {
                    if (isDel[E[v][k]]) continue;
                    count += num[E[v][k]];
                }
                if (count < limitL) {
                    needDel = 1;
                    break;
                }
            }
        }
        if (needDel) {
            nex.isDel[u] = 1;
            for (int j = 0; j < E[u].size(); ++j) {
                --nex.siz[E[u][j]];
            } 
        }
    }

    //puts("T2");


    nex.CL.clear();
    nex.CR.clear();
    for (int i = 0; i < CL.size(); ++i)
        if (!nex.isDel[CL[i]]) nex.CL.push_back(CL[i]);
    for (int i = 0; i < CR.size(); ++i)
        if (!nex.isDel[CR[i]]) nex.CR.push_back(CR[i]);
    nex.totL = nex.CL.size() + nex.SL.size();
    nex.totR = nex.CR.size() + nex.SR.size();
    return nex;
}

MVQBPack MVQBPack::Prune_Force(int LBL, int LBR, vector<vector<int> >& E) {
    int limitR = 2 * (int)(ceil(alpha * SR.size() * 100) / 100) - SR.size() - 1;
    int limitL = 2 * (int)(ceil(beta * SL.size() * 100) / 100) - SL.size() - 1;
    MVQBPack nex = *this;
    int wR = ceil(alpha * LBR);
    for (int i = 0; i < CL.size(); ++i) {
        int u = CL[i]; 
        bool needDel = (siz[u] < wR);
        if (needDel) {
            nex.isDel[u] = 1;
            for (int j = 0; j < E[u].size(); ++j) {
                --nex.siz[E[u][j]];
            } 
        }
    }

    int wL = ceil(beta * LBL);
    for (int i = 0; i < CR.size(); ++i) {
        int u = CR[i]; 
        bool needDel = (siz[u] < wL);
        if (needDel) {
            nex.isDel[u] = 1;
            for (int j = 0; j < E[u].size(); ++j) {
                --nex.siz[E[u][j]];
            } 
        }
    }

    //puts("T2");


    nex.CL.clear();
    nex.CR.clear();
    for (int i = 0; i < CL.size(); ++i)
        if (!nex.isDel[CL[i]]) nex.CL.push_back(CL[i]);
    for (int i = 0; i < CR.size(); ++i)
        if (!nex.isDel[CR[i]]) nex.CR.push_back(CR[i]);
    nex.totL = nex.CL.size() + nex.SL.size();
    nex.totR = nex.CR.size() + nex.SR.size();
    return nex;
}

MVQBPack MVQBPack::JudgeC(int LBL, int LBR, vector<vector<int> >& E) {
    int limitR = 2 * (int)(ceil(alpha * SR.size() * 100) / 100) - SR.size() - 1;
    int limitL = 2 * (int)(ceil(beta * SL.size() * 100) / 100) - SL.size() - 1;
    MVQBPack nex = *this;
    int wR = ceil(alpha * LBR);
    for (int i = 0; i < CL.size(); ++i) {
        int u = CL[i]; 
        bool needDel = (siz[u] < wR);
        if (!needDel) {
            vector<int> num(graphSize, 0);
            for (int j = 0; j < E[u].size(); ++j) {
                if (!isDel[E[u][j]]) ++num[E[u][j]];
            }
            for (int j = 0; j < SL.size(); ++j) {
                int v = SL[j]; 
                int count = 0;
                for(int k = 0; k < E[v].size(); ++k) {
                    if (isDel[E[v][k]]) continue;
                    count += num[E[v][k]];
                }
                if (count < limitR) {
                    needDel = 1;
                    break;
                }
            }
        }
        if (needDel) {
            nex.isDel[u] = 1;
            for (int j = 0; j < E[u].size(); ++j) {
                --nex.siz[E[u][j]];
            } 
        }
    }

    int wL = ceil(beta * LBL);
    for (int i = 0; i < CR.size(); ++i) {
        int u = CR[i]; 
        bool needDel = (siz[u] < wL);
        if (!needDel) {
            vector<int> num(graphSize, 0);
            for (int j = 0; j < E[u].size(); ++j) {
                if (!isDel[E[u][j]]) ++num[E[u][j]];
            }
            for (int j = 0; j < SR.size(); ++j) {
                int v = SR[j]; 
                int count = 0;
                for(int k = 0; k < E[v].size(); ++k) {
                    if (isDel[E[v][k]]) continue;
                    count += num[E[v][k]];
                }
                if (count < limitL) {
                    needDel = 1;
                    break;
                }
            }
        }
        if (needDel) {
            nex.isDel[u] = 1;
            for (int j = 0; j < E[u].size(); ++j) {
                --nex.siz[E[u][j]];
            } 
        }
    }

    //puts("T2");


    nex.CL.clear();
    nex.CR.clear();
    for (int i = 0; i < CL.size(); ++i)
        if (!nex.isDel[CL[i]]) nex.CL.push_back(CL[i]);
    for (int i = 0; i < CR.size(); ++i)
        if (!nex.isDel[CR[i]]) nex.CR.push_back(CR[i]);
    nex.totL = nex.CL.size() + nex.SL.size();
    nex.totR = nex.CR.size() + nex.SR.size();
    return nex;
}

MVQBPack MVQBPack::Prune7_2(int LB_L, int LB_R, vector<vector<int> >& E) {
    MVQBPack nex = *this;

    int limitR = 2 * (int)(ceil(alpha * LB_R * 100) / 100) - LB_R - 1;
    int wR = ceil(alpha * LB_R);
    for (int i = 0; i < CL.size(); ++i) {
        int u = CL[i]; 
        bool needDel = (siz[u] < wR);
        if (!needDel) {
            vector<int> num(graphSize, 0);
            for (int j = 0; j < E[u].size(); ++j) {
                if (!isDel[E[u][j]]) ++num[E[u][j]];
            }
            for (int j = 0; j < SL.size(); ++j) {
                int v = SL[j]; 
                int count = 0;
                for(int k = 0; k < E[v].size(); ++k) {
                    if (isDel[E[v][k]]) continue;
                    count += num[E[v][k]];
                }
                if (count < limitR) {
                    needDel = 1;
                    break;
                }
            }
        }
        if (needDel) {
            nex.isDel[u] = 1;
            for (int j = 0; j < E[u].size(); ++j) {
                --nex.siz[E[u][j]];
            } 
        }
    }
    
    //puts("T1");
    int limitL = 2 * (int)(ceil(beta * LB_L * 100) / 100) - LB_L - 1;
    int wL = ceil(beta * LB_L);
    for (int i = 0; i < CR.size(); ++i) {
        int u = CR[i]; 
        bool needDel = (siz[u] < wL);
        if (!needDel) {
            vector<int> num(graphSize, 0);
            for (int j = 0; j < E[u].size(); ++j) {
                if (!isDel[E[u][j]]) ++num[E[u][j]];
            }
            for (int j = 0; j < SR.size(); ++j) {
                int v = SR[j]; 
                int count = 0;
                for(int k = 0; k < E[v].size(); ++k) {
                    if (isDel[E[v][k]]) continue;
                    count += num[E[v][k]];
                }
                if (count < limitL) {
                    needDel = 1;
                    break;
                }
            }
        }
        if (needDel) {
            nex.isDel[u] = 1;
            for (int j = 0; j < E[u].size(); ++j) {
                --nex.siz[E[u][j]];
            } 
        }
    }

    //puts("T2");


    nex.CL.clear();
    nex.CR.clear();
    for (int i = 0; i < CL.size(); ++i)
        if (!nex.isDel[CL[i]]) nex.CL.push_back(CL[i]);
    for (int i = 0; i < CR.size(); ++i)
        if (!nex.isDel[CR[i]]) nex.CR.push_back(CR[i]);
    nex.totL = nex.CL.size() + nex.SL.size();
    nex.totR = nex.CR.size() + nex.SR.size();
    return nex;
}

void MVQBPack::PreBranch1L(vector<int>& vec, int &p, int &q, int w, int UBR, vector<vector<int> >& E) {  
    p = q = 0; 
    vec.clear();
    vector<bool> isNeighbor(graphSize, 0);
    //printf("E[w] : %d\n", E[w].size());
    for (int i = 0; i < E[w].size(); ++i) {
        int u = E[w][i];
        if (isDel[u]) continue;
        //printf("%d\n", u);
        isNeighbor[u] = 1;
        if (!isSR[u]) {
            ++p;
        }
    }

    //printf("p : %d\n", p);

    //p = graphSize - bipartite - SR.size() - p;
    for (int i = bipartite; i < graphSize; ++i) {
        if (isDel[i]) continue;
        if (!isSR[i] && !isNeighbor[i]) {
            vec.push_back(i);
        }
        else if (isSR[i] && !isNeighbor[i]) {
            ++q;
        }
    }

    p = vec.size();
    //printf("%d\n", q);
    q = (int)floor((100 - 100 * alpha) * UBR) / 100 - q;
} 

void MVQBPack::PreBranch1R(vector<int>& vec, int &p, int &q, int w, int UBL, vector<vector<int> >& E) {  
    p = q = 0; 
    vec.clear();
    vector<bool> isNeighbor(graphSize, 0);
    for (int i = 0; i < E[w].size(); ++i) {
        int u = E[w][i];
        if (isDel[u]) continue;
        isNeighbor[u] = 1;
        if (!isSL[u]) {
            ++p;
        }
    }

    //p = graphSize - bipartite - SR.size() - p;
    for (int i = 0; i < bipartite; ++i) {
        if (isDel[i]) continue;
        if (!isSL[i] && !isNeighbor[i]) {
            vec.push_back(i);
        }
        else if (isSL[i] && !isNeighbor[i]) {
            ++q;
        }
    }

    p = vec.size();
    q = (int)floor(((100 - beta * 100) * UBL)) / 100 - q;
} 

MVQBPack MVQBPack::Branch1(vector<int>& vec, int w, int p, int q, int pos, bool isL, vector<vector<int> >& E) {
    MVQBPack nex = *this;
    //printf("pos q : %d %d\n", pos, q);
    //delete
    if (pos != q) {
        //printf("vec[pos] : %d\n", vec[pos]);
        nex.isDel[vec[pos]] = 1;
        int u = vec[pos];
        for (int j = 0; j < E[u].size(); ++j) {
            int v = E[u][j];
            if (nex.isDel[v]) continue;
            --nex.siz[v];
        }
    } else {
        for (int i = pos; i < p; ++i) {
            nex.isDel[vec[i]] = 1;
            //assert(vec[i] < graphSize);
        }

        for (int i = pos; i < p; ++i) {
            int u = vec[i];
            for (int j = 0; j < E[u].size(); ++j) {
                int v = E[u][j];
                if (nex.isDel[v]) continue;
                --nex.siz[v];
            }
        }
    }

    //insert to S
    //nex.sizS[w] += pos;
    if (isL) {
        for (int i = 0; i < pos; ++i) {
            nex.isSR[vec[i]] = 1;
            int u = vec[i];
            for (int j = 0; j < E[u].size(); ++j) {
                int v = E[u][j];
                if (nex.isDel[v]) continue;
                ++nex.sizS[v];
            }
        }
    } else {
        for (int i = 0; i < pos; ++i) {
            nex.isSL[vec[i]] = 1;
            int u = vec[i];
            for (int j = 0; j < E[u].size(); ++j) {
                int v = E[u][j];
                if (nex.isDel[v]) continue;
                ++nex.sizS[v];
            }
        }
    }

    //fix C and S
    nex.SL.clear();
    nex.CL.clear();
    nex.SR.clear();
    nex.CR.clear();
    for (int i = 0; i < bipartite; ++i) {
        if (nex.isDel[i]) continue;
        if (nex.isSL[i]) nex.SL.push_back(i);
        else nex.CL.push_back(i); 
    }

    for (int i = bipartite; i < graphSize; ++i) {
        if (nex.isDel[i]) continue;
        if (nex.isSR[i]) nex.SR.push_back(i);
        else nex.CR.push_back(i);
    }

    nex.totL = nex.CL.size() + nex.SL.size();
    nex.totR = nex.CR.size() + nex.SR.size();
    return nex;
}

MVQBPack MVQBPack::Branch2Add(int w, bool isL, vector<vector<int> >& E) {
    //insert to S
    MVQBPack nex = *this;
    if (isL) {
        nex.isSL[w] = 1;
    } else nex.isSR[w] = 1;

    for (int i = 0; i < E[w].size(); ++i) {
        int u = E[w][i];
        if (nex.isDel[u]) continue;
        ++nex.sizS[u];    
    }

     //fix C and S
    nex.SL.clear();
    nex.CL.clear();
    nex.SR.clear();
    nex.CR.clear();
    for (int i = 0; i < bipartite; ++i) {
        if (nex.isDel[i]) continue;
        if (nex.isSL[i]) nex.SL.push_back(i);
        else nex.CL.push_back(i); 
    }

    for (int i = bipartite; i < graphSize; ++i) {
        if (nex.isDel[i]) continue;
        if (nex.isSR[i]) nex.SR.push_back(i);
        else nex.CR.push_back(i);
    }

    nex.totL = nex.CL.size() + nex.SL.size();
    nex.totR = nex.CR.size() + nex.SR.size();
    return nex;
}

MVQBPack MVQBPack::Branch2Del(int w, bool isL, vector<vector<int> >& E) {
    MVQBPack nex = *this;
    //delete
    nex.isDel[w] = 1;
    for (int i = 0; i < E[w].size(); ++i) {
        int u = E[w][i];
        if (nex.isDel[u]) continue;
        --nex.siz[u];
    }
    //fix C and S
    nex.SL.clear();
    nex.CL.clear();
    nex.SR.clear();
    nex.CR.clear();
    for (int i = 0; i < bipartite; ++i) {
        if (nex.isDel[i]) continue;
        if (nex.isSL[i]) nex.SL.push_back(i);
        else nex.CL.push_back(i); 
    }

    for (int i = bipartite; i < graphSize; ++i) {
        if (nex.isDel[i]) continue;
        if (nex.isSR[i]) nex.SR.push_back(i);
        else nex.CR.push_back(i);
    }
    nex.totL = nex.CL.size() + nex.SL.size();
    nex.totR = nex.CR.size() + nex.SR.size();
    return nex;
}

void MVQBPack::Debug() {
    //freopen("a1.txt", "w", stdout);
    printf("n bip: %d %d\n", graphSize, bipartite);
    printf("a : %.2lf b : %.2lf\n", alpha, beta);
    printf("SL: ");
    for (int i = 0; i < SL.size(); ++i)
        printf("%d ", SL[i]);
    puts("");
    printf("SR: ");
    for (int i = 0; i < SR.size(); ++i)
        printf("%d ", SR[i]);
    puts("");
    printf("CL: ");
    for (int i = 0; i < CL.size(); ++i)
        printf("%d ", CL[i]);
    puts("");
    printf("CR: ");
    for (int i = 0; i < CR.size(); ++i)
        printf("%d ", CR[i]);
    puts("");
    // // for (int i = 0; i < graphSize; ++i)
    // //     printf("siz[%d] : %d ", i, siz[i]);
    // // puts("");
    // // printf("totL : %d totR : %d\n", totL, totR);
    // // for (int i = 0; i < graphSize; ++i)
    // //     if (!isDel[i]) printf("%d ", i);
    // // puts("");
    // printf("n bip: %d %d\n", graphSize, bipartite);
    // cout << totL + totR << endl;
}

void MVQBPack::Debug(CorePack T) {
    printf("n bip: %d %d\n", graphSize, bipartite);
    printf("%d\n", (int)ceil(0.9 * 7));
    printf("%d\n", totL + totR);
    printf("SL: \n");
    for (int i = 0; i < SL.size(); ++i) {
        // printf("%d\n", T.trueID[SL[i]]);
        printf("%d\n", SL[i]);
        int t = 0;
        int u = SL[i];
        for (int j = 0; j < T.E[u].size(); ++j) {
            int v = T.E[u][j];
            if (isDel[v]) continue;
            ++t;
            printf("%d ", T.trueID[v]);
        } puts("");
        // assert(t == siz[u]);
    }
    puts("");
    printf("SR: \n");
    for (int i = 0; i < SR.size(); ++i) {
        printf("%d\n", T.trueID[SR[i]]);
        int u = SR[i];
        int t = 0;
        for (int j = 0; j < T.E[u].size(); ++j) {
            int v = T.E[u][j];
            if (isDel[v]) continue;
            ++t;
            printf("%d ", T.trueID[v]);
        } puts("");
        // assert(t == siz[u]);
    }
    puts("");
    printf("CL: \n");
    for (int i = 0; i < CL.size(); ++i) {
        printf("%d\n", T.trueID[CL[i]]);
        int u = CL[i];
        int t = 0;
        for (int j = 0; j < T.E[u].size(); ++j) {
            int v = T.E[u][j];
            if (isDel[v]) continue;
             ++t;
            printf("%d ", T.trueID[v]);
        } puts("");
        // assert(t == siz[u]);
    }
    puts("");
    printf("CR: \n");
    for (int i = 0; i < CR.size(); ++i) {
        printf("%d\n", T.trueID[CR[i]]);
        int u = CR[i];
        int t = 0;
        for (int j = 0; j < T.E[u].size(); ++j) {
            int v = T.E[u][j];
            if (isDel[v]) continue;
            ++t;
            printf("%d ", T.trueID[v]);
        } puts("");
        // assert(t == siz[u]);
    }
    puts("");
}
/* 

*/