#include <vector>
#include <algorithm>
#include <cstdio>
#include <queue>
#include "ColorProcess.h"
#include "CoreProcess.h"
#include <cmath>
#define ull unsigned long long 
using namespace std;


class PrePack {
    int thetaL;
    int thetaR;
    public:
        PrePack(){}
        PrePack(int thetaL, int thetaR) : thetaL(thetaL), thetaR(thetaR){}
        CorePack PreProcess(CorePack cur);
};

CorePack PrePack::PreProcess(CorePack cur) {
    while(true) {
        CorePack nex = cur.FindBiCore(cur.lowL, cur.lowR);
        // puts("?");
        if (nex.graphSize == cur.graphSize) break;
        cur = nex;
        
        //puts("test");
        //for (int i = 0; i < cur.graphSize; ++i) printf("%d ", cur.trueID[i]);

        //color

        int bipartite = cur.bipartite;
        int graphSize = cur.graphSize;

        vector<vector<int> > colorEL, colorER;
        int limit = 2 * cur.lowL - thetaR - 1; 

        colorEL.resize(bipartite);
        // printf("%d %d\n", limit, bipartite);
        for (int i = 0; i < bipartite; ++i) {
            vector<int> deg(bipartite, 0);
            //printf("? %d\n", cur.E[i].size());
            for (int j = 0; j < cur.E[i].size(); ++j) {
                int u = cur.E[i][j];
                //printf("?? %d\n", u);
                for (int z = 0; z < cur.E[u].size(); ++z) {
                    // if (cur.E[u][z] >= bipartite) {
                    //     printf("%d %d\n", cur.trueID[i], cur.trueID[u]);
                    //     puts("err");
                    // }
                    ++deg[cur.E[u][z]];
                }
            }
            for (int j = i + 1; j < bipartite; ++j) {
                if (deg[j] >= limit) {
                    //printf("%d %d %d\n", i, j, cur.E[i].size());
                    colorEL[i].push_back(j);
                    colorEL[j].push_back(i);
                }
            }
        }

        // puts("??");
        
        vector<bool> isValidL, isValidR, isValid(graphSize);
        ColorPack colorpackL(colorEL, bipartite, thetaL - 1);
        colorpackL.DrawColor(isValidL);
        // puts("???");
        limit = 2 * cur.lowR - thetaL - 1;
        //ColorE.clear();

        colorER.resize(graphSize - bipartite);

        for (int i = bipartite; i < graphSize; ++i) {
            vector<int> deg(graphSize, 0);
            //printf("? %d %d\n", i, cur.E[i].size());
            for (int j = 0; j < cur.E[i].size(); ++j) {
                int u = cur.E[i][j];
                //printf("?? %d\n", u);
                for (int z = 0; z < cur.E[u].size(); ++z) 
                    ++deg[cur.E[u][z]];
            }
            for (int j = i + 1; j < graphSize; ++j) {
                if (deg[j] >= limit) {
                    //printf("%d %d\n", i, j);
                    colorER[i - bipartite].push_back(j - bipartite);
                    colorER[j - bipartite].push_back(i - bipartite);
                }
            }
        }
        //puts("???");
        //printf("%d\n", graphSize - bipartite);
        ColorPack colorpackR(colorER, graphSize - bipartite, thetaR - 1);
        colorpackR.DrawColor(isValidR);
        //puts("???");
        for (int i = 0; i < bipartite; ++i) isValid[i] = isValidL[i];
        for (int i = bipartite; i < graphSize; ++i) isValid[i] = isValidR[i - bipartite];

        int validCount = 0;
        nex.trueID.clear();
        vector<int> pos(graphSize);
        for (int i = 0; i < bipartite; ++i) {
            if (isValid[i]) {
                nex.trueID.push_back(cur.trueID[i]);
                pos[i] = validCount++;
            }
        }

        nex.bipartite = validCount;

        for (int i = bipartite; i < graphSize; ++i) {
            if (isValid[i]) {
                nex.trueID.push_back(cur.trueID[i]);
                pos[i] = validCount++;
            }
        }

        nex.graphSize = validCount;
        if (nex.graphSize == cur.graphSize) break;
        nex.E.resize(validCount);
        nex.lowL = cur.lowL;
        nex.lowR = cur.lowR;   

        for (int i = 0; i < graphSize; ++i) {
            if (isValid[i]) {
                nex.E[pos[i]].clear();
                for (int j = 0; j < cur.E[i].size(); ++j) {
                    if (isValid[cur.E[i][j]]) {
                        nex.E[pos[i]].push_back(pos[cur.E[i][j]]);
                    }
                }
            }
        }

        cur = nex;
    }
    //puts("test");
    //for (int i = 0; i < cur.graphSize; ++i) printf("%d ", cur.trueID[i]);
    return cur;
}

