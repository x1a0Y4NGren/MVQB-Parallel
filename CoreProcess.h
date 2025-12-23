#include <vector>
#include <algorithm>
#include <cstdio>
#include <queue>

using namespace std;

class CorePack {
    public:
        vector<vector<int> > E;
        int graphSize;
        int bipartite;
        int lowL;
        int lowR;
        vector<int> trueID;

        CorePack(){}
        CorePack(int **&new_graph, int *&degree, int graphSize, int bipartite, int thetaL, int thetaR); //init
        CorePack FindBiCore(int thetaL, int thetaR);  
};

CorePack::CorePack(int **&new_graph, int *&degree, int graphSize, int bipartite, int lowL, int lowR) {
    this->graphSize = graphSize;
    this->bipartite = bipartite;
    this->lowL = lowL;
    this->lowR = lowR;

    E.resize(graphSize);
    trueID.resize(graphSize);
    for (int i = 0; i < graphSize; ++i) {
        trueID[i] = i;
        E[i].resize(degree[i]);
        for (int j = 0; j < degree[i]; ++j) {
            E[i][j] = new_graph[i][j];
        }
        sort(E[i].begin(), E[i].end());
    }
}

CorePack CorePack::FindBiCore(int lowL, int lowR) {//find bicore
    vector<bool> isValid(graphSize, 1);
    vector<int> pos(graphSize, -1);
    vector<int> deg(graphSize);
    queue<int> q;
    int validCount = 0;
    CorePack nex;

    for (int i = 0; i < bipartite; ++i) deg[i] = E[i].size() - lowL;
    for (int i = bipartite; i < graphSize; ++i) deg[i] = E[i].size() - lowR;
    for (int i = 0; i < graphSize; ++i) {
        if (deg[i] < 0) {
            isValid[i] = 0;
            q.push(i);
        }
    }

    while (!q.empty()) {
        int u = q.front(); q.pop();
        //printf("u: %d\n", u);
        //printf("%d\n", cur.E[u].size());
        for (int i = 0; i < E[u].size(); ++i) {
            int v = E[u][i];
            if (!isValid[v]) continue;
            --deg[v];
            //printf("%d %d\n", v, deg[v]);
            if (deg[v] < 0) {
                isValid[v] = 0;
                q.push(v);
            }
            //puts("???");
        }
    }

    //puts("???");
    nex.trueID.clear();

    for (int i = 0; i < bipartite; ++i) {
        if (isValid[i]) {
            //nex.trueID[validCount] = cur.trueID[i];
            nex.trueID.push_back(trueID[i]);
            pos[i] = validCount++;
        }
    }

    nex.bipartite = validCount;

    for (int i = bipartite; i < graphSize; ++i) {
        if (isValid[i]) {
            nex.trueID.push_back(trueID[i]);
            //nex.trueID[validCount] = cur.trueID[i];
            pos[i] = validCount++;
        }
    }
    //puts("???");

    nex.graphSize = validCount;
    nex.E.resize(validCount);
    nex.lowL = lowL;
    nex.lowR = lowR;

    for (int i = 0; i < graphSize; ++i) {
        if (isValid[i]) {
            nex.E[pos[i]].clear();
            for (int j = 0; j < E[i].size(); ++j) {
                if (isValid[E[i][j]]) {
                    nex.E[pos[i]].push_back(pos[E[i][j]]);
                }
            }
        }
    }

    return nex;
}