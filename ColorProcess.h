#include <vector>
#include <algorithm>
#include <cstdio>
#include <queue>

using namespace std;

class ColorPack {
    vector<vector<int> > E;
    int graphSize;
    int limit;
    public:
        ColorPack(){}
        ColorPack(vector<vector<int> > E, int graphSize, int limit);
        void DrawColor(vector<bool>& isValid);
};

ColorPack::ColorPack(vector<vector<int> > E, int graphSize, int limit) {
    this->E = E;
    this->graphSize = graphSize;
    this->limit = limit;
}

void ColorPack::DrawColor(vector<bool>& isValid) {
    vector<int> color(graphSize, -1), neColor(graphSize, -1); // neighborColor
    isValid.resize(graphSize);
    //for (int i = 0; i < graphSize; ++i) isValid[i] = 1;
    int clr = 0;
    vector<int> unDraw;
    for (int i = 0; i < graphSize; ++i) unDraw.push_back(i);
    while (!unDraw.empty()) {
        vector<int> nexUnDraw; nexUnDraw.clear();
        for (int i = 0; i < unDraw.size(); ++i) {
            int u = unDraw[i];
            //printf("? : %d\n ", u);
            if (neColor[u] == clr) {
                nexUnDraw.push_back(u);
                continue;
            }       
            color[u] = clr;
            for (int j = 0; j < E[u].size(); ++j) {
                int v = E[u][j];
                if (color[v] == -1) {
                    neColor[v] = clr;
                }
            }
        } ++clr;
        //printf("clr : %d ", clr);
        unDraw = nexUnDraw;
    }

    //puts("end");

    for (int i = 0; i < graphSize; ++i) {
        vector<int> tempColor(clr, 0);
        //printf("i:%d %d\n", i, E[i].size());
        for (int j = 0; j < E[i].size(); ++j) {
            tempColor[color[E[i][j]]]++;
            //printf("%d\n", j);
        }
        //puts("End");
        int colorCount = 0;
        for (int j = 0; j < clr; ++j) 
            colorCount += tempColor[j] > 0;
        //puts("End");
        isValid[i] = colorCount >= limit;
    }

    //puts("end");
}