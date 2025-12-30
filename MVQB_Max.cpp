#include <iostream>
#include <string>
#include <time.h>
#include <list>
#include <vector>
#include <cstring>
#include <cmath>
#include <iostream>
#include <thread>
#include <chrono>
#include <atomic>
#include <cstdlib>
#include <mutex>
#include <omp.h> // 引入 OpenMP 库

#include "args.hxx"
#include "Util.h"
#include "preProcess.h"
#include "MVQBP.h"

#define FILELEN 1024

std::atomic<int> curBest; 
MVQBPack debug;
std::mutex debugMutex;

const int PARALLEL_THRESHOLD = 16; 

void MVQBP_SUB(MVQBPack cur, vector<vector<int> >& E, int depth) {
    if (curBest.load(std::memory_order_relaxed) >= cur.validSize()) return; //Proposition 6
    if (cur.SL.size() + cur.CL.size() < cur.thetaL || cur.SR.size() + cur.CR.size() < cur.thetaR) return;
    
    if (cur.JudgeBiclique()) {
        int currentVal = cur.validSize();
        if (currentVal > curBest.load(std::memory_order_relaxed)) {
            std::lock_guard<std::mutex> lock(debugMutex);
            if (currentVal > curBest.load(std::memory_order_relaxed)) {
                curBest.store(currentVal);
                debug = cur;
            }
        }
        return;
    }

    int UBL, UBR, LBL, LBR, minSL, minSR;
    cur.CalcLimit(UBL, UBR, LBL, LBR, minSL, minSR);
    if (UBL + UBR <= curBest.load(std::memory_order_relaxed) || minSL < (int)ceil(cur.alpha * LBR) || minSR < (int)ceil(cur.beta * LBL)) return;
    if (UBL < LBL || UBR < LBR) return;


    cur = cur.Prune7_2(LBL, LBR, E);
    if (curBest.load(std::memory_order_relaxed) >= cur.validSize()) return;
    if (cur.SL.size() + cur.CL.size() < cur.thetaL || cur.SR.size() + cur.CR.size() < cur.thetaR) return;
    if (cur.JudgeBiclique()) {
        int currentVal = cur.validSize();
        if (currentVal > curBest.load(std::memory_order_relaxed)) {
            std::lock_guard<std::mutex> lock(debugMutex);
            if (currentVal > curBest.load(std::memory_order_relaxed)) {
                curBest.store(currentVal);
                debug = cur;
            }
        }
        return;
    }

    cur.CalcLimit(UBL, UBR, LBL, LBR, minSL, minSR);
    if (UBL + UBR <= curBest.load(std::memory_order_relaxed) || minSL < (int)ceil(cur.alpha * LBR) || minSR < (int)ceil(cur.beta * LBL)) return;
    if (UBL < LBL || UBR < LBR) return; //Proposition 5

    int wSR, wSL, wCR, wCL;
    cur.FindPivot(wSL, wSR, wCL, wCR);

    bool spawn_tasks = (depth < PARALLEL_THRESHOLD);

    if (wSL == -1 && wSR == -1) {
        if (wCR == -1 || (wCL != -1 && cur.siz[wCL] < cur.siz[wCR])) {
            if (spawn_tasks) {
                #pragma omp task shared(E) firstprivate(cur)
                {
                    MVQBPack nex = cur.Branch2Add(wCL, 1, E);
                    MVQBP_SUB(std::move(nex), E, depth + 1);
                }
                #pragma omp task shared(E) firstprivate(cur)
                {
                    MVQBPack nex = cur.Branch2Del(wCL, 1, E);
                    MVQBP_SUB(std::move(nex), E, depth + 1);
                }
                #pragma omp taskwait
            } else {
                MVQBPack nex = cur.Branch2Add(wCL, 1, E);
                MVQBP_SUB(std::move(nex), E, depth + 1);
                nex = cur.Branch2Del(wCL, 1, E);
                MVQBP_SUB(std::move(nex), E, depth + 1);
            }
        } else {
            if (spawn_tasks) {
                #pragma omp task shared(E) firstprivate(cur)
                {
                    MVQBPack nex = cur.Branch2Add(wCR, 0, E);
                    MVQBP_SUB(std::move(nex), E, depth + 1);
                }
                #pragma omp task shared(E) firstprivate(cur)
                {
                    MVQBPack nex = cur.Branch2Del(wCR, 0, E);
                    MVQBP_SUB(std::move(nex), E, depth + 1);
                }
                #pragma omp taskwait
            } else {
                MVQBPack nex = cur.Branch2Add(wCR, 0, E);
                MVQBP_SUB(std::move(nex), E, depth + 1);
                nex = cur.Branch2Del(wCR, 0, E);
                MVQBP_SUB(std::move(nex), E, depth + 1);
            }
        }
    } else {
        if (wSR == -1 || (wSL != -1 && cur.siz[wSL] < cur.siz[wSR])) {
            vector<int> vec; int p, q;
            cur.PreBranch1L(vec, p, q, wSL, UBR, E);
            
            for (int i = 0; i <= q; ++i) {
                if (spawn_tasks) {
                    #pragma omp task shared(E) firstprivate(i, vec, wSL, p, q, cur)
                    {
                        MVQBPack nex = cur.Branch1(vec, wSL, p, q, i, 1, E);
                        MVQBP_SUB(std::move(nex), E, depth + 1);
                    }
                } else {
                    MVQBPack nex = cur.Branch1(vec, wSL, p, q, i, 1, E);
                    MVQBP_SUB(std::move(nex), E, depth + 1);
                }
            }
            if (spawn_tasks) {
                #pragma omp taskwait
            }
        } else {
            vector<int> vec; int p, q;
            cur.PreBranch1R(vec, p, q, wSR, UBL, E);
            
            for (int i = 0; i <= q; ++i) {
                if (spawn_tasks) {
                    #pragma omp task shared(E) firstprivate(i, vec, wSR, p, q, cur)
                    {
                        MVQBPack nex = cur.Branch1(vec, wSR, p, q, i, 0, E);
                        MVQBP_SUB(std::move(nex), E, depth + 1);
                    }
                } else {
                    MVQBPack nex = cur.Branch1(vec, wSR, p, q, i, 0, E);
                    MVQBP_SUB(std::move(nex), E, depth + 1);
                }
            }
            if (spawn_tasks) {
                #pragma omp taskwait
            }
        }
    }
}

const int MAX_RUNTIME_SECONDS = 30000;

void timerFunction()
{
    std::this_thread::sleep_for(std::chrono::seconds(MAX_RUNTIME_SECONDS));
    std::cout << "Program exiting due to timeout." << std::endl;
    exit(0);
}


int main(int argc, char** argv) {
    char filepath[1024] = ".........";
    double alpha = 0.0, beta = 0.0;

    args::ArgumentParser parser(
        "An algorithm for computing the maximum quasi-biclique\n");

    args::HelpFlag help(parser, "help", "Display this help menu",
                        {'h', "help"});
    args::Group required(parser, "", args::Group::Validators::All);

    args::ValueFlag<std::string> benchmark_file(
        parser, "benchmark", "Path to benchmark", {'f', "file"},
        "");

    args::ValueFlag<int> Threshold_l(parser, "Lower Bound", "The lower bound of the size of quasi-biclique", {'u', "u"}, 1);

    args::ValueFlag<int> Threshold_r(parser, "Lower Bound", "The lower bound of the size of quasi-biclique", {'v', "v"}, 1);

    args::ValueFlag<double> Alpha(parser, "para alpha", "The parameter alpha", {'a', "alpha"}, 1);
    
    args::ValueFlag<double> Beta(parser, "para beta", "The parameter beta", {'b', "beta"}, 1);

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 0;
    } catch (args::ValidationError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 0;
    }

    strncpy(filepath, args::get(benchmark_file).c_str(), FILELEN);
    ofstream file("./mvqb.txt", ios::app);
    streambuf *coutbuf = cout.rdbuf(); // 保存 cout 的缓冲区指针
    cout.rdbuf(file.rdbuf()); // 将 cout 的缓冲区指针指向文件流的缓冲区

    alpha = args::get(Alpha);
    beta = args::get(Beta);
    int theta_l = args::get(Threshold_l);
    int theta_r = args::get(Threshold_r);
    // cout<<"---------------"<<filepath<<" a,b,lb_L,lb_R,r: "<<alpha<<" "<<beta<<" "<<theta_l<<" "<<theta_r<<" "<<isTwoHopReduce<<"--------------"<<endl;
    if(alpha<=0.5 || beta<=0.5){
        cout<<"alpha and beta should be larger than 0.5"<<endl;
        exit(-1);
    }

    Util util;
    int bi=0;
    int *degree=NULL;
    int **Graph=NULL;
    int graph_size=util.ReadGraph(filepath,Graph,degree,bi);
    
    std::thread timerThread(timerFunction);
    timerThread.detach(); 

    double s1 = omp_get_wtime();
    CorePack corepack(Graph, degree, graph_size,bi,ceil(theta_r*alpha),ceil(theta_l*beta));
    PrePack prepack(theta_l, theta_r);
    CorePack cur = prepack.PreProcess(corepack);
    
    time_t s3 = clock();
    curBest.store(0);
    MVQBPack st(cur.E, alpha, beta, cur.bipartite, cur.graphSize, theta_l, theta_r);
    
    #pragma omp parallel
    {
        #pragma omp single
        {
            MVQBP_SUB(st, cur.E, 0);
            std::cout << ">>> [DEBUG] OpenMP Threads detected: " << omp_get_num_threads() << " <<<" << std::endl;
        }
    }
    
    double s2 = omp_get_wtime();
    cout<<"---------------"<<filepath<<" a,b,lb_L,lb_R: "<<alpha<<" "<<beta<<" "<<theta_l<<" "<<theta_r<<" "<<"--------------"<<endl 
    << debug.validSize() << endl
    <<"Running Time: "<<s2-s1<<" sec"<<endl;

    cout.rdbuf(coutbuf);

    // 关闭文件流
    file.close();
    return 0;
}