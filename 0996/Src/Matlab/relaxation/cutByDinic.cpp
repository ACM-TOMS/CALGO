#include <vector>
#include <limits>
#include <queue>
#include <string.h>
// #include <iostream>
#include "mex.h"
#include "matrix.h"

const int INF = std::numeric_limits<int>::max();

struct edge { int to; double cap; int rev; }; // rev-id

typedef std::vector<edge> edges;
typedef std::vector<edges> graph;

void add_edge(graph &G, int from, int to, double cap) {
    G[from].push_back((edge){to, cap, (int)G[to].size()});
    G[to].push_back((edge){from, 0, (int)G[from].size() - 1 });
}

void compute_min_dist_by_bfs(graph &G, int *level, int s) {
    memset(level, -1, sizeof(int) * G.size());
    std::queue<int> que;
    level[s] = 0;
    que.push(s);
    while (!que.empty()) {
        int v = que.front(); que.pop();
        for (int i = 0; i < G[v].size(); i++) {
            edge &e = G[v][i];
            if (e.cap > 0 && level[e.to] < 0) {
                level[e.to] = level[v] + 1;
                que.push(e.to);
            }
        }
    }
}

double find_augment_path_by_dfs(graph &G, int *level, int *iter, int v, int t, double f) {
    if (v == t) return f;
    for (int &i = iter[v]; i < G[v].size(); i++){
        edge &e = G[v][i];
        if (e.cap > 0 && level[v] < level[e.to]) {
            double d = find_augment_path_by_dfs(G, level, iter, e.to, t, std::min(f, e.cap));
            if (d > 0) {
                e.cap -= d;
                G[e.to][e.rev].cap += d;
                return d;
            }
        }
    }
    return 0;
}

double max_flow(graph &G, int s, int t) {
    double flow = 0;
    int n = G.size();
    int *level = new int[n];
    int *iter = new int[n];
    for (;;) {
        compute_min_dist_by_bfs(G, level, s);
        if (level[t] < 0) {
        	delete[] level;
            delete[] iter;
        	return flow;
        }
        memset(iter, 0, sizeof(int) * n);
        double f;
        while ((f = find_augment_path_by_dfs(G, level, iter, s, t, INF)) > 0) {
            flow += f;
        }
    }
}

void check_reachability_from(graph &G, bool *reachable, int s) {
    memset(reachable, false, sizeof(bool) * G.size());
    std::queue<int> que;
    reachable[s] = true;
    que.push(s);
    while (!que.empty()) {
        int v = que.front(); que.pop();
        for (int i = 0; i < G[v].size(); i++) {
            edge &e = G[v][i];
            if (e.cap > 0 && !reachable[e.to]) {
                reachable[e.to] = true;
                que.push(e.to);
            }
        }
    }
}

void check_reachability_to(graph &G, bool *reachable, int t) {
    memset(reachable, false, sizeof(bool) * G.size());
    std::queue<int> que;
    reachable[t] = true;
    que.push(t);
    while (!que.empty()) {
        int v = que.front(); que.pop();
        for (int i = 0; i < G[v].size(); i++) {
            edge &e = G[v][i];
            if (G[e.to][e.rev].cap > 0 && !reachable[e.to]) {
                reachable[e.to] = true;
                que.push(e.to);
            }
        }
    }
}

double min_cut(graph &G, bool *reachable, int s, int t) {
    double f; 
    f = max_flow(G, s, t);
    check_reachability_to(G, reachable, t);
    return f;
}

// int main(){
// 	int MAX_V = 5;
// 	graph G(MAX_V);
// 	int s = 0;
// 	int t = 4;
// 	add_edge(G, 0, 1, 10);
// 	add_edge(G, 0, 2, 2);
// 	add_edge(G, 1, 2, 6);
// 	add_edge(G, 1, 3, 6);	
// 	add_edge(G, 3, 2, 3);	
// 	add_edge(G, 2, 4, 5);	
// 	add_edge(G, 3, 4, 8);	
//
// 	bool *reachable = new bool[MAX_V];
//
// 	std::cout << min_cut(G, reachable, s, t) << std::endl;
// 	for (int i = 0; i < MAX_V; i++) {
// 		std::cout << reachable[i];
// 	}
// 	std::cout << std::endl;
// 	delete[] reachable;
//     return 0;
// }

void mexFunction(int nlhs,   mxArray  *plhs[], 
                 int nrhs,   const mxArray  *prhs[] )
{
/* CHECK FOR PROPER NUMBER OF ARGUMENTS */

    if (nrhs != 3) { mexErrMsgTxt("[val, Z] = cutByDinic(CapacityMat, s, t) requires 3 input arguments."); }
    if (nlhs != 2) { mexErrMsgTxt("[val, Z] = cutByDinic(CapacityMat, s, t) requires 2 output arguments."); }

/* Check the proper input type */

    if (!mxIsSparse(prhs[0])) { mexErrMsgTxt("CoapacityMat must be sparse"); }
    // if (!mxIsScalar(prhs[1]) || !mxIsScalar(prhs[2])) { mexErrMsgTxt("s and t must be scalar."); }
    if (!(mxGetM(prhs[1]) == 1) || !(mxGetM(prhs[2]) == 1)
        || !(mxGetN(prhs[1]) == 1) || !(mxGetN(prhs[2]) == 1)) 
        { mexErrMsgTxt("s and t must be scalar."); }

/* Read inputs */
    mwSize m = mxGetM(prhs[0]); 
    mwSize n = mxGetN(prhs[0]); 
    if (m != n) {mexErrMsgTxt("CapacityMat must be square.");}
    //double *CapacityMat = mxGetPr(prhs[0]);

    double *pr = mxGetPr(prhs[0]);
    mwIndex *jc = mxGetJc(prhs[0]);
    mwIndex *ir = mxGetIr(prhs[0]);
    graph G(n);
    for (mwIndex j = 0; j < n; j++) {
        for (mwIndex k = jc[j]; k < jc[j+1]; k++) {
            mwIndex i = ir[k];
            double capacity = pr[k];
            //g[i].push_back(Edge(i , j, capacity));
            add_edge(G, i, j, capacity);
        }
    }

    mwIndex s = (mwIndex)mxGetScalar(prhs[1])-1;
    mwIndex t = (mwIndex)mxGetScalar(prhs[2])-1;

/* Allocate space for output vector */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateLogicalMatrix(1, n);
    double *val = mxGetPr(plhs[0]);
    bool *Z = mxGetLogicals(plhs[1]);

/* solve */
    *val = min_cut(G, Z, s, t);
}


