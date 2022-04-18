//
// cppr.cpp
//
// Main driver for max flow solver
//
// (http://lpsolve.sourceforge.net/5.5/DIMACS_maxf.htm)
//

#include <fstream>

#include "dimacs.hpp"
#include "push_relabel.hpp"

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::queue;
using std::streambuf;
using std::string;

void
usage(const string& name)
{
    cerr << "\n"
         << "  usage: " << name << " [-g GRAPH] [-a ALGORITHM] [INP] [-o OUT] [-m] [-c] [-f]\n"
         << "\n"
         << "where GRAPH is one of the followings:\n"
         << "\n"
         << "  list           : Adjacency list graph representation\n"
         << "  star (default) : Star graph representation\n"
         << "\n"
         << "and ALGORITHM is one of the followings:\n"
         << "\n"
         << "  pr             : Sequential FIFO Push Relabel algorithm\n"
         << "  cppr (default) : Colored Parallel Push Relabel algorithm\n"
         << "                   (Number of threads is $OMP_NUM_THREADS)\n"
         << "\n"
         << "and INP/OUT (defaults to stdin/stdout) are DIMACS maximum flow files\n"
         << "and -m flag calculates only the minimum cut (phase 1)\n"
         << "and -c flag checks the solution for feasibility\n"
         << "and -f flag includes flow values on arcs in the output\n"
         << endl;
}

int
main(int argc, char* argv[])
{
    bool flows = false;
    bool mincut = false;
    bool check = false;

    enum class Alg { PR, CPPR } algorithm = Alg::CPPR;
    enum class Graph { LIST, STAR } graph_type = Graph::STAR;

    string inp_name;
    string out_name;

    for (int i = 1; i < argc; ++i) {
        string opt = argv[i];
        if (opt == "-a") {
            ++i;
            if (i == argc) {
                cerr << "error: expected algorithm" << endl;
                usage(argv[0]);
                return 1;
            }
            string val = argv[i];
            if (val == "pr") {
                algorithm = Alg::PR;
            } else if (val == "cppr") {
                algorithm = Alg::CPPR;
            } else {
                cerr << "error: unknown algorithm: " << val << endl;
                usage(argv[0]);
                return 1;
            }
        } else if (opt == "-g") {
            ++i;
            if (i == argc) {
                cerr << "error: expected graph type" << endl;
                usage(argv[0]);
                return 1;
            }
            string val = argv[i];
            if (val == "list") {
                graph_type = Graph::LIST;
            } else if (val == "star") {
                graph_type = Graph::STAR;
            } else {
                cerr << "error: unknown graph type: " << val << endl;
                usage(argv[0]);
                return 1;
            }
        } else if (opt == "-o") {
            ++i;
            if (i == argc) {
                cerr << "error: expected output file name" << endl;
                usage(argv[0]);
                return 1;
            }
            out_name = argv[i];
        } else if (opt == "-m") {
            mincut = true;
        } else if (opt == "-c") {
            check = true;
        } else if (opt == "-f") {
            flows = true;
        } else if (opt.front() == '-') {
            cerr << "error: unknown option: " << opt << endl;
            usage(argv[0]);
            return 1;
        } else {
            if (!inp_name.empty()) {
                cerr << "error: unexpected argument " << argv[i] << endl;
                usage(argv[0]);
                return 1;
            }
            inp_name = argv[i];
        }
    }

    ifstream inp;
    streambuf *cinbuf = cin.rdbuf();
    if (!inp_name.empty()) {
        inp.open(inp_name);
        if (!inp) {
            cerr << "error: can not open file: " << inp_name << endl;
            return -1;
        }

        cin.rdbuf(inp.rdbuf());
    }

    ofstream out;
    streambuf *coutbuf = cout.rdbuf();
    if (!out_name.empty()) {
        out.open(out_name);
        if (!out) {
            cerr << "error: can not open file: " << out_name << endl;
            return -1;
        }

        cout.rdbuf(out.rdbuf());
    }

    ulval_t flow = 0;

    cout << "c\nc cppr v1.0.0" << endl;

    switch (graph_type) {
    case Graph::LIST: {
        auto graph = ListGraph(read_dimacs(cin));

        cout << "c\n"
             << "c # nodes : " << graph.num_nodes() << "\n"
             << "c # arcs  : " << graph.num_arcs() << endl;

        cout << "c\nc graph     : list" << endl;

        if (mincut) {
            cout << "c problem   : minimum cut\n";
        } else {
            cout << "c problem   : maximum flow\n";
        }

        switch (algorithm) {
        case Alg::PR:
            cout << "c algorithm : push-relabel (sequential)" << endl;

            flow = push_relabel<ListGraph, queue<node_t>, queue<node_t>>(graph, mincut);

            STAT_OUT(cout, pr_stat);

            if (check && !mincut) {
                check_solution_max(cout, graph);
            } else {
                cout << "c\nc solution is not checked\n";
            }

            cout << "c\ns " << flow << endl;

            if (flows) {
                print_dimacs(graph);
            }

            TIME_OUT(cout, pr_time);

            break;
        case Alg::CPPR:
            cout << "c algorithm : colored parallel push-relabel ("
                 << omp_get_max_threads()
                 << " threads)"
                 << endl;

            flow = push_relabel<ListGraph, ColorQueue, MultiQueue>(graph, mincut);

            STAT_OUT(cout, cppr_stat);

            if (check && !mincut) {
                check_solution_max(cout, graph);
            } else {
                cout << "c\nc solution is not checked\n";
            }

            cout << "c\ns " << flow << endl;

            if (flows) {
                print_dimacs(graph);
            }

            TIME_OUT(cout, cppr_time);

            break;
        }

        break;
    }
    case Graph::STAR: {
        auto graph = StarGraph(read_dimacs(cin));

        cout << "c\n"
             << "c # nodes : " << graph.num_nodes() << "\n"
             << "c # arcs  : " << graph.num_arcs() << endl;

        cout << "c\nc graph     : star" << endl;

        if (mincut) {
            cout << "c problem   : minimum cut\n";
        } else {
            cout << "c problem   : maximum flow\n";
        }

        switch (algorithm) {
        case Alg::PR:
            cout << "c algorithm : push-relabel (sequential)" << endl;

            flow = push_relabel<StarGraph, queue<node_t>, queue<node_t>>(graph, mincut);

            STAT_OUT(cout, pr_stat);

            if (check && !mincut) {
                check_solution_max(cout, graph);
            } else {
                cout << "c\nc solution is not checked\n";
            }

            cout << "c\ns " << flow << endl;

            if (flows) {
                print_dimacs(graph);
            }

            TIME_OUT(cout, pr_time);

            break;
        case Alg::CPPR:
            cout << "c algorithm : colored parallel push-relabel ("
                 << omp_get_max_threads() << " threads)" << endl;

            flow = push_relabel<StarGraph, ColorQueue, MultiQueue>(graph, mincut);

            STAT_OUT(cout, cppr_stat);

            if (check && !mincut) {
                check_solution_max(cout, graph);
            } else {
                cout << "c\nc solution is not checked\n";
            }

            cout << "c\ns " << flow << endl;

            if (flows) {
                print_dimacs(graph);
            }

            TIME_OUT(cout, cppr_time);

            break;
        }

        break;
    }
    }

    cin.rdbuf(cinbuf);
    cout.rdbuf(coutbuf);

    return 0;
}
