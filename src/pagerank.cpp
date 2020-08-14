/*
 * CS 106X Final Project. Fall 2019.
 * Rishi Desai
 *
 * The "proteins" folder should go in the "res" folder
 */
#include "graphtypes.h"
#include <iostream>
#include <string>
#include <fstream>
#include "tokenscanner.h"
#include "hashset.h"
#include "simpio.h"
#include "strlib.h"
#include "console.h"
#include "filelib.h"

using namespace std;

/*
 * Constants for PageRank calculation
 */
const static int MAX_ITER = 50;             // max #iterations - set low at 50 for student demo
const static double ALPHA = 1;              // dampener
const static double EPSILON = .000001;      // convergence constant

/*
 * Matrix-Vector multiplication with scaling factor of ALPHA
 */
static Vector<double> matrixVectorMult(Vector<Vector<double>> &matrix, Vector<double> &vector) {
    if(matrix.size() != vector.size()) { error("Invalid matrix and vector size!"); }
    for (int i = 0; i < matrix.size(); i++) {
        if (matrix.size() != matrix[i].size()) cout << "FAILED";
    }

    Vector<double> ret(vector.size(), 0.0);
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            ret[i] += matrix[i][j] * vector[j];
        }
        ret[i] *= ALPHA;
    }
    return ret;
}

/*
 * Vector addition
 */
static Vector<double> vectorAdd(Vector<double> &vector1, Vector<double> &vector2) {
    if (vector1.size() != vector2.size()) { error("Cannot Add Vectors!"); }

    Vector<double> ret(vector1.size(), 0.0);
    for (int i = 0; i < ret.size(); i++) { ret[i] = vector1[i] + vector2[2]; }
    return ret;
}

/*
 * Tests for matrix convergence by randomly sampling ~10 nodes and averaging the
 * difference - if avg is less than EPSILON, matrix has converged
 */
static bool isConverged(Vector<double> &curr, Vector<double> &next) {
    if (curr.size() != next.size()) { error("Invalid vectors for convergence!"); }

//    double error = 0 ;
//    for (int i = 0; i < curr.size(); i++) {
//        error += abs(curr[i] - next[i]);
//    }
//    error /=  curr.size();
//    cout << error << endl;

    HashSet<int> randIndices;
    for (int i = 0; i < 10; i++) { randIndices.add(randomInteger(0, curr.size() - 1)); }

    double avg_diff = 0;
    for (int i : randIndices) {
        avg_diff += abs(curr[i] - next[i]);
    }
    avg_diff /= randIndices.size();
    return avg_diff <= EPSILON;
}

/*
 * PageRank iterative algorithm:
 * R(t+1) = alpha*M*R(t) + (1-alpha)/N * 1V
 * where 1V = 1-vector, M = PR matrix, R(t) = PR values, alpha = dampener value
 */
static void computePageRank(Graph *graph) {
    const int size = graph->nodes.size();

    Vector<Vector<double>> matrix(size);            // set the Markov matrice
    for (int i = 0; i < matrix.size(); i++) {
        Node *ni = graph->IDtoNode.get(i);

        matrix[i]= Vector<double>(size, 0.0);
        for (int j = 0; j < matrix[i].size(); j++) {    // M[i][j] = 1/L(n_i)
            Node *nj = graph->IDtoNode.get(j);          // M is based on #links pointing to neighbors

            matrix[i][j] = nj->outEdgeNodeNames.contains(ni->name) ?
                        1.0 / ((double) nj->outEdges.size()) : 0.0;
        }
    }

    Vector<double> dampener(size, (1.0 - ALPHA) /((double) size));    // dampener has no effect for alpha = 1
    Vector<double> currPR(size, 1.0 / ((double) size));               // initialize each PR to 1/size(V)

//    cout << endl;

    int iter;
    for (iter = 1; iter <= MAX_ITER; iter++) {   // iterate the PR values until convergence or maximum reached
        Vector<double> product = matrixVectorMult(matrix, currPR);
        Vector<double> nextPR = vectorAdd(product, dampener);

        if (isConverged(currPR, nextPR)) {      // test for convergence
            currPR = nextPR;
            break;
        }
        currPR = nextPR;
        cout << "+" << flush;
    }

    cout << flush;
    cout << " (Converged in " << iter << " iterations!)" << endl;

    for (int i = 0; i < size; i++) {           // assign values to nodes
        graph->IDtoNode[i]->rank = currPR[i];
    }
}

/*
 * Create alignment by mapping the best PR-values nodes from each network to each other
 */
static Alignment *createAlignment(Graph *graphA, Graph *graphB) {
    sort(graphA->nodes.begin(), graphA->nodes.end(), Node::compareNodes);
    sort(graphB->nodes.begin(), graphB->nodes.end(), Node::compareNodes);

    Graph *graph1, *graph2;
    if (graphA->nodes.size() <= graphB->nodes.size()) {    // Graph1 must have <= nodes than Graph2
        graph1 = graphA;
        graph2 = graphB;
    } else {
        graph1 = graphB;
        graph2 = graphA;
    }

    Map<Node *, Node *> mapG1toG2, mapG2toG1;
    for (int i = 0; i < graph1->nodes.size(); i++) {       // populate the alignment
        Node *nodeG1 = graph1->nodes[i];
        Node *nodeG2 = graph2->nodes[i];
        mapG1toG2[nodeG1] = nodeG2;
        mapG2toG1[nodeG2] = nodeG1;
    }
    return new Alignment(graph1, graph2, mapG1toG2, mapG2toG1);
}

int main() {
    cout << "Welcome to the PageRank Network Aligner!" << endl;
    cout << "A network alignment allows you to compare two protein networks." << endl;
    cout << "Edge Coverage and Symmetric Substructure Score are two established measures " <<
            "used to calculated the topological similarity between the two protein networks." << endl;
    cout << "This algorithm utilizes PageRank values to generate the network alignment." << endl << endl;
    cout << "Choose two protein networks (by number) you would like to align!" << endl;

    Vector<string> proteinFiles{"CElegans", "SPombe", "MMusculus", "AThaliana",
                                "syeast0", "syeast05", "syeast10", "syeast20"};
    Vector<string> proteinNames{"C.Elegans", "S.Pombe", "M.Musculus", "A.Thaliana",
                                "S.Cerevisiae 00", "S.Cerevisiae 05",
                                "S.Cerevisiae 10", "S.Cerevisiae 20"};

    Vector<string> commonNames{"(roundworm)", "(fission yeast)", "(mouse)", "(plant)",
                               "(baker's yeast)", "(yeast + 5% noise)", "(yeast + 10% noise)",
                               "(yeast + 20% noise)"};

    for (int i = 1; i <= proteinNames.size(); i += 2) {
        cout << i << ") " + proteinNames[i - 1] << " " << commonNames[i - 1]  << "\t\t";
        cout << i + 1 << ") " + proteinNames[i] << " " << commonNames[i] << endl;
    }

    cout << endl;

    while (true) {
        string resp = getLine("Press enter to continue, type \"quit\" to stop:");
        if (toLowerCase(resp) == "quit") { break; }

        int index1 = getIntegerBetween("Choice for Graph 1: ", 1, proteinFiles.size());
        int index2 = getIntegerBetween("Choice for Graph 2: ", 1, proteinFiles.size());

        cout << "Loading networks..." << endl;

        string protein1 = proteinFiles[index1 - 1];
        string protein2 = proteinFiles[index2 - 1];

        string filename1 = "proteins/" + protein1 + ".gw";
        string filename2 = "proteins/" + protein2 + ".gw";
                                                              // error shouldn't happen
        if (!fileExists(filename1) || !fileExists(filename2)) { error("Invalid file(s) chosen!");}

        ifstream file1(filename1);
        ifstream file2(filename2);

        Graph *graph1 = new Graph(file1);
        Graph *graph2 = new Graph(file2);

        cout << "Calculating PageRanks for " << proteinNames[index1 - 1] << ": " << flush;
        computePageRank(graph1);

        cout << "Calculating PageRanks for " << proteinNames[index2 - 1] << ": " << flush;
        computePageRank(graph2);

        cout << "Topological Similarity: " << endl;

        Alignment *align = createAlignment(graph1, graph2);
        align->calcTopologicalMeasures();

        cout << "Edge Coverage: " << align->edgeCoverage << endl;
        cout << "Symmetric Substructure Score: " << align->substructureScore << endl << endl;

        delete align;                   // free the resources
    }

    cout << "Thank you for aligning!" << endl;
    return 0;
}
