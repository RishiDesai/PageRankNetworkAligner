#ifndef GRAPHTYPES_H
#define GRAPHTYPES_H

/*
 * CS 106X Final Project. Fall 2019.
 * PageRank Network Aligner
 * Rishi Desai
 */

#include <iostream>
#include <string>
#include "set.h"
#include "hashmap.h"
#include "hashset.h"
using namespace std;

struct Edge;
struct Node;
struct Graph;
struct Alignment;

struct Node {

    string name;
    Set<Edge *> outEdges, inEdges;
    HashSet<string> outEdgeNodeNames;

    int ID;
    double rank;

    Node(string &name, int ID);
    ~Node();

    bool operator <(Node *other);

    string toString();

    /*
     * Sorts nodes by rank greatest to least
     */
    static bool compareNodes(Node *node1, Node *node2);
};

struct Edge {

    Node *start, *finish;
    string name;

    Edge(Node *start, Node *finish);
    ~Edge();

    bool operator <(Edge *other);

    string toString();

    /*
     * String representation of edge (A->B and B->A give same string)
     */
    static string makeEdgeString(Edge *edge);

    static string makeEdgeString(Node *node1, Node *node2);

    /*
     * Extract nodes from string representation of edge
     */
    static pair<Node *, Node *> extractNodes(Graph *graph,  string &edge);
};

struct Graph {

    Graph(ifstream &file);
    ~Graph();

    HashMap<string, Node *> nodeIndex;  // node name to node object
    HashMap<string, Edge *> edgeIndex;

    HashMap<int , Node *> IDtoNode;     // nodes given ID from 0...n-1
    Vector<Node *> nodes;
    HashSet<string> edges;              // set of string reps of each edge

    string toString();

private:
    void readGraphFile(ifstream &infile);
};

/*
 * An alignment is a mapping of each node from G1=(V1,E1) to a unique node in G2=(V2,E2),
 * where |V1| <= |V2|
 */
struct Alignment {

    Map<Node *, Node *> alignG1toG2, alignG2toG1;     // the mapping from V1 to V2, and map reversed for V2 to V1
    Graph *graph1, *graph2;
    double edgeCoverage = -1, substructureScore = -1; // topological measures
    double nodeCorrectness = -1;                      // correctness measure (not displayed to user)

    Alignment(Graph *graph1, Graph *graph2, Map<Node *, Node *> &alignG1toG2,
                                            Map<Node *, Node *> &alignG2toG1);

    Alignment(ifstream &file, Graph *graph1, Graph *graph2);

    ~Alignment();

    /*
     * Calculate Edge Coverage and Symmetric Substructure Score (aka S3)
     */
    void calcTopologicalMeasures();

    /*
     * Calculate Node Correctness given the correct alignment
     */
    void calcNodeCorrectness(Alignment *correct);

    /*
     * Write out alignment to a file
     */
    void writeAlignFile(ofstream &file);

private:
    void loadAlignFile(ifstream &file);
};

#endif // GRAPHTYPES_H
