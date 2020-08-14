/*
 * CS 106X Final Project. Fall 2019.
 * Rishi Desai
 */

#include <iostream>
#include <fstream>
#include <string>
#include "graphtypes.h"
#include "strlib.h"
#include "filelib.h"
using namespace std;

/************************************
 * Node Method Implementations
 */

Node::Node(string &name, int ID) {
    this->name = name;
    this->ID = ID;
}

Node::~Node() {
    for (Edge *edge : outEdges) { delete edge; }
}

bool Node::operator <(Node *other) { return this->name < other->name; }

string Node::toString() {
    string ret = name + ": ";
    for (Edge *edge : outEdges) { ret += edge->finish->name + ", "; }
    return ret;
}

bool Node::compareNodes(Node *node1, Node *node2) {
    return node1->rank > node2->rank;
}

/************************************
 * Edge Method Implementations
 */

Edge::Edge(Node *start, Node *finish) {
    if (start == nullptr) cout << "FAILED";
    if (finish == nullptr) cout << "FAILED2";
    this->start = start;
    this->finish = finish;
    this->name = start->name + finish->name;
}

Edge::~Edge() {}

bool Edge::operator <(Edge *other) { return this->name < other->name; }

string Edge::toString() { return start->name + ":" + finish->name;}

string Edge::makeEdgeString(Edge *edge) {
    return makeEdgeString(edge->start, edge->finish);
}

string Edge::makeEdgeString(Node *node1, Node *node2) {
    Vector<string> names{node1->name, node2->name};
    sort(names.begin(), names.end());
    return names[0] + "_" + names[1];
}

pair<Node *, Node *> Edge::extractNodes(Graph *graph, string &edge) {
    int pos = edge.find('_');
    string nameA = edge.substr(0, pos), nameB = edge.substr(pos + 1);

    Node *nodeA = graph->nodeIndex[nameA], *nodeB = graph->nodeIndex[nameB];
    return make_pair(nodeA, nodeB);
}

/************************************
 * Graph Method Implementations
 */

Graph::Graph(ifstream &file) {
    readGraphFile(file);
}

Graph::~Graph() {
    for (Node *node : nodes) { delete node; }
}

void Graph::readGraphFile(ifstream &infile) {
    string line;                      // skip version info
    for (int i = 0; i < 4; i++ ) { getline(infile, line); }

    getline(infile, line);            // read the nodes
    const int numNodes = stringToInteger(line);
    for (int i = 0; i < numNodes; i++) {
        getline(infile, line);
        line = trim(line);
        string nodeName = line.substr(2, line.size() - 4);

        Node *node = new Node(nodeName, i);
        nodes.add(node);
        nodeIndex[node->name]= node;
        IDtoNode[node->ID] = node;
    }

    getline(infile, line);             // read the edges
    const int numEdges = stringToInteger(line);
    for (int i = 0; i < numEdges; i++) {
        getline(infile, line);

        string index1, index2;         // .gw files use 1-based indexing when listing edges
        stringstream ss(line);
        getline(ss, index1, ' ');
        getline(ss, index2, ' ');

        int ID1 = stringToInteger(index1) - 1, ID2 = stringToInteger(index2) - 1;

        Node *node1 = IDtoNode[ID1], *node2 = IDtoNode[ID2];
        Edge *edge12 = new Edge(node1, node2), *edge21 = new Edge(node2, node1);

        node1->outEdges.add(edge12);
        node1->inEdges.add(edge21);
        node1->outEdgeNodeNames.add(node2->name);

        node2->outEdges.add(edge21);
        node2->inEdges.add(edge12);
        node2->outEdgeNodeNames.add(node1->name);

        edgeIndex.put(edge12->name, edge12);
        edgeIndex.put(edge21->name, edge21);

        this->edges.add(Edge::makeEdgeString(edge12));
    }
    infile.close();
}

string Graph::toString() {
    string ret = "";
    for (Node *node : this->nodes) { ret += node->toString() + "\n"; }
    return ret;
}

/************************************
 * Alignment Method Implementations
 */

Alignment::Alignment(Graph *graph1, Graph *graph2, Map<Node *, Node *> &alignG1toG2,
                     Map<Node *, Node *> &alignG2toG1) {
    this->alignG1toG2 = alignG1toG2;
    this->alignG2toG1 = alignG2toG1;
    this->graph1 = graph1;
    this->graph2 = graph2;
}

Alignment::Alignment(ifstream &file, Graph *graph1, Graph *graph2) {
    this->graph1 = graph1;
    this->graph2 = graph2;
    loadAlignFile(file);
}

Alignment::~Alignment() {
    delete this->graph1;
    delete this->graph2;
}

void Alignment::loadAlignFile(ifstream &file) {
    string line;
    while (getline(file, line)) {
        string nameG1, nameG2;
        stringstream ss(line);

        getline(ss, nameG1, '\t');
        getline(ss, nameG2, '\t');

        Node *nodeG1 = graph1->nodeIndex[nameG1];
        Node *nodeG2 = graph2->nodeIndex[nameG2];

        this->alignG1toG2[nodeG1] = nodeG2;
        this->alignG2toG1[nodeG2] = nodeG1;
    }
    file.close();
}

void Alignment::calcTopologicalMeasures() {
    int coveredEdges = 0;
    for (string edge1 : this->graph1->edges) {
        pair<Node *, Node *> edgeNodes = Edge::extractNodes(graph1, edge1);

        Node *graph2A = alignG1toG2[edgeNodes.first], *graph2B = alignG1toG2[edgeNodes.second];

        string graph2Edge = Edge::makeEdgeString(graph2A, graph2B);
        if (this->graph2->edges.contains(graph2Edge)) { coveredEdges++; }
    }

    int inducedEdges = 0;                       // # edges between aligned nodes in G2
    for (string edge2 : this->graph2->edges) {
        pair<Node *, Node *> edgeNodes = Edge::extractNodes(graph2, edge2);

        if (alignG2toG1.containsKey(edgeNodes.first) &&
            alignG2toG1.containsKey(edgeNodes.second)) { inducedEdges++; }
    }

    const int edgesG1size = this->graph1->edges.size();
    edgeCoverage = coveredEdges / ((double) edgesG1size);
    substructureScore = coveredEdges / (double (edgesG1size + inducedEdges - coveredEdges));
}

void Alignment::calcNodeCorrectness(Alignment *correct) {
    int numCorrect = 0;
    for (Node *nodeG1 : this->alignG1toG2.keys()) {
        if (this->alignG1toG2[nodeG1]->name == correct->alignG1toG2[nodeG1]->name) {
            numCorrect++;
        }
    }
    nodeCorrectness = numCorrect / ((double) this->alignG1toG2.size());
}

void Alignment::writeAlignFile(ofstream &file) {
    for (Node *nodeG1 : this->alignG1toG2.keys()) {
        file << nodeG1->name << "\t" << this->alignG1toG2[nodeG1]->name << endl;
    }
    file.close();
}
