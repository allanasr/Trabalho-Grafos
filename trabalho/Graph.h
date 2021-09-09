/**************************************************************************************************
 * Implementation of the TAD Graph
**************************************************************************************************/

#ifndef GRAPH_H_INCLUDED
#define GRAPH_H_INCLUDED
#include "Node.h"
#include "Edge.h"
#include "Cluster.h"
#include <fstream>
#include <stack>
#include <list>
#include <vector>

using namespace std;

class Graph
{

    //Atributes
private:
    int order;
    int number_edges;
    bool directed;
    bool weighted_edge;
    bool weighted_node;
    bool has_clusters;
    int node_cont;
    Node *first_node;
    Node *last_node;
    list<int> *adjacencia; /// lista de adjacencia
    Cluster *first_cluster;
    int number_clusters;

public:
    //Constructor
    Graph(int order, bool directed, bool weighted_edge, bool weighted_node);
    Graph(int order, bool directed, bool weighted_edge, bool weighted_node, bool has_clusters);
    //Destructor
    ~Graph();
    //Getters
    int getOrder();
    int getNumberEdges();
    bool getDirected();
    bool getWeightedEdge();
    bool getWeightedNode();
    Node *getFirstNode();
    Node *getLastNode();
    Cluster *getCluster(int id);
    //Other methods
    void insertNode(int id);
    void insertNode(int id, int cluster);
    void insertEdge(int id, int target_id, float weight);
    Cluster *insertCluster(int id);
    void removeNode(int id);
    bool searchNode(int id);
    Node *getNode(int id);

    //methods phase1
    void topologicalSorting(Graph *graph);
    void auxTopologicalSorting(int index, vector<bool> &nosVisitados, stack<int> &Pilha);
    void breadthFirstSearch(ofstream &output_file);
    Graph *getVertexInduced(bool *vertices, int x, ofstream &output_file);
    Graph *agmKuskal(Graph *graph,ofstream &output_file);
    Graph *agmPrim();
    float floydMarshall(int idSource, int idTarget, ofstream &output_file);
    float dijkstra(int idSource, int idTarget, ofstream &output_file);

    //methods phase 2
    Graph *greed();
    Graph *greedRandom();

private:
};

#endif // GRAPH_H_INCLUDED
