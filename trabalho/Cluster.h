#ifndef CLUSTER_H
#define CLUSTER_H

#include "Node.h"
#include "vector"

class Cluster
{
private:
    int id;
    vector<Node *> elements;
    Cluster *next_cluster;

public:
    Cluster(int id);
    ~Cluster();

    int getId();
    int getSize();
    Node *getElement(int i);
    Cluster *getNextCluster();
    void setNextCluster(Cluster *c);

    void insertElement(Node *n);
};
#endif

