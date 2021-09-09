#include "Cluster.h"
#include "Node.h"
#include <iostream>

Cluster::Cluster(int id)
{
    this->id = id;
}

Cluster::~Cluster()
{
    if (next_cluster != nullptr)
    {
        delete next_cluster;
        next_cluster = nullptr;
    }
}

int Cluster::getId()
{
    return id;
}

int Cluster::getSize()
{
    return elements.size();
}

Node *Cluster::getElement(int i)
{
    return elements[i];
}

Cluster *Cluster::getNextCluster()
{
    return next_cluster;
}

void Cluster::setNextCluster(Cluster *c)
{
    next_cluster = c;
}

void Cluster::insertElement(Node *n)
{
    elements.push_back(n);
}
