#include "Graph.h"
#include "Node.h"
#include "Edge.h"
#include "MinHeapData.h"
#include <iostream>
#include <fstream>
#include <stack>
#include <queue>
#include <list>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <float.h>
#include <iomanip>
#include <algorithm>
#include <string.h>
#include <vector>
#include <iomanip>
#include <climits>

using namespace std;

using namespace std;

/**************************************************************************************************
 * Defining the Graph's methods
**************************************************************************************************/


Graph::Graph(int order, bool directed, bool weighted_edge, bool weighted_node)
{

    this->order = order;
    this->directed = directed;
    this->weighted_edge = weighted_edge;
    this->weighted_node = weighted_node;
    this->first_node = this->last_node = nullptr;
    this->number_edges = 0;
    this->node_cont = 0;
    adjacencia = new list<int>[order + 1];
}

Graph::Graph(int order, bool directed, bool weighted_edge, bool weighted_node, bool has_clusters)
{

    this->order = order;
    this->directed = directed;
    this->weighted_edge = weighted_edge;
    this->weighted_node = weighted_node;
    this->first_node = this->last_node = nullptr;
    this->has_clusters = has_clusters;
    this->first_cluster = nullptr;
    this->number_edges = 0;
    this->node_cont = 0;
    adjacencia = new list<int>;
}


Graph::~Graph()
{

    Node *next_node = this->first_node;

    while (next_node != nullptr)
    {

        next_node->removeAllEdges();
        Node *aux_node = next_node->getNextNode();
        delete next_node;
        next_node = aux_node;
    }
}

vector<Edge> edges;

int Graph::getOrder()
{
    return this->order;
}
int Graph::getNumberEdges()
{

    return this->number_edges;
}

bool Graph::getDirected()
{

    return this->directed;
}

bool Graph::getWeightedEdge()
{

    return this->weighted_edge;
}

bool Graph::getWeightedNode()
{

    return this->weighted_node;
}


Node *Graph::getFirstNode()
{

    return this->first_node;
}

Node *Graph::getLastNode()
{

    return this->last_node;
}


void Graph::insertNode(int id)
{
    Node *node = new Node(id);

    if (first_node == nullptr)
    {
        first_node = last_node = node;
    }
    else
    {
        last_node->setNextNode(node);
        last_node = node;
    }

    node_cont++;

    if (node_cont > order)
    {
        order++;
    }
}

Cluster *Graph::getCluster(int id)
{
    Cluster *c = first_cluster;

    while (c != nullptr && c->getId() != id)
    {
        c = c->getNextCluster();
    }

    return c;
}

void Graph::insertNode(int id, int clusterId)
{
    Node *node = new Node(id, clusterId);
    Cluster *c = getCluster(clusterId);

    if (c == nullptr)
    {
        c = insertCluster(clusterId);
    }

    c->insertElement(node);

    if (first_node == nullptr)
    {
        first_node = last_node = node;
    }
    else
    {
        last_node->setNextNode(node);
        last_node = node;
    }

    node_cont++;

    if (node_cont > order)
    {
        order++;
    }
}


void Graph::insertEdge(int id, int target_id, float weight)
{
    Edge edge(id, target_id, weight);
    edges.push_back(edge);

    if (!has_clusters)
        {
            // a lista de adjacencia e montada adicionando o vertice alvo no array referente ao vertice origem
            adjacencia[id].push_back(target_id);
        }

    Node *node, *target_node;
    node = getNode(id);

    if (node != nullptr)
    {
        target_node = getNode(target_id);


        if (target_node != nullptr)
        {
            node->insertEdge(target_id, weight);

            if (directed)
            {
                node->incrementOutDegree();
                target_node->incrementInDegree();
            }
            else
            {
                node->incrementInDegree();
                target_node->insertEdge(id, weight);
                target_node->incrementOutDegree();
            }
        }
    }
}

Cluster *Graph::insertCluster(int id)
{
    Cluster *c = new Cluster(id);
    if (first_cluster == nullptr)
    {
        first_cluster = c;
    }
    else
    {
        c->setNextCluster(first_cluster);
        first_cluster = c;
    }
    number_clusters++;

    return c;
}


void Graph::removeNode(int id)
{

}

bool Graph::searchNode(int id)
{
    Node *node = first_node;

    while (node != nullptr)
    {
        if (node->getId() == id)
        {
            return true;
        }

        node = node->getNextNode();
    }

    return false;
}


Node *Graph::getNode(int id)
{
    Node *node = first_node;

    while (node != nullptr && node->getId() != id)
    {
        node = node->getNextNode();
    }

    return node;
}

void Graph::breadthFirstSearch(ofstream &output_file)
{
    int tam = this->getOrder();
    bool *visitados = new bool[tam];
    vector<Node *> nos(tam);
    queue<Node *> fila;
    vector<Edge *> tree;
    Node *auxNode = this->getFirstNode();
    Edge *auxEdge;
    for (int i = 0; i < tam; i++)
    {
        *(visitados + i) = false;
        nos[auxNode->getId()] = auxNode;
        auxNode = auxNode->getNextNode();
    }
    fila.push(nos[0]);
    *(visitados + fila.front()->getId()) = true;
    while (!fila.empty() && tree.size() != (tam - 1))
    {
        auxNode = fila.front();
        fila.pop();
        auxEdge = auxNode->getFirstEdge();
        int cont = 0;
        while ((cont < (auxNode->getOutDegree() + auxNode->getInDegree())) && auxEdge != nullptr)
        {

            if (*(visitados + auxEdge->getTargetId()) == false)
            {

                *(visitados + auxEdge->getTargetId()) = true;
                fila.push(nos[auxEdge->getTargetId()]);
                Edge *galho = new Edge(auxNode->getId(), auxEdge->getTargetId(), auxEdge->getWeight());
                tree.push_back(galho);
                if (tree.size() == (tam - 1))
                    break;
            }
            auxEdge = auxEdge->getNextEdge();
            cont++;
        }
    }
    float weightResult = 0;
    cout << endl
         << "Busca em Largura" << endl
         << endl;
    if (this->getDirected())
    {
        output_file << "digraph busca{" << endl;
        for (int i = 0; i < tree.size(); i++)
        {
            cout << "(" << tree[i]->getOriginId() << ", " << tree[i]->getTargetId() << ") - peso = " << tree[i]->getWeight() << endl;
            weightResult += tree[i]->getWeight();
            output_file << "\t" << tree[i]->getOriginId() << " -> " << tree[i]->getTargetId() << ";" << endl;
        }
    }
    else
    {
        output_file << "graph busca{" << endl;
        for (int i = 0; i < tree.size(); i++)
        {
            cout << "(" << tree[i]->getOriginId() << ", " << tree[i]->getTargetId() << ") - peso = " << tree[i]->getWeight() << endl;
            weightResult += tree[i]->getWeight();
            output_file << "\t" << tree[i]->getOriginId() << " -- " << tree[i]->getTargetId() << ";" << endl;
        }
    }
    output_file << "}" << endl;
    cout << endl
         << "Peso total da arvore: " << weightResult << endl
         << endl;
}

float Graph::floydMarshall(int idSource, int idTarget, ofstream &output_file)
{
    Node *node;
    Edge *edge;

    bool isSource = false, isTarget = false;
    for (node = first_node; node != nullptr; node = node->getNextNode())
    {
        if (node->getId() == idSource)
            isSource = true;
        if (node->getId() == idTarget)
            isTarget = true;
    }

    if (!isSource || !isTarget)
    {
        cout << "Entrada invalida!" << endl;
        return INT_MAX;
    }

    int i, j;

    float matDistancias[order][order];

    for (int i = 0; i < order; i++)
    {
        for (int j = 0; j < order; j++)
        {
            matDistancias[i][j] = INT_MAX;
        }
    }


    Node *x;
    Edge *y;
    int z;
    x = first_node;
    y = x->getFirstEdge();
    z = y->getTargetId();

    for (x = first_node; x != nullptr; x = x->getNextNode())
    {
        for (y = x->getFirstEdge(); y != x->getLastEdge(); y = y->getNextEdge())
        {
            z = y->getTargetId();
            matDistancias[x->getId() - 1][z - 1] = y->getWeight();
        }
        y = x->getLastEdge();
        z = y->getTargetId();
        matDistancias[x->getId() - 1][z - 1] = y->getWeight();
    }

    for (int i = 0; i < order; i++)
    {
        for (int j = 0; j < order; j++)
        {
            if (i == j)
                matDistancias[i][j] = 0;
        }
    }

    int matR[order][order];
    for (int i = 0; i < order; i++)
    {
        for (int j = 0; j < order; j++)
        {
            matR[i][j] = j;
        }
    }

    for (int k = 0; k < order; k++)
    {
        for (i = 0; i < order; i++)
        {
            for (j = 0; j < order; j++)
            {
                if (matDistancias[i][j] > matDistancias[i][k] + matDistancias[k][j])
                {
                    matDistancias[i][j] = matDistancias[i][k] + matDistancias[k][j];
                    matR[i][j] = matR[i][k];
                }
            }
        }
    }

    int s = idSource - 1;
    int e = idTarget - 1;
    int cont = 0;

    int vet[order];

    for (int i = 0; i < order; i++)
    {
        vet[i] = -1;
    }

    while (s != e)
    {
        vet[cont] = s + 1;
        s = matR[s][e];
        cont++;
    }

    vet[cont] = e + 1;

    output_file << "graph caminho_minimo{" << endl;
    for (int i = 0; vet[i + 1] != -1; i++)
    {
        output_file << vet[i] << " -- " << vet[i + 1] << endl;
    }
    output_file << "}";

    return cont;
}


float Graph::dijkstra(int idSource, int idTarget, ofstream &output_file)
{
    float pi[order];
    int prev[order];
    MinHeap s_barra(order);
    MinHeapNode *minNode;
    Node *node;
    Edge *edge;

    bool isSource = false, isTarget = false;
    for (node = first_node; node != nullptr; node = node->getNextNode())
    {
        if (node->getId() == idSource)
            isSource = true;
        if (node->getId() == idTarget)
            isTarget = true;
    }
    if (!isSource || !isTarget)
    {
        std::cout << "Entrada inválida!" << endl;
        return INT_MAX;
    }

    int infinity = INT_MAX / 2;

    node = first_node;
    for (int i = 0; i < order; i++)
    {
        pi[node->getId()] = infinity;
        s_barra.insertKey(new MinHeapNode(node->getId(), infinity));
        node = node->getNextNode();
    }
    pi[idSource] = 0;
    s_barra.decreaseKey(idSource, 0);

    int prevId = -1;
    prev[idSource] = prevId;

    while (!s_barra.isEmpty())
    {
        minNode = s_barra.extractMin();
        prevId = minNode->getId();

        for (node = first_node; node->getId() != minNode->getId(); node = node->getNextNode())
            ;

        edge = node->getFirstEdge();
        while (edge != nullptr)
        {
            int pi_estrela = pi[minNode->getId()] + edge->getWeight();

            int edgeTargetId = edge->getTargetId();

            if (pi_estrela < pi[edgeTargetId])
            {
                pi[edgeTargetId] = pi_estrela;
                prev[edgeTargetId] = prevId;

                int idx = s_barra.getIndexOf(edgeTargetId);
                if (idx == -1)
                {
                    s_barra.insertKey(new MinHeapNode(edgeTargetId, pi[edgeTargetId]));
                }
                else
                {
                    s_barra.decreaseKey(idx, pi[edgeTargetId]);
                }
            }

            edge = edge->getNextEdge();
        }

        delete minNode;
    }

    int prevNode = idTarget;
    int cont = 0;
    output_file << "graph caminho_minimo{" << endl;

    for (int node = prev[idTarget]; node != -1; node = prev[node])
    {
        output_file << prevNode << " -- " << node << endl;
        prevNode = node;
        cont++;
    }
    output_file << "}";

    return cont;
}

void Graph::topologicalSorting(Graph *graph)
{
    stack<int> pilhaTopologica;
    int tamGrafo = graph->getOrder() + 1;
    vector<bool> nosVisitados(tamGrafo, false);

    for (int i = 0; i < tamGrafo; i++)
    {
        if (nosVisitados[i] == false)
        {
            auxTopologicalSorting(i, nosVisitados, pilhaTopologica);
        }
    }

    cout << "\nOrdenacao topologica:" << endl
         << "< ";

    while (pilhaTopologica.empty() == false && pilhaTopologica.top() != 0)
    {
        if (pilhaTopologica.size() > 2)
        {
            cout << pilhaTopologica.top() << ", ";
            pilhaTopologica.pop();
        }
        else
        {
            cout << pilhaTopologica.top() << " ";
            pilhaTopologica.pop();
        }
    }
    cout << ">" << endl
         << endl;
}

void Graph::auxTopologicalSorting(int index, vector<bool> &nosVisitados, stack<int> &pilhaTopologica)
{
    nosVisitados[index] = true;


    list<int>::iterator i;
    for (i = adjacencia[index].begin(); i != adjacencia[index].end(); ++i)
    {
        if (!nosVisitados[*i])
        {
            auxTopologicalSorting(*i, nosVisitados, pilhaTopologica);
        }
    }

    pilhaTopologica.push(index);
}

void breadthFirstSearch(ofstream &output_file)
{
}

Graph *Graph::getVertexInduced(bool *vertices, int x, ofstream &output_file)
{
    vector<Edge> arestas;
    Graph *g1 = new Graph(x, this->getDirected(), this->getWeightedEdge(), this->getWeightedNode());

    Node *node = this->getFirstNode();

    Edge *edge;


    for (int i = 0; i < this->getOrder(); i++)
    {
        if (vertices[node->getId()])
        {
            for (edge = node->getFirstEdge(); edge != node->getLastEdge(); edge = edge->getNextEdge())
            {
                if (vertices[edge->getTargetId()])
                {
                    arestas.push_back(Edge(node->getId(), edge->getTargetId(), edge->getWeight()));
                }
            }
            if (vertices[edge->getTargetId()])
            {
                arestas.push_back(Edge(node->getId(), edge->getTargetId(), edge->getWeight()));
            }
            vertices[node->getId()] = false;
            g1->insertNode(node->getId());
        }
        node = node->getNextNode();
    }


    cout << endl
         << "Subgrafo induzido por um conjunto de vertices" << endl
         << endl;

    for (int j = 0; j < arestas.size();  j++)
    {


        g1->insertEdge(arestas[j].getOriginId(), arestas[j].getTargetId(), arestas[j].getWeight());
    }

    return g1;
}

int searchForSubset(int subset[], int i)
{
    if (subset[i] == -1)
        return i;
    return searchForSubset(subset, subset[i]);
}
void join(int subset[], int v1, int v2)
{
    int v1_set = searchForSubset(subset, v1);
    int v2_set = searchForSubset(subset, v2);
    subset[v1_set] = v2_set;
}

Graph *Graph::agmKuskal(Graph *graph, ofstream &output_file)
{
    vector<Edge> tree;
    int size_edges = edges.size();

    sort(edges.begin(), edges.end());

    int V = graph->getOrder();
    int *subset = new int[V + 1];

    memset(subset, -1, sizeof(int) * V);

    for (int i = 0; i < size_edges; i++)
    {
        int v1 = searchForSubset(subset, edges[i].getOriginId());
        int v2 = searchForSubset(subset, edges[i].getTargetId());

        if (v1 != v2)
        {
            tree.push_back(edges[i]);
            join(subset, v1, v2);
        }
    }

    int size_tree = tree.size();

    cout << endl;
    cout << "Arvore Geradora Minima usando algoritmo de Kruskal" << endl;
    float weightResult = 0;
    for (int i = 0; i < size_tree; i++)
    {
        int v1 = tree[i].getOriginId();
        int v2 = tree[i].getTargetId();
        int w = tree[i].getWeight();
        weightResult = w + weightResult;
        cout << "(" << v1 << ", " << v2 << ") - peso = " << w << endl;
    }
    cout << "Peso total do arvore: " << weightResult << endl;
    cout << endl;

    output_file << "strict graph kruskal{" << endl;
    for (int i = 0; i < size_tree; i++)
    {
        output_file << "\t" << tree[i].getOriginId() << " -- " << tree[i].getTargetId() << ";" << endl;
    }
    output_file << "}";
}

Graph *Graph::agmPrim()
{
    if (!this->getDirected())
    {
        int tam = this->getOrder();
        int prox[tam];
        Graph *tree = new Graph(this->getOrder(), this->getDirected(), this->getWeightedEdge(), this->getWeightedNode());
        vector<Edge> custo;
        vector<Node *> nos(tam);
        Node *auxNode = this->getFirstNode();
        int primeiro = auxNode->getId();
        Edge *auxEdge;
        for (int i = 0; i < tam; i++)
        {
            prox[i] = primeiro;
            custo.push_back(Edge(i, primeiro, INT_MAX));
            nos[auxNode->getId()] = auxNode;
            tree->insertNode(auxNode->getId());
            auxNode = auxNode->getNextNode();
        }
        int i, j, k = 0;
        j = (primeiro);
        while (k < tam)
        {
            prox[j] = -1;
            auxNode = nos[j];
            auxEdge = auxNode->getFirstEdge();
            int cont = 0;
            while ((cont < (auxNode->getOutDegree() + auxNode->getInDegree())) && auxEdge != nullptr)
            {
                if (prox[auxEdge->getTargetId()] != -1 && (custo[auxEdge->getTargetId()].getWeight() > auxEdge->getWeight()))
                {
                    custo[auxEdge->getTargetId()] = Edge(auxNode->getId(), auxEdge->getTargetId(), auxEdge->getWeight());
                    prox[auxEdge->getTargetId()] = auxNode->getId();
                }
                auxEdge = auxEdge->getNextEdge();
                cont++;
            }

            for (i = 0; i < tam; i++)
                if (prox[i] != -1)
                {
                    j = i;
                    break;
                }
            for (; i < tam; i++)
                if (prox[i] != -1 && custo[i].getWeight() < custo[j].getWeight())
                    j = i;
            if (custo[j].getWeight() == INT_MAX)
            {
                cout << "Arvore Geradora Minima usando algoritmo de Prim" << endl
                     << endl;
                cout << "Grafo desconexo" << endl
                     << endl;
                return nullptr;
            }
            k++;
        }
        cout << "zé" << endl;
        sort(custo.begin(), custo.end());
        cout << "Arvore Geradora Minima usando algoritmo de Prim" << endl
             << endl;
        float weightResult = 0;
        for (int i = 0; i < tam - 1; i++)
        {

            cout << "(" << custo[i].getOriginId() << ", " << custo[i].getTargetId() << ") - peso = " << custo[i].getWeight() << endl;
            weightResult += custo[i].getWeight();
            tree->insertEdge(custo[i].getOriginId(), custo[i].getTargetId(), custo[i].getWeight());
        }
        cout << endl
             << "Peso total da arvore: " << weightResult << endl
             << endl;
        return tree;
    }
    else
    {
        cout << "Arvore Geradora Minima usando algoritmo de Prim" << endl
             << endl;
        cout << "Grafo direcionado" << endl
             << endl;
        return nullptr;
    }
}

void addEdges(vector<Edge *> *vet, Node *sourceNode, Graph *g, bool visitedClusters[])
{
    Edge *edge = g->getNode(sourceNode->getId())->getFirstEdge();

    while (edge != nullptr)
    {
        if (!visitedClusters[g->getNode(edge->getTargetId())->getCluster() - 1])
            vet->push_back(edge);
        edge = edge->getNextEdge();
    }
}

// Algoritmo guloso para o problema de otmizaçao da árvore geradora mínima generalizada
// Adaptaçao de PRIM
Graph *Graph::greed()
{
    if (!this->has_clusters)
    {
        cout << "Erro: Grafo nao tem grupos para se realizar arvore minima generalizada. Tente o algoritmo de Prim ou de Kruskal." << endl;
        return nullptr;
    }

    Graph *minimalTree = nullptr;
    Graph *tree;
    int minCost = INT_MAX;
    int currentCost;
    Cluster *cluster;
    vector<Edge *> k(number_edges);
    bool visitedClusters[number_clusters];
    Node *node;
    Edge *edge;

    for (cluster = first_cluster; cluster != nullptr; cluster = cluster->getNextCluster())
    {
        tree = new Graph(0, directed, weighted_edge, weighted_node, has_clusters);

        for (int i = 0; i < number_clusters; i++)
        {
            visitedClusters[i] = false;
        }
        currentCost = 0;
        visitedClusters[cluster->getId() - 1] = true;

        for (int i = 0; i < cluster->getSize(); i++)
        {
            addEdges(&k, cluster->getElement(i), this, visitedClusters);
        }

        for (int i = 1; i < number_clusters; i++)
        {
            for (node = tree->getFirstNode(); node != nullptr; node = node->getNextNode())
            {
                addEdges(&k, node, this, visitedClusters);
            }

            edge = k[0];
            for (int i = 1; i < k.size(); i++)
            {
                if (k[i]->getWeight() < edge->getWeight())
                {
                    edge = k[i];
                }
            }

            if (!tree->searchNode(edge->getOriginId()))
            {
                tree->insertNode(edge->getOriginId());
            }
            if (!tree->searchNode(edge->getTargetId()))
            {
                tree->insertNode(edge->getTargetId());
            }
            tree->insertEdge(edge->getOriginId(), edge->getTargetId(), edge->getWeight());
            currentCost += edge->getWeight();
            visitedClusters[this->getNode(edge->getTargetId())->getCluster() - 1] = true;

            k.clear();
        }

        if (currentCost < minCost)
        {
            if (minimalTree != nullptr)
            {
                delete minimalTree;
            }
            minimalTree = tree;
            minCost = currentCost;
        }
        else
        {
            delete tree;
        }
    }

    std::cout << "Custo total da árvore: " << minCost << endl;
    return minimalTree;
}

bool compareEdgesWeight(Edge *a, Edge *b)
{
    return a->getWeight() < b->getWeight();
}

Graph *Graph::greedRandom()
{
    if (!this->has_clusters)
    {
        std::cout << "Erro: Grafo nao tem grupos para se realizar arvore minima generalizada. Tente o algoritmo de Prim ou de Kruskal." << endl;
        return nullptr;
    }

    Graph *minimalTree = nullptr;
    Graph *tree;
    int minCost = INT_MAX;
    int currentCost;
    Cluster *cluster;
    vector<Edge *> k(number_edges);
    float alpha = 0.05;
    float max;
    bool visitedClusters[number_clusters];
    Node *node;
    Edge *edge;

    srand(time(NULL));

    for (cluster = first_cluster; cluster != nullptr; cluster = cluster->getNextCluster())
    {
        tree = new Graph(0, directed, weighted_edge, weighted_node, has_clusters);

        for (int i = 0; i < number_clusters; i++)
        {
            visitedClusters[i] = false;
        }
        currentCost = 0;
        visitedClusters[cluster->getId() - 1] = true;

        for (int i = 0; i < cluster->getSize(); i++)
        {
            addEdges(&k, cluster->getElement(i), this, visitedClusters);
        }

        for (int i = 1; i < number_clusters; i++)
        {
            for (node = tree->getFirstNode(); node != nullptr; node = node->getNextNode())
            {
                addEdges(&k, node, this, visitedClusters);
            }
            sort(k.begin(), k.end(), compareEdgesWeight);
            max = k.front()->getWeight() + alpha * (k.back()->getWeight() - k.front()->getWeight());

            for (int i = 0; i < k.size(); i++)
            {
                if (k[i]->getWeight() == max)
                {
                    edge = k[rand() % (i + 1)];
                    break;
                }
                else if (k[i]->getWeight() > max)
                {
                    edge = k[rand() % i];
                    break;
                }
            }

            if (!tree->searchNode(edge->getOriginId()))
            {
                tree->insertNode(edge->getOriginId());
            }
            if (!tree->searchNode(edge->getTargetId()))
            {
                tree->insertNode(edge->getTargetId());
            }
            tree->insertEdge(edge->getOriginId(), edge->getTargetId(), edge->getWeight());
            currentCost += edge->getWeight();
            visitedClusters[this->getNode(edge->getTargetId())->getCluster() - 1] = true;

            k.clear();
        }


        if (currentCost < minCost)
        {
            if (minimalTree != nullptr)
            {
                delete minimalTree;
            }
            minimalTree = tree;
            minCost = currentCost;
        }
        else
        {
            delete tree;
        }
    }

    std::cout << "Custo total da árvore: " << minCost << endl;
    return minimalTree;
}


