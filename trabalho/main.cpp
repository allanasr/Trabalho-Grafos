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

Graph *leitura(ifstream &input_file, int directed, int weightedEdge, int weightedNode)
{

    int idNodeSource;
    int idNodeTarget;
    int order;

    input_file >> order;

    Graph *graph = new Graph(order, directed, weightedEdge, weightedNode);

    if (!graph->getWeightedEdge() && !graph->getWeightedNode())
    {

        while (input_file >> idNodeSource >> idNodeTarget)
        {
            if (!graph->searchNode(idNodeSource))
                graph->insertNode(idNodeSource);
            if (!graph->searchNode(idNodeTarget))
                graph->insertNode(idNodeTarget);
            graph->insertEdge(idNodeSource, idNodeTarget, 0);
        }
    }
    else if (graph->getWeightedEdge() && !graph->getWeightedNode())
    {

        float edgeWeight;

        while (input_file >> idNodeSource >> idNodeTarget >> edgeWeight)
        {

            if (!graph->searchNode(idNodeSource))
                graph->insertNode(idNodeSource);
            if (!graph->searchNode(idNodeTarget))
                graph->insertNode(idNodeTarget);
            graph->insertEdge(idNodeSource, idNodeTarget, edgeWeight);
        }
    }
    else if (graph->getWeightedNode() && !graph->getWeightedEdge())
    {

        float nodeSourceWeight, nodeTargetWeight;

        while (input_file >> idNodeSource >> nodeSourceWeight >> idNodeTarget >> nodeTargetWeight)
        {
            if (!graph->searchNode(idNodeSource))
                graph->insertNode(idNodeSource);
            if (!graph->searchNode(idNodeTarget))
                graph->insertNode(idNodeTarget);
            graph->insertEdge(idNodeSource, idNodeTarget, 0);
            graph->getNode(idNodeSource)->setWeight(nodeSourceWeight);
            graph->getNode(idNodeTarget)->setWeight(nodeTargetWeight);
        }
    }
    else if (graph->getWeightedNode() && graph->getWeightedEdge())
    {

        float nodeSourceWeight, nodeTargetWeight, edgeWeight;

        while (input_file >> idNodeSource >> nodeSourceWeight >> idNodeTarget >> nodeTargetWeight >> edgeWeight)
        {
            if (!graph->searchNode(idNodeSource))
                graph->insertNode(idNodeSource);
            if (!graph->searchNode(idNodeTarget))
                graph->insertNode(idNodeTarget);
            graph->insertEdge(idNodeSource, idNodeTarget, edgeWeight);
            graph->getNode(idNodeSource)->setWeight(nodeSourceWeight);
            graph->getNode(idNodeTarget)->setWeight(nodeTargetWeight);
        }
    }

    return graph;
}

Graph* leituraInstancia(ifstream& input_file, int directed, int weightedEdge, int weightedNode)
{

 //Variáveis para auxiliar na criação dos nós no Grafo
    int idNodeSource;
    int idNodeTarget;
    int numEdges;
    int id = 0;
    string clusterId;
    Cluster *cluster;

    // Criando objeto grafo
    Graph *graph = new Graph(0, directed, weightedEdge, weightedNode, true);

    //Leitura de arquivo
    while (getline(input_file, clusterId))
    {
        if (clusterId == "" || clusterId[0] == '\r')
        {
            break;
        }
        graph->insertNode(id, stoi(clusterId));
        id++;
    }

    if (!graph->getWeightedEdge() && !graph->getWeightedNode())
    {

        while (input_file >> idNodeSource >> idNodeTarget)
        {
            if (!graph->searchNode(idNodeSource) || !graph->searchNode(idNodeTarget))
            {
                cout << "Error: trying to add edge to a non-existent vetex" << endl;
                exit(1);
            }
            graph->insertEdge(idNodeSource, idNodeTarget, 0);
        }
    }
    else if (graph->getWeightedEdge() && !graph->getWeightedNode())
    {

        float edgeWeight;

        while (input_file >> idNodeSource >> idNodeTarget >> edgeWeight)
        {

            if (!graph->searchNode(idNodeSource) || !graph->searchNode(idNodeTarget))
            {
                cout << "Error: trying to add edge to a non-existent vetex" << endl;
                exit(1);
            }
            graph->insertEdge(idNodeSource, idNodeTarget, edgeWeight);
        }
    }
    else if (graph->getWeightedNode() && !graph->getWeightedEdge())
    {

        float nodeSourceWeight, nodeTargetWeight;

        while (input_file >> idNodeSource >> nodeSourceWeight >> idNodeTarget >> nodeTargetWeight)
        {
            if (!graph->searchNode(idNodeSource) || !graph->searchNode(idNodeTarget))
            {
                cout << "Error: trying to add edge to a non-existent vetex" << endl;
                exit(1);
            }
            graph->insertEdge(idNodeSource, idNodeTarget, 0);
            graph->getNode(idNodeSource)->setWeight(nodeSourceWeight);
            graph->getNode(idNodeTarget)->setWeight(nodeTargetWeight);
        }
    }
    else if (graph->getWeightedNode() && graph->getWeightedEdge())
    {

        float nodeSourceWeight, nodeTargetWeight, edgeWeight;

        while (input_file >> idNodeSource >> nodeSourceWeight >> idNodeTarget >> nodeTargetWeight >> edgeWeight)
        {
            if (!graph->searchNode(idNodeSource) || !graph->searchNode(idNodeTarget))
            {
                cout << "Error: trying to add edge to a non-existent vetex" << endl;
                exit(1);
            }
            graph->insertEdge(idNodeSource, idNodeTarget, edgeWeight);
            graph->getNode(idNodeSource)->setWeight(nodeSourceWeight);
            graph->getNode(idNodeTarget)->setWeight(nodeTargetWeight);
        }
    }

    return graph;
}

void escrita(Graph *graph, ofstream &output_file)
{
    if (graph != nullptr)
    {
        Node *auxNode = graph->getFirstNode();
        Edge *auxEdge;
        int tam = graph->getOrder();
        if (graph->getDirected())
        {
            output_file << "strict digraph grafo{" << endl;
            for (int i = 0; i < tam; i++)
            {
                auxEdge = auxNode->getFirstEdge();
                int cont = 0;
                while ((cont < (auxNode->getOutDegree() + auxNode->getInDegree())) && auxEdge != nullptr)
                {
                    output_file << "\t" << auxNode->getId() << " -> " << auxEdge->getTargetId() << ";" << endl;
                    auxEdge = auxEdge->getNextEdge();
                    cont++;
                }
                auxNode = auxNode->getNextNode();
            }
            output_file << "}" << endl;
        }
        else
        {
            output_file << "strict graph grafo{" << endl;
            for (int i = 0; i < tam; i++)
            {
                auxEdge = auxNode->getFirstEdge();
                int cont = 0;
                while ((cont < (auxNode->getOutDegree() + auxNode->getInDegree())) && auxEdge != nullptr)
                {
                    output_file << "\t" << auxNode->getId() << " -- " << auxEdge->getTargetId() << ";" << endl;
                    auxEdge = auxEdge->getNextEdge();
                    cont++;
                }
                auxNode = auxNode->getNextNode();
            }
            output_file << "}" << endl;
        }
    }
}

int menu()
{

    int selecao;

    cout << "MENU" << endl;
    cout << "----" << endl;
    cout << "[1] Subgrafo induzido por conjunto de vertices" << endl;
    cout << "[2] Caminho Minimo entre dois vertices - Dijkstra" << endl;
    cout << "[3] Caminho Minimo entre dois vertices - Floyd" << endl;
    cout << "[4] Arvore Geradora Minima de Prim" << endl;
    cout << "[5] Arvore Geradora Minima de Kruskal" << endl;
    cout << "[6] Imprimir caminhamento em largura" << endl;
    cout << "[7] Imprimir ordenacao topologica" << endl;
    cout << "[8] Algoritmo Guloso" << endl;
    cout << "[9] Algoritmo Guloso Randomizado " << endl;
    cout << "[0] Sair" << endl;

    cin >> selecao;

    return selecao;

}


void selecionar(int selecao, Graph *graph, ofstream &output_file)
{

    switch (selecao)
    {

    //Subgrafo induzido por um conjunto de vértices X;
    case 1:
    {
        bool *vertices = new bool[graph->getOrder()];

        for (int i = 0; i < graph->getOrder(); i++)
        {
            vertices[i] = false;
        }
        int x;
        cout << "quantos vertices formam o grafo induzido?" << endl;
        cin >> x;

        int nos;
        cout << "digite os vertices que fazem parte desse grafo induzido: " << endl;
        for (int j = 0; j < x; j++)
        {
            cin >> nos;
            vertices[nos] = true;
        }

        graph = graph->getVertexInduced(vertices, x, output_file);
        escrita(graph, output_file);

        break;
    }
        //Caminho mínimo entre dois vértices usando Dijkstra;
    case 2:
    {
        int source, target;
        cout << "Caminho mínimo entre dois vértices usando Dijkstra;" << endl;
        cout << "Digite o id do vértice de partida: ";
        cin >> source;
        cout << "Digite o id do vértice de destino: ";
        cin >> target;
        cout << "Comprimento do caminho minimo: " << graph->dijkstra(source, target, output_file) << endl;
        break;
    }

        //Caminho mínimo entre dois vértices usando Floyd;
    case 3:
    {
        int source, target;
        cout << "Caminho minimo entre dois vertices usando Floyd;" << endl;
        cout << "Digite o id do vertice de partida: ";
        cin >> source;
        cout << "Digite o id do vertice de destino: ";
        cin >> target;
        cout << "Comprimento do caminho minimo entre os vertices escolhidos: ";
        cout << graph->floydMarshall(source, target, output_file) << endl;
        break;
    }

    case 4:
    {

        Graph *tree = graph->agmPrim();
        escrita(tree, output_file);
        break;
    }

        //AGM Prim;
    case 5:
    {
        graph->agmKuskal(graph, output_file);
        break;
    }

        //Busca em largura;
    case 6:
    {
        graph->breadthFirstSearch(output_file);
        break;
    }
        //Ordenação Topologica;
    case 7:
    {
        graph->topologicalSorting(graph);
        break;
    };
    case 8:
    {
        time_t start, end;
        Graph *g;

        time(&start);
        g = graph->greed();
        time(&end);

        double time_taken = double(end - start);
        cout << "Execution time : " << fixed
             << time_taken << setprecision(3);
        cout << " sec " << endl;
        escrita(g, output_file);
        break;
    }
    case 9:
    {
        time_t start, end;
        Graph *g;

        time(&start);
        g = graph->greedRandom();
        time(&end);

        double time_taken = double(end - start);
        cout << "Execution time : " << fixed
             << time_taken << setprecision(3);
        cout << " sec " << endl;
        escrita(g, output_file);
        break;
    }
    case 0:
    {
        cout << "Execution finished" << endl;
        break;
    }
    default:
    {
        cout << " Error!!! invalid option!!" << endl;
    }
    }
}


int mainMenu(ofstream& output_file, Graph* graph)
{
    int selecao = 1;

    while(selecao != 0)
    {

        selecao = menu();

        if(output_file.is_open())
        {
            selecionar(selecao, graph, output_file);

        }

        else
            cout << "Unable to open the output_file" << endl;

        output_file << endl;

    }

    return 0;
}



int main(int argc, char const *argv[])
{

    if (argc != 6)
    {

        cout << "ERROR: Expecting: ./<program_name> <input_file> <output_file> <directed> <weighted_edge> <weighted_node> " << endl;
        return 1;

    }

    string program_name(argv[0]);
    string input_file_name(argv[1]);

    string instance;
    if(input_file_name.find("v") <= input_file_name.size())
    {
        string instance = input_file_name.substr(input_file_name.find("v"));
        cout << "Running " << program_name << " with instance " << instance << " ... " << endl;
    }

    ifstream input_file;
    input_file.open(argv[1], ios::in);

    std::ofstream output_file(argv[2], std::ofstream::out);


    Graph* graph;

    if(input_file.is_open())
    {

        graph = leitura(input_file, atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));

    }
    else
        cout << "Unable to open " << argv[1];


    mainMenu(output_file, graph);



    input_file.close();

    output_file.close();

    return 0;
}

