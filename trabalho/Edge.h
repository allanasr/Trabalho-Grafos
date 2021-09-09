/**************************************************************************************************
 * Implementation of the TAD Edge
**************************************************************************************************/

#ifndef EDGE_H_INCLUDED
#define EDGE_H_INCLUDED

using namespace std;

// Definition of the Edge class
class Edge
{

    // Attributes
private:
    int target_id;
    int origin_id;
    Edge* next_edge;
    float weight;

public:
    // Constructor
    Edge(int target_id);
    Edge(int id, int target_id, float weight);
    // Destructor
    ~Edge();
    // Getters
    int getTargetId();
    int getOriginId();
    Edge* getNextEdge();
    float getWeight();
    // Setters
    void setNextEdge(Edge* edge);
    void setWeight(float weight);
    bool operator<(const Edge &edge2) const
    {
        return (weight < edge2.weight);
    }

};

#endif // EDGE_H_INCLUDED
