#ifndef MIN_HEAP_NODE_H
#define MIN_HEAP_NODE_H

class MinHeapNode
{
private:
    int id;
    int targetId;
    float weight;
    int targetCluster;

public:
    MinHeapNode() {};
    MinHeapNode(int id, float weight)
    {
        this->id = id;
        this->weight = weight;
        this->targetCluster = -1;
    };

    MinHeapNode(int id, int targetId, float weight, int cluster)
    {
        this->id = id;
        this->targetId = targetId;
        this->weight = weight;
        this->targetCluster = cluster;
    };
    ~MinHeapNode() {};
    int getId()
    {
        return id;
    };
    int getTargetId()
    {
        return targetId;
    };
    int getTargetCluster()
    {
        return targetCluster;
    };
    float getWeight()
    {
        return weight;
    };
    void setId(int id)
    {
        this->id = id;
    }
    void setTargetId(int targetId)
    {
        this->targetId = targetId;
    };
    void setTargetCluster(int cluster)
    {
        this->targetCluster = cluster;
    }
    void setWeight(int weight)
    {
        this->weight = weight;
    }
};
#endif
