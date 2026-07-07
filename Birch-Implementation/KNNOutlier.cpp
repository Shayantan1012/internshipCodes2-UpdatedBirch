#include "KNNOutlier.h"

#include <queue>

namespace
{
struct KDNode
{
    size_t point=0, axis=0;
    int left=-1, right=-1;
};

class KDTree
{
    const std::vector<Point>& points;
    std::vector<KDNode> nodes;
    int root=-1;

    int build(std::vector<size_t>& indices,size_t begin,size_t end,size_t depth)
    {
        if(begin==end) return -1;
        const size_t axis=depth%points.front().size();
        const size_t middle=begin+(end-begin)/2;
        std::nth_element(indices.begin()+begin,indices.begin()+middle,indices.begin()+end,
            [&](size_t a,size_t b){return points[a][axis]<points[b][axis];});
        const int current=static_cast<int>(nodes.size());
        nodes.push_back({indices[middle],axis,-1,-1});
        nodes[current].left=build(indices,begin,middle,depth+1);
        nodes[current].right=build(indices,middle+1,end,depth+1);
        return current;
    }

    void nearest(int nodeIndex,size_t query,size_t k,
                 std::priority_queue<std::pair<double,size_t>>& heap) const
    {
        if(nodeIndex<0) return;
        const KDNode& node=nodes[nodeIndex];
        const double delta=points[query][node.axis]-points[node.point][node.axis];
        const int nearChild=delta<=0.0?node.left:node.right;
        const int farChild=delta<=0.0?node.right:node.left;
        nearest(nearChild,query,k,heap);
        if(node.point!=query)
        {
            const double distance=distanceSquared(points[query],points[node.point]);
            if(heap.size()<k) heap.push({distance,node.point});
            else if(distance<heap.top().first){heap.pop();heap.push({distance,node.point});}
        }
        if(heap.size()<k||delta*delta<=heap.top().first) nearest(farChild,query,k,heap);
    }

public:
    explicit KDTree(const std::vector<Point>& input):points(input)
    {
        std::vector<size_t> indices(points.size());
        std::iota(indices.begin(),indices.end(),size_t{0});
        nodes.reserve(points.size());
        root=build(indices,0,indices.size(),0);
    }

    double kDistance(size_t query,size_t k) const
    {
        std::priority_queue<std::pair<double,size_t>> heap;
        nearest(root,query,k,heap);
        return std::sqrt(heap.top().first);
    }
};
}

KNNOutlierResult detectKNNOutliers(const std::vector<Point>& points,
                                   size_t k,double outlierFraction)
{
    if(points.size()<2||k<1||k>=points.size())
        throw std::invalid_argument("Invalid k for kNN outlier detection");
    if(outlierFraction<0.0||outlierFraction>=1.0)
        throw std::invalid_argument("kNN outlier fraction must be in [0,1)");
    KDTree tree(points);
    KNNOutlierResult result;
    result.kDistances.resize(points.size());
    result.outlier.assign(points.size(),false);
    for(size_t i=0;i<points.size();++i) result.kDistances[i]=tree.kDistance(i,k);
    result.count=static_cast<size_t>(std::ceil(outlierFraction*points.size()));
    std::vector<size_t> order(points.size());
    std::iota(order.begin(),order.end(),size_t{0});
    std::sort(order.begin(),order.end(),[&](size_t a,size_t b){
        if(result.kDistances[a]!=result.kDistances[b]) return result.kDistances[a]>result.kDistances[b];
        return a<b;
    });
    for(size_t rank=0;rank<result.count;++rank) result.outlier[order[rank]]=true;
    result.cutoff=result.count?result.kDistances[order[result.count-1]]:
                               std::numeric_limits<double>::infinity();
    return result;
}
