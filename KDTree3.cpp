/*
KDTree3.cpp
Written by Matthew Fisher

A 3D KD-tree that supports fast points-within-rectangle queries
*/
#include "KDTree3.h"

#pragma comment(lib, "ANN.lib")

KDTree3::KDTree3()
{
    nnIdx = NULL;
    dists = NULL;
    queryPt = NULL;
    dataPts = NULL;
    kdTree = NULL;
}

KDTree3::~KDTree3()
{
    FreeMemory();
}

void KDTree3::FreeMemory()
{
    if(nnIdx)
    {
        delete[] nnIdx;
        nnIdx = NULL;
    }
    if(dists)
    {
        delete[] dists;
        dists = NULL;
    }
    if(kdTree)
    {
        delete kdTree;
        kdTree = NULL;
    }
    if(queryPt)
    {
        annDeallocPt(queryPt);
        queryPt = NULL;
    }
    if(dataPts)
    {
        annDeallocPts(dataPts);
        dataPts = NULL;
    }
}

void KDTree3::BuildTree(const std::vector<Eigen::Vector3f> &Points)
{
    FreeMemory();
    unsigned int PointCount = Points.size();
   
    queryPt = annAllocPt(3); // allocate query point
    dataPts = annAllocPts(PointCount, 3); // allocate data points
    nnIdx = new ANNidx[KDTree3MaxK];  // allocate near neigh indices
    dists = new ANNdist[KDTree3MaxK]; // allocate near neighbor dists
    for(unsigned int i = 0; i < PointCount; i++)
    {
        for(unsigned int ElementIndex = 0; ElementIndex < 3; ElementIndex++)
        {
            dataPts[i][ElementIndex] = Points[i][ElementIndex];
        }
    }

    kdTree = new ANNkd_tree( // build search structure
        dataPts,    // the data points
        PointCount, // number of points
        3);         // dimension of space   
}

unsigned int KDTree3::Nearest(const Eigen::Vector3f &Pos)
{
    std::vector<unsigned int> Result(1);
    KNearest(Pos, 1, Result, 0.0f);
    return Result[0];
}

void KDTree3::KNearest(const Eigen::Vector3f &Pos, unsigned int k, std::vector<unsigned int> &Result, float Epsilon)
{
	assert(k <= KDTree3MaxK);//, "k too large");
    for(unsigned int ElementIndex = 0; ElementIndex < 3; ElementIndex++)
    {
        queryPt[ElementIndex] = Pos[ElementIndex];
    }
    kdTree->annkSearch( // search
        queryPt,        // query point
        k,                // number of near neighbors
        nnIdx,            // nearest neighbors (returned)
        dists,            // distance (returned)
        Epsilon);        // error bound

    if(Result.size() < k)
    {
        Result.resize(k);
    }
    for(unsigned int i = 0; i < k; i++)
    {
        Result[i] = nnIdx[i];
    }
}

void KDTree3::KNearest(const Eigen::Vector3f &Pos, unsigned int k, unsigned int *Result, float Epsilon)
{
	assert(k <= KDTree3MaxK);// , "k too large");
    for(unsigned int ElementIndex = 0; ElementIndex < 3; ElementIndex++)
    {
        queryPt[ElementIndex] = Pos[ElementIndex];
    }
    kdTree->annkSearch( // search
        queryPt,        // query point
        k,                // number of near neighbors
        nnIdx,            // nearest neighbors (returned)
        dists,            // distance (returned)
        Epsilon);        // error bound
    for(unsigned int i = 0; i < k; i++)
    {
        Result[i] = nnIdx[i];
    }
}

void KDTree3::WithinDistance(const Eigen::Vector3f &Pos, float Radius, std::vector<unsigned int> &Result) const
{
    for(unsigned int ElementIndex = 0; ElementIndex < 3; ElementIndex++)
    {
        queryPt[ElementIndex] = Pos[ElementIndex];
    }

    int NeighborCount = kdTree->annkFRSearch(
        queryPt,
        Radius * Radius,
        KDTree3MaxK,
        nnIdx,
        dists,
        0.0f);

	NeighborCount = (NeighborCount > KDTree3MaxK) ? KDTree3MaxK : NeighborCount;
	Result.resize(NeighborCount);
    for(unsigned int i = 0; int(i) < Result.size(); i++)
    {
        Result[i] = nnIdx[i];
    }
}
