/*
KDTree3.h
Written by Matthew Fisher

A 3D KD-tree that supports fast points-within-rectangle queries
*/
#pragma once

#include <ANN/ANN.h>
#include <Eigen/Core>
#include <vector>

const unsigned int KDTree3MaxK = 5000;

class KDTree3
{
public:
    KDTree3();
    ~KDTree3();
    void FreeMemory();

    void BuildTree(const std::vector<Eigen::Vector3f> &Points);
    void KNearest(const Eigen::Vector3f &Pos, unsigned int k, std::vector<unsigned int> &Result, float Epsilon);
    void KNearest(const Eigen::Vector3f &Pos, unsigned int k, unsigned int *Result, float Epsilon);
    void WithinDistance(const Eigen::Vector3f &Pos, float Radius, std::vector<unsigned int> &Result) const;
    unsigned int Nearest(const Eigen::Vector3f &Pos);
    __forceinline Eigen::Vector3f GetPoint(unsigned int Index)
    {
        Eigen::Vector3f Result;
        for(unsigned int ElementIndex = 0; ElementIndex < 3; ElementIndex++)
        {
            Result(ElementIndex) = float(dataPts[Index][ElementIndex]);
        }
        return Result;
    }

private:
    __forceinline Eigen::Vector3f GetDataPoint(unsigned int PointIndex)
    {
        Eigen::Vector3f Result;
        for(unsigned int ElementIndex = 0; ElementIndex < 3; ElementIndex++)
        {
            Result(ElementIndex) = float(dataPts[PointIndex][ElementIndex]);
        }
        return Result;
    }
    ANNpointArray    dataPts; // data points
    ANNpoint         queryPt; // query point
    ANNidxArray      nnIdx;   // near neighbor indices
    ANNdistArray     dists;   // near neighbor distances
    ANNkd_tree*      kdTree;  // search structure
    //Mutex          _Lock;
};
