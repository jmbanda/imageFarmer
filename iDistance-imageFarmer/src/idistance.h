/*////////////////////////////////////////////////////////////////////
idistance.h

Copyright (C) 2012  Michael Schuh, Timothy Wylie, Rafal Angryk
  Data Mining Laboratory at Montana State University 
  Principle author contact: michael.schuh@cs.montana.edu

This file is part of iDistance*.

iDistance* is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; only version 2
of the License. See the COPYING file for more information.
////////////////////////////////////////////////////////////////////*/

#ifndef _IDISTANCE_H_
#define _IDISTANCE_H_

#include "definitions.h"


////////////////////////////////////////////////////
//the main class implementing the idistance method
////////////////////////////////////////////////////
class iDistance
{
public:
    iDistance();
    iDistance(double* data, int n, int d);
    ~iDistance();
    
    //main query calls
    vector<int> QueryPoint(double q[]);
    vector<int> QueryRange(double low_range[], double high_range[]);
    vector<int> QuerySphere(double q[], double rad);
    vector<int> QueryKNN(double q[], int num);
    
protected:
    vector<int> Range_SIMPLE(double low_range[], double high_range[]);
    vector<int> Range_ADV1(double low_range[], double high_range[]);
    vector<int> KNN_KEEP(double q[], int num);
    vector<int> KNN_RESTART(double q[], int num);
    vector<int> Point_SEQSEARCH(double q[]);
    vector<int> Range_SEQSEARCH(double low_range[], double high_range[]);
    vector<int> KNN_SEQSEARCH(double q[], int num);
        
    
    //knn helper functions
    void SearchO_RESTART(double q[], double r, bool restrict_num);
    void SearchInward_RESTART(btree_mm::iterator node, double ivalue, double q[], int partition_i, bool restrict_num); //left
    void SearchOutward_RESTART(btree_mm::iterator node, double ivalue, double q[], int partition_i, bool restrict_num); //right
    void AddNeighbor_RESTART(int indexID, double distnode, bool sort, bool restrict_num);
    
    void SearchO_KEEP(double q[], double r, bool restrict_num);
    btree_mm::iterator SearchInward_KEEP(btree_mm::iterator node, double ivalue, double q[], int partition_i, bool restrict_num); //left
    btree_mm::iterator SearchOutward_KEEP(btree_mm::iterator node, double ivalue, double q[], int partition_i, bool restrict_num); //right
    void AddNeighbor_KEEP(int indexID, double distnode, bool sort);
    
    
public:
    //class helper functions
    void setData(double* data, int n, int d);
    void setDataWithIDs(double* data, int n, int d);
    void buildTree();
    void buildIndexValues();
    void setDataIDs(double* ids);
    void buildReferencePoints(int num);
    void buildReferencePoints();
    void buildSplits();
    
protected:    
    void buildRefs_HALF_POINTS();
    void buildRefs_HALF_POINTS_CENTER();
    void buildRefs_HALF_POINTS_OUTSIDE();
    void buildRefs_RANDOM();
    void buildRefs_CLUSTER();
    void buildRefs_UNIFORM();
    void buildIndex_CLOSEST();
    void buildIndex_ASSIGN();
    void buildSplits_TOP_CENTERS();
    
public:
    void setReferencePoints(double* refs, int num);
    void saveTree(const char* fStr);
    void loadTree(const char* fStr);
    void initOptions();
    void setSequentialSearch(int s);
    void getStats(stats* theStats);
    void outputReferencePoints();

    void setSplitPoints(double* splits, int num);
    void setPointAssignment(double* assigns, int num);
    
    //options
    void setOptions(configOptions co);
    vector<int> getIDs(vector<int> indices);
    vector<double> getKNN_dists();
    vector<int> getKNN_nodes();
    double getKNN_avgDist();
    int calculateC(int d, int e);
    
    double calcIndex(int p, double c, double d);
    vector<int> getKNN_candidates();
    
    
    
    //Stats stuff
    template<typename T>
    void printVector(const std::vector<T>& vec, ostream &out, bool commas);
    double AverageDist(vector<double>& distances);
    void PrintNodeRange(ostream &out);
    void PrintNodeRange(const char* fp);
    void setQueryTime(double qTime);
    stats theStats;
    void resetNodeCount();
    int getNodeCount();
    void resetCandidateCount();
    int getCandidateCount();
    void resetPartitionCount();
    int getPartitionCount();  
    int getNumberOfPoints();
    int getNumberOfDimensions();
    
protected:
    double dist(double p1[], double p2[]);
    double roundOff(double in);
    double* findMedians(double* data, int p, int d, bool hasID);
    
    //the b+ tree
    btree_mm btree;
    
    //the nxd array storing our datapoints
    double* datapoints;
    
    //the reference points for each partition
    double* reference_points;
    
    //the reference point to assign each data point
    //only used in certain configs
    double* point_assignments;
    
    //the max distance for each partition
    double* partition_dist_max;
    
    //the array index of the point that's the max dist for that partition
    unsigned int* partition_dist_max_index;
    
    //the ids for each data point
    double* datapoint_ids;   //INT?????????????????????????????????
    
    //the y index value for each data point
    double* datapoint_index;
    
    unsigned int number_partitions;
    unsigned int number_dimensions;
    unsigned int number_points;
    double constant_c;
    double radius_r;
    configOptions settings; //algorithm settings
    int doSeqSearch; // 0 = no, 1 = yes, 2 = verify
    int candidates;
    int partitionsChecked;
    
    //KNN Global Vars
    vector<int> knn_S;
    vector<double> knn_S_dists;
    vector<int> knn_S_nodes;
    vector<int> knn_S_candidates;
    
    unsigned int knn_K;
    bool knn_stopFlag;
    
    //Additional KNN_KEEP Global Vars
    double knn_pfarthest;
    int* knn_partition_checked;
    btree_mm::iterator* knn_left_pointers;
    btree_mm::iterator* knn_right_pointers;
};

#endif //_IDISTANCE_H_

