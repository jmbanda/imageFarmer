/*////////////////////////////////////////////////////////////////////
idistance.cpp

Copyright (C) 2012  Michael Schuh, Timothy Wylie, Rafal Angryk
  Data Mining Laboratory at Montana State University 
  Principle author contact: michael.schuh@cs.montana.edu

This file is part of iDistance*.

iDistance* is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; only version 2
of the License. See the COPYING file for more information.
////////////////////////////////////////////////////////////////////*/

#include <iostream>
#include <fstream>
#include <sys/time.h>
#include "string.h"
#include <math.h>
#include <numeric>

#include "idistance.h"


   
//*******************************************************************//
// Basic Constructor
//*******************************************************************//
iDistance::iDistance()
{
    datapoint_index = NULL;
    datapoint_ids = NULL;
    partition_dist_max_index = NULL;
    partition_dist_max = NULL;
    reference_points = NULL;
    constant_c = 1; // basic default, not guaranteed separation
    doSeqSearch = 0; // basic default, don't do sequential search
}

//*******************************************************************//
// Constructor passing in data
//*******************************************************************//
iDistance::iDistance(double* data, int n, int d)
{
    iDistance();
    setData(data,n,d);
}

//*******************************************************************//
//*******************************************************************//

iDistance::~iDistance()
{
    //delete allocated data
    if(datapoint_index != NULL)
    {
        delete [] datapoint_index;
        datapoint_index = NULL;
    }
    if(datapoint_ids != NULL)
    {
        delete [] datapoint_ids; //causes seg fault
        datapoint_ids = NULL;
    }
    if(partition_dist_max_index != NULL)
    {
        delete [] partition_dist_max_index;
        partition_dist_max_index = NULL;
    }
    if(partition_dist_max != NULL)
    {
        delete [] partition_dist_max;
        partition_dist_max = NULL;
    }
    if(reference_points != NULL)
    {
        delete [] reference_points;
        reference_points = NULL;
    }
    if(datapoints != NULL)
    {
        delete [] datapoints; //causes seg fault
        datapoints = NULL;
    }
    //erase tree
    btree.clear();


}


//*******************************************************************//
//Returns a vector of indices of all the points matching q
//*******************************************************************//
vector<int> iDistance::QueryPoint(double q[])
{
    /*
    cout << "query point = ";
    
    for(int i = 0; i < number_dimensions; i++)
    {
        cout << q[i] << ", ";
    }
    cout << endl;
    */
    
    if(doSeqSearch == 1) //do sequential search instead!
    {
        return Point_SEQSEARCH(q); // return list of final set
    }

    vector<int> results;
    btree_mm::iterator bi;
    double ref_dist;
    double y;
    int y_rnd;
    //cout<<"#Partitions= "<<number_partitions;
    //find which partition(s) it's in
    //stx::CountQueryNodes("Start"); //reset query node count
    for(int i=0 ; i < number_partitions ; i++)
    {
//        cout<<"Check Partition "<<i<<endl;
        ref_dist = dist(q,&reference_points[i*number_dimensions]);
        //cout << "  DISTMAX - REFDIST = " << partition_dist_max[i] - ref_dist << endl;
        if(ref_dist <= partition_dist_max[i]) //querydist = r = 0
        {
            partitionsChecked++;
            //cout<<": True"<<endl;
            //if in this partition, calculate key index value
            y = roundOff(i*constant_c + ref_dist);
            
            bi = btree.find( y ); 
            
//            cout << "BI | Y = " << bi.key() << " | " << y << endl;
            while(bi.key() == y) //for every match of y
            {
                candidates++; //one more possible candidate!
                double d = dist(&datapoints[bi.data()*number_dimensions],q);
                //cout << "QP DIST = " << d << endl;
               //check the candidate to see if actual point in full dimensions
                if(d == 0.0)
                {
                    //pushing back our internal id found in b+ tree
                    results.push_back(bi.data());
                }
                bi++;
            }
        }
    }
    
    
    return results;
}


vector<int> iDistance::Point_SEQSEARCH(double q[])
{
    vector<int> results; //final result set (internal IDs)
    
    bool candidate = true; //point is still a possibility
    
    for(int i = 0; i < number_points; i++) //for each point
    {
        candidate = true;
        for(int j = 0; j < number_dimensions; j++) //for each dim
        {
            if(datapoints[i*number_dimensions + j] != q[j])
            {
                //point does not equal query point in this dim
                candidate = false;
                j = number_dimensions; // no need to check the other dims
            }
        }
        
        if(candidate) //add points to results that are still candidates 
        {
            results.push_back(i);
        }
        
    }
    
    return results;
}


//*******************************************************************//
// Master range search function which delegates the range search to the
// method specified in the current configuration options.
//*******************************************************************//
vector<int> iDistance::QueryRange(double low_range[], double high_range[])
{
    
    if(doSeqSearch == 1) //do sequential search instead!
    {
        return Range_SEQSEARCH(low_range, high_range); // return list of final set
    }
        
    vector<int> results;
    candidates = 0; //reset candidates counter
    
    switch(settings.range_method)
    {
      case SIMPLE:
//        cout << "SIMPLE!" << endl;
        results = Range_SIMPLE(low_range, high_range);
        break;
      case ADV1:
        //do this
        break;
    }
 
    return results;    
  
}

vector<int> iDistance::Range_SEQSEARCH(double low_range[], double high_range[])
{
    vector<int> results; //final result set (internal IDs)
    
    bool candidate = true; //point is still a possibility
    
    for(int i = 0; i < number_points; i++) //for each point
    {
        candidate = true;
        for(int j = 0; j < number_dimensions; j++) //for each dim
        {
            
            if( (datapoints[i*number_dimensions + j] < low_range[j] ) || 
                (datapoints[i*number_dimensions + j] > high_range[j] ) )
            {
                //point is not inside query range in this dim
                candidate = false;
                j = number_dimensions; // no need to check the other dims
            }
        }
        
        if(candidate) //add points to results that are still candidates 
        {
            results.push_back(i);
        }
        
    }
    
    return results;  
}


/*********************************************************
 * Simple range search function which uses an enveloping spherical query
 * to retrieve candidates for the query hypercube. All candidates must then
 * be brute forced checked in each dimension for hypercube containment.
 ********************************/ 
vector<int> iDistance::Range_SIMPLE(double low_range[], double high_range[])
{
    //set radius of search sphere to reach the corners of hypercube
    double r = dist(low_range, high_range)/2.0;
    cout<<"Range radius: "<<r<<endl;
    
    //allocate array for center pt
    double* center_pt = new double[number_dimensions];
    
    cout << "center = " << endl;
    //calculate center point of hypercube
    for(int dim=0; dim < number_dimensions; dim++)
    {
        //average of each dimension
        center_pt[dim] = (high_range[dim] + low_range[dim])/2.0;
        
        cout << center_pt[dim] << ", ";
    }
    
    cout << endl;
    
    //do a radius-based search for candidates (filter)
    vector<int> possible_pts = QuerySphere(center_pt, r);
    
    cout << "possible pts = " << possible_pts.size() << endl;
    
    candidates = possible_pts.size();
    
    //now check each point (refine)
    int index;
    for(int point = possible_pts.size()-1; point >= 0; point--)
    {
        //cout<<"ID: "<<datapoint_ids[possible_pts[point]]<<endl;
        
        //get the point's ID into raw data array -- assumes the internal 
        //data ids which are sequential
        
        //THIS SHOUDL BE FIXED 
        index = possible_pts[point] * number_dimensions;
        for(int dim=0; dim < number_dimensions; dim++) //check each dim
        {
            //cout<<datapoints[index + dim]<<endl;
            //if the point is outside the range
            if(datapoints[index + dim] < low_range[dim] 
                || datapoints[index + dim] > high_range[dim])
            {
                //delete element and exit for:dims loop
                possible_pts.erase(possible_pts.begin()+point);
                break;
            }
        }
    }

    return possible_pts;
  
}


vector<int> iDistance::Range_ADV1(double low_range[], double high_range[])
{
    vector<int> results;
    return results;
}

/********************************************************************
 * KNN Master function to delegate control based on the type of KNN
 * the config file stated.
 * 
 *********************************************************************/ 
vector<int> iDistance::QueryKNN(double q[], int num)
{
    /*
    cout << "query point = ";
    
    for(int i = 0; i < number_dimensions; i++)
    {
        cout << q[i] << ", ";
    }
    cout << endl;
    */
    
    if(doSeqSearch == 1) //do sequential search instead!
    {
        return KNN_SEQSEARCH(q, num); // return list of final set
    }
    
    // else proceed as normal
    
    vector<int> results;

    switch(settings.knn_method)
    {
        case KEEP:
        {
            results = KNN_KEEP(q, num);
        }
        break;
        case RESTART:
        {
            results = KNN_RESTART(q, num);
        }
        break;
    }

    return results;
    
}

/****************************************************************
 * KNN function used for range queries. Returns ALL points found
 * within a given radius. Uses the RESTART methods for simplicity
 * but there is no incremental searching.
 * 
 ***************************************************************/ 
vector<int> iDistance::QuerySphere(double q[], double rad)
{
    knn_K = number_points;
    knn_stopFlag = false;
    knn_S.clear();
    knn_S_dists.clear();
    knn_S_nodes.clear();
    knn_S_candidates.clear();
    radius_r = rad;
    
    SearchO_RESTART(q, rad, false);
    /*
    switch(settings.knn_method)
    {
        case KEEP:
        {
            knn_pfarthest = 0;
            knn_partition_checked = new int[number_partitions];
            for(int i=0;i<number_partitions;i++)
                knn_partition_checked[i] = 0;
            //array of iterators
            knn_left_pointers = new btree_mm::iterator[number_partitions];
            knn_right_pointers = new btree_mm::iterator[number_partitions];
    
            SearchO_KEEP(q, rad, false);
            
            delete [] knn_partition_checked;
            delete [] knn_left_pointers;
            delete [] knn_right_pointers;
        }
        break;
        case RESTART:
        {
            SearchO_RESTART(q, rad, false);
        }
        break;
    }
    */ 
    return knn_S;
}


//////////////////////////////////////////////////////////////////
//////////// specific KNN functions  ////////////////////////////
/////////////////////////////////////////////////////////////////
//*******************************************************************//
//************cleaned up version of knn
//*******************************************************************//
vector<int> iDistance::KNN_RESTART(double q[], int num)
{
    radius_r = settings.r_init;
    knn_stopFlag = false;
    knn_K = num;
    knn_S.clear();
    knn_S_dists.clear();
    knn_S_nodes.clear();
    knn_S_candidates.clear();
    
    //knn_pfarthest = 0;
        
    stx::CountQueryNodes("Start"); //reset query node count
    while(!knn_stopFlag)  //knn_ret.size() < num)//should probably have safety net (suppose num = 4, but only 3 points exist)
    {
        //increase radius (this could be a class parameter)
        radius_r += settings.r_delta;
        //check cross with each partition

        //cout << "searching radius : " << radius_r << endl;
        SearchO_RESTART(q, radius_r, true);
        
        //cout<<"********************************************"<<endl;
        //cout<<"With radius: " << radius_r << " KNN="<<knn_S.size()<<endl;
        if(radius_r > constant_c)   //assuming that constant_c is greater than max distance of space
            knn_stopFlag = true;
    }

    return knn_S;   
}

/********************************************************
 * Searches through all partitions for points within radius r of a 
 * query point q. If for true KNN, we are restricted on the number of
 * data points being retrieved.
 **************************/ 
void iDistance::SearchO_RESTART(double q[], double r, bool restrict_num)
{
    ///////////////////set stop flag////////////////////////
    if(knn_S.size() == knn_K)  //also changed if our radius gets too big
        knn_stopFlag = true;
    
    //variable to store the distance from partition to query pt
    double distp;
    //the index value of the query point
    double q_index;
    
    partitionsChecked = 0;
    //check for each partition
    for(int i=0; i < number_partitions; i++)
    {   
        //cout << "PARTITION: " << i << endl;
        
        //calc distance from query point to partition i
        distp = dist(&reference_points[i*number_dimensions],q);
        
        //calculate the q_index
        q_index = i*constant_c + distp;
        
        //check if query is interacting with partition (filter):
        //dist(O_i, q) - querydist(q)="r" <= dist_max_i
        if(distp - r <= partition_dist_max[i])
        {
            partitionsChecked++;
            //if query sphere is inside this partition
            if(distp <= partition_dist_max[i])
            {   
                //cout<<"Inside Partition "<<i<<endl;
                //find query pt and search in/out
                btree_mm::iterator lnode = btree.upper_bound(q_index);  
                btree_mm::iterator rnode;// = btree.lower_bound(q_index);
                rnode = lnode;
                SearchInward_RESTART(lnode, q_index - r, q, i, restrict_num);
                SearchOutward_RESTART(rnode--, q_index + r, q, i, restrict_num);
            }
            else //it intersects
            {
                //cout<<"Intersect Partition "<<i<<endl;
                //get the index value(y) of the pt of dist max
                btree_mm::iterator dist_max_i = btree.find(datapoint_index[partition_dist_max_index[i]]);
                SearchInward_RESTART(dist_max_i, q_index - r, q, i, restrict_num);
            }
        }
    }
}

/***************************************************
 * Search the B+ tree from the given node pointer inward (lower values)
 * until reaching the stopping value (ivalue)
 ***************************************/ 
void iDistance::SearchInward_RESTART(btree_mm::iterator node, double ivalue, double q[], int partition_i, bool restrict_num)
{ 
    //variable for distance
    double distnode=constant_c;
    double partition_offset = roundOff(partition_i*constant_c);
    int past_data = -1;
    //i think these will be faster
    double node_key = node.key();
    int node_data = node.data();
    //cout<<"--SearchInward--"<<endl;
    //cout<<"Part_Offset: "<<partition_offset<<" iValue: "<<ivalue<<endl;
    //must get each node of that value  
    //while not to stopping value and still inside partition
    while(node_key >= ivalue && node_key >= partition_offset && node_data != past_data)
    {
        candidates++;
        distnode = dist(&datapoints[node_data*number_dimensions], q);
        //cout<<"INWARD KEY: "<<node.key()<<" DATA: "<<node.data()<<" DIST: "<<distnode<<endl;
        if(distnode <= radius_r)
        {
            AddNeighbor_RESTART(node_data, distnode, restrict_num, restrict_num);    //i'm assuming you only sort when you need to
        }
        past_data = node_data;
        node--;
        node_key = node.key();
        node_data = node.data();
    }
}

/***************************************************
 * Search the B+ tree from the given node pointer outward (higher values)
 * until reaching the stopping value (ivalue)
 ***************************************/ 
void iDistance::SearchOutward_RESTART(btree_mm::iterator node, double ivalue, double q[], int partition_i, bool restrict_num)
{
    //variable for distance
    double distnode=constant_c;
    double dist_max = roundOff(partition_i*constant_c + partition_dist_max[partition_i]);
    int past_data = -1;
    //i think these will be faster
    double node_key = node.key();
    int node_data = node.data();
    //cout<<"--SearchOutward--"<<endl;
    //cout<<"Dist_max: "<<dist_max<<" iValue: "<<ivalue<<endl;
    //must get each node of that value    
    //while not to stopping value and still inside partition
    while(node_key <= ivalue && node_key <= dist_max && node_data != past_data)
    {
        candidates++;
        distnode = dist(&datapoints[node_data*number_dimensions], q);
        //cout<<"OUTWARD KEY: "<<node.key()<<" DATA: "<<node.data()<<" DIST: "<<distnode<<endl;
        if(distnode <= radius_r)
        {
            AddNeighbor_RESTART(node_data, distnode, restrict_num, restrict_num);    //i'm assuming you only sort when you need to
        }
        past_data = node_data;
        node++;
        node_key = node.key();
        node_data = node.data();
    }
}

//implicit index
//the distance from q
//if the data needs to be sorted
//if we're restricted to a certain num

//you can sort and not restrict, but you can't restrict without sort
void iDistance::AddNeighbor_RESTART(int indexID, double distnode, bool sort, bool restrict_num)
{
    //cout<<"NEIGHBOR: "<<indexID<<" Dist: "<<distnode<<endl;
    int k;
    
    //first we check that the data point has not already been added
    bool okay=true;
    for(int i=0;i<knn_S.size();i++)
    {
        if(knn_S.at(i) == indexID)
            okay=false;
    }
    if(okay) //if it hasn't yet...
    {
        if(!knn_S.empty() && sort) //restrict num?
        {   
            //if it needs to be inserted before last element
            if(distnode < knn_S_dists.back())
            {
                //finding the right insert spot
                for(k=knn_S_dists.size()-1; distnode < knn_S_dists.at(k) && k>0; k--)
                {}
                //if k=0 we need to know if we insert at 0 or 1
                if(k==0 && distnode < knn_S_dists.at(0))
                {
                    knn_S.insert(knn_S.begin(), indexID);
                    knn_S_dists.insert(knn_S_dists.begin(), distnode);
                    knn_S_nodes.insert(knn_S_nodes.begin(), stx::CountQueryNodes("Get"));
                    knn_S_candidates.insert(knn_S_candidates.begin(), candidates);
                }
                else
                {
                    //insert at one greater than k
                    knn_S.insert(knn_S.begin() + k + 1, indexID);
                    knn_S_dists.insert(knn_S_dists.begin() + k + 1, distnode);
                    knn_S_nodes.insert(knn_S_nodes.begin() + k + 1, stx::CountQueryNodes("Get"));       
                    knn_S_candidates.insert(knn_S_candidates.begin() + k + 1 , candidates);
                }
                //if we're restricting the number of elements
                if(restrict_num)
                {
                    //if we now have num+1 neighbors
                    if(knn_S.size() > knn_K)
                    {
                        knn_S.pop_back();
                        knn_S_dists.pop_back();
                        knn_S_nodes.pop_back();
                        knn_S_candidates.pop_back();
                    }
                }
            }
            else
            {
                //if we're restricting the number of elements
                if(restrict_num)
                {   //if we're still less than how many neighbors we want,
                    //we can just add it to the end
                    if(knn_S.size() < knn_K)  //class variables
                    {
                        //add to end of list
                        knn_S.push_back(indexID);
                        knn_S_dists.push_back(distnode);
                        knn_S_nodes.push_back(stx::CountQueryNodes("Get"));
                        knn_S_candidates.push_back(candidates);
                    }
                }
                else
                {
                    //add to end of list
                    knn_S.push_back(indexID);
                    knn_S_dists.push_back(distnode);
                    knn_S_nodes.push_back(stx::CountQueryNodes("Get"));
                    knn_S_candidates.push_back(candidates);
                }
            }
        }
        else
        {   //add to end of list
            knn_S.push_back(indexID);
            knn_S_dists.push_back(distnode);
            knn_S_nodes.push_back(stx::CountQueryNodes("Get"));
            knn_S_candidates.push_back(candidates);
        }
    }
}


vector<int> iDistance::KNN_SEQSEARCH(double q[], int num)
{
    knn_K = num;
    knn_S.clear();
    knn_S_dists.clear();
    knn_S_nodes.clear();
    knn_S_candidates.clear();
    
    //vector<int> results; //final result set (internal IDs)
    //vector<double> distances; // distances vector, same order as above
    
    if(num < 1)
    {
        return knn_S;
    }
    
    for(int i = 0; i < num; i++) //just fill the node counter with max
    {
        knn_S_nodes.push_back(theStats.treeNodes);
        knn_S_candidates.push_back(theStats.treeNodes);
    }
    
    int k = 0; // total k out of num saved
    double maxdist = 0.0; // max distance out of the k
    double tempdist;
    
    //prime the loop
    maxdist = dist(&datapoints[0], q);
    knn_S.push_back(0);
    knn_S_dists.push_back(maxdist);
    k = 1;
    
    //cout << "primed loop with 0 : " << maxdist << endl;
    
    //now continue on without base case checks
    for(int i = 1; i < number_points; i++) //for each point
    {
        tempdist = dist(&datapoints[(i*number_dimensions)],q);
        //cout << "checking " << i << " : " << tempdist << " vs. " << maxdist << endl;
        if(tempdist < maxdist || k < num)
        {
            //cout << "\t keep " << i << endl;
            //we want to keep it, but need to find where it goes in sorted list
            
            int j = 0;
            vector<int>::iterator r = knn_S.begin();
            vector<double>::iterator d = knn_S_dists.begin();
            while(j < k && tempdist >= knn_S_dists[j])
            {
                d++;
                r++;
                j++;
            }
            
            //cout << "inserting at pos: " << j << endl;
            
            if(j >= k) // end of list
            {
                //cout << "j == k!!!" << endl;
                //insert at position j, and push the rest back
                knn_S_dists.push_back(tempdist);
                knn_S.push_back(i);
                maxdist = tempdist;
                //cout << "inserted at end!" << endl;
            }
            else
            {
                knn_S_dists.insert(d, tempdist);
                knn_S.insert(r, i);
                //cout << "inserted!" << endl;
            }
            
            if(k == num)
            {
                //cout << "k = num!" << endl;
                //already full, bumped one
                //drop last one
                r = knn_S.end();
                r--;
                knn_S.erase(r);
                maxdist = knn_S_dists[k-1];
                d = knn_S_dists.end();
                d--;
                knn_S_dists.erase(d);
            }
            else
            {
                //cout << "inc k " << endl;
                k++;
            }
            
            
        }
        //cout << "so far: " << k << " of " << num << endl;
        
        //for(int z = 0; z < k; z++)
        //{
            //cout << knn_S[z] << ":" << knn_S_dists[z] << ", ";
        //}
        //cout << endl;
            
    }
    
    return knn_S;  
}

//*******************************************************************//
//does search with radius r from pt q and returns all results
//*******************************************************************//
vector<int> iDistance::KNN_KEEP(double q[], int num)
{
    radius_r = settings.r_init;
    knn_stopFlag = false;
    knn_K = num;
    knn_S.clear();
    knn_S_dists.clear();
    knn_S_nodes.clear();
    knn_S_candidates.clear();
    
    //variable to mark for partitions checked
    knn_partition_checked = new int[number_partitions];
    for(int i=0;i<number_partitions;i++)
        knn_partition_checked[i] = 0;
    knn_pfarthest = 0;
    //array of iterators
    knn_left_pointers = new btree_mm::iterator[number_partitions];
    knn_right_pointers = new btree_mm::iterator[number_partitions];
    
    //stx::CountQueryNodes("Start"); //reset query node count
    while(!knn_stopFlag)  //knn_ret.size() < num)//should probably have safety net (suppose num = 4, but only 3 points exist)
    {
        //increase radius (this could be a class parameter)
        radius_r += settings.r_delta;
        //check cross with each partition

        SearchO_KEEP(q, radius_r, true);
        //cout<<"********************************************"<<endl;
        //cout<<"With radius: " << radius_r << " KNN="<<knn_S.size()<<endl;
        if(radius_r > constant_c)   //assuming that constant_c is greater than max distance of space
            knn_stopFlag = true;

    }

    //look at partitions checked before finishing up
    partitionsChecked = 0;
    for(int i=0; i < number_partitions; i++)
    {
        if(knn_partition_checked[i] == 1)
        {
          partitionsChecked++;
        }
    }


    delete [] knn_partition_checked;
    delete [] knn_left_pointers;
    delete [] knn_right_pointers;
    //theStats.ind = knn_S;
    //theStats.knn = true;
    
    return knn_S;
}


//searches partitions
void iDistance::SearchO_KEEP(double q[], double r, bool sort)
{
    //pfarthest in code is checking for the farthest neighbor in the answer
    //set S.  If |S| == #neighbors wanted and pfarthest is < r from q we
    //have all of our values
    ///////////////////set stop flag////////////////////////
    if(knn_S.size() == knn_K && knn_pfarthest <= r)
        knn_stopFlag = true;
    
    //variable to store the distance from partition to query pt
    double distp;
    //the index value of the query point
    double q_index;
    //check for each partition
    for(int i=0; i < number_partitions; i++)
    {   
        //cout<<"PARTITION: "<<i<<endl;
        
        //calc distance from query point to partition i
        distp = dist(&reference_points[i*number_dimensions],q);
        //calculate the q_index - can be moved around
        q_index = roundOff(i*constant_c + distp);
        //cout<<"Q: "<<q_index<<endl;                             
        //have we checked this partition before?
        if(knn_partition_checked[i] != 1) 
        {
            //filter dist(O_i, q) - querydist(q)="r" <= dist_max_i
            if(distp - r <= partition_dist_max[i])
            {
                knn_partition_checked[i] = 1;  //mark it as checked
                //if query sphere is inside this partition
                if(distp <= partition_dist_max[i])
                {   
                    //cout<<"Inside Partition "<<i<< ", q_index = " << q_index << endl;
                    
                    //find query pt and search in/out
                    btree_mm::iterator lnode = btree.upper_bound(q_index);
                    btree_mm::iterator rnode;// = btree.lower_bound(q_index);
                    rnode = lnode;
                    knn_left_pointers[i] = SearchInward_KEEP(lnode, q_index - r, q, i, sort);
                    knn_right_pointers[i] = SearchOutward_KEEP(rnode--, q_index + r, q, i, sort);
                }
                else //it intersects
                {
                    //cout<<"Intersect Partition "<<i<<", q_index = " << q_index <<endl;
                    //get the index value(y) of the pt of dist max
                    btree_mm::iterator dist_max_i = btree.find(datapoint_index[partition_dist_max_index[i]]);
                    knn_left_pointers[i] = SearchInward_KEEP(dist_max_i, q_index - r, q, i, sort);
                }
            }
        }
        else //we've checked it before
        {
            btree_mm::iterator null_iterator;
            if(knn_left_pointers[i] != null_iterator) //can't actually check if it's null
            {
                knn_left_pointers[i] = SearchInward_KEEP(knn_left_pointers[i]--, q_index - r, q, i, sort);
            }
            if(knn_right_pointers[i] != null_iterator)
            {
                knn_right_pointers[i] = SearchOutward_KEEP(knn_right_pointers[i]++, q_index + r, q, i, sort);
            }
        }
    }
}
//node is a pointer to the node to start from
//ivalue is the stopping value
btree_mm::iterator iDistance::SearchInward_KEEP(btree_mm::iterator node, double ivalue, double q[], int partition_i, bool sort)
{ 
    //variable for distance
    double distnode=constant_c;
    double partition_offset = roundOff(partition_i*constant_c);
    int past_data = -1;
    //i think these will be faster
    double node_key = node.key();
    int node_data = node.data();
    //cout<<"--SearchInward--"<<endl;
    //cout<<"Part_Offset: "<<partition_offset<<" iValue: "<<ivalue<<endl;
    //cout<<"Init Node Inward: "<<node_data<<" value: "<<node_key<<endl;

    //must get each node of that value  
    //while not to stopping value and still inside partition
    
    while(node_key >= ivalue && node_key >= partition_offset && node_data != past_data)
    {        
        //candidates++;
        distnode = dist(&datapoints[node_data*number_dimensions], q);

        //cout<<"INWARD KEY: "<<node.key()<<" DATA: "<<node.data()<<" DIST: "<<distnode<<endl;
        AddNeighbor_KEEP(node_data, distnode, sort);
        past_data = node_data;
        node--;
        node_key = node.key();
        node_data = node.data();
    }
    //past partition, set node to null pointer  
    if(node_key < partition_offset || node_data == past_data) 
    {
        btree_mm::iterator node1;
        return node1;
    }
    return node;    //maybe i don't need to ++
}

btree_mm::iterator iDistance::SearchOutward_KEEP(btree_mm::iterator node, double ivalue, double q[], int partition_i, bool sort)
{
    //variable for distance
    double distnode=constant_c;
    double dist_max = roundOff(partition_i*constant_c + partition_dist_max[partition_i]);
    int past_data = -1;
    //i think these will be faster
    double node_key = node.key();
    int node_data = node.data();
    //cout<<"--SearchOutward--"<<endl;
    //cout<<"Dist_max: "<<dist_max<<" iValue: "<<ivalue<<endl;
    //cout<<"Init Node Outward: "<<node_key<<endl;

    //must get each node of that value    
    //while not to stopping value and still inside partition
    while(node_key <= ivalue && node_key <= dist_max && node_data != past_data)
    {        
        //candidates++;
        distnode = dist(&datapoints[node_data*number_dimensions], q);
        //cout<<"OUTWARD KEY: "<<node.key()<<" DATA: "<<node.data()<<" DIST: "<<distnode<<endl;
        AddNeighbor_KEEP(node_data, distnode, sort);

        past_data = node_data;
        node++;
        node_key = node.key();
        node_data = node.data();
    }
    //past partition, set node to null pointer  
    if(node_key > dist_max || node_data == past_data) 
    {
        btree_mm::iterator node1;
        return node1;
    }
    return node;    //maybe i don't need to --
}


//you can sort and not restrict, but you can't restrict without sort
void iDistance::AddNeighbor_KEEP(int indexID, double distnode, bool sort)
{
    //cout<<"NEIGHBOR: "<<indexID<<" Dist: "<<distnode<<endl;
    int k;
    bool okay=true;
    for(int i=0;i<knn_S.size();i++)
    {
        if(knn_S.at(i) == indexID)
            okay=false;
    }
    if(okay) //if it hasn't yet...
    {
        candidates++;
        if(!knn_S.empty() && sort) //restrict num?
        {   
            //if it needs to be inserted before last element
            if(distnode < knn_S_dists.back())
            {
                //finding the right insert spot
                for(k=knn_S_dists.size()-1; distnode < knn_S_dists.at(k) && k>0; k--)
                {}
                //if k=0 we need to know if we insert at 0 or 1
                if(k==0 && distnode < knn_S_dists.at(0))
                {
                    knn_S.insert(knn_S.begin(), indexID);
                    knn_S_dists.insert(knn_S_dists.begin(), distnode);
                    knn_S_nodes.insert(knn_S_nodes.begin(), stx::CountQueryNodes("Get"));
                    knn_S_candidates.insert(knn_S_candidates.begin(), getCandidateCount());
                }
                else
                {
                    //insert at one greater than k
                    knn_S.insert(knn_S.begin() + k + 1, indexID);
                    knn_S_dists.insert(knn_S_dists.begin() + k + 1, distnode);
                    knn_S_nodes.insert(knn_S_nodes.begin() + k + 1, stx::CountQueryNodes("Get"));               
                    knn_S_candidates.insert(knn_S_candidates.begin() + k + 1, getCandidateCount());
                }
            }
            else
            {
                //if we're still less than how many neighbors we want,
                //we can just add it to the end
                if(knn_S.size() < knn_K)  //class variables
                {
                    //add to end of list
                    knn_S.push_back(indexID);
                    knn_S_dists.push_back(distnode);
                    knn_S_nodes.push_back(stx::CountQueryNodes("Get"));
                    knn_S_candidates.push_back(getCandidateCount());
                }
                else
                {
                    //add to end of list
                    knn_S.push_back(indexID);
                    knn_S_dists.push_back(distnode);
                    knn_S_nodes.push_back(stx::CountQueryNodes("Get"));
                    knn_S_candidates.push_back(getCandidateCount());
                }
            }
            
            
            if(knn_S.size() > knn_K)  //could also change to while
            {   
                knn_S.pop_back();
                knn_S_dists.pop_back();
                knn_S_nodes.pop_back();
                knn_S_candidates.pop_back();
                knn_pfarthest = knn_S_dists[knn_K-1];
            }
            else
            {
                knn_pfarthest = knn_S_dists.back();
            }
        }
        else
        {   
            if(distnode <= radius_r || sort) //or sort???
            {
                //add to end of list
                knn_S.push_back(indexID);
                knn_S_dists.push_back(distnode);
                knn_S_nodes.push_back(stx::CountQueryNodes("Get"));
                knn_S_candidates.push_back(getCandidateCount());
                
          
                if(distnode > knn_pfarthest)
                    knn_pfarthest = distnode;
            }
        }
    }
}


//*******************************************************************//
//*******************************************************************//
void iDistance::setData(double* data, int n, int d)
{
    datapoints = data;
    number_points = n;
    number_dimensions = d;
}
//*******************************************************************//
//This function takes in a pointer to a data array and the number
//of points and dimensions contained in the array. It sets this
//information globally and splits the data points from their IDs into
//two class global arrays.
//*******************************************************************//
void iDistance::setDataWithIDs(double* data, int n, int d)
{
    number_points = n;
    number_dimensions = d;
    datapoint_ids = new double[number_points];
    datapoints = new double[number_points*number_dimensions];
    
    for(int i=0; i < number_points; i++) //for each data point
    {
        //dims+1 because the ID is still in the 1D array
        for(int j=0; j < number_dimensions+1; j++)
        {
            if(j==0) //looking at first data value, its the ID
            {
                datapoint_ids[i] = data[i*(number_dimensions+1)];
            }
            else
            {
                datapoints[i*number_dimensions + (j-1)] = data[i*(number_dimensions+1) + j];
            }
        }
    }
}

/******************************
 * These options have already been established by the config file
 * but some require knowing information about the data. this function
 * is called after data is set.
 ************************************/ 
void iDistance::initOptions()
{
  
  if(settings.c_type == CALCULATE)
  {
    constant_c = calculateC(number_dimensions, settings.refs_dist);
  }
  else if(settings.c_type == MANUAL)
  {
    constant_c = settings.c_val;
  }
  
  if(settings.refs_build != REF_FILE)
  {
    buildReferencePoints();
  }
  
  
  if(settings.algo_version == EXT1)
  {
	  if(settings.splits_build != SPLIT_FILE)
	  {
		buildSplits();
	  }
  }
  
  //cout << "constant_c = " << constant_c << endl;
  //outputReferencePoints();
  
}

void iDistance::setSequentialSearch(int s)
{
    doSeqSearch = s;
}

int iDistance::getNumberOfPoints()
{
    return number_points;
}

int iDistance::getNumberOfDimensions()
{
    return number_dimensions;
}

//*******************************************************************//
//This function takes the list of reference points and uses them to
//to assign a partition to all of the data points.  Then it calculates
//the index value and stores it in an array.
//*******************************************************************//
void iDistance::buildIndexValues()
{
  
    switch(settings.refs_assign)
    {
      case CLOSEST:
        buildIndex_CLOSEST();
        break;
      case ASSIGN_FILE:
        buildIndex_ASSIGN();
        break;
    }
    
    
}

void iDistance::buildIndex_CLOSEST()
{
    
    //initialize the arrays
    datapoint_index = new double[number_points];
    partition_dist_max = new double[number_partitions];
    partition_dist_max_index = new unsigned int[number_partitions];
    for (int k=0; k < number_partitions; k++) {
        partition_dist_max[k] = 0.0;    // Initialize all elements to zero.
        partition_dist_max_index[k] = -1; // if -1, partition contains NO points!
    }
    
    
    
    double ref_dist, temp_dist;
    int current_partition=0;
    int offset;
    int temp_dist_rnd;
    
    for(int i=0 ; i < number_points ; i++) //for each point
    {
        
        //cout << "Point " << i << ": " << endl;
        
        //use idp for inspecting data points
        //bool idp = false;
        //if(datapoint_ids[i] == 700 || datapoint_ids[i] == 872 || datapoint_ids[i] == 927) { idp = true; }
        
        offset = i*number_dimensions;   //loop invariant
        
        current_partition = 0;
        ref_dist = dist(&datapoints[offset],&reference_points[0]);
        
        //if(idp) { cout << "checking point " << i << endl; }
        
        //if(idp) { cout << "p 0: " << ref_dist << endl; }
         
        //check each partition (baseline from partition 0)
        for(int part=1; part < number_partitions; part++)
        {
            //get distance to this partition
            temp_dist = dist(&datapoints[offset],&reference_points[part*number_dimensions]);
        
            //if(idp) { cout << "p " << part << ": " << temp_dist << endl; }
            
            if(temp_dist < ref_dist) //if closer to this partition
            {
                //if(idp) { cout << "  - updating from " << current_partition << endl; }
                
                current_partition = part; //update
                ref_dist = temp_dist;
                
            }
        }
        
        
        //if(idp) { cout << "final part = " << current_partition << " : " << ref_dist << endl; }
        
        //calculate key index value
        //datapoint_index[i] = current_partition*constant_c + ref_dist;
        
        //updated to try to fix double precision issues
        temp_dist = roundOff(current_partition*constant_c + ref_dist);
        
        //if(idp) { cout << "  index_dist = " << temp_dist << endl; }
        
        datapoint_index[i] = temp_dist;
        
        //cout << "  DATA INDEX: " << datapoint_index[i] << " P: " <<
        //  current_partition << " D: " << ref_dist << endl;
        
        //also need to update partition dist_max
        if(ref_dist > partition_dist_max[current_partition])
        {
            partition_dist_max[current_partition] = ref_dist;
            partition_dist_max_index[current_partition] = i;
        }
    }
}


void iDistance::buildIndex_ASSIGN()
{
    //given that we already have the assignments loaded up from file
    if(point_assignments == NULL)
    {
      cout << "No assignment file loaded up!" << endl;
      return;
    }
    
    //initialize the arrays
    datapoint_index = new double[number_points];
    partition_dist_max = new double[number_partitions];
    partition_dist_max_index = new unsigned int[number_partitions];
    for (int k=0; k < number_partitions; k++) {
        partition_dist_max[k] = 0.0;    // Initialize all elements to zero.
        partition_dist_max_index[k] = -1; // if -1, partition contains NO points!
    }
    
    
    double ref_dist, temp_dist;
    int part, offset;
    
    for(int i=0 ; i < number_points ; i++) //for each point
    {
        
        cout << "  Point " << i << " assigned to: " << point_assignments[i] << endl;
        
        offset = i*number_dimensions;   //loop invariant
        
        part = point_assignments[i];
        ref_dist = dist(&datapoints[offset],&reference_points[part*number_dimensions]);

        //updated to try to fix double precision issues
        temp_dist = roundOff(part*constant_c + ref_dist);
        
        //if(idp) { cout << "  index_dist = " << temp_dist << endl; }
        
        datapoint_index[i] = temp_dist;
        
        //cout << "  DATA INDEX: " << datapoint_index[i] << " P: " <<
        //  current_partition << " D: " << ref_dist << endl;
        
        //also need to update partition dist_max
        if(ref_dist > partition_dist_max[part])
        {
            partition_dist_max[part] = ref_dist;
            partition_dist_max_index[part] = i;
        }
        
        
    }
}


//*******************************************************************//
//*******************************************************************//
void iDistance::setDataIDs(double* ids)
{
    datapoint_ids = ids;
}


//*******************************************************************//
// This outputs the reference points in order of partition.
//*******************************************************************//
void iDistance::outputReferencePoints()
{
  int offset;
  for(int i=0; i < number_partitions; i++) //for each partition
  {
      cout << "Partition " << i << ": ";
      offset = i*number_dimensions; //partition offset into 1D array
      for(int j=0; j < number_dimensions; j++) //for each dimension
      {
          if(j != 0)
          {
            cout << ", ";
          }
          cout << reference_points[offset+j];
      }
      cout << endl;
  }
  
  
}


void iDistance::buildSplits()
{
	    switch(settings.refs_build)
    {
      case HALF_POINTS:
        buildSplits_TOP_CENTERS();
        break;
      case REF_FILE:
        //already done before we init to here
        break;
    }
	
}

void iDistance::buildSplits_TOP_CENTERS()
{
	
	//use findmedians (copied below from Helper class)
	//pick top X splits near the center of the data space (0.5)

	
	double* medianVals = findMedians(datapoints, number_points, number_dimensions, false);
	
	cout << "MEDIANS: " << endl;
	for(int i = 0; i < number_dimensions; i++)
	{
		cout << medianVals[i] << "," << endl;
	}
	
}




/*****************************************************
 * Brute force method to find the median of each
 * dimension and return them as an array. Note the medians
 * are in the same order as the original point dimensions.
 * If the data still has an ID attribute (as the first dim)
 * then the input flag will remove it and return just
 * the (d-1) median values
 * 
 * Input: data, number of points, number of dimensions, Id dimension flag
 * Output: array of median values as doubles
 */ 
 double* iDistance::findMedians(double* data, int p, int d, bool hasID)
{
  
  double* meds = new double[d]; //found median array
  
  //set up for odd number of points
  bool medSplit = false;  
  int m = (p / 2);
  int m2 = -1; //not used for odd
  
  if(p % 2 == 0) //if even number of points
  {
    medSplit = true;
    m = (p / 2) - 1;
    m2 = (p / 2);
  }
  
//  cout << "m | m2 = " << m << " | " << m2 << endl;

  int i = 0;
  if(hasID)
  {
    i = 1; //skip first dimension if its the ID attribute
  }
  while(i < d)
  {
    //double* points = new double[p];
    vector<double> points;
    
//    cout << i << ": ";
    
    for(int j = 0; j < p; j++)
    {
      //add each point's ith value to the vector
      points.push_back(  data[(j*d) + i] );
//      cout << points[j] << ",";
    }
//    cout << endl;
    
    //sort points vector
    sort(points.begin(), points.end());
    
    //get and save median
//    cout << "sorted: ";
//    for(int j = 0; j < p; j++)
//    {
//      cout << points[j] << ",";
//    }
//    cout << endl;
  
    if(medSplit)
    {
      //average of the two middle points
      meds[i] = (points[m] + points[m2]) / 2.0;
    }
    else
    {
      meds[i] = points[m];
    }
    
    i++; //go to next dimension
  }
  
  return meds;
}



//*******************************************************************//
// Master Reference point building function which calls the specialized
// methods. These methds have the common naming convention
// "buildRefs_" + enum name. 
// Also note, FROM_FILE build method occurs before this because of the file load
//*******************************************************************//
void iDistance::buildReferencePoints()
{
    //cout << "building with refs: " << settings.refs_build << endl;
    switch(settings.refs_build)
    {
      case HALF_POINTS:
        buildRefs_HALF_POINTS();
        break;
      case HALF_POINTS_CENTER:
        buildRefs_HALF_POINTS_CENTER();
        break;
      case HALF_POINTS_OUTSIDE:
        //cout << "HERE!!!" << endl;
        buildRefs_HALF_POINTS_OUTSIDE();
        break;
      case RANDOM:
        buildRefs_RANDOM();
        break;
      case CLUSTER:
        buildRefs_CLUSTER();
        break;
      case UNIFORM:
        buildRefs_UNIFORM();
        break;
      case REF_FILE:
        //already done.
        break;
    }

}



void iDistance::buildRefs_HALF_POINTS()
{
  //for every dimension i, two points:
  //(0.5, 0.5, ..., i = 0.0, ..., 0.5)
  //(0.5, 0.5, ..., i = 1.0, ..., 0.5)
  
  number_partitions = number_dimensions*2;
  reference_points = new double[number_partitions*number_dimensions];
  
  //initialize everything to 0.5
  for(int i = 0; i < number_partitions * number_dimensions; i++)
  {
    reference_points[i] = 0.5;
  }
   
  //set 0, 1 values
  //half_offset is to the start of the second half of ref points
  //the dim^2 is because i want dim*(.5*partss) and parts = 2*dim
  int half_offset = number_dimensions*number_dimensions;
  for(int i=0; i < number_dimensions; i++) //for each dimension
  {
      //set only its points to 0 and 1, 
      //walking diagonally through array in both spots at once
      reference_points[(i*number_dimensions) + i] = 0;
      reference_points[(half_offset) + (i*number_dimensions) + i] = 1;
  }

}

void iDistance::buildRefs_HALF_POINTS_CENTER()
{
  
    //plus one more ref point in the space center,
    //so this point remains all 0.5 values.
    number_partitions = number_dimensions*2 + 1;
    reference_points = new double[number_partitions*number_dimensions];
    
    //initialize everything to 0.5
    for(int i = 0; i < number_partitions * number_dimensions; i++)
    {
      reference_points[i] = 0.5;
    }
    
    //set 0, 1 values
    //half_offset is to the start of the second half of ref points
    //the dim^2 is because i want dim*(.5*partss) and parts = 2*dim
    int half_offset = number_dimensions*number_dimensions;
    for(int i=0; i < number_dimensions; i++)
    {
        //set only its points to 0 and 1, 
        //walking diagonally through array in both spots at once
        reference_points[i*number_dimensions + i] = 0;
        reference_points[(half_offset) + (i*number_dimensions) + i] = 1;
    }
    //it will leave the last partition alone with all .5s.

}


void iDistance::buildRefs_HALF_POINTS_OUTSIDE()
{
  //for every dimension i, two points:
  //(0.5, 0.5, ..., i = MIN, ..., 0.5)
  //(0.5, 0.5, ..., i = MAX, ..., 0.5)
  
  //where MIN and MAX are the arbitrarily far outside points from the dataspace
  // example, lets set them to -10, +11
  
  
//  cout << "HERE! outside" << endl;
  
  int outside_min = 0 - settings.refs_dist;
  int outside_max = 1 + settings.refs_dist;
  
  number_partitions = number_dimensions*2;
  reference_points = new double[number_partitions*number_dimensions];
  
  //initialize everything to 0.5
  for(int i = 0; i < number_partitions * number_dimensions; i++)
  {
    reference_points[i] = 0.5;
  }
   
  //set 0, 1 values
  //half_offset is to the start of the second half of ref points
  //the dim^2 is because i want dim*(.5*partss) and parts = 2*dim
  int half_offset = number_dimensions*number_dimensions;
  for(int i=0; i < number_dimensions; i++) //for each dimension
  {
      //set only its points to 0 and 1, 
      //walking diagonally through array in both spots at once
      reference_points[(i*number_dimensions) + i] = outside_min;
      reference_points[(half_offset) + (i*number_dimensions) + i] = outside_max;
  }

}


void iDistance::buildRefs_RANDOM()
{
  
  number_partitions = number_dimensions*2;
  reference_points = new double[number_partitions*number_dimensions];

  //cout.precision(10);
  
  for(int i = 0; i < number_partitions * number_dimensions; i++)
  {
    double temp = rand();
    temp = temp / 100000.0;
    temp = roundOff(temp - int(temp));
    //cout << "temp3 = " << temp << endl;

    reference_points[i] = temp;
  }
  
}


void iDistance::buildRefs_CLUSTER()
{
    
    
}

void iDistance::buildRefs_UNIFORM()
{
  
  //DOES NOT WORK YET!!
  number_partitions = number_dimensions*2;
  reference_points = new double[number_partitions*number_dimensions];

  //cout.precision(10);
  
  int splits_per_dim = number_partitions / number_dimensions;
  double split_length = 1.0 / (splits_per_dim + 1.0);
  
  cout << "splits per dimension = " << endl;
  for(int i = 1; i <= splits_per_dim; i++)
  {
    cout << split_length * i << endl;
  }
  
  for(int i = 0; i < number_partitions; i++)
  {
      double temp = 0.0;

      reference_points[i] = temp;  
  }
    
}


//*******************************************************************//
// if given the split points: splits 1D array, num is the number of splits
// first line is the dims indexes
// second line is the split vals for each dim
//*******************************************************************//
void iDistance::setSplitPoints(double* splits, int num)
{
    //cout << "HERE!" << endl;
    
	for(int i = 0; i < num; i++)
	{
		cout << "split " << i << " on dim " << splits[i] 
			 << " at val " << splits[num + i] << endl; 
	}
	
	
}


//*******************************************************************//
// if given the reference points
//*******************************************************************//
void iDistance::setReferencePoints(double* refs, int num)
{
    reference_points = refs;
    number_partitions = num;
}


void iDistance::setPointAssignment(double* assigns, int num)
{
  point_assignments = assigns;
  
}

//*******************************************************************//
// This function inserts all of the ids into the b+tree.  This is
// simply so that we can time how long it takes.
//*******************************************************************//
void iDistance::buildTree()
{
	//cout<<"BUILD TREE"<<endl;
    for(int i = 0; i < number_points; i++)
    {
        //cout << i << "=" << datapoint_index[i]<<endl;
        
        btree.insert2(datapoint_index[i],i);
    }
        
     //___________STAT COLLECTION___________//
	btree_mm :: tree_stats ts;
	ts = btree.get_stats();
	theStats.treeItems = ts.itemcount;
	theStats.treeNodes = ts.nodes(); //fulfills node count
	theStats.treeLeaves = ts.leaves;
	theStats.treeInner = ts.innernodes;
	theStats.treeLevels = ts.levels;
	theStats.treeAver = ts.avgfill_leaves();

	//ofstream fragFile;
	//fragFile.open ("treeFrag.txt");
	//btree.print_tree_segments(fragFile);
	//fragFile.close();
     //___________STAT COLLECTION___________//
}

//*******************************************************************//
// writes tree to serialized file
// 
// Takes in a filename without extension
// and writes out .tree and .data binary files. The .tree is a serialized
// output included with the standard stx b+ tree. The .data file contains
// a sequence of important data required in addition to the tree data
//*******************************************************************//
void iDistance::saveTree(const char* fStr)
{
    char* ext= new char[5];
    strcpy(ext, ".tree");
    char* fname = new char[strlen(fStr)+6]; 
    strcpy(fname, fStr);
    
    ofstream fout;
    fout.open(strcat(fname,ext), ios::out);
    btree.dump(fout);
    fout.close();
    
    //save other data structures
    strcpy(ext, ".data");
    strcpy(fname, fStr);
    
    fout.open(strcat(fname,ext), ios::out | ios::binary);
    
    //write needed metadata
    fout.write((char *)(&number_points), sizeof(number_points));
    fout.write((char *)(&number_dimensions), sizeof(number_dimensions));
    fout.write((char *)(&number_partitions), sizeof(number_partitions));
    fout.write((char *)(&constant_c), sizeof(constant_c));
    //write datapoints
    fout.write(reinterpret_cast<const char*> (datapoints), sizeof(double)*number_points*number_dimensions);
    //write reference points
    fout.write(reinterpret_cast<const char*> (reference_points), sizeof(double)*number_partitions*number_dimensions);
    //write partition dist max
    fout.write(reinterpret_cast<const char*> (partition_dist_max), sizeof(double)*number_partitions);
    //write partition_dist_max_index
    fout.write(reinterpret_cast<const char*> (partition_dist_max_index), sizeof(int)*number_partitions);
    //write datapoint ids
    fout.write(reinterpret_cast<const char*> (datapoint_ids), sizeof(double)*number_points);
    //write datapoint index
    fout.write(reinterpret_cast<const char*> (datapoint_index), sizeof(double)*number_points);
    
    fout.close();
}

//*******************************************************************//
// loads tree from serialized file
//
// this takes in a .tree file name and reads in the
// binary tree file and associated .data files with the same name.
//*******************************************************************//
void iDistance::loadTree(const char* fp)
{
    
    string dataStr = fp;
    dataStr = dataStr.substr(0, strlen(fp)-5);
    dataStr = dataStr + ".data";
    
    ifstream fin;
    fin.open(fp, ios::in);
    btree.restore(fin);
    fin.close();
    
    //load other data
    fin.open(dataStr.c_str(), ios::in | ios::binary);
    
    //read needed metadata
    fin.read((char *)(&number_points), sizeof(number_points));
    fin.read((char *)(&number_dimensions), sizeof(number_dimensions));
    fin.read((char *)(&number_partitions), sizeof(number_partitions));
    fin.read((char *)(&constant_c), sizeof(constant_c));
    
    
    //read datapoints
    datapoints = new double[number_points*number_dimensions];
    fin.read(reinterpret_cast<char *> (datapoints), sizeof(double)*number_points*number_dimensions);
    //read reference points
    reference_points = new double[number_partitions*number_dimensions];
    fin.read(reinterpret_cast<char *> (reference_points), sizeof(double)*number_partitions*number_dimensions);
    //read partition dist max
    partition_dist_max = new double[number_partitions];
    fin.read(reinterpret_cast<char *> (partition_dist_max), sizeof(double)*number_partitions);
    //read partition_dist_max_index
    partition_dist_max_index = new unsigned int[number_partitions];
    fin.read(reinterpret_cast<char *> (partition_dist_max_index), sizeof(int)*number_partitions);
    //read datapoint ids
    datapoint_ids = new double[number_points];
    fin.read(reinterpret_cast<char *> (datapoint_ids), sizeof(double)*number_points);
    //read datapoint index
    datapoint_index = new double[number_points];
    fin.read(reinterpret_cast<char *> (datapoint_index), sizeof(double)*number_points);
    
    fin.close();
}
    
//*******************************************************************//
//Euclidean distance function.
//*******************************************************************//
double iDistance::dist(double p1[], double p2[])
{
    double sum = 0.0;
    double current_dist, diff;
    
    for(int j = 0; j < number_dimensions; j++) //for each dimension
    {
        diff = p1[j] - p2[j]; //subtract
        sum += diff*diff; //and square, keep the total sum
        
        //cout<<" DIFF: "<<diff;
        //cout<<" SUM: "<<sum<<endl;
    }
    
    current_dist = roundOff(sqrt(sum));
    
    //cout<<"DIST ROOT: "<<current_dist<<endl;
    return current_dist;
}
//*******************************************************************//
// options are set via a struct with the possibilities
//*******************************************************************//
void iDistance::setOptions(configOptions co)
{
    settings = co;
}

//*******************************************************************//
// Given a vector of indices, return a vector of point IDs
//*******************************************************************//
vector<int> iDistance::getIDs(vector<int> indices)
{
    vector<int> ret;
    for(int i=0;i<indices.size();i++)
    {
        ret.push_back(datapoint_ids[indices[i]]);
    }
    return ret;
}


void iDistance::resetCandidateCount()
{
  
  candidates = 0;
  
}

int iDistance::getCandidateCount()
{
  if(doSeqSearch == 1) //do sequential search instead!
  {
      return number_points; // return list of final set
  }
  return candidates;
  
}



void iDistance::resetPartitionCount()
{
  
  partitionsChecked = 0;
  
}

int iDistance::getPartitionCount()
{
  if(doSeqSearch == 1) //do sequential search instead!
  {
      return number_partitions; // return list of final set
  }
  return partitionsChecked;
  
}
//*******************************************************************//
// Return the vector of knn distances 
//*******************************************************************//
vector<double> iDistance::getKNN_dists()
{
    return knn_S_dists;
}
//*******************************************************************//
// Return the vector of knn num nodes visited
//*******************************************************************//
vector <int> iDistance::getKNN_nodes()
{
	return knn_S_nodes;
}

vector <int> iDistance::getKNN_candidates()
{
	return knn_S_candidates;
}

double iDistance::getKNN_avgDist()
{
	return AverageDist(knn_S_dists);
}

//*******************************************************************//
// Externally reset the node visit counter
//*******************************************************************//
void iDistance::resetNodeCount()
{
  stx::CountQueryNodes("Start");
}

int iDistance::getNodeCount()
{
  if(doSeqSearch == 1) //do sequential search instead!
  {
      return theStats.treeNodes; // return total nodes in tree
  }
  
  return stx::CountQueryNodes("Get");
}

//*******************************************************************//
// Prints the query statistics, including
// Query answers, the number of nodes visited for each node in the answer, 
// the distance from the query point to the nodes, the average of those 
// distances, and the total time in milliseconds.   
//*******************************************************************//
/*
void iDistance::PrintQueryStats(ostream &out)
{
	//check to make sure there are results to match
     if(theStats.ind.empty())
	out <<"There are no matches for this query." << endl;
     else{
	    out<<"  RESULTS:"<<endl;
	    out << "    ";
	    printVector(getIDs(theStats.ind), out, true);

	    out<<"  Nodes Visited: " <<endl;
	    out << "    ";
	    printVector(theStats.numNodes, out, true);

	    if (theStats.knn)
	    {
		out<<"  kNN Distances: " <<endl;
		printVector(knn_S_dists, out, false);
	    
		out<<"  kNN Average Distance: " << AverageDist(knn_S_dists) << endl;
	    }
      }
    out<<"  Query Search Time: " << theStats.queryTime << " Milliseconds" <<endl;
    out<<endl;
}
*/
//*******************************************************************//
/*
void iDistance::PrintQueryStats(const char* fp)
{
    char* ext= new char[8];
    strcpy(ext, ".results");
    char* fname = new char[strlen(fp)+9]; 
    strcpy(fname, fp);
    
    //Print to File 
    ofstream resFile;
    resFile.open(strcat(fname,ext));
    if ( !resFile.is_open() )
    {
      cout << "Unable to open file!" << endl;
      return;
    }
    PrintQueryStats(resFile);
    resFile.flush();
    resFile.close();   
}

*/

//*******************************************************************//
void iDistance::getStats(stats* s)
{
  *s = theStats;
}

void iDistance::PrintNodeRange(const char* fp)
{
  //Print to File 
    char* ext= new char[9];
    strcpy(ext, ".segments");
    char* fname = new char[strlen(fp)+10]; 
    strcpy(fname, fp);
 
    ofstream segFile;
    segFile.open(strcat(fname,ext));
    if( !segFile.is_open() )
    {
      cout << "Unable to open file!" << endl;
      return;
    }
    
    PrintNodeRange(segFile);
    
    segFile.flush();
    segFile.close();
    
}

//*******************************************************************//
//print node ranges
//*******************************************************************//
void iDistance::PrintNodeRange(ostream &out)
{
    int leaves = theStats.treeLeaves;
//    out<<"Number of leaves: "<<leaves<<endl;
    //print the ranges in each node
    double low,high;
    //put iterator at very first key entry
    btree_mm::iterator node=btree.lower_bound(0.0);
    //nodecount starts at 0;
    int nodecount=stx::CountQueryNodes("Start");
    
    low = node.key();
    while(nodecount<leaves-1)   
    {
       //we've moved to a new node
       if(nodecount < stx::CountQueryNodes("Get"))
       {
 //          out<<"Node "<<nodecount+1<<" Low: "<<low<<" High: "<<high<<endl;
           out << low << " " << high << endl;
 
           low = node.key();
           high = 0;
           //nodecount = stx::CountQueryNodes("Get");
           nodecount++;
       }
       else
       {
           high = node.key();
           node++;       
       }
    }
    
    double old_high = -1;
    while(high != old_high)
    {
        old_high = high;
        node++;
        if(node.key() > 0)
            high=node.key();
    }
 //   out<<"Node "<<nodecount+1<<" Low: "<<low<<" High: "<<high<<endl;
    out << low << " " << high << endl;
    out.flush();
}


//*******************************************************************//
// Return the average distance of a vector
//   input: a vector of doubles 
//   output: the average of the input
//*******************************************************************//
double iDistance::AverageDist(vector<double>& distances)
{
	double aver;	
	for (int i = 0; i < distances.size(); i++){
		aver = aver + distances[i];
	}
	aver = aver / distances.size();

	return roundOff(aver);
}

//*******************************************************************//
//Print the elements of a vector 
//  input: Vector vec to be printed
//         Boolean commas, true -separate elements with a comma
//*******************************************************************//
template<typename T>
void iDistance::printVector(const std::vector<T>& vec, ostream &out, bool commas)
{
    if(commas)
    {
        for (int i = 0; i < vec.size()-1; i++)
        {
            out<< vec[i] << ", ";
        }
        out<<vec[vec.size()-1]<<endl;
    }
    else
    {
        for (int i = 0; i < vec.size(); i++)
        {
            out<< "    " << vec[i] << endl;
        }
    }
}

/***************************************************
 * Calculates the theoretical upper boundary distance of two points
 * in d dimensions of [0,1] space. This can be set as an non-optimal
 * c value that will gaurantee separation in the b+ tree.
 ****************************************************/
int iDistance::calculateC(int d, int e)
{
  
  //cout << "calculate c given " << d << ",  " << e << endl;
  
  //added so we can dynamically set c based on what config file refs-dist is
  //space given extra e (default is e=0, no affect)
  double edist = (2 * e) + 1;
  
  //since each dim is max length of 1
  //multiply by 2 for safety 
  double diag = 2 * sqrt( (edist * d));
  
  //old
//  double diag = sqrt(d);
  
  //added ONLY for Tim's EM algorithm 4/25/2012
  //diag = diag * 3;
  

  //round up
  int c = (int)(diag + 0.5);

  //cout << "final c = " << c << endl;
  
  return c;

}


double iDistance::roundOff(double in)
{
    int temp = (int)((in * 1000) + 0.5);
    return (temp / 1000.0);
}


double iDistance::calcIndex(int p, double c, double d)
{
    return roundOff( (p*c + d) );
}

//*******************************************************************//
//this does the knn search.
//*******************************************************************//


