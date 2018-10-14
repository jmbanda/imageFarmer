/*////////////////////////////////////////////////////////////////////
definitions.h

Copyright (C) 2012  Michael Schuh, Timothy Wylie, Rafal Angryk
  Data Mining Laboratory at Montana State University 
  Principle author contact: michael.schuh@cs.montana.edu

This file is part of iDistance*.

iDistance* is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; only version 2
of the License. See the COPYING file for more information.
////////////////////////////////////////////////////////////////////*/


#ifndef _DEFINITIONS_H_
#define _DEFINITIONS_H_

#include<vector>

//debugging preprocessor macro in order to use B+Tree print/etc
//#define BTREE_DEBUG
#ifdef BTREE_DEBUG
const bool DEBUGGING = true;
#else
const bool DEBUGGING = false;
#endif

#include "stx/btree_multimap.h"

using namespace std;

//this is the default traits struct with debug option
template <typename _Key, typename _Data, bool _Debug>
struct btree_map_traits
{
    static const bool   selfverify = false; 
    static const bool   debug = _Debug; 
    static const int    leafslots = BTREE_MAX( 8, 256 / (sizeof(_Key) + sizeof(_Data)) ); 
    static const int    innerslots = BTREE_MAX( 8, 256 / (sizeof(_Key) + sizeof(void*)) );
};


//the b+tree key/value pair to store
typedef double key_type; 
typedef int value_type;  //could have int[] for two ints?

//the template btree with correct types for easy use later
typedef stx::btree_multimap<key_type, value_type, std::less<key_type>,
        btree_map_traits<key_type,value_type,DEBUGGING>, std::allocator<std::pair<key_type, value_type> > > btree_mm;


/////////////////////////////////////////////////////
//algorithm configuration options

//#define stringify( name ) # name


static int num_Options = 14;
static const char* def_Options[] = {"VERSION", "RANGE", "KNN", 
    "C-TYPE", "C-VAL", "R-TYPE", "R-INIT", "R-DELTA", 
    "REFS-BUILD", "REFS-NUM", "REFS-DIST", "REFS-ASSIGN", 
	  "SPLITS-BUILD", "OUTPUT-MODE"};
enum {VERSION, RANGE, KNN, 
    C_TYPE, C_VAL, R_TYPE, R_INIT, R_DELTA,
    REFS_BUILD, REFS_NUM, REFS_DIST, REFS_ASSIGN,
	SPLITS_BUILD, OUTPUT_MODE };

static int num_Versions = 3;
static const char* def_Versions[] = {"ORIGINAL", "EXT0", "EXT1"};
enum { ORIGINAL, EXT0, EXT1 };
    
static int num_Ranges = 2;
static const char* def_Ranges[] = {"SIMPLE", "ADV1"};
enum { SIMPLE, ADV1 };

static int num_KNNs = 2;
static const char* def_KNNs[] = {"KEEP","RESTART"};
enum { KEEP, RESTART };

static int num_Types = 2;
static const char* def_Types[] = {"CALCULATE", "MANUAL"};
enum { CALCULATE, MANUAL };

static int num_RefBuilds = 7;
static const char* def_RefBuilds[] = { "HALF_POINTS", 
    "HALF_POINTS_CENTER", "HALF_POINTS_OUTSIDE", "RANDOM", "CLUSTER", "UNIFORM", "REF_FILE" };
enum { HALF_POINTS, HALF_POINTS_CENTER, HALF_POINTS_OUTSIDE, RANDOM, CLUSTER, UNIFORM, REF_FILE };

static int num_RefAssigns = 3;
static const char* def_RefAssigns[] = {"CLOSEST", "FARTHEST", "ASSIGN_FILE"};
enum { CLOSEST, FARTHEST, ASSIGN_FILE };

static int num_SplitBuilds = 2;
static const char* def_SplitBuilds[] = { "TOP_CENTERS", "SPLIT_FILE" };
enum { TOP_CENTERS, SPLIT_FILE  };

static int num_OutputModes = 2;
static const char* def_OutputModes[] = { "NORMAL", "OUTPUT_DIR" };
enum { NORMAL, OUTPUT_DIR };


//struct to hold all config options for the idistance algorithm
struct configOptions
{
    configOptions()
    {
        //initialize with defaults
        algo_version = ORIGINAL;
        range_method = SIMPLE;
        knn_method = KEEP;
        
        c_type = CALCULATE;
        c_val = -1; //calculated
        
        r_type = MANUAL;
        r_init = 0.01;
        r_delta = 0.01;
        
        refs_build = HALF_POINTS;
        refs_num = -1; //defined by half_points
        refs_dist = 0; //max distance beyond the dataspace refs can go
          //used for half points outside, as well as cluster center pushes
        
        refs_assign = CLOSEST;
        refs_file = ""; //not needed for half_points
        assigns_file = "";
        
        splits_build = TOP_CENTERS;
        splits_file = "";
        
        output_mode = NORMAL;
        output_dir = "";
    }
    
    int algo_version, range_method, knn_method;
    
    int c_type, r_type;
    double c_val, r_init, r_delta;
    
    int refs_build, refs_num, refs_dist, refs_assign;
    const char* refs_file;
    const char* assigns_file;
    
    int splits_build;
    const char* splits_file;
    
    int output_mode;
    const char* output_dir;
};

//iDistance statistics collected during program runs
struct stats {
	
  int points, dims; //dataset stats
  
  //tree stats
	int treeItems, treeNodes, treeLeaves, treeInner, treeLevels;
  double treeAver, buildTime, indexTime;

};



#endif //_DEFINITIONS_H_
