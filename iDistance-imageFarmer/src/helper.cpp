/*////////////////////////////////////////////////////////////////////
helper.cpp

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
#include <stdio.h>
#include <stdlib.h>
#include "string.h"
#include <vector>
#include <time.h>
#include <algorithm>

#include <sys/time.h>

#include "helper.h"

using namespace std;

/********************************************************
 * Trims leading and trailing whitespace characters,
 *   including space, tab, and new lines
 * 
 * Input: string
 * Output: trimmed string
 * 
 */
inline std::string trimStr(std::string& str)
{
  str.erase(0, str.find_first_not_of(" \t\r\n"));  
  str.erase(str.find_last_not_of(" \t\r\n")+1);         
  
  return str;
}

std::string trimToName(std::string str)
{
  str.erase(0, str.find_last_of("/\\")+1);  
  str.erase(str.find_last_of("."));         
  
  return str;
}


/*********************************************
 * Simple helper function to print to a given filename and 
 * possibly cout as well. If the file exists, it will append to 
 * the file. By using filename instead of stream, it is 
 * more encapsulated and removed from user, also we can print 
 * to 0,1,2 streams with one external line of code.
 * 
 * A good tip is to use a ostringstream to manually buffer up output
 * and then call this function even fewer times.
 * 
 * Input: filename, string to print, console flag
 * 
 */ 
void printOut(const char* fp, string printStr, bool toConsole)
{
  
  if(toConsole) //do this first!
  {
    cout << printStr;
  }
  
  if(fp == NULL) //technically a way to get out without printing
  {
    return;
  }
  
  ofstream fileOut;
  fileOut.open(fp, ios_base::app);
  
  if( !fileOut.is_open() )
  {
    cout << "Unable to open file for printing!" << endl;
    return;
  }
  
  fileOut << printStr; // don't add a newline here -- up to input
  
  fileOut.flush();
  fileOut.close();
  
}

/*********************************************
 * Reads a given csv file and saves the data found
 * in a 2D array of [points][dimensions]
 * 
 * Input: data file name, number of points, number of dimensions
 * Output: 2D double array of points and their dimensions
 * 
 */ 
double** readDataFile2D(const char* fp, int p, int d)
{
  /* Data file line:
   * id, dim1, dim2, ..., dimN
   * 
   * Query file line:
   * dim1, dim2, ..., dimN
   * 
   * example:
   * uniform: 0,0.0173,0.3729,0.1081,0.0466,0.5339,0.3966,0.3892,0.6968,0.8678,0.5521,0.3693,0.8507,0.0219,0.6543,0.5903,0.8401
   * cluster: 0,0.7578,0.8032,0.7382,0.309,0.4025,0.3287,0.1351,0.6352,0.2762,0.4595,0.523,0.6555,0.159,0.8006,0.8809,0.2546
   * real: 155,0.3132,0.4807,0.4896,0.3274,0.3117,0.5368,0.6084,0.4306,0.2527,0.4752,0.3269,0.3477,0.2661,0.2839,0.2907,0.2188
   * 
   * Query file: 0.0173,0.3729,0.1081,0.0466,0.5339,0.3966,0.3892,0.6968,0.8678,0.5521,0.3693,0.8507,0.0219,0.6543,0.5903,0.8401
   * 
   */
  
  double** data = NULL; //this will hold the entire 2D data array
    
  string line;
  int lineI = 0; 
  
  ifstream fileIn(fp);
  if ( !fileIn.is_open() )
  {
    cout << "Unable to open file!" << endl;
    return NULL;
  }
  
  data = new double*[p]; // initialize data array of points
  
  while ( fileIn.good() )
  {
    getline(fileIn,line);
    
    line = trimStr(line);
    
    if(line.length() != 0)
    {

      data[lineI] = new double[d];
    
      int tokI = 0;
      char* tok = NULL;
      char* cStr = new char[line.length()+1]; //+1 for auto added null terminator
      strcpy(cStr, line.c_str());
      
      tok = strtok(cStr, ",");
      while (tok != NULL)
      {
        data[lineI][tokI] = strtod(tok, NULL);
    
        tok = strtok(NULL, ",");
        tokI++;
      }
      
      lineI++;
      
      delete [] cStr;
    }
  }
  
  fileIn.close();
  
  return data;
  
}

/************************************************
 * Read a given csv file and save the data as a
 * 1D array of all the points' dimension values 
 * sequentially. Requires more complicated index math
 * to find specific data point, but might be faster 
 * with memory accesses.
 * 
 * Input: file name, number of points, number of dimensions
 * Output: 1D array of all points' dimensions
 */ 
double* readDataFile(const char* fp, int p, int d)
{
  /* Data file line:
   * id, dim1, dim2, ..., dimN
   * 
   * Query file line:
   * dim1, dim2, ..., dimN
   * 
   * example:
   * uniform: 0,0.0173,0.3729,0.1081,0.0466,0.5339,0.3966,0.3892,0.6968,0.8678,0.5521,0.3693,0.8507,0.0219,0.6543,0.5903,0.8401
   * cluster: 0,0.7578,0.8032,0.7382,0.309,0.4025,0.3287,0.1351,0.6352,0.2762,0.4595,0.523,0.6555,0.159,0.8006,0.8809,0.2546
   * real: 155,0.3132,0.4807,0.4896,0.3274,0.3117,0.5368,0.6084,0.4306,0.2527,0.4752,0.3269,0.3477,0.2661,0.2839,0.2907,0.2188
   * 
   * Query file: 0.0173,0.3729,0.1081,0.0466,0.5339,0.3966,0.3892,0.6968,0.8678,0.5521,0.3693,0.8507,0.0219,0.6543,0.5903,0.8401
   * 
   */
  
  double* data = NULL; //this will hold the entire 2D data array
    
  string line;
  int lineI = 0; 
  
  ifstream fileIn(fp);
  if( !fileIn.is_open() )
  {
    cout << "Unable to open file" << endl;
    return NULL;
  }
  
  data = new double[p*d]; // initialize data array
  
  while ( fileIn.good() )
  {
    getline(fileIn,line);
    
    //cout << line << endl;
     //lets just trim the line anyways, won't hurt
    line = trimStr(line);
    
    if(line.length() != 0) //nonempty line
    {
    
      int tokI = 0;
      char* tok = NULL;
      char* cStr = new char[line.length()+1]; //+1 for auto added null terminator
      strcpy(cStr, line.c_str());
      
      tok = strtok(cStr, ",");
      
      //use idp flag to investigate a specific point if needed
      //bool idp = false;
      //if(strcmp(tok, "6131") == 0){idp = true;}
      
      while (tok != NULL)
      {
        data[ (lineI * d) + tokI ] = strtod(tok, NULL);
    
        //if(idp){cout << tok << ", "; }
        
        tok = strtok(NULL, ",");
        tokI++;
      }
      
      //if(idp){cout << endl;}
      
      lineI++;
      
      delete [] cStr;
    }   
  }
    
  fileIn.close();
  
  return data;
  
}

/********************************************
 * Reads a given csv file and calculates the 
 * number of give non-empty lines and the number
 * of values in the first non-empty line
 * 
 * Input: File name
 * Output: Int Array of: line count, 
 * 
 */ 
int* getFileInfo(const char* fp)
{
  string line;
  int lineC = 0; 
  int valC = 0;
  bool checkVals = true;
  ifstream fileIn(fp);
  
  if( !fileIn.is_open() ) //check good open
  {
    cout << "Unable to open file!" << endl;
    return NULL;
  }
  
  while( fileIn.good() ) //check for EOF
  {
    getline(fileIn,line); //get next line as string
    
    /* Line character debugging
     *
    
    cout << fileIn.good() << endl;
    
    cout << "line='" << line << "'" << endl;

    //doesn't find any of these, even though you'd expect so..
    if(strcmp(line.c_str(), "\r") == 0)
    {
      cout << "r!" << endl;
    }
    if(strcmp(line.c_str(), "\n") == 0)
    {
      cout << "n!" << endl;
    }
    if(strcmp(line.c_str(), "\r\n") == 0)
    {
      cout << "rn!" << endl;
    }
    
    cout << "len = " << line.length() << endl;
    
    *
    */ 
    
   
    //lets just trim the line anyways, won't hurt
    line = trimStr(line);
    
    //cout << line << endl;
    
    if(line.length() != 0) //nonempty line
    {
   
      lineC++; //increment to current number of lines
      
      if(checkVals) //get dims too, only needed once
      {
        char* tok = NULL;
        char* cStr = new char[line.length()+1]; //+1 for auto added null terminator
        strcpy(cStr, line.c_str());
   
        //don't actually care about token value, just that they exist
        tok = strtok(cStr, ","); 
        while (tok != NULL)
        {
          tok = strtok(NULL, ",");
          valC++;
        }
        
        checkVals = false;       
        
        delete [] cStr;
      }
    }    
  }
  fileIn.close();
  int* info = new int[2];
  info[0] = lineC;
  info[1] = valC;
  
  
  return info;
  
}

/* This function read an INI formatted text file and creates an
 * iDistance configuration struct to update algorithmic settings
 */
int readConfigFile(const char* fp, configOptions* settings)
{
  /* Sample file:
   * 
   * ; comment lines
   * [section header]
   * param1 = value1
   * param2 = value2
   * 
   */
  string line;

  ifstream fileIn(fp);
  if( !fileIn.is_open() )
  {
    cout << "Unable to open file!" << endl;
    return -1;
  }
    
  while ( fileIn.good() )
  {
    getline(fileIn,line);
    line = trimStr(line);
    
    char* tok1 = NULL; //option name (key)
    char* tok2 = NULL; //option value
    
    char* cStr = new char[line.length()+1]; //+1 for auto added null terminator
    strcpy(cStr, line.c_str());
    
    //check for proper lines
    if(strlen(cStr) == 0)
    {
      //empty line
    }
    else if(cStr[0] == ';' || cStr[0] == '[')
    {
      //comments and headers, skip
    }
    else //actual param lines
    {
      
      tok1 = strtok(cStr, "=");
      string tokStr1 = tok1; 
      tokStr1 = trimStr(tokStr1); //clean up whitespace
  //    cout << "Option: '" << tokStr1 << "'" << endl;
        
      tok2 = strtok(NULL, "=");
      string tokStr2 = tok2; 
      tokStr2 = trimStr(tokStr2); //clean up whitespace
  //    cout << "Value: '" << tokStr2 << "'" << endl;
      
      setConfigOption(tokStr1.c_str(), tokStr2.c_str(), settings);
      
    }
    delete [] cStr;
    
  }//end while line
    
  fileIn.close();
  
  return 0;
}

const char* stringToUpper(const char* cStr)
{
  string str = cStr;
  int i=0;
  char c;
  while (cStr[i])
  {
    c=cStr[i];
    str[i] = (toupper(c));
    i++;
  }
  return str.c_str();
}

int stringToEnum(const char* cStr, const char* opts[], int optCount)
{
  for(int i = 0; i < optCount; i++)
  {
    if( strcmp(cStr, opts[i]) == 0 )
    {
      return i; //assumes enum value is array index
    }
  }
  return -1; //didn't match anything, so return -1 as error
}


void setConfigOption(const char * key, const char * val, configOptions * settings)
{
    int enumKey = stringToEnum(stringToUpper(key), def_Options, num_Options);
    
    switch(enumKey)
    {
        case VERSION:
        {
            int enumVal = stringToEnum(stringToUpper(val), def_Versions, num_Versions);
            settings->algo_version = (enumVal != -1) ? enumVal : 0;
        }
        break;
        case RANGE:
        {
            int enumVal = stringToEnum(stringToUpper(val), def_Ranges, num_Ranges);
            settings->range_method = (enumVal != -1) ? enumVal : 0;
        }
        break;
        case KNN:
        {
            int enumVal = stringToEnum(stringToUpper(val), def_KNNs, num_KNNs);
            //cout << "ENUM KNN = " << enumVal << endl;
            settings->knn_method = (enumVal != -1) ? enumVal : 0;
        }
        break;
        case C_TYPE:
        {
            int enumVal = stringToEnum(stringToUpper(val), def_Types, num_Types);
            settings->c_type = (enumVal != -1) ? enumVal : 0;
        }
        break;
        case C_VAL:
        {
            settings->c_val = atoi(val);
        }
        break;
        case R_TYPE:
        {
            int enumVal = stringToEnum(stringToUpper(val), def_Types, num_Types);
            settings->r_type = (enumVal != -1) ? enumVal : 0;
        }
        break;
        case R_INIT:
        {
            settings->r_init = atof(val);
        }
        break;
        case R_DELTA:
        {
            settings->r_delta = atof(val);
        }
        break;
        case REFS_BUILD:
        {
            int enumVal = stringToEnum(stringToUpper(val), def_RefBuilds, num_RefBuilds);
            //cout << " refsbuild enum = " << enumVal << endl;
            settings->refs_build = (enumVal != -1) ? enumVal : 0;
        }
        break;
        case REFS_NUM:
        {
            settings->refs_num = atoi(val);
        }
        break;
        case REFS_DIST:
        {
            settings->refs_dist = atoi(val);
        }
        break;
        case REFS_ASSIGN:
        {
            int enumVal = stringToEnum(stringToUpper(val), def_RefAssigns, num_RefAssigns);
            settings->refs_assign = (enumVal != -1) ? enumVal : 0;
        }
        break;
        case SPLITS_BUILD:
        {
          int enumVal = stringToEnum(stringToUpper(val), def_SplitBuilds, num_SplitBuilds);
          settings->splits_build = (enumVal != -1) ? enumVal : 0;
        }
        break;
        case OUTPUT_MODE:
        {
          int enumVal = stringToEnum(stringToUpper(val), def_OutputModes, num_OutputModes);
          settings->output_mode = (enumVal != -1) ? enumVal : 0;
        }
        break;
    	
        
        
    } //end switch
  
  
  
  
  
}


void outputConfigOptions(configOptions* settings, ostream& out)
{
    
  out << "Configuration file details: " << endl;
  
  out << "  algo_version: " << settings->algo_version << " : " <<
    def_Versions[settings->algo_version] << endl;
  
  out << "  range_method: " << settings->range_method << " : " <<
    def_Ranges[settings->range_method] << endl;
  
  out << "  knn_method: " << settings->knn_method << " : " <<
    def_KNNs[settings->knn_method] << endl;
  
  out << "  c_type: " << settings->c_type << " : " <<
    def_Types[settings->c_type] << endl;
  
  out << "  c_value: " << settings->c_val << endl;
  
  out << "  r_type: " << settings->r_type << " : " <<
    def_Types[settings->r_type] << endl;
  
  out << "  r_init: " << settings->r_init << endl;
  
  out << "  r_delta: " << settings->r_delta << endl;
  
  out << "  refs_build: " << settings-> refs_build << " : " <<
    def_RefBuilds[settings->refs_build] << endl;
  
  out << "  refs_num: " << settings->refs_num << endl;

  out << "  refs_dist: " << settings->refs_dist << endl;
  
  out << "  refs_assign: " << settings->refs_assign << " : " <<
    def_RefAssigns[settings->refs_assign] << endl;
  
  out << "  refs_file: " << settings->refs_file << endl; 
  
  out << "  splits_build: " << settings->splits_build << " : " <<
    def_SplitBuilds[settings->splits_build] << endl;
  
  out << "  splits_file: " << settings->splits_file << endl; 

  out << "  output_mode: " << settings->output_mode << " : " <<
    def_OutputModes[settings->output_mode] << endl;
  
  
}

/***************************************
 * Simple start timer function to essentially 
 * increase code readability elsewhere, and 
 * make it easier to modify ALL timers later
 * 
 * Input: nothing
 * Output: a clock_t variable with the current time information 
 * 
 */ 
timeval startTimer()
{
  timeval start;
  gettimeofday(&start, NULL);
  return start;
}
/*
clock_t startTimer()
{
  return clock();
}
*/ 



void quickStats(vector<double> data, double* stats)
{
 
  double min = 9999, max = 0, mean = 0, temp = 0;
  
  for(int i = 0; i < data.size(); i++)
  {
    temp = data[i];
    if(temp < min)
    { min = temp; }
    if(temp > max)
    { max = temp; }
    
    mean = mean + temp;
  }
  
  stats[0] = min;
  stats[1] = max;
  stats[2] = mean / data.size();
  
  
}


void quickStats(vector<int> data, double* stats)
{
  double min = 9999, max = 0, mean = 0, temp = 0;
  
  for(int i = 0; i < data.size(); i++)
  {
    temp = data[i];
    if(temp < min)
    { min = temp; }
    if(temp > max)
    { max = temp; }
    
    mean = mean + temp;
  }
  
  stats[0] = min;
  stats[1] = max;
  stats[2] = mean / data.size();
  
  
}


/********************************************
 * A simple stop timer that subtracts the
 * current time info from a supplied time variable
 * and then converts the value to the number 
 * of milliseconds (as a double)
 * 
 * Input: clock_t variable as the start time
 * Output: double value of milliseconds elapsed until now
 * 
 */ 
double stopTimer(timeval start)
{
  timeval stop;
  gettimeofday(&stop, NULL);
  long seconds  = stop.tv_sec - start.tv_sec;
  long useconds = stop.tv_usec - start.tv_usec;
  double mtime = ((seconds) * 1000 + useconds/1000.0);
  
  return mtime;
}
/*
double stopTimer(clock_t start)
{
  clock_t end = clock() - start;
  return ((double)end / (double)CLOCKS_PER_SEC);
}
*/
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
double* findMedians(double* data, int p, int d, bool hasID)
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
  	meds[0] = -1;
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
void printTreeStats_CSV(ostream &out, stats* theStats)
{
    out << theStats->dims << ", " << theStats->points << ", ";
    out << theStats->treeNodes << ", " << theStats->treeLeaves << ", ";
    out << theStats->treeInner << ", " << theStats->treeLevels << ", ";
    out << theStats->treeAver << ", ";
    out << theStats->buildTime << ", " << theStats->indexTime << endl;
    
}

//*******************************************************************//
void printTreeStats(ostream &out, stats* theStats)
{
    out << "Generated Tree Results (.results file)" << endl;
    out << "  see .segments file for tree segments information\n" << endl;
    
    out << "Total Data Points: " << theStats->treeItems << endl;
    out << "Total Nodes: " << theStats->treeNodes << endl;
    out << "  Leaves: " << theStats->treeLeaves << endl;
    out << "  InnerNodes: " << theStats->treeInner << endl;   
    out << "  Levels: " << theStats->treeLevels << endl;   
    out << "  Average Fill: " << theStats->treeAver << endl;
    
    out << endl;
    out << "Tree Construction Time: " << theStats->buildTime << " milliseconds" << endl;
    out << "Index Construction Time: " << theStats->indexTime << " milliseconds" << endl;   

}

//*******************************************************************//
// Print the statistics about the tree, including:
// number of items, nodes, leaves, inner leaves, levels, the average fill
// and the time to construct the tree and the index

//*******************************************************************//
void printTreeStats(const char* fp, stats* treeStats)
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
      cout << "Unable to open file for tree stats!" << endl;
      return;
    }

    printTreeStats(resFile, treeStats);
    
    resFile.flush();
    resFile.close();   
}





void printResultsQP(ostream& out, vector<int>& results, int numNodes, double queryT)
{
    out << "  Found " << results.size() << " points: ";
    for(int j=0; j < results.size(); j++)
    {
        if(j != 0)
        {
          out << ", ";
        }
        if(j % 10 == 0)
        {
          out << endl << "    ";
        }
        out << results[j];
    }
    out << endl;
    out << "  Number of Nodes Visited: ";
    out << numNodes << endl;
    out << "  Query search time: " << queryT << " milliseconds." << endl << endl;
  
}


void printResultsQP_CSV(ostream& out, vector<double> qstats, vector<int> inds)
{
    for(int i = 0; i < qstats.size(); i++)
    {
      if(i != 0)
      {
        out << ", ";
      }
      out << qstats[i];
    }
    
    out << endl;
    
    for(int i = 0; i < inds.size(); i++)
    {
      if(i != 0)
      {
        out << ", ";
      }
      out << inds[i];
    }
    
    out << endl;
}


void printResultsQR(ostream& out, vector<int>& results, int numNodes, double queryT)
{
    out << "  Found " << results.size() << " points: ";
    for(int j=0; j < results.size(); j++)
    {
        if(j != 0)
        {
          out << ", ";
        }
        if(j % 10 == 0)
        {
          out << endl << "    ";
        }
        out << results[j];
    }
    out << endl;
    out << "  Number of Nodes Visited: ";
    out << numNodes << endl;
    out << "  Query search time: " << queryT << " milliseconds." << endl << endl;
  
}

void printResultsQR_CSV(ostream& out, vector<double> qstats, vector<int> inds)
{
    for(int i = 0; i < qstats.size(); i++)
    {
      if(i != 0)
      {
        out << ", ";
      }
      out << qstats[i];
    }
    
    out << endl;
    
    for(int i = 0; i < inds.size(); i++)
    {
      if(i != 0)
      {
        out << ", ";
      }
      out << inds[i];
    }
    
    out << endl;
}

void printResultsQN_CSV(ostream& out, vector<double> qstats, vector<int> results, vector<int> nodes, vector<int> candidates, vector<double> dists)
{
    for(int i = 0; i < qstats.size(); i++)
    {
      if(i != 0)
      {
        out << ", ";
      }
      out << qstats[i];
    }
    
    out << endl;
    
    for(int i = 0; i < results.size(); i++)
    {
      if(i != 0)
      {
        out << ", ";
      }
      out << results[i];
    }
    
    out << endl;
  
    for(int i = 0; i < candidates.size(); i++)
    {
      if(i != 0)
      {
        out << ", ";
      }
      out << candidates[i];
    }
    
    out << endl;
    
    for(int i = 0; i < nodes.size(); i++)
    {
      if(i != 0)
      {
        out << ", ";
      }
      out << nodes[i];
    }
    
    out << endl;
    
    for(int i = 0; i < dists.size(); i++)
    {
      if(i != 0)
      {
        out << ", ";
      }
      out << dists[i];
    }
    
    out << endl;
    
  
}


void printResultsQN_CSV_IF(ostream& out, vector<double> qstats, vector<int> results, vector<int> nodes, vector<int> candidates, vector<double> dists)
{
	/* Special output for imageFarmer
	 * 
	 * 3 lines per query:
	 * query point ID
	 * returned neighbor IDs
	 * returned neighbor distances
	 * returned neighbor nodes accessed
	 * total query time
	 */
	
	out << qstats[0] << endl;

	
    for(int i = 0; i < results.size(); i++)
    {
      if(i != 0)
      {
        out << ", ";
      }
      out << results[i];
    }
    
    out << endl;
  
    
    for(int i = 0; i < dists.size(); i++)
    {
      if(i != 0)
      {
        out << ", ";
      }
      out << dists[i];
    }
    
    out << endl;
    
    
    for(int i = 0; i < nodes.size(); i++)
    {
      if(i != 0)
      {
        out << ", ";
      }
      out << nodes[i];
    }
    
    out << endl;
    
    out << qstats[5] << endl;
	
  
}



/*
void printKNNQueryResults_RAW(ostream& out, const char* filename, int dims, int points, const char* qtype, vector<double> dists, vector<int> numNodes, vector<int> numCandidates, double queryT, int k)
{
    
    double* distStats = new double[3];
    double* nodeStats = new double[3];
    double* candidateStats = new double[3];
    
    quickStats(dists, distStats);
    quickStats(numNodes, nodeStats);
    quickStats(numCandidates, candidateStats);
                    
    //avg dist, max dist, avg nodes, max nodes, avg candidates, max candidates                
    out << filename[0] << ", " << dims << ", " << points << ", " << qtype[2] << ", ";
    out << k << ", " <<  distStats[2] << ", " <<  distStats[1] << ", ";
    out << nodeStats[2] << ", " << nodeStats[1] << ", ";
    out << candidateStats[2] << ", " << candidateStats[1] << ", " << queryT << endl;
}
*/


void printResultsQN(ostream& out, vector<int>& results, vector<double> dists, vector<int> numNodes, double avgDist, double queryT)
{
    out << "  Found closest " << results.size() << " points: ";
    for(int j=0; j < results.size(); j++)
    {
        if(j != 0)
        {
          out << ", ";
        }
        if(j % 10 == 0)
        {
          out << endl << "    ";
        }
        out << results[j];
    }
    out << endl;
    
    out << "  kNN Distances: ";
    for(int j=0; j < dists.size(); j++)
    {
        if(j != 0)
        {
          out << ", ";
        }
        if(j % 10 == 0)
        {
          out << endl << "    ";
        }
        out << dists[j];
    }
    out << endl;
    
    
    out << "  Number of Nodes Visited: ";
    for(int j=0; j < numNodes.size(); j++)
    {
        if(j != 0)
        {
          out << ", ";
        }
        if(j % 10 == 0)
        {
          out << endl << "    ";
        }
        out << numNodes[j];
    }
    out << endl;
    
    out << "  The average distance is: " << avgDist << endl; 
    out << "  Query search time: " << queryT << " milliseconds." << endl << endl;
  
}


/********************************************
 * A method to test possible dataset sizes,
 * and see how much memory it approximately 
 * takes in real time.
 * 
 * Input: Number of data points, Number of dimensions
 * Ouput: nothing
 * 
 */ 
void testMemoryData(int n, int d)
{
  //float p = 0.01234567;
  double p = 0.01234567;
  
  cout << "Testing n = " << n << ", and d = " << d << endl;
  
  double* pp = new double[n*d];
  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < d; j++)
    {
      pp[(i*d) + j] = p;
    }
  }
  
  cout << "All Data now in array!" << endl;
  
  //crude pause loop so you have time to check your system memory
  int t = 0;
  double pause = 200000000;
  for(double i = 0; i < pause; i++)
  {
    t = t * i; 
    //something to force loop despite compiler optimiztion
  }
  delete [] pp;
}



