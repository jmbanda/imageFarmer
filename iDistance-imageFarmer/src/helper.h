/*////////////////////////////////////////////////////////////////////
helper.h

Copyright (C) 2012  Michael Schuh, Timothy Wylie, Rafal Angryk
  Data Mining Laboratory at Montana State University 
  Principle author contact: michael.schuh@cs.montana.edu

This file is part of iDistance*.

iDistance* is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; only version 2
of the License. See the COPYING file for more information.
////////////////////////////////////////////////////////////////////*/


#ifndef _HELPER_H_
#define _HELPER_H_

#include "definitions.h"

const char* stringToUpper(const char* cStr);
std::string trimToName(std::string str);

inline std::string trimStr(std::string& str);

int stringToEnum(const char* cStr, const char* opts[], int optCount);

void printTreeStats_CSV(ostream &out, stats* theStats);
void printTreeStats(ostream &out, stats* theStats);
void printTreeStats(const char* fp, stats* treeStats);

void printResultsQP(ostream& out, vector<int>& results, int numNodes, double queryT);
void printResultsQP_CSV(ostream& out, vector<double> qstats, vector<int> inds);
void printResultsQR(ostream& out, vector<int>& results, int numNodes, double queryT);
void printResultsQR_CSV(ostream& out, vector<double> qstats, vector<int> inds);
void printResultsQN(ostream& out, vector<int>& results, vector<double> dists, vector<int> numNodes, double avgDist, double queryT);
void printResultsQN_CSV(ostream& out, vector<double> qstats, vector<int> results, vector<int> numNodes, vector<int> numCandidates, vector<double> dists);
void printResultsQN_CSV_IF(ostream& out, vector<double> qstats, vector<int> results, vector<int> numNodes, vector<int> numCandidates, vector<double> dists);


void printOut(const char* fp, string printStr, bool toConsole);

configOptions readConfigFile(const char* fp);

double* readDataFile(const char* fp, int p, int d);

double** readDataFile2D(const char* fp, int p, int d);

int* getFileInfo(const char* fp);

int readConfigFile(const char* fp, configOptions* settings);

void setConfigOption(const char * key, const char * val, configOptions * settings);
    
void outputConfigOptions(configOptions* settings, ostream& out);

double* findMedians(double* data, int p, int d, bool hasID);

void quickStats(vector<int> data, double* stats);
void quickStats(vector<double> data, double* stats);


timeval startTimer();
double stopTimer(timeval start);

void testMemoryData(int n, int d);

#endif
