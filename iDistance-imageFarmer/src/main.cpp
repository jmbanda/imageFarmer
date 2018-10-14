/*////////////////////////////////////////////////////////////////////
Main.cpp : The main driver file for iDistance*.

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
#include <sstream>
#include "string.h"

#include "helper.h"
#include "idistance.h"

#define CONSOLE_PRINT false //print file output to console as well
#define SS_CHECK false //use this to verify iDist results
   // this circumvents normal stat outputs for SS and instead outputs when
   // there is a mismatch between the results of both methods
   // DO NOT USE -SS with this check because BOTH results will be SS then!

#define IFOUTPUT true //special output file for imageFarmer

//command line args defined
int num_Args = 7;
const char* def_Args[] = {"-C", "-SS", "-B", "-L", "-QP", "-QR", "-QN"};
enum {LOADCONFIG, SEQSEARCH, BUILDTREE, LOADTREE, 
    QUERYPOINT, QUERYRANGE, QUERYKNN};

using namespace std;

//print the correct program arguments to the console
void failedArgs()
{
    cout<<endl;
    cout<<"Error loading program arguments. Correct format is:"<<endl;
    cout<<"./iDist {-C \"configfile\" {additional files}} {-SS}"<<endl; 
    cout<<"  {-B \"datasetfile\" | -L \"treefile\"}"<<endl; 
    cout<<"  {-QP \"qfile\" | -QR \"qfile\" | -QN \"qfile\" kval}"<<endl;
    cout<<endl;
    cout<<"See readme.txt for more information."<<endl;
    cout<<endl;
}

//prints the command line arguments to the console
void dumpArgs(int argc, char* argv[])
{
    cout << "argc = " << argc << endl; 
    for(int i = 0; i < argc; i++)
    {
        cout << "argv[" << i << "] = '" << argv[i] << "'" << endl;
    }
  
}

/*************************************************
 * Runs through the given command line args and performs
 * the required functions. Consider this the master function
 * which dictates what happens. Note that it requires 
 * the correct order: {config} {tree} {query} args
 */
int runArgs(int argc, char* argv[])
{ 
    bool treeDone = false; //ensures tree before query

   	iDistance iDist;
    
    configOptions settings; //declare our default configuration
    settings.output_dir = "results/"; //and default results folder

    string statusFile = "status.txt"; //append ALL results output
    string treeStatsFile = "stats_trees.txt"; //raw stats of ALL results
    string queryStatsFile = "stats_queries.txt"; //raw stats of ALL results
    string queryStatsFile_IF = "stats_IF.txt"; //raw stats of ALL results for imageFarmer
    
        
    int num_points; //in dataset
    int num_dims; //in dataset, these used by query stats output
    int argI = 1; //init arg index, 0 = program name
    int arg = -1; //enum value for arg switch
    ostringstream ss; //string stream used for buffering up output
               
    ss << "\niDistance* v0.9\n\n";
    
    ss << "Program call:\n";
    for(int i = 0; i < argc; i++)
    {
      ss << argv[i] << " ";
    }
    ss << "\n\n";

    //useful to see when scripting many runs
    cout << ss.str() << endl;
               
    string argStr;
    while(argI < argc) //loop over all command line args
    {
    
        //get the command arg enum value
        argStr = argv[argI];
        arg = stringToEnum(stringToUpper(argStr.c_str()), def_Args, num_Args);
      
        //ensure we got a valid one, and enough remaining args for it
        if( (arg == -1) ||
          ((arg == QUERYKNN) && (argc - argI < 2)) ||
          ((arg != QUERYKNN) && (arg != SEQSEARCH) && (argc - argI < 1)) )
        {
            failedArgs();
            return 1;
        }
        
        //if query arg, make sure tree ready
        if((!treeDone) && 
            ((arg == QUERYPOINT) || (arg == QUERYRANGE) || (arg == QUERYKNN)))
        {
            cout << "No Tree has been established yet!" << endl;
            failedArgs();
            return 1;
        }
        
        switch(arg)
        {
        ////////////////////////////////////////////////////////////////
            case LOADCONFIG:
            {
                argI++; //inc and get next arg
                const char* argCFile = argv[argI];
                ss << "Loading config file: '" << argCFile << "'" << endl;
            
                //load config file and get settings
                int configLoad = readConfigFile(argCFile, &settings);
                
                if(configLoad == 0)
                {
                    outputConfigOptions(&settings, ss);
                  
                    ss << endl;
                    
                    if(settings.refs_build == REF_FILE)
                    {
                    	//instead, lets read the very next command arg for the filename
                    	argI++;
                    	const char* argRefFile = argv[argI];
                    	
                    	ss << "  Loading reference points from file: " << argRefFile << endl;
                        int* refFileInfo = getFileInfo(argRefFile);
                        
                        double* data = readDataFile(argRefFile, refFileInfo[0], refFileInfo[1]);
                        
                        if(data != NULL)
                        {
                            iDist.setReferencePoints(data, refFileInfo[0]);
                        }
                        else
                        {
                        	cout << "Error loading ref file!" << endl;
                        	return 1;
                        }
                    	
                    }
                    
                    
                    if(settings.refs_assign == ASSIGN_FILE)
                    {
                    	argI++;
                    	const char* argAssignFile = argv[argI];
                    	
                    	ss << "  Loading point assignments from file: " << argAssignFile << endl;
                        int* assignFileInfo = getFileInfo(argAssignFile);
                        
                        double* data = readDataFile(argAssignFile, assignFileInfo[0], assignFileInfo[1]);
                        
                        if(data != NULL)
                        {
                            iDist.setPointAssignment(data, assignFileInfo[0]);
                        }
                        else
                        {
                        	cout << "Error loading assign file!" << endl;
                        	return 1;
                        }
                    	
                    }
                       
                      /* 
                                        
                    if(settings.algo_version == EXT1 && 
                    	settings.splits_build == SPLIT_FILE)
                    {
                    	//instead, lets read the very next command arg for the filename
                    	argI++;
                    	const char* argSplitFile = argv[argI];
                      
                    	ss << "  Loading medians from file: " << argSplitFile << endl;
                        int* splitFileInfo = getFileInfo(argSplitFile);
                        
                        double* data = readDataFile(argSplitFile, splitFileInfo[0], splitFileInfo[1]);
                
                        if(data != NULL)
                        {
                            iDist.setSplitPoints(data, splitFileInfo[1]);
                        }
                        else
                        {
                        	cout << "Error loading splits file!" << endl;
                        	return 1;
                        }

                    }
                    
                    */
                    
                    if(settings.output_mode == OUTPUT_DIR)
                    {
                    	argI++;
                    	const char* argOutputDir = argv[argI];
                    	
                    	ss << "  Saving program output to directory: " << argOutputDir << endl;
                        
                      settings.output_dir = argOutputDir;
                      
                      
                      statusFile = settings.output_dir + statusFile;
                      treeStatsFile = settings.output_dir + treeStatsFile;
                      queryStatsFile = settings.output_dir + queryStatsFile;
                      queryStatsFile_IF = settings.output_dir + queryStatsFile_IF;
                       
              
                    }

                    iDist.setOptions(settings);
                }
                else
                {
                    ss << "  Error in configuration file!" << endl;
                }
                
                ss << endl << endl;
                
                string mkdir = "mkdir -p ";
                mkdir = mkdir + settings.output_dir;
                system(mkdir.c_str());
                
            }
            break;
        ////////////////////////////////////////////////////////////////
            case BUILDTREE: 
            {
                argI++; //inc and get next arg
                const char* argTFile = argv[argI];
     
                ss << "Building index on data file: '" << argTFile << "'" << endl;
                
                //NOTE: that the data contains the id as the first value in each point
                int* fileInfo = getFileInfo(argTFile);
                ss << "  total points = " << fileInfo[0] << endl;
                ss << "  total dims = " << fileInfo[1]-1 << endl; 
                
                num_points = fileInfo[0];
                num_dims = fileInfo[1]-1;
                
                double* data = readDataFile(argTFile, fileInfo[0], fileInfo[1]);

                if(data == NULL) //error loading data
                {
                    cout << "  Error loading data!" << endl;
                    return 1;
                }
                
                //load data into iDistance (and data IDs)
                iDist.setDataWithIDs(data, fileInfo[0], fileInfo[1]-1);
                
                iDist.initOptions(); //after data loaded
                
                timeval timer = startTimer(); //start timer
                iDist.buildIndexValues();
                double calcTime = stopTimer(timer); //time taken to calculate iDist values
        
                iDist.buildTree();
                double indexTime = stopTimer(timer); //calcTime + time to load B+ tree
                
                double buildTime = (indexTime - calcTime); //just the time taken to load B+ tree
                
                stats treeStats;
                iDist.getStats(&treeStats);
                treeStats.buildTime = buildTime;
                treeStats.indexTime = indexTime;
                treeStats.points = fileInfo[0];
                treeStats.dims = fileInfo[1]-1;
                
                //get the base file name for printing out all new files
                //this strips the extension and all folder structures
                string baseFileStr = argTFile;
                baseFileStr = trimToName(baseFileStr);
                string resultsFileStr = settings.output_dir + baseFileStr;                
                const char* baseFile = baseFileStr.c_str();
                const char* resultsFile = resultsFileStr.c_str();
                
                //print out tree segments as .segments file
                iDist.PrintNodeRange(resultsFile);
                //iDist.PrintNodeRange(cout);
                
                //save the build as .tree and .data files
                iDist.saveTree(resultsFile);
                
                printTreeStats(resultsFile, &treeStats);
                printTreeStats(ss, &treeStats);
                ss << endl << endl;
                
                ostringstream statsStr;
                printTreeStats_CSV(statsStr, &treeStats);
                printOut(treeStatsFile.c_str(), statsStr.str(), false);
                
                treeDone = true;
          
                delete [] data;
                delete [] fileInfo;
            }
            break;
        ////////////////////////////////////////////////////////////////
            case LOADTREE:
            {
                argI++; //inc and get next arg
                const char * argTFile = argv[argI];  
                ss << "Loading index from: '" << argTFile << "'" << endl;
                
                iDist.loadTree(argTFile);
                
                num_points = iDist.getNumberOfPoints();
                num_dims = iDist.getNumberOfDimensions();
                
                ss << "  Index and data files loaded!" << endl;
                ss << endl;
                //printOut(statusFile, ss.str(), CONSOLE_PRINT);
                treeDone = true;
                
            }
            break;
        ////////////////////////////////////////////////////////////////
            case SEQSEARCH:
            {
                //set the sequential search flag in iDistance class
                
                iDist.setSequentialSearch(1);
                ss << "\n---- Performing Sequential Scan ----\n" << endl;
                

            }
            break;
        ////////////////////////////////////////////////////////////////
            case QUERYPOINT:
            {
                argI++; //inc and get next arg
                const char * argQFile = argv[argI];
                ss << "Loading query point(s) from file: '" << argQFile << "'" << endl;
                
                string baseFileStr = argQFile;
                baseFileStr = trimToName(baseFileStr);
                string resultsFileStr = settings.output_dir + baseFileStr + ".results";
                const char* baseFile = baseFileStr.c_str();
                const char* resultsFile = resultsFileStr.c_str();
                
                int* fileInfo = getFileInfo(argQFile);
                double* data = readDataFile(argQFile, fileInfo[0], fileInfo[1]);
                
                if(data == NULL)
                {
                    cout << "  Error loading query file!" << endl;
                    return 1;
                }
                
                ostringstream statsStr; //separate ss for csv stats
                
                for(int i = 0; i < fileInfo[0]; i++)
                {
                    ss << "  Query Point " << i << " (ID "<< data[i*fileInfo[1]] << ")" << endl;

                    //RESET stats before each query
                    iDist.resetNodeCount();
                    iDist.resetCandidateCount();
                    iDist.resetPartitionCount();
                    timeval qtimer = startTimer(); //and start timer 
                    
                    vector<int> results_ind = iDist.QueryPoint(&data[i*fileInfo[1]+1]); //move to start of data
                    
                    double queryT = stopTimer(qtimer); //return elapsed time as a double in ms	

                    int numNodes = iDist.getNodeCount();  //get node count
                    int numCandidates = iDist.getCandidateCount();
                    int numPartitions = iDist.getPartitionCount(); 
                    vector<int> results = iDist.getIDs(results_ind); //get real IDs
                    
                    vector<double> qstats;
                    qstats.push_back(data[i*fileInfo[1]]);
                    qstats.push_back(numPartitions);
                    qstats.push_back(numCandidates);
                    qstats.push_back(results.size());
                    qstats.push_back(numNodes);
                    qstats.push_back(queryT);
                    
                    
                    if(SS_CHECK)
                    {
                        iDist.setSequentialSearch(1);
                        vector<int> results_ind_ss = iDist.QueryPoint(&data[i*fileInfo[1]+1]); //move to start of data
                        iDist.setSequentialSearch(0);
                        
                        vector<int> results_ss = iDist.getIDs(results_ind_ss);
                        
                        for(int i = 0; i < results_ss.size(); i++)
                        {
                            if (results_ss[i] != results[i])
                            {
                                ss << "  MISMATCH!\n\nSS:\n";
                                printResultsQP(ss, results_ss, -1, -1);
                                i = results_ss.size();
                                ss << "iDist:\n";
                            }
                                
                        }
                    
                    }

                    printResultsQP(ss, results, numNodes, queryT);
                    printResultsQP_CSV(statsStr, qstats, results);

                } 
             
                ss << endl << endl;
                
                printOut(resultsFile, ss.str(), false);
                printOut(queryStatsFile.c_str(), statsStr.str(), false);    

                delete [] data;
                delete [] fileInfo;
            }
            break;
        ////////////////////////////////////////////////////////////////
            case QUERYRANGE:
            {
              
                argI++; //inc and get next arg
                const char * argQFile = argv[argI];
                ss << "Loading query range(s) from: '" << argQFile << "'" << endl;

                string baseFileStr = argQFile;
                baseFileStr = trimToName(baseFileStr);
                string resultsFileStr = settings.output_dir + baseFileStr + ".results";
                const char* baseFile = baseFileStr.c_str();
                const char* resultsFile = resultsFileStr.c_str();

                int* fileInfo = getFileInfo(argQFile);
                double* data = readDataFile(argQFile, fileInfo[0], fileInfo[1]);
                
                if(data == NULL)
                {
                    cout << "  Error loading query file!" << endl;
                    return 1;
                }
            
                ostringstream statsStr;
                int range = 1;
                for(int i = 0; i < fileInfo[0]; i=i+2) //increment two lines at a time
                {

                    ss << "  Range Query " << range << endl;

                    //RESET stats before each query
                    iDist.resetNodeCount();
                    iDist.resetCandidateCount();
                    iDist.resetPartitionCount();
                    timeval qtimer = startTimer(); // and start timer
                    
                    vector<int> results_ind = iDist.QueryRange(&data[i*fileInfo[1]], &data[(i+1)*fileInfo[1]]); //move to start of data
                    
                    double queryT = stopTimer(qtimer); 
                    
                    int numNodes = iDist.getNodeCount();
                    int numCandidates = iDist.getCandidateCount();
                    int numPartitions = iDist.getPartitionCount();
                    
                    vector<int> results = iDist.getIDs(results_ind);
                                        
                    vector<double> qstats;
                    qstats.push_back(range); //no query ID, just range
                    qstats.push_back(numPartitions);
                    qstats.push_back(numCandidates);
                    qstats.push_back(results.size());
                    qstats.push_back(numNodes);
                    qstats.push_back(queryT);
                    
                    //cout << "HERE: " << results.size() << " : " << numNodes << endl;
                    
                    if(SS_CHECK)
                    {
                        iDist.setSequentialSearch(1);
                        vector<int> results_ind_ss = iDist.QueryRange(&data[i*fileInfo[1]], &data[(i+1)*fileInfo[1]]); //move to start of data
                        iDist.setSequentialSearch(0);
                        
                        vector<int> results_ss = iDist.getIDs(results_ind_ss);
                        
                        for(int i = 0; i < results_ss.size(); i++)
                        {
                            if (results_ss[i] != results[i])
                            {
                                ss << "  MISMATCH!\nSS:\n";
                                printResultsQR(ss, results_ss, -1, -1);
                                i = results_ss.size();
                                ss << "iDist:\n";
                            }
                                   
                        }
                    
                    }
                    
                    printResultsQR(ss, results, numNodes, queryT);
                    printResultsQR_CSV(statsStr, qstats, results);
                    
                    range++;
                }

                ss << endl;
                printOut(resultsFile, ss.str(), false);
                printOut(queryStatsFile.c_str(), statsStr.str(), false);  
                
                delete [] data;
                delete [] fileInfo;
            
            }
            break;
        ////////////////////////////////////////////////////////////////
            case QUERYKNN:
            {
                argI++; //inc and get next arg
                const char * argQFile = argv[argI];
                argI++; //inc and get next arg again
                int k = atoi(argv[argI]);           
                
                ss << "Loading knn query point(s) from file: '" << argQFile << "'" << endl;
                ss << "  searching for " << k << " nearest neighbors.\n" << endl;
                
                string baseFileStr = argQFile;
                baseFileStr = trimToName(baseFileStr);
                string resultsFileStr = settings.output_dir + baseFileStr + "_" + argv[argI] + ".results";
                const char* baseFile = baseFileStr.c_str();
                const char* resultsFile = resultsFileStr.c_str();
                
                int* fileInfo = getFileInfo(argQFile);
                double* data = readDataFile(argQFile, fileInfo[0], fileInfo[1]);
                
                if(data == NULL)
                {
                    cout << "  Error loading query file!" << endl;
                    return 1;
                }

                ostringstream statsStr;
                for(int i = 0; i < fileInfo[0]; i++)
                {
                    ss << "  KNN Query " << i << " (ID " << data[i*fileInfo[1]] << ")" << endl;

                    //RESET stats before each query
                    iDist.resetNodeCount();
                    iDist.resetCandidateCount();
                    iDist.resetPartitionCount();
                    timeval qtimer = startTimer(); //and start timer
                    
                    vector<int> results_ind = iDist.QueryKNN(&data[i*fileInfo[1]+1], k); //move to start of data
                    
                    double queryT = stopTimer(qtimer); //return elapsed time as a double
                    
                    vector<int> results = iDist.getIDs(results_ind);
                    vector<double> dists = iDist.getKNN_dists();
                    vector<int> nodes = iDist.getKNN_nodes();
                    vector<int> candidates = iDist.getKNN_candidates();
                    
                    int numPartitions = iDist.getPartitionCount();
                    int numNodes = iDist.getNodeCount();  //get node count
                    int numCandidates = iDist.getCandidateCount();
                    double avgDist = iDist.AverageDist(dists);

                    vector<double> qstats;
                    qstats.push_back(data[i*fileInfo[1]]);
                    qstats.push_back(numPartitions);
                    qstats.push_back(numCandidates);
                    qstats.push_back(results.size());
                    qstats.push_back(numNodes);
                    qstats.push_back(queryT);
                    
                    if(SS_CHECK)
                    {
                        iDist.setSequentialSearch(1);
                        vector<int> results_ind_ss = iDist.QueryKNN(&data[i*fileInfo[1]+1], k); //move to start of data
                        iDist.setSequentialSearch(0);
                        
                        vector<int> results_ss = iDist.getIDs(results_ind_ss);
                        vector<double> dists_ss;
                        vector<int> nodes_ss;
                        
                        for(int i = 0; i < results_ss.size(); i++)
                        {
                            if (results_ss[i] != results[i])
                            {                      
                                ss << "  MISMATCH!\nSS:\n";
                                printResultsQN(ss, results_ss, dists_ss, nodes_ss, -1, -1);
                                i = results_ss.size();
                                ss << "iDist:\n";
                            }     
                                
                        }
                    
                    }
                    
                    printResultsQN(ss, results, dists, nodes, avgDist, queryT);
                    printResultsQN_CSV(statsStr, qstats, results, nodes, candidates, dists);
                     
   					if(IFOUTPUT)
					{
						ostringstream ifStr;
						ifStr.str("");
						printResultsQN_CSV_IF(ifStr, qstats, results, nodes, candidates, dists);
						printOut(queryStatsFile_IF.c_str(), ifStr.str(), false);
                    
					} 
                      
                }
            
                ss << endl;
                //printOut(statusFile, ss.str(), CONSOLE_PRINT);
                printOut(resultsFile, ss.str(), false);
                printOut(queryStatsFile.c_str(), statsStr.str(), false);
                    
                delete [] data;
                delete [] fileInfo;
            }
            break; 
         
        } //end switch statement on args
                
        //write out status string to file (and maybe console)
        printOut(statusFile.c_str(), ss.str(), CONSOLE_PRINT);
        ss.str(""); //clear string
              
        argI++; //increment arg pointer, and print separation in status file

    } //end while loop on command args
    
    printOut(statusFile.c_str(), 
            "\n----------------------------------------\n\n", false);
    
    return 0;
}







 









/////////////////////////////////////////////////////////////////////
//////////  MAIN PROGRAM RUN

int main(int argc, char* argv[])
{
    ////////////////////////////////
    //Testing area -- simply circumvent rest of program
    /////////////////////////////////
    
    
    // ./extractMedians utility code
    // inputs a data filename, outputs a median file
    
    /*
    ostringstream ss;
    const char* dataFile = argv[1]; //first and only arg
    cout << "file: " << dataFile << endl;
    //NOTE: that the data contains the id as the first value in each point
	int* fileInfo = getFileInfo(dataFile);
	cout << "  total points = " << fileInfo[0] << endl;
	cout << "  total dims = " << fileInfo[1]-1 << endl; 
    
    double* data = readDataFile(dataFile, fileInfo[0], fileInfo[1]);

    double* medians = findMedians(data, fileInfo[0], fileInfo[1], true);
    
    for(int i = 1; i < fileInfo[1]; i++)
	{
		ss << medians[i]  << endl;
	}
	
	string baseFileStr = dataFile;
    baseFileStr = trimToName(baseFileStr);
    string medianFileStr = baseFileStr + "_" + "medians.txt";
    const char* medianFile = medianFileStr.c_str();
                
	printOut(medianFile, ss.str(), false);
    
	return 0;
	
	*/
    ///////////////////////////////////
    // END TESTING AREA
    ///////////////////////////////
  
    if(argc <= 2) //has to be at least three args (prog name first)
    {
        failedArgs();
        return 1;
    }
    
    timeval start = startTimer();
    int retStatus = runArgs(argc, argv);
    double runTime = stopTimer(start);
    runTime = runTime / 1000.0; // in seconds
    
    cout << endl << "Program" << 
      (retStatus == 0 ? 
        " finished successfully!" : " terminated with errors!") << endl << endl;
    
    cout << "Total elapsed time: " << runTime << " seconds." << endl << endl;
       
    return 0;

}

