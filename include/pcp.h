/**
 * Author: Kartik Lakhotia
           Sourav Pati
 * Email id: klakhoti@usc.edu
             spati@usc.edu
 * Date: 27-Feb-2018
 *
 * This code implements work optimized propagation blocking with
 * transposed bin graph to reduce cache misses in scatter
 */


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <assert.h>
#include <string.h>
#include <immintrin.h>
#include "../include/gas.h"

intV binWidth = (256*1024)/sizeof(float); //512kB
unsigned int binOffsetBits = (unsigned int)std::log2((double)binWidth); 
intV NUM_BINS = 10000000/binWidth;






#define DEBUG
#undef DEBUG

//////////////////////////////////////////
// performance monitoring via PCM
//////////////////////////////////////////
#define PERF_MON
#undef PERF_MON


//////////////////////////////////////////
// level 2 debugging - asserts enabled
//////////////////////////////////////////
#define DEBUGL2
#undef DEBUGL2


#define ITERTIME
#undef ITERTIME

int NUM_THREADS = std::max(omp_get_max_threads(), 1);
unsigned int MAX_ITER = 10;



template<class graph>
void initialize(graph* G, int argc, char** argv)
{
    G->start = 1;
    G->rounds= 3;
    for (int i = 1; i < argc; i++)
    {  
        if (i + 1 != argc)
        {
            if (strcmp(argv[i], "-s") == 0) // seed node 
            {                 
                G->start = (intV)atoi(argv[i + 1]);    
                i++;    
            }
            if (strcmp(argv[i], "-t") == 0) // num threads 
            {                 
                NUM_THREADS = (unsigned int)atoi(argv[i + 1]);    
                i++;    
            }
            if (strcmp(argv[i], "-iter") == 0)  // num iterations
            {                 
                MAX_ITER = (unsigned int)atoi(argv[i + 1]);    
                i++;    
            }
            if (strcmp(argv[i], "-rounds") == 0) 
            {                 
                G->rounds = (unsigned int)atoi(argv[i + 1]);   
                i++;   
            }
        }
    }
    if (argc < 2)
    {
        printf("Usage : %s <filename> -s <start node(not needed for pr)> -t <numThreads(optional)> -iter <#iterations(optional) -rounds <#rounds(default 3)> \n", argv[0]);
        exit(1);
    }

    omp_set_num_threads(NUM_THREADS);

    //////////////////////////////////////////
    // read csr file
    //////////////////////////////////////////
    if (read_csr(argv[1], G)==-1)
    {
        printf("couldn't read %s\n", argv[1]);
        exit(1);
    }
    
    intV numVerticesPerBin= (G->numVertex/(NUM_THREADS*4));
    numVerticesPerBin = (numVerticesPerBin < binWidth) ? numVerticesPerBin : binWidth;
    intV pow2=1;
    while(pow2<=numVerticesPerBin)
        pow2*=2;
    pow2/=2;
    if(pow2==0) binWidth=4;
    else binWidth = pow2;
    NUM_BINS = (G->numVertex-1)/binWidth + 1;
    G->numBins = NUM_BINS;
    printf("number of partitions %d, size of partitions %d\n", NUM_BINS, binWidth);
    binOffsetBits = (unsigned int)std::log2((double)binWidth);
    //////////////////////////////////////////
    //initialize graph frontier, degree etc.//
    //////////////////////////////////////////
    initGraph (G);
}





template<class type, class graph>
void initBin(graph* G)
{
    //////////////////////////////////////////
    //static work allocation to threads
    //equal no. of edges to all bins
    //////////////////////////////////////////
    G->TD = (partitionData*) malloc (sizeof(partitionData)*NUM_BINS);
    partition(G->TD, G);

    printf("partitioning successful\n");


    //////////////////////////////////////////////////
    //compute storage space required for each bin and
    //offsets for storage in bins for a partition
    //1 column -> 1 gather bin; 1 row -> 1 scatter bin
    //bin[i][j] -> stores what i sends to j
    //////////////////////////////////////////////////
    G->updateBinAddrSize = allocateBinMat<intE>(NUM_BINS, NUM_BINS);
    G->destIdBinAddrSize = allocateBinMat<intE>(NUM_BINS, NUM_BINS);
    G->binFlag = allocateBinMat<bool>(NUM_BINS, NUM_BINS);
    G->activeBins = allocateBinMat<intV>(NUM_BINS, NUM_BINS);

//    struct timespec preStart, preEnd; 
//    float preTime;
//    if( clock_gettime(CLOCK_REALTIME, &preStart) == -1) { perror("clock gettime");}

    //////////////////////////////////////////
    //// transpose and compute offsets ///////
    //////////////////////////////////////////
    #pragma omp parallel for schedule (dynamic, 1)
    for (intV i=0; i<NUM_BINS; i++)
        transposePartition(G, &(G->TD[i]), G->updateBinAddrSize[i], G->destIdBinAddrSize[i]);

    printf("PNG construction successful\n");



//    if( clock_gettime( CLOCK_REALTIME, &preEnd) == -1 ) { perror("clock gettime");}      
//    preTime = (preEnd.tv_sec - preStart.tv_sec)+ (int)(preEnd.tv_nsec - preStart.tv_nsec)/1e9;
//    printf("%s, preprocessing time - %lf\n", argv[1], preTime);

//////////////////////////////////////////
//////////////// BINNING ////////////////
//////////////////////////////////////////

    //////////////////////////////////////////
    ////individual bins to->fro each partition //////
    //////////////////////////////////////////
    G->indUpdateBins = allocateBinMatPtr<type>(NUM_BINS, NUM_BINS);
    G->indDestIdBins = allocateBinMatPtr<intV>(NUM_BINS, NUM_BINS);
    G->sparseDestIdBins = allocateBinMatPtr<intV>(NUM_BINS, NUM_BINS);
#ifdef WEIGHTED
    G->indWeightBins = allocateBinMatPtr<unsigned int>(NUM_BINS, NUM_BINS);
#endif
    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic, 1)
    for (intV i=0; i<NUM_BINS; i++)
    {
        for (intV j=0; j<NUM_BINS; j++)
        {
            G->indUpdateBins[i][j] = new type [G->destIdBinAddrSize[i][j]];
            G->indDestIdBins[i][j] = new intV [G->destIdBinAddrSize[i][j]];
            G->sparseDestIdBins[i][j] = new intV [G->destIdBinAddrSize[i][j]];
#ifdef WEIGHTED
            G->indWeightBins[i][j] = new unsigned int [G->destIdBinAddrSize[i][j]];
#endif
        }
    }

    //pointers for each (i,j) bin for later use //
    G->updateBinPointers = allocateBinMat<intE>(NUM_BINS, NUM_BINS);
    G->destIdBinPointers = allocateBinMat<intE>(NUM_BINS, NUM_BINS);


    #pragma omp parallel for num_threads(NUM_THREADS) schedule (dynamic, 1)
    for (intV i=0; i<NUM_BINS; i++)
    {
#ifdef WEIGHTED
        writeDestIds(G, &G->TD[i], G->indDestIdBins[i], G->indWeightBins[i], G->destIdBinPointers[i]);
#else
        writeDestIds(G, &G->TD[i], G->indDestIdBins[i], G->destIdBinPointers[i]);
#endif
    }

#ifdef DEBUG
    printf("binning complete\n");
#endif 

//////////////////////////////////////////
//////////// BINNING COMPLETE ////////////
//////////////////////////////////////////

}



template<class type, class graph, class userArg>
void scatter_and_gather(graph* G, userArg UA)
{
#ifdef ITERTIME
    float time;
	struct timespec scatterStart, scatterEnd,gatherStart,gatherEnd,start,end;
	float scatterTime = 0.0,gatherTime=0.0;
    struct timespec iterStart,iterEnd;
    float iterTime;
    if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}
#endif
    intV numActiveBins;
///////////////////////////////////////
////Set FLAG For Scatter and Gather////
///////////////////////////////////////


#ifndef DENSE
        G->frontierSize = 0;
#endif
//    printf("\n"); 
//	for(intV i=0;i<NUM_BINS; i++)
//    {
//        for (intV j=0; j<G->TD[i].frontierSize; j++)
//            printf("%d, ", G->TD[i].frontier[j]);
//    }
//    printf("\n"); 
#ifdef ITERTIME
    	if( clock_gettime(CLOCK_REALTIME, &scatterStart) == -1) { perror("clock gettime");}
#endif

        numActiveBins = G->partListPtr;

#ifndef DENSE
        G->partListPtr = 0;
#ifndef ASYNCH
        #pragma omp parallel for num_threads(NUM_THREADS) schedule (dynamic, 1)
        for (intV i=0; i<numActiveBins; i++)
            densityCheck(&G->TD[G->activeScatter[i]]);
#endif
#endif

        #pragma omp parallel for schedule(dynamic,1) num_threads(NUM_THREADS)
        for (intV ptr=0; ptr<numActiveBins; ptr++)
        {
            intV i = G->activeScatter[ptr];
#ifdef ASYNCH
            sgMix(G, &G->TD[i], G->indUpdateBins, G->indDestIdBins, G->sparseDestIdBins, G->TD, G->destIdBinAddrSize, G->destIdBinPointers, G->updateBinPointers, G->scatterDone, UA);
#else
            scatter<type>(G, &G->TD[i], G->indUpdateBins[i], G->sparseDestIdBins[i], G->updateBinPointers[i], G->destIdBinPointers[i], UA);
#endif
        }
        
        #pragma omp parallel for
        for (intV ptr=0; ptr<numActiveBins; ptr++)
            G->scatterDone[G->activeScatter[ptr]] = false; //reset scatter done status of partitions

#ifdef ITERTIME

    	if( clock_gettime( CLOCK_REALTIME, &scatterEnd) == -1 ) { perror("clock gettime");}        


   		scatterTime += (scatterEnd.tv_sec - scatterStart.tv_sec) + (float)(scatterEnd.tv_nsec - scatterStart.tv_nsec)/1e9; 

     	if( clock_gettime(CLOCK_REALTIME, &gatherStart) == -1) { perror("clock gettime");}
#endif   
        #pragma omp parallel for schedule(dynamic,1) num_threads(NUM_THREADS)
        for (intV ptr=0; ptr<G->partListPtr; ptr++)
        {
            intV i=G->activeGather[ptr];
            gather<type>(G, &G->TD[i], G->indUpdateBins, G->indDestIdBins, G->sparseDestIdBins, G->TD, G->destIdBinAddrSize, G->destIdBinPointers, G->updateBinPointers, UA);
            G->activeScatter[ptr] = i;
        }
#ifdef ITERTIME
   	if( clock_gettime( CLOCK_REALTIME, &gatherEnd) == -1 ) { perror("clock gettime");}        
    	gatherTime += (gatherEnd.tv_sec - gatherStart.tv_sec) + (float)(gatherEnd.tv_nsec - gatherStart.tv_nsec)/1e9; 
#endif

#ifdef ITERTIME

        if( clock_gettime(CLOCK_REALTIME, &iterEnd) == -1) { perror("clock gettime");}
        iterTime = (iterEnd.tv_sec - iterStart.tv_sec)+ (int)(iterEnd.tv_nsec - iterStart.tv_nsec)/1e9;
        printf("scatter time= %lf gather time = %lf, total time = %lf\n", scatterTime,gatherTime,scatterTime+gatherTime);
#endif



    

    // free allocated memory//
//    freeMem(&G);
//    free(TD);
//    freeMat<intE>(updateBinAddrSize, NUM_BINS+1);
//    freeMat<intE>(destIdBinAddrSize, NUM_BINS+1);
//    freeMat<intE>(updateBinPointers, NUM_BINS);
//    freeMat<intE>(destIdBinPointers, NUM_BINS);
//    freeMatPtr<type>(indUpdateBins, NUM_BINS,NUM_BINS);
//    freeMatPtr<intV>(indDestIdBins, NUM_BINS,NUM_BINS);
//#ifdef WEIGHTED
//    freeMatPtr<unsigned int>(indWeightBins, NUM_BINS, NUM_BINS);
//#endif
    
}



