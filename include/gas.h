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
#include <omp.h>
#include <assert.h>
#include <vector>
//#include <atomic>
#include "../include/partition.h"

using namespace std;


template<class graph>
void getFrontier(graph* G)
{
    unsigned int* prefix = new unsigned int [(G->partListPtr)+1]();
    prefix[0] = 0;
    for (unsigned int i=0; i<G->partListPtr; i++)
    {
        prefix[i+1] = prefix[i] + G->TD[G->activeScatter[i]].frontierSize;
    } 
    #pragma omp parallel for
    for (unsigned int i=0; i<G->partListPtr; i++)
    {
        for (unsigned int j=0; j<G->TD[G->activeScatter[i]].frontierSize; j++)
            G->frontier[prefix[i]+j] = G->TD[G->activeScatter[i]].frontier[j];
    }
    delete[] prefix;
}


template<class graph>
void resetFrontier(graph* G)
{
    #pragma omp parallel for
    for (unsigned int i=0; i<G->partListPtr; i++)
    {
        for (unsigned int j=0; j<G->TD[G->activeScatter[i]].frontierSize; j++)
            G->inFrontier[G->TD[G->activeScatter[i]].frontier[j]] = false;
        G->TD[G->activeScatter[i]].frontierSize = 0;
    }
    G->frontierSize = 0;
}

template<class graph>
void loadFrontier(graph* G, unsigned int* initFrontier, unsigned int initFrontierSize)
{
    partitionData* allTD = G->TD;
    unsigned int vertexId, pId;
    G->partListPtr = 0;
    for (unsigned int i=0; i<initFrontierSize; i++)
    {
        vertexId = initFrontier[i];
        G->inFrontier[vertexId] = true;
        pId = (vertexId >> binOffsetBits);
        allTD[pId].frontier[allTD[pId].frontierSize++] = vertexId;
        allTD[pId].activeEdges += G->outDeg[vertexId];
        if (!G->flag[pId])
        {
            G->flag[pId] = true;
            G->activeScatter[G->partListPtr++] = pId;
        }
    }
    for (unsigned int i=0; i<initFrontierSize; i++)
    {
        pId = (initFrontier[i] >> binOffsetBits);
        if (G->flag[pId])
            G->flag[pId] = false;
    }
    G->frontierSize = initFrontierSize;
}

template<class graph>
void loadFrontierPar(graph* G, unsigned int* initFrontier, unsigned int initFrontierSize)
{
    partitionData* allTD = G->TD;
    G->partListPtr = 0;
    #pragma omp parallel for
    for (unsigned int i=0; i<initFrontierSize; i++)
    {
        unsigned int vertexId = initFrontier[i];
        G->inFrontier[vertexId] = true;
        unsigned int pId = (vertexId >> binOffsetBits);
        unsigned int ptr = __sync_fetch_and_add(&allTD[pId].frontierSize, 1);
        allTD[pId].frontier[ptr] = vertexId;
        __sync_fetch_and_add(&allTD[pId].activeEdges, G->outDeg[vertexId]);
        if (__sync_bool_compare_and_swap(&G->flag[pId], false, true))
        {
            ptr = G->partListPtr.fetch_add(1);
            G->activeScatter[ptr] = pId;
        }
    }
    #pragma omp parallel for
    for (unsigned int i=0; i<initFrontierSize; i++)
    {
        unsigned int pId = (initFrontier[i] >> binOffsetBits);
        if (G->flag[pId])
            G->flag[pId] = false;
    }
    G->frontierSize = initFrontierSize;
}


template <class graph, class userArg>
void reInitializeSparseFrontier(graph* G, partitionData* TD, userArg UA)
{
    unsigned int trueSize = 0;
    for (unsigned int i=0; i<TD->frontierSize; i++)
    {
        G->inFrontier[TD->frontier[i]] = UA.reInit(TD->frontier[i]);
        if (G->inFrontier[TD->frontier[i]])
        {
            TD->frontier[trueSize++] = TD->frontier[i];
        }
    }
    if ((trueSize > 0) && (__sync_bool_compare_and_swap(&G->flag[TD->tid], false, true)))
    {
       unsigned int listPtr = G->partListPtr.fetch_add(1);
       G->activeGather[listPtr] = TD->tid; 
    }
    TD->frontierSize = trueSize;
}

template <class graph, class userArg>
void reInitializeDenseFrontier(graph* G, partitionData* TD, userArg UA)
{
    for (unsigned int i=TD->startVertex; i<TD->endVertex; i++)
    {
        UA.reInit(i);
    }
}

template <class graph, class userArg>
void filterFrontier(graph*G, partitionData* TD, userArg UA)
{
    unsigned int trueSize = 0;
    for (unsigned int i=0; i<TD->frontierSize; i++)
    {
        G->inFrontier[TD->frontier[i]] = UA.apply(TD->frontier[i]);
        if (G->inFrontier[TD->frontier[i]])
        {
            TD->activeEdges += G->outDeg[TD->frontier[i]];
            TD->frontier[trueSize++] = TD->frontier[i];
        }
    }
    TD->frontierSize = trueSize;
    if (TD->frontierSize > 0)
        G->frontierSize.fetch_add(TD->frontierSize);
}


void densityCheck(partitionData* TD)
{
    TD->isDense = ((28.0 * (float)TD->activeEdges) > ((float)(TD->PNG->numEdges)*10.5 + 2.67*(float)TD->totalEdges + 4.0*(float)NUM_BINS)); 
}


template <class type,class graph, class userArg>
void scatterVC(graph* G, partitionData* TD, type** updateBins, unsigned int** destIdBins, unsigned int* updateBinPointers, unsigned int* destIdBinPointers, userArg UA)
{
    unsigned int destId = 0;
    unsigned int destBin = 0;
    unsigned int vertexId = 0;
    unsigned int cond = 0;
    unsigned int prevBin = 0;
    unsigned int listPtr = 0;
    type userReturn;

#ifdef WEIGHTED
    type weightedVal;
#endif

    for (unsigned int i=0; i<TD->frontierSize; i++)
    {
        vertexId = TD->frontier[i];
        prevBin = NUM_BINS;
        userReturn = UA.scatterFunc(vertexId); 
        for (unsigned int j=G->VI[vertexId]; j<G->VI[vertexId+1]; j++)
        {
            destId = G->EI[j];
            destBin = (destId >> binOffsetBits);

#ifdef WEIGHTED
            weightedVal = UA.applyWeight(userReturn, G->EW[j]); 
            updateBins[destBin][destIdBinPointers[destBin]] = weightedVal; 
            destIdBins[destBin][destIdBinPointers[destBin]++] = destId;
#else
            ///////////////////////////////////
            ///// branch avoiding approach ////
            ///////////////////////////////////
            cond = (destBin != prevBin);
            updateBins[destBin][updateBinPointers[destBin]] = userReturn;
            updateBinPointers[destBin] += cond;
            destId |= (cond << 31);
            destIdBins[destBin][destIdBinPointers[destBin]++] = destId;
            prevBin = destBin;

            ///////////////////////////////////
            //////// branched approach ////////
            ///////////////////////////////////
            //if (destBin!=prevBin)
            //{
            //    updateBins[destBin][updateBinPointers[destBin]++] = userReturn;
            //    destId |= (1 << 31);
            //    prevBin = destBin;
            //}
            //destIdBins[destBin][destIdBinPointers[destBin]++] = destId;

           
#endif
            if (!G->binFlag[TD->tid][destBin])
            {
                G->binFlag[TD->tid][destBin] = true;
                listPtr = G->TD[destBin].binListPtr.fetch_add(1);
                G->activeBins[destBin][listPtr] = TD->tid; 
                if (__sync_bool_compare_and_swap(&G->flag[destBin], false, true))
                {
                   listPtr = G->partListPtr.fetch_add(1);
                   G->activeGather[listPtr] = destBin; 
                }
            }
        }
    }
}


template <class type,class graph, class userArg>
void scatterPC(graph*G, partitionData* TD, type** updateBins, userArg UA)
{
    partitionGraph* PNG = TD->PNG;
    unsigned int pointer;
    unsigned int listPtr;
    type userReturn;
#ifndef DENSE
    for (unsigned int i=0; i<NUM_BINS; i++)
    {
        G->binFlag[TD->tid][i] = true;
        listPtr = G->TD[i].binListPtr.fetch_add(1);
        G->activeBins[i][listPtr] = TD->tid;
        if (__sync_bool_compare_and_swap(&G->flag[i], false, true))
        {
            listPtr = G->partListPtr.fetch_add(1);
            G->activeGather[listPtr] = i; 
        }
    }
#endif
    for (unsigned int i=0; i<PNG->numVertex; i++)
    {
        pointer = 0;
        for (unsigned int j=PNG->VI[i]; j<PNG->VI[i+1]; j++){
            userReturn = UA.scatterFunc(PNG->EI[j]);
            updateBins[i][pointer++] = userReturn;
        }
    } 
}

template <class type, class graph, class userArg>
void scatter(graph*G, partitionData* TD, type** updateBins, unsigned int** destIdBins, unsigned int* updateBinPointers, unsigned int* destIdBinPointers, userArg UA)
{
#ifndef DENSE
    if (TD->isDense)
        scatterPC<type>(G, TD, updateBins, UA);
    else
        scatterVC<type>(G, TD, updateBins, destIdBins, updateBinPointers, destIdBinPointers, UA);
    reInitializeSparseFrontier(G, TD, UA);
#else
    scatterPC<type>(G, TD, updateBins, UA);
    reInitializeDenseFrontier(G,TD,UA);
#endif
}


//////////////////////////////////
/////////// PC GATHER ////////////
//////////////////////////////////
template <class type, class graph, class userArg>
#ifdef WEIGHTED
void gatherPC(graph* G, partitionData* TD, type* updateBin, unsigned int* destIdBin, unsigned int* weightBin, unsigned int binSize, userArg UA)
#else
void gatherPC(graph* G, partitionData* TD, type* updateBin, unsigned int* destIdBin, unsigned int binSize, userArg UA)
#endif
{
    unsigned int destId = 0;
    unsigned int updateBinPointer = MAX_UINT;
    type updateVal;
    bool cond;
    for (unsigned int j=0; j<binSize; j++)
    {
        destId = destIdBin[j];
        updateBinPointer += (destId >> 31);
        destId = destId & MAX_POS;
#ifdef WEIGHTED
        updateVal = UA.applyWeight(updateBin[updateBinPointer], weightBin[j]);
#else
        updateVal = updateBin[updateBinPointer];
#endif
        cond = UA.gatherFunc(updateVal, destId);
        if (!G->inFrontier[destId] && cond)
        {
            TD->frontier[TD->frontierSize++] = destId;
            G->inFrontier[destId] = true;
//            G->inFrontier[destId] = 1;
        }
   } 
}

///////////////////////////////////
///////// DENSE GATHER ///////////
//////////////////////////////////
template <class type, class graph, class userArg>
#ifdef WEIGHTED
void gatherDense(graph* G, partitionData* TD, type* updateBin, unsigned int* destIdBin, unsigned int* weightBin, unsigned int binSize, userArg UA)
#else
void gatherDense(graph* G, partitionData* TD, type* updateBin, unsigned int* destIdBin, unsigned int binSize, userArg UA)
#endif
{
    
    unsigned int destId = 0; 
    unsigned int updateBinPointer = MAX_UINT;
    type updateVal;
    for (unsigned int j=0; j<binSize; j++)
    {
        destId = destIdBin[j];
        updateBinPointer += (destId >> 31);
        destId = destId & MAX_POS;
#ifdef WEIGHTED
        updateVal = UA.applyWeight(updateBin[updateBinPointer], weightBin[j]);
#else
        updateVal = updateBin[updateBinPointer];  
#endif
        UA.gatherFunc(updateVal, destId); 
    } 
}

#ifdef WEIGHTED
///////////////////////////////////
///////// VC GATHER ///////////
//////////////////////////////////
template <class type, class graph, class userArg>
void gatherVC(graph* G, partitionData* TD, type* updateBin, unsigned int* destIdBin, unsigned int binSize, userArg UA)
{
    
    unsigned int destId = 0; 
    unsigned int updateBinPointer = MAX_UINT;
    bool cond;
    type updateVal;
    for (unsigned int j=0; j<binSize; j++)
    {
        destId = destIdBin[j];
        updateVal = updateBin[j];  
        cond = UA.gatherFunc(updateVal, destId);
        if (!G->inFrontier[destId] && cond)
        {
            TD->frontier[TD->frontierSize++] = destId;
            G->inFrontier[destId] = true;
//            G->inFrontier[destId] = 1;
        }
    } 
}
#endif



template <class type, class graph, class userArg>
void gather(graph* G, partitionData* TD, type*** updateBins, unsigned int*** denseDestIdBins, unsigned int*** sparseDestIdBins, partitionData* allTD, unsigned int** destIdBinAddrSize, unsigned int** destIdBinPointers, unsigned int** updateBinPointers, userArg UA)
{
    TD->activeEdges = 0;
#ifdef WEIGHTED
    unsigned int*** weightBin = G->indWeightBins;
#endif
#ifndef DENSE
    G->flag[TD->tid] = false;
#endif

    for (unsigned int ptr=0; ptr<TD->binListPtr; ptr++)
    {
#ifndef DENSE
        unsigned int i = G->activeBins[TD->tid][ptr];
        G->binFlag[i][TD->tid] = false;
#else
        unsigned int i = ptr;
#endif
#ifdef WEIGHTED
#ifndef DENSE
        if (allTD[i].isDense)
            gatherPC<type>(G, TD, updateBins[i][TD->tid], denseDestIdBins[i][TD->tid], weightBin[i][TD->tid], destIdBinAddrSize[i][TD->tid], UA);
        else
            gatherVC<type>(G, TD, updateBins[i][TD->tid], sparseDestIdBins[i][TD->tid], destIdBinPointers[i][TD->tid], UA); 
#else
        gatherDense(G, TD, updateBins[i][TD->tid], denseDestIdBins[i][TD->tid], weightBin[i][TD->tid], destIdBinAddrSize[i][TD->tid], UA);
#endif
#else
#ifndef DENSE
        if (allTD[i].isDense)
            gatherPC<type>(G, TD, updateBins[i][TD->tid], denseDestIdBins[i][TD->tid], destIdBinAddrSize[i][TD->tid], UA);
        else
            gatherPC<type>(G, TD, updateBins[i][TD->tid], sparseDestIdBins[i][TD->tid], destIdBinPointers[i][TD->tid], UA); 
#else
        gatherDense(G, TD, updateBins[i][TD->tid], denseDestIdBins[i][TD->tid], destIdBinAddrSize[i][TD->tid], UA);
#endif
#endif
        destIdBinPointers[i][TD->tid] = 0;
        updateBinPointers[i][TD->tid] = 0;
    }
#ifndef DENSE
    TD->binListPtr = 0;
    filterFrontier(G, TD, UA);
#else
    for (unsigned int i=TD->startVertex; i<TD->endVertex; i++)
        UA.apply(i);
#endif
} 


