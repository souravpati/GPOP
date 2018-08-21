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
#include <math.h>
#include "../include/graph.h"

#define MAX_NEG 0x80000000
#define MAX_POS 0x7fffffff
#define MAX_UINT 0xffffffff


//////////////////////////////////////////
//partition centric programming variables
//////////////////////////////////////////
extern unsigned int binWidth;
extern unsigned int binOffsetBits; 
extern unsigned int NUM_BINS;


//////////////////////////////////////////
//////// pre-processing functions ////////
//////////////////////////////////////////

//compute equal sized (vertex/edge) partitions
template<class graph>
void partition(partitionData* TD, graph* G)
{
    unsigned int numVertexPerBin = binWidth;
    int vcount = 0;
    G->activeScatter = new unsigned int [G->numBins]();
    G->activeGather = new unsigned int [G->numBins]();
    G->partListPtr = 0;
    #pragma omp parallel for
    for (int i=0; i<G->numBins; i++)
    {
        TD[i].tid = i;
        TD[i].PNG = NULL;
        TD[i].isDense = false;
        TD[i].startVertex = i*numVertexPerBin;
        TD[i].endVertex = (i+1)*numVertexPerBin;
        TD[i].endVertex = (TD[i].endVertex > G->numVertex) ? G->numVertex : TD[i].endVertex;
        TD[i].frontier = new unsigned int [TD[i].endVertex - TD[i].startVertex];
        TD[i].frontierSize = 0;
        TD[i].activeEdges = 0;
        TD[i].totalEdges = G->VI[TD[i].endVertex] - G->VI[TD[i].startVertex];
        TD[i].binListPtr = 0;
#ifdef DENSE
        G->activeScatter[i] = i;
        G->activeGather[i] = i;
        G->partListPtr = G->numBins;
        TD[i].binListPtr = G->numBins;
#endif
    }

}

//transpose the graph to sort on destination
template<class graph>
void transposePartition(graph* G, partitionData* TD, unsigned int* updateBinAddrOffset, unsigned int* destIdBinAddrOffset)
{
//    graph* G = TD->G;
    partitionGraph* GSort = new partitionGraph [1];
    unsigned int currBin, prevBin;

    GSort->numVertex = G->numBins;
    GSort->VI = new unsigned int [GSort->numVertex+1]();

    for (unsigned int i=TD->startVertex; i<TD->endVertex; i++)
    {
        prevBin = G->numBins+1;
        for (unsigned int j=G->VI[i]; j<G->VI[i+1]; j++)
        {
            currBin = (G->EI[j] >> binOffsetBits);
            destIdBinAddrOffset[currBin]++;
            if (currBin == prevBin)
                continue;
            GSort->VI[currBin+1]++;
            prevBin = currBin;
        }
    }

    for (unsigned int i=0; i<G->numBins; i++)
        updateBinAddrOffset[i] = GSort->VI[i+1];

    for (unsigned int i=0; i<GSort->numVertex; i++)
        GSort->VI[i+1] += GSort->VI[i];

    GSort->numEdges = GSort->VI[GSort->numVertex];
    GSort->EI = new unsigned int [GSort->numEdges];

    unsigned int* binOffset = new unsigned int [G->numBins]();
    for (unsigned int i=TD->startVertex; i<TD->endVertex; i++)
    {
        prevBin = G->numBins;
        for (unsigned int j=G->VI[i]; j<G->VI[i+1]; j++)
        {
            currBin = (G->EI[j] >> binOffsetBits);
            if (currBin == prevBin)
                continue; 
            
            GSort->EI[GSort->VI[currBin] + (binOffset[currBin]++)] = i;
            prevBin = currBin;
        }
    }
    
    TD->PNG = GSort;


    delete[] binOffset;

}

template<class graph>
#ifdef WEIGHTED
void writeDestIds(graph* G, partitionData* TD, unsigned int** destIdBins, unsigned int** weightBins, unsigned int* destIdBinPointers)
#else
void writeDestIds(graph* G, partitionData* TD, unsigned int** destIdBins, unsigned int* destIdBinPointers)
#endif
{
    unsigned int destId = 0;
    unsigned int destBin = 0;
    unsigned int prevBin = 0;
#ifdef DEBUGL2
    for (unsigned int i=0; i<G->numBins; i++)
        assert(destIdBinPointers[i] == 0);
#endif
    for (unsigned int i=TD->startVertex; i<TD->endVertex; i++)
    {
        prevBin = G->numBins;
        for (unsigned int j=G->VI[i]; j<G->VI[i+1]; j++)
        {
            destId = G->EI[j];
            destBin = (destId >> binOffsetBits);
            if (destBin != prevBin)
            {
                destId |= MAX_NEG;
                prevBin = destBin;
            } 
#ifdef WEIGHTED
            weightBins[destBin][destIdBinPointers[destBin]] = G->EW[j];
#endif
            destIdBins[destBin][destIdBinPointers[destBin]++] = destId;
        }
    }
    for (unsigned int i=0; i<G->numBins; i++)
        destIdBinPointers[i] = 0;
}







//////////////////////////////////////////
//////////// memory allocation ///////////
//////////////////////////////////////////
//allocate BIN x BIN space for offsets, pointers
template <class T> T** allocateBinMat (unsigned int numRows, unsigned int numCols)
{
    T** pointerMat;
    pointerMat = new T* [numRows];
    for (unsigned int i=0; i<numRows; i++)
        pointerMat[i] = new T [numCols]();
    return pointerMat;
}

//allocate BIN x BIN space for pointers
template <class T> T*** allocateBinMatPtr (unsigned int numRows, unsigned int numCols)
{
    T*** pointerMat;
    pointerMat = new T** [numRows];
    for (unsigned int i=0; i<numRows; i++)
        pointerMat[i] = new T* [numCols];
    return pointerMat;
}

//////////////////////////////////////////
////////////free the memory //////////////
//////////////////////////////////////////
template <class T> void freeMat (T** mat, unsigned int numRows)
{
    for (unsigned int i=0; i<numRows; i++)
        delete[] mat[i];
    delete[] mat; 
}

template <class T> void freeMatPtr (T*** mat, unsigned int numRows, unsigned int numCols)
{
    for (unsigned int i=0; i<numRows; i++)
        for (unsigned int j=0; j<numCols; j++)
            delete[] mat[i][j];
    for (unsigned int i=0; i<numRows; i++)
        delete[] mat[i];
    delete[] mat; 
}




