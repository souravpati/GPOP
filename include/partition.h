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

//#define MAX_UINT 0xffffffffffffffff
#if defined (HUGE_EDGE) || defined (HUGE_VERTEX)
#define MAX_UINT 0xffffffffffffffff
#else
#define MAX_UINT 0xffffffff
#endif

//#define MAX_NEG 0x80000000
//#define MAX_POS 0x7fffffff
//#define MSB_ROT 31
#ifdef HUGE_VERTEX
#define MAX_NEG 0x8000000000000000
#define MAX_POS 0x7fffffffffffffff
#define MSB_ROT 63
#else
#define MAX_NEG 0x80000000
#define MAX_POS 0x7fffffff
#define MSB_ROT 31
#endif


//////////////////////////////////////////
//partition centric programming variables
//////////////////////////////////////////
extern intV binWidth;
extern unsigned int binOffsetBits; 
extern intV NUM_BINS;


//////////////////////////////////////////
//////// pre-processing functions ////////
//////////////////////////////////////////

//compute equal sized (vertex/edge) partitions
template<class graph>
void partition(partitionData* TD, graph* G)
{
    intV numVertexPerBin = binWidth;
    G->activeScatter = new intV [G->numBins]();
    G->activeGather = new intV [G->numBins]();
    G->partListPtr = 0;
    #pragma omp parallel for
    for (intV i=0; i<G->numBins; i++)
    {
        TD[i].tid = i;
        TD[i].PNG = NULL;
        TD[i].IPG = NULL; //for asynch processing
        TD[i].isDense = false;
        TD[i].startVertex = i*numVertexPerBin;
        TD[i].endVertex = (i+1)*numVertexPerBin;
        TD[i].endVertex = (TD[i].endVertex > G->numVertex) ? G->numVertex : TD[i].endVertex;
        TD[i].frontier = new intV [TD[i].endVertex - TD[i].startVertex];
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
void transposePartition(graph* G, partitionData* TD, intE* updateBinAddrOffset, intE* destIdBinAddrOffset)
{
//    graph* G = TD->G;
    partitionGraph* GSort = new partitionGraph [1];
    intV currBin, prevBin;

    GSort->numVertex = G->numBins;
    GSort->VI = new intE [GSort->numVertex+1]();

    //for intra partition asynch processing
    partitionGraph* IGSort = new partitionGraph [1]; 
    IGSort->numVertex = TD->endVertex - TD->startVertex;
    IGSort->VI = new intE [IGSort->numVertex+1]();

    for (intV i=TD->startVertex; i<TD->endVertex; i++)
    {
        prevBin = G->numBins+1;
        intV intraPartId = i - TD->startVertex;
        for (intE j=G->VI[i]; j<G->VI[i+1]; j++)
        {
            currBin = (G->EI[j] >> binOffsetBits);
            destIdBinAddrOffset[currBin]++;
            if (currBin == prevBin)
                continue;
            GSort->VI[currBin+1]++;
            prevBin = currBin;
            IGSort->VI[intraPartId+1] += (currBin == TD->tid); //for intra partition asynch processing
        }
    }

    for (intV i=0; i<G->numBins; i++)
        updateBinAddrOffset[i] = GSort->VI[i+1];

    for (intV i=0; i<GSort->numVertex; i++)
        GSort->VI[i+1] += GSort->VI[i];

    GSort->numEdges = GSort->VI[GSort->numVertex];
    GSort->EI = new intV [GSort->numEdges];
    intE* binOffset = new intE [G->numBins]();

    //for intra partition asynch processing
    for (intV i=0; i<IGSort->numVertex; i++)
        IGSort->VI[i+1] += IGSort->VI[i];
    IGSort->numEdges = IGSort->VI[IGSort->numVertex];
    IGSort->EI = new intV [IGSort->numEdges];  
#ifdef WEIGHTED
    IGSort->EW = new intV [IGSort->numEdges]; 
#endif
    intE IGoffset = 0;

    for (intV i=TD->startVertex; i<TD->endVertex; i++)
    {
        prevBin = G->numBins;
        for (intE j=G->VI[i]; j<G->VI[i+1]; j++)
        {
            currBin = (G->EI[j] >> binOffsetBits);
            if (currBin == prevBin)
                continue; 
            GSort->EI[GSort->VI[currBin] + (binOffset[currBin]++)] = i;

            if (currBin == TD->tid) //for intra partition asynch processing
            {   
#ifdef WEIGHTED
                IGSort->EW[IGoffset] = G->EW[j];
#endif
                IGSort->EI[IGoffset++] = G->EI[j];
            }

            prevBin = currBin;
        }
    }
    
    TD->PNG = GSort;
    TD->IPG = IGSort; //for intra partition asynch processing


    delete[] binOffset;

}

template<class graph>
#ifdef WEIGHTED
void writeDestIds(graph* G, partitionData* TD, intV** destIdBins, unsigned int** weightBins, intE* destIdBinPointers)
#else
void writeDestIds(graph* G, partitionData* TD, intV** destIdBins, intE* destIdBinPointers)
#endif
{
    intV destId = 0;
    intV destBin = 0;
    intV prevBin = 0;
#ifdef DEBUGL2
    for (intV i=0; i<G->numBins; i++)
        assert(destIdBinPointers[i] == 0);
#endif
    for (intV i=TD->startVertex; i<TD->endVertex; i++)
    {
        prevBin = G->numBins;
        for (intE j=G->VI[i]; j<G->VI[i+1]; j++)
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
    for (intV i=0; i<G->numBins; i++)
        destIdBinPointers[i] = 0;
}







//////////////////////////////////////////
//////////// memory allocation ///////////
//////////////////////////////////////////
//allocate BIN x BIN space for offsets, pointers
template <class T> T** allocateBinMat (intV numRows, intV numCols)
{
    T** pointerMat;
    pointerMat = new T* [numRows];
    for (intV i=0; i<numRows; i++)
        pointerMat[i] = new T [numCols]();
    return pointerMat;
}

//allocate BIN x BIN space for pointers
template <class T> T*** allocateBinMatPtr (intV numRows, intV numCols)
{
    T*** pointerMat;
    pointerMat = new T** [numRows];
    for (intV i=0; i<numRows; i++)
        pointerMat[i] = new T* [numCols];
    return pointerMat;
}

//////////////////////////////////////////
////////////free the memory //////////////
//////////////////////////////////////////
template <class T> void freeMat (T** mat, intV numRows)
{
    for (intV i=0; i<numRows; i++)
        delete[] mat[i];
    delete[] mat; 
}

template <class T> void freeMatPtr (T*** mat, intV numRows, intV numCols)
{
    for (intV i=0; i<numRows; i++)
        for (intV j=0; j<numCols; j++)
            delete[] mat[i][j];
    for (intV i=0; i<numRows; i++)
        delete[] mat[i];
    delete[] mat; 
}




