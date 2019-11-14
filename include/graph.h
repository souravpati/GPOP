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
#include <atomic>
#include <boost/dynamic_bitset.hpp>

#if defined (HUGE_EDGE) || defined (HUGE_VERTEX)
typedef unsigned long long int intE;
#else
typedef unsigned int intE;
#endif
//typedef unsigned long long int intE;
//typedef unsigned int intV;

#ifdef HUGE_VERTEX
typedef unsigned long long int intV;
#else
typedef unsigned int intV;
#endif

//////////////////////////////////////////
//partition centric programming data types
//////////////////////////////////////////
typedef struct partitionGraph 
{
    intV numVertex;
    intE numEdges;
    intE* VI;
    intV* EI;
    unsigned int* EW;
    intV* outDeg;
} partitionGraph;

typedef struct partitionData
{
    intV tid; //partition ID
    intV startVertex; 
    intV endVertex;
    partitionGraph* PNG; //for dense scatter
    partitionGraph* IPG; //for intra-partition asynch processing
    intV* frontier; //active vertex list within the partition
    intV frontierSize; 
    bool isDense; //true if partition should be scattered in dense mode
    intE totalEdges;
    intE activeEdges;
    std::atomic<intV> binListPtr; //number of bins with messages to be gathered??
}partitionData;

template<class type>
struct graph 
{
    intV numBins; // # partitions
    intV numVertex;
    intE numEdges;
    intE* VI; //vertex offset in CSR
    intV* EI; //edge array in CSR
    unsigned int* EW; //edge weights
    intV* frontier; //
    intV start; //starting vertex??
    unsigned int rounds;
    bool* inFrontier; //array to indicate if a vertex is in frontier
    std::atomic<intV> frontierSize;
    partitionData* TD; //array of partitions
    intV* outDeg;
    intV* inDeg;
    bool* scatterDone; // to indicate if a partition is done scattering
    bool* flag; //to indicate if a partition has a message or not
    bool** binFlag; //to indicate if a bin has any message or not
    intE** updateBinAddrSize; //2D array for no. of updates 
    intE** destIdBinAddrSize; //2D array for no. of dest IDs
    intE** updateBinPointers; //running address (counters) pointing to next write location
    intE** destIdBinPointers; 
    unsigned int*** indWeightBins; //weight bins
    intV*** indDestIdBins; //dest ID bins (dense, fixed order)
    intV*** sparseDestIdBins; //sparse dest ID bins (empty, to be filled on scatter)
    type*** indUpdateBins;   //update bins (same for dense and sparse)
    intV* activeScatter; //list of partitions with active vertices
    intV* activeGather; //list of partitions who have received a message
    std::atomic<intV> partListPtr; //# partitions to be processed (scatter or gather)
    intV** activeBins; //array of bins with messages
};


template<class graph>
int read_csr (char* filename, graph* G)
{
    FILE* graphFile = fopen(filename, "rb");
    if (graphFile == NULL)
    {
        fputs("file error", stderr);
        return -1;
    }
    fread (&(G->numVertex), sizeof(intV), 1, graphFile);
    std::cout << G->numVertex << std::endl;
    
    fread (&(G->numEdges), sizeof(intE), 1, graphFile);
    std::cout << G->numEdges << std::endl;

    G->VI = new intE[G->numVertex+1];
    fread (G->VI, sizeof(intE), G->numVertex, graphFile);
    std::cout << "read offsets" << std::endl;
    if (feof(graphFile))
    {
        delete[] G->VI;
        printf("unexpected end of file while reading vertices\n");
        return -1;
    }
    else if (ferror(graphFile))
    {
        delete[] G->VI;
        printf("error reading file\n");
        return -1;
    }
    G->VI[G->numVertex] = G->numEdges;

    G->EI = new intV[G->numEdges];
    fread (G->EI, sizeof(intV), G->numEdges, graphFile);
    std::cout << "read edges" << std::endl;
    if (feof(graphFile))
    {
        delete[] G->EI;
        delete[] G->VI;
        printf("unexpected end of file while reading edges\n");
        return -1;
    }
    else if (ferror(graphFile))
    {
        delete[] G->EI;
        delete[] G->VI;
        printf("error reading file\n");
        return -1;
    }

#ifdef WEIGHTED
    G->EW = new unsigned int[G->numEdges];
    fread (G->EW, sizeof(unsigned int), G->numEdges, graphFile);
    if (feof(graphFile))
    {
        delete[] G->EW;
        delete[] G->EI;
        delete[] G->VI;
        printf("unexpected end of file while reading edge weights\n");
        return -1;
    }
    else if (ferror(graphFile))
    {
        delete[] G->EW;
        delete[] G->EI;
        delete[] G->VI;
        printf("error reading file\n");
        return -1;
    }
    std::cout << "read weights" << std::endl;
#endif

    fclose(graphFile);

    return 1;
}

template<class graph>
void write_csr (char* filename, graph* G)
{
    FILE* fp = fopen(filename, "wb");
    if (fp == NULL)
    {
        fputs("file error", stderr);
        return;
    }
    fwrite(&G->numVertex, sizeof(intV), 1, fp); 
    fwrite(&G->numEdges, sizeof(intE), 1, fp); 
    fwrite(G->VI, sizeof(intE), G->numVertex, fp); 
    fwrite(G->EI, sizeof(intV), G->numEdges, fp); 
    fclose(fp); 
}

template<class graph>
void printGraph(graph* G)
{
    printf("num vertices = %d\n numEdges = %d\n", G->numVertex, G->numEdges);
    for (intV i=0; i<=G->numVertex; i++)
    {
        for (intE j=G->VI[i]; j<G->VI[i+1]; j++)
            printf("%d, %d\n", i, G->EI[j]);
    }
}

template<class graph>
void transposeCSR(graph* G1)
{
    intE* newVI = new intE[G1->numVertex+1]();
    intV* newEI = new intV[G1->numEdges]; 

    for (intE i=0; i<G1->numEdges; i++)
    {
        newVI[G1->EI[i]+1]++;
    }
    for (intV i=0; i<G1->numVertex; i++)
        newVI[i+1] += newVI[i];

    intV* tempId = new intV [G1->numVertex]();
    for (intV i=0; i<G1->numVertex; i++)
    {
        for (intE j=G1->VI[i]; j<G1->VI[i+1]; j++)
        {
            newEI[newVI[G1->EI[j]] + tempId[G1->EI[j]]] = i;
            tempId[G1->EI[j]]++;
        } 
    }
    delete[] G1->VI;
    delete[] G1->EI;
    delete[] tempId;
    G1->VI = newVI;
    G1->EI = newEI;
}


template<class graph>
void findOutDeg(graph* G)
{
    #pragma omp parallel for 
    for (intV i=0; i<G->numVertex; i++)
    {
        intE outDeg = G->VI[i+1] - G->VI[i];
        G->outDeg[i] = outDeg;
        for (intE j=G->VI[i]; j<G->VI[i+1]; j++)
        {
            #pragma omp atomic
            G->inDeg[G->EI[j]]++;
        }
    }
    return;
}

template<class graph>
void initGraph (graph* G)
{
    G->inFrontier = new bool [G->numVertex](); 
    G->outDeg = new intV [G->numVertex]();
    G->inDeg = new intV [G->numVertex]();
    findOutDeg(G);
    G->frontierSize = 0;
    G->frontier = new intV [G->numVertex]; 
    G->flag = new bool [G->numBins]();
    G->scatterDone = new bool [G->numBins]();
    return;
}

template<class graph>
void freeMem (graph* G)
{
    if (G->VI != NULL)
    {
        delete[] G->VI;
        G->VI = NULL;
    }
    if (G->EI != NULL)
    {
        delete[] G->EI;
        G->EI = NULL;
    }
#ifdef WEIGHTED
    if (G->EW != NULL)
    {
        delete[] G->EW;
        G->EW = NULL;
    }
#endif
    if (G->outDeg != NULL)
    {
        delete[] G->outDeg;
        G->outDeg = NULL;
    }
}

template<class graph>
intV findFrontierSize(graph* G)
{
    return G->frontierSize.load();
}


