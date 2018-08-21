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
#include "../include/sort.h"

//////////////////////////////////////////
//partition centric programming data types
//////////////////////////////////////////

typedef struct partitionGraph 
{
    unsigned int numVertex;
    unsigned int numEdges;
    unsigned int* VI;
    unsigned int* EI;
    unsigned int* outDeg;
} partitionGraph;

typedef struct partitionData
{
    unsigned int tid;
    unsigned int startVertex;
    unsigned int endVertex;
    partitionGraph* PNG;
    unsigned int* frontier;
    unsigned int frontierSize;
    bool isDense;
    unsigned int totalEdges;
    unsigned int activeEdges;
    std::atomic<unsigned int> binListPtr;
}partitionData;

template<class type>
struct graph 
{
    unsigned int numBins;
    unsigned int numVertex;
    unsigned int numEdges;
    unsigned int* VI;
    unsigned int* EI;
    unsigned int* EW;
    unsigned int* frontier;
    unsigned int start;
    unsigned int rounds;
//    boost::dynamic_bitset<> inFrontier; 
    bool* inFrontier;
    std::atomic<unsigned int> frontierSize;
    partitionData* TD;
    unsigned int* outDeg;
    unsigned int* inDeg;
    bool* flag;
    bool** binFlag;
    unsigned int** updateBinAddrSize;
    unsigned int** destIdBinAddrSize;
    unsigned int** updateBinPointers;
    unsigned int** destIdBinPointers;  
    unsigned int*** indWeightBins;
    unsigned int*** indDestIdBins;
    unsigned int*** sparseDestIdBins;
    type*** indUpdateBins;   
    unsigned int* activeScatter;
    unsigned int* activeGather;
    std::atomic<unsigned int> partListPtr;
    unsigned int** activeBins;
};


//int read_csr (char* filename, graph* G);
//void write_csr (char* filename, graph* G);
//void printGraph(graph* G);
//void transposeCSR(graph* G1);
//void sortEdges(graph* G);
//void initGraph (graph* G);
//void freeMem (graph* G);





template<class graph>
int read_csr (char* filename, graph* G)
{
    FILE* graphFile = fopen(filename, "rb");
    if (graphFile == NULL)
    {
        fputs("file error", stderr);
        return -1;
    }
    fread (&(G->numVertex), sizeof(unsigned int), 1, graphFile);
    
    fread (&(G->numEdges), sizeof(unsigned int), 1, graphFile);


    G->VI = new unsigned int[G->numVertex+1];
    fread (G->VI, sizeof(unsigned int), G->numVertex, graphFile);
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

    G->EI = new unsigned int[G->numEdges];
    fread (G->EI, sizeof(unsigned int), G->numEdges, graphFile);
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
    fwrite(&G->numVertex, sizeof(unsigned int), 1, fp); 
    fwrite(&G->numEdges, sizeof(unsigned int), 1, fp); 
    fwrite(G->VI, sizeof(unsigned int), G->numVertex, fp); 
    fwrite(G->EI, sizeof(unsigned int), G->numEdges, fp); 
    fclose(fp); 
}

template<class graph>
void printGraph(graph* G)
{
    printf("num vertices = %d\n numEdges = %d\n", G->numVertex, G->numEdges);
    for (unsigned int i=0; i<=G->numVertex; i++)
    {
        for (unsigned int j=G->VI[i]; j<G->VI[i+1]; j++)
            printf("%d, %d\n", i, G->EI[j]);
    }
}

template<class graph>
void transposeCSR(graph* G1)
{
    unsigned int* newVI = new unsigned int[G1->numVertex+1]();
    unsigned int* newEI = new unsigned int[G1->numEdges]; 

    for (unsigned int i=0; i<G1->numEdges; i++)
    {
        newVI[G1->EI[i]+1]++;
    }
    for (unsigned int i=0; i<G1->numVertex; i++)
        newVI[i+1] += newVI[i];

    unsigned int* tempId = new unsigned int [G1->numVertex]();
    for (unsigned int i=0; i<G1->numVertex; i++)
    {
        for (unsigned int j=G1->VI[i]; j<G1->VI[i+1]; j++)
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
void sortEdges(graph* G)
{
    #pragma omp parallel for
    for (unsigned int i=0; i<G->numVertex; i++)
    {
        if (G->VI[i+1] > (G->VI[i]+1))
            mergeSortWOkey<unsigned int>(G->EI, G->VI[i], G->VI[i+1]-1);
    }
    return;
}

template<class graph>
void findOutDeg(graph* G)
{
    #pragma omp parallel for
    for (unsigned int i=0; i<G->numVertex; i++)
    {
        G->outDeg[i] = G->VI[i+1] - G->VI[i];
        for (unsigned int j=G->VI[i]; j<G->VI[i+1]; j++)
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
    G->outDeg = new unsigned int [G->numVertex]();
    G->inDeg = new unsigned int [G->numVertex]();
    findOutDeg(G);
    G->frontierSize = 0;
    G->frontier = new unsigned int [G->numVertex]; 
    G->flag = new bool [G->numBins]();
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
unsigned int findFrontierSize(graph* G)
{
    return G->frontierSize.load();
}


