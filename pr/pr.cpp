/**
 * Author: Kartik Lakhotia
           Sourav Pati
 * Email id: klakhoti@usc.edu
             spati@usc.edu
 * Date: 27-Feb-2018
 *
 * Compute PageRank using Partition-centric graph processing
 *
 */

#define DENSE
//#undef DENSE

unsigned int numIter = 0;

#include "../include/pcp.h"


float damping =0.15;

struct PR_F{
    float* pageRank;
    intV* deg;
    PR_F(float* _pcurr, intV* _outDeg):pageRank(_pcurr), deg(_outDeg){}
    inline float scatterFunc (intV node)
    {
        return pageRank[node];
    }
#ifndef DENSE
    inline bool initFunc(intV node)
    {
        pageRank[node]=0;
        return true;
    }
    inline bool gatherFunc (float updateVal, intV destId)
    {
        pageRank[destId] += updateVal;
        return true;
    }
    inline bool filterFunc(intV node)
    {
        pageRank[node] = ((damping) + (1-damping)*pageRank[node]);
        if (deg[node]>0)
            pageRank[node] = pageRank[node]/deg[node];
        return true;
    } 
#else
    inline void initFunc(intV node)
    {
        pageRank[node]=0;
    }
    inline void gatherFunc (float updateVal, intV destId)
    {
        pageRank[destId] += updateVal;
    }
    inline void filterFunc(intV node)
    {
        pageRank[node] = ((damping) + (1-damping)*pageRank[node]);
        if (deg[node]>0)
            pageRank[node] = pageRank[node]/deg[node];
    }
#endif
};




int main(int argc, char** argv)
{
    graph<float> G;
    initialize(&G, argc, argv);    
    initBin<float>(&G);    
    float* pcurr = new float [G.numVertex]();
    for(int i=0;i<G.numVertex;i++){
        if (G.outDeg[i] > 0)
            pcurr[i] = 1.0/G.outDeg[i];
        else
            pcurr[i] = 1.0;
    }
#ifndef DENSE
    intV* initFrontier = new intV [G.numVertex];
#pragma omp parallel for
    for (intV i=0; i<G.numVertex; i++)
        initFrontier[i] = i;
    loadFrontier (&G, initFrontier, G.numVertex);
#endif


    struct timespec start, end;
    float time;

    int ctr =0;
    while(ctr < G.rounds)
    {
        numIter = 0;
        for(int i=0;i<G.numVertex;i++){
            if (G.outDeg[i] > 0)
                pcurr[i] = 1.0/G.outDeg[i];
            else
                pcurr[i] = 1.0;
        }

        if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}

        while(numIter < MAX_ITER)
        {
            scatter_and_gather<float>(&G, PR_F(pcurr, G.outDeg));
            numIter++;
        }

        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");}
        time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
        printf("pr_dense, %d, %s, %lf\n",NUM_THREADS, argv[1], time);
        ctr++;
    }
    printf("\n");

    return 0;
}



