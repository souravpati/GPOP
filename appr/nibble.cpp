/**
 * Author: Kartik Lakhotia
           Sourav Pati
 * Email id: klakhoti@usc.edu
             spati@usc.edu
 * Date: 27-Feb-2018
 *
 * This code implements nibble algorithm
 * for computing probability distribution of 
 * a seeded random walk
 */


#define DUMP
#undef DUMP

unsigned int numIter = 0;


float threshold = 0.000000001;

#include "../include/pcp.h"



struct PR_F{
    float* pageRank;
    float* pageRankScat;
    intV* deg;
    PR_F(float* _pcurr, float* _pscat, intV* _outDeg):pageRank(_pcurr), pageRankScat(_pscat), deg(_outDeg){}
    inline float scatterFunc (intV node)
    {
        return pageRankScat[node];
    }
    inline bool initFunc(intV node)
    {
        pageRank[node]=pageRank[node]/2;
        pageRankScat[node] = 0;
        return (pageRank[node] >= threshold*deg[node]);
    }
    inline bool gatherFunc (float updateVal, intV destId)
    {
        pageRank[destId] += updateVal;
        return (updateVal > 0);
    }
    inline bool filterFunc(intV node)
    {
        bool cond = (pageRank[node] >= threshold*deg[node]);
        if (!cond)
            pageRank[node] = 0;
        if (cond && (deg[node]>0))
            pageRankScat[node] = pageRank[node]/(2*deg[node]);
        return cond; 
    } 
};




int main(int argc, char** argv)
{
    graph<float> G;
    initialize(&G, argc, argv);    
    initBin<float>(&G);    
    float* pcurr = new float [G.numVertex]();
    float* pscat = new float [G.numVertex]();
    for (intV i=0; i<G.numVertex; i++)
    {
        pcurr[i] = 0;
        pscat[i] = 0;
    }
    intV initFrontierSize = 1;
    intV* initFrontier = new intV [initFrontierSize];
    intV actVertex;
 
    for (intV i=0; i<initFrontierSize; i++)
        initFrontier[i] = G.start;

    struct timespec start, end;
    float time;

    int ctr =0;
    while(ctr < G.rounds)
    {
        resetFrontier(&G);

        for (intV i=0; i<initFrontierSize; i++)
        {
            actVertex = initFrontier[i];
            pcurr[actVertex] = 1;
            if (G.outDeg[actVertex] > 0)
                pscat[actVertex] = pcurr[actVertex]/(2*G.outDeg[actVertex]);
        }
        loadFrontier (&G, initFrontier, initFrontierSize);
        numIter = 0;

        if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}

        while(numIter < MAX_ITER)
        {
            scatter_and_gather<float>(&G, PR_F(pcurr, pscat, G.outDeg));   
            numIter++;
        }   

        getFrontier(&G);
        #pragma omp parallel for
        for (intV i=0; i<G.frontierSize; i++)
        {
            pcurr[G.frontier[i]]=0;
            pscat[G.frontier[i]]=0;
        }
 
        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");}
        time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
        printf("appr, %d, %s, %lf\n",NUM_THREADS, argv[1], time);
        ctr++;
    }
    printf("\n");

    return 0;
}
