/**
 * Author: Kartik Lakhotia
           Sourav Pati
 * Email id: klakhoti@usc.edu
             spati@usc.edu
 * Date: 27-Feb-2018
 *
 * Weakly connected components
 * 
 */

unsigned int numIter = 0;

#include "../include/pcp.h"



struct CC_F{
    unsigned int* label;
    CC_F(unsigned int* _label):label(_label){}

    inline unsigned int scatterFunc (unsigned int node)
    {
        return label[node];
    }

    inline bool reInit(unsigned int node)
    {
        return false;
    }

    inline bool gatherFunc (unsigned int updateVal, unsigned int destId)
    {
        bool cond = (updateVal < label[destId]);
        if (cond)
            label[destId] = updateVal;
        return cond;
    }  
    
    inline bool apply(unsigned int node)
    {
        return true;
    } 

};


int main(int argc, char** argv)
{
    graph<unsigned int> G;
    initialize(&G, argc, argv);
    initBin<unsigned int>(&G);
    unsigned int n = G.numVertex;
    unsigned int* label = new unsigned int [n]();
    unsigned int initFrontierSize = n;
    unsigned int* initFrontier = new unsigned int [initFrontierSize];
    for (unsigned int i=0; i<initFrontierSize; i++)
        initFrontier[i] = i;

    struct timespec start, end, half;
    float time;

    while((G.frontierSize > 0))
    {
         pcpm<unsigned int>(&G, CC_F(label));
         numIter++;
    }
    numIter = 0;
    int ctr =0;
    while(ctr < G.rounds){
        for(int i=0;i<n;i++)
            label[i] = i;
        loadFrontier(&G, initFrontier, initFrontierSize);

        if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}

        while((G.frontierSize > 0))
        {
            pcpm<unsigned int>(&G, CC_F(label));
            numIter++;
        }

        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");}
        time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
        printf("cc, %d, %s, %lf\n", NUM_THREADS, argv[1], time);
        ctr++;
    }
    printf("\n");

    return 0;
}



