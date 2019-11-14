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

//for asynchronous update propagation//
//converges faster//
#define ASYNCH

#include "../include/pcp.h"



struct CC_F{
    intV* label;
    CC_F(intV* _label):label(_label){}

    inline intV scatterFunc (intV node)
    {
        return label[node];
    }

    inline bool initFunc(intV node)
    {
        return false;
    }

    inline bool gatherFunc (intV updateVal, intV destId)
    {
        bool cond = (updateVal < label[destId]);
        if (cond)
            label[destId] = updateVal;
        return cond;
    }  
    
    inline bool filterFunc(intV node)
    {
        return true;
    } 

};


int main(int argc, char** argv)
{
    graph<intV> G;
    initialize(&G, argc, argv);
    initBin<intV>(&G);
    intV n = G.numVertex;
    intV* label = new intV [n]();
    intV initFrontierSize = n;
    intV* initFrontier = new intV [initFrontierSize];
    for (intV i=0; i<initFrontierSize; i++)
        initFrontier[i] = i;

    struct timespec start, end, half;
    float time;

    while((G.frontierSize > 0))
    {
         scatter_and_gather<intV>(&G, CC_F(label));
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
            scatter_and_gather<intV>(&G, CC_F(label));
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



