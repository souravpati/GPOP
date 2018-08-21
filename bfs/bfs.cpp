/**
 * Author: Kartik Lakhotia
           Sourav Pati
 * Email id: klakhoti@usc.edu
             spati@usc.edu
 * Date: 27-Feb-2018
 *
 * This code implements work efficient BFS
 * 
 */

unsigned int numIter = 0;


#include "../include/pcp.h"

#define MAX_NEG 0x80000000


struct BFS_F{
    unsigned int* parent;
    bool* visited;
    BFS_F(unsigned int* _parent, bool* _visited):parent(_parent), visited(_visited){}
    inline unsigned int scatterFunc (unsigned int node)
    {
        return (((!visited[node])<<31) | node);
    }

    inline bool reInit(unsigned int node)
    {
        return false;
    }

    inline bool gatherFunc (unsigned int updateVal, unsigned int destId)
    {
//        if((!visited[destId]) && (!(updateVal>>31)))
//        {
//            parent[destId] = updateVal;
//            visited[destId] = true;
//            return true;
//        }
        if (!visited[destId])
        {
            parent[destId] = updateVal;
            visited[destId] = (!(updateVal>>31));
            return visited[destId];
        }
        return false;
    }  
    
    inline bool apply(unsigned int node)
    {
        return true;
    } 

};

struct BFS_F2{
    unsigned int* parent;
    BFS_F2(unsigned int* _parent):parent(_parent){}
    inline unsigned int scatterFunc (unsigned int node)
    {
        return ((parent[node] & MAX_NEG)| node);
    }

    inline bool reInit(unsigned int node)
    {
        return false;
    }

    inline bool gatherFunc (unsigned int updateVal, unsigned int destId)
    {
        if((parent[destId] & MAX_NEG) && (!(updateVal>>31)))
        {
            parent[destId] = updateVal;
            return true;
        }
        return false;
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
    unsigned int* parent = new unsigned int [n]();
    bool* visited = new bool [n]();
    unsigned int initFrontierSize = 1;
    unsigned int* initFrontier = new unsigned int [initFrontierSize];
    for (unsigned int i=0; i<initFrontierSize; i++)
        initFrontier[i] = G.start;  

	struct timespec start, end;
    float time;

    int ctr = 0;
    while(ctr < G.rounds){
   
        for(int i=0;i<n;i++)
        {
            parent[i] = 1<<31;
            visited[i] = false;
        }
        visited[G.start] = true;
        parent[G.start] = G.start;

        if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}
    
//        loadFrontier(&G, initFrontier, initFrontierSize);
        loadFrontierPar(&G, initFrontier, initFrontierSize);

        numIter=0;
        while((G.frontierSize > 0))
        {
            pcpm<unsigned int>(&G, BFS_F(parent, visited));
//            pcpm<unsigned int>(&G, BFS_F2(parent));
            numIter++;
         }   

        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");}
        time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
        printf("bfs, %d, %s, %lf\n", NUM_THREADS, argv[1], time);
        ctr++; 
    }


    printf("\n");
    return 0;
}



