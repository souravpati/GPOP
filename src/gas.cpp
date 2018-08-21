//#include "../include/gas.h"

//void loadFrontier(partitionData* allTD, unsigned int* initFrontier, unsigned int initFrontierSize)
//{
//    graph*G = allTD[0].G;
//    unsigned int vertexId, pId;
//    for (unsigned int i=0; i<initFrontierSize; i++)
//    {
//        vertexId = initFrontier[i];
//        G->inFrontier[vertexId] = true;
//        G->visited[vertexId] = true;
//        pId = (vertexId >> binOffsetBits);
//        allTD[pId].frontier[allTD[pId].frontierSize++] = vertexId;
//    }
//    G->frontierSize = initFrontierSize;
//}
//
//void resetFrontier(partitionData* TD)
//{
//    graph* G = TD->G;
//    for (unsigned int i=0; i<TD->frontierSize; i++)
//        G->inFrontier[TD->frontier[i]] = false;
//    TD->frontierSize = 0;
//}

//inline void buildFrontier (graph* G, partitionData* TD, unsigned int destID)
//{
//    if (!G->inFrontier[destId])
//    {
//        TD->frontier[TD->frontierSize++] = destId;
//        G->inFrontier[destId] = true;
//    }
//}
//
// filter the frontier //
//void filterFrontier(partitionData* TD)
//{
//    graph* G = TD->G;
//    unsigned int trueSize = 0;
//    for (unsigned int i=0; i<TD->frontierSize; i++)
//    {
//        // BFS condition. Can be some other condition as well //
//        if (G->visited[TD->frontier[i]] == false)
//        {
//            TD->activeEdges += G->outDeg[TD->frontier[i]];
//            G->visited[TD->frontier[i]] = true;
//            TD->frontier[trueSize++] = TD->frontier[i];
//        }
//    }
//    TD->frontierSize = trueSize;
//    if (TD->frontierSize > 0)
//        G->frontierSize += TD->frontierSize;
//    
////    std::atomic_fetch_add(&G->frontierSize, trueSize);    
//}

//void densityCheck(partitionData* TD)
//{
//    unsigned int e = TD->PNG->numEdges;
//    unsigned int m = TD->activeEdges;
//    TD->isDense = (17.0 * (float)m > (float(e)*10.5 + 0.67*(float)m));
//    
//}
//
//void scatterVC(partitionData* TD, bool** updateBins, unsigned int** destIdBins, unsigned int* updateBinPointers, unsigned int* destIdBinPointers)
//{
//    graph* G = TD->G;
//    unsigned int destId = 0;
//    unsigned int destBin = 0;
//    unsigned int vertexId = 0;
//    unsigned int cond = 0;
//    unsigned int prevBin = 0;
//#ifdef DEBUGL2
//    for (unsigned int i=0; i<NUM_BINS; i++)
//    {
//        assert(updateBinPointers[i]==0);
//        assert(destIdBinPointers[i]==0);
//    }
//#endif
//    for (unsigned int i=0; i<TD->frontierSize; i++)
//    {
//        vertexId = TD->frontier[i];
//#ifdef DEBUGL2
//        assert(vertexId >= TD->startVertex);
//        assert(vertexId < TD->endVertex);
//#endif
//        G->inFrontier[vertexId] = false;
//        prevBin = NUM_BINS;
//        for (unsigned int j=G->VI[vertexId]; j<G->VI[vertexId+1]; j++)
//        {
//            destId = G->EI[j];
//            destBin = (destId >> binOffsetBits);
//            ///////////////////////////////////
//            ///// branch avoiding approach ////
//            ///////////////////////////////////
//            //cond = (destBin != prevBin);
//            //updateBins[destBin][updateBinPointers[destBin]] = true;
//            //updateBinPointers[destBin] += cond;
//            //destId |= (cond << 31);
//            //destIdBins[destBin][destIdBinPointers[destBin]++] = destId;
//            //prevBin = destBin;
//
//            ///////////////////////////////////
//            //////// branched approach ////////
//            ///////////////////////////////////
//            //if (destBin!=prevBin)
//            //{
//            //    updateBins[destBin][updateBinPointers[destBin]++] = prVal;
//            //    destId |= (MAX_NEG);
//            //    prevBin = destBin;
//            //}
//            //destIdBins[destBin][destIdBinPointers[destBin]++] = destId;
//
//            ///////////////////////////////////
//            ////// VC update propagation //////
//            ///////////////////////////////////
//            updateBins[destBin][destIdBinPointers[destBin]] = G->visited[vertexId];
//            destIdBins[destBin][destIdBinPointers[destBin]++] = destId;
//            if (!G->flag[destBin])
//	            G->flag[destBin] = true;    
//        }
//    }
//    TD->frontierSize = 0;
//}
//
//void scatterPC(partitionData* TD, bool** updateBins)
//{
//    graph* PNG = TD->PNG;
//    graph* G = TD->G;
//    unsigned int pointer;
//    for (unsigned int i=0; i<PNG->numVertex; i++)
//    {
//        pointer = 0;
//        if (!G->flag[i])
//            G->flag[i] = ((PNG->VI[i+1]-PNG->VI[i])>0);
//        for (unsigned int j=PNG->VI[i]; j<PNG->VI[i+1]; j++)
//            updateBins[i][pointer++] = G->visited[PNG->EI[j]];
//    } 
//    resetFrontier(TD); 
//}
//
//
//void scatter(partitionData* TD, bool** updateBins, unsigned int** destIdBins, unsigned int* updateBinPointers, unsigned int* destIdBinPointers)
//{
//    if (TD->isDense)
//        scatterPC(TD, updateBins);
//    else
//        scatterVC(TD, updateBins, destIdBins, updateBinPointers, destIdBinPointers);
//}
//
//////////////////////////////////
///////// DENSE GATHER ///////////
//////////////////////////////////
//template <class userArg>
//void gatherPC(partitionData* TD, bool* updateBin, unsigned int* destIdBin, unsigned int binSize, userArg UA)
//{
//    graph*G = TD->G;
//    unsigned int destId = 0;
//    unsigned int updateBinPointer = MAX_UINT;
//    bool updateVal;
//    bool userReturn;
//    for (unsigned int j=0; j<binSize; j++)
//    {
//        destId = destIdBin[j];
//        updateBinPointer += (destId >> 31);
//        destId = destId & MAX_POS;
//        updateVal = updateBin[updateBinPointer];
//        
////        if (updateVal)
////        {
////            userReturn = true;
////        }
//        userReturn = UA.gatherFunc(G, updateVal, destId);
//        if (!G->inFrontier[destId] && userReturn)
//        {
//            TD->frontier[TD->frontierSize++] = destId;
//            G->inFrontier[destId] = true;
//        }
//    } 
//    //////////////////////////////////
//    //////////////////////////////////
//	
//}
//
////////////////////////////////////
/////////// SPARSE GATHER ///////////
////////////////////////////////////
//template <class userArg>
//void gatherVC(partitionData* TD, bool* updateBin, unsigned int* destIdBin, unsigned int binSize, userArg UA)
//{
//    graph* G = TD->G;
//    unsigned int destId = 0; 
//    bool updateVal;
//    bool userReturn; 
//    for (unsigned int j=0; j<binSize; j++)
//    {
//        destId = destIdBin[j];
//        updateVal = updateBin[j];  
////        if (updateVal)
////        {
////            userReturn = true;
////        }
//        userReturn = UA.gatherFunc(G, updateVal, destId); 
//        if (!G->inFrontier[destId] && userReturn)
//        {
//            TD->frontier[TD->frontierSize++] = destId;
//            G->inFrontier[destId] = true;
//        }
//    } 
//}

//void gather(partitionData* TD, bool*** updateBins, unsigned int*** denseDestIdBins, unsigned int*** sparseDestIdBins, partitionData* allTD, unsigned int** destIdBinAddrSize, unsigned int** destIdBinPointers, unsigned int** updateBinPointers, userArg UA)
//{
//    graph*G = TD->G;
//    TD->activeEdges = 0;
//    if(G->flag[TD->tid] == false)
//        return;
//
//    G->flag[TD->tid] = false;
//
//    for (unsigned int i=0; i<NUM_BINS; i++)
//    {
//        if (allTD[i].isDense)
//            gatherPC(TD, updateBins[i][TD->tid], denseDestIdBins[i][TD->tid], destIdBinAddrSize[i][TD->tid], UA);
//        else
//            gatherVC(TD, updateBins[i][TD->tid], sparseDestIdBins[i][TD->tid], destIdBinPointers[i][TD->tid], UA); 
//        destIdBinPointers[i][TD->tid] = 0;
//        updateBinPointers[i][TD->tid] = 0;
//    }
//    filterFrontier(TD);
//} 
//
//
