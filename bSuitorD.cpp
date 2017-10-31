
// ***********************************************************************
//
//            Vite: A C++ library for shared-memory graph b-Matching 
//                  using OpenMP
// 
//               Arif Khan (ariful.khan@pnnl.gov)
//               Mahantesh Halappanavar (hala@pnnl.gov)
//               Pacific Northwest National Laboratory       
//
//               Alex Pothen (apothen@purdue.edu)
//               Purdue University
//
// ***********************************************************************
//
//       Copyright (2017) Battelle Memorial Institute
//                      All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************
 
#include "bMatching.h"
#include <cstring>
//#include <parallel/algorithm>
using namespace std;

void bSuitorBPQD(CSR* G, int *b, int *nlocks, Node* S, int* start, int* end, char* mark, 
        int type, int stepM, bool verbose) 
{

    int numThreads;
    #pragma omp parallel
    numThreads=omp_get_num_threads();
    
    Edge* verInd=G->verInd;
    int n=G->nVer;
    int m=G->nEdge;
    int* ver=G->verPtr;
    int* Q1=new int[n];
    int* Q2=new int[n];
    int* Qb1=new int[n];
    int* Qb2=new int[n];
    
    #ifdef LBUF
    int** TQ=new int*[numThreads];
    int* tindx=new int[numThreads];

    for(int i=0;i<numThreads;i++)
    {
        tindx[i]=0;
        TQ[i]=new int[BSIZE];
    }
    #endif
    
    int* Qtemp;
    int Qsize=n,Qindx=0,iter=1;
    double t1,t2,t3,t4,t5;
    long long edgeC=0,heapC=0,kickC=0;
    

    if(verbose)
        t1=omp_get_wtime();

    #ifdef DYN
        #pragma omp parallel for schedule(dynamic,CHUNK)
    #else
        #pragma omp parallel for schedule(static, CHUNK)
    #endif
    for(int i=0;i<n;i++)  
    {
        
        int tid=omp_get_thread_num();
        nlocks[i]=0;            // Initialize locks
        S[i].curSize=0;         // current matching=0
        S[i].maxSize=b[i];         // maximum matching=b
        S[i].minEntry.id=-1;
        S[i].minEntry.weight=0.0;
        

        if(type!=1)
        {
            start[i]=ver[i];    // Adj List start pointer
            if(b[i]>0)
                end[i]=custom_sort(verInd,ver[i],ver[i+1],stepM*b[i],NULL,type);
            else
            {
                end[i]=-1;
                S[i].minEntry.id=G->nEdge+1;
                S[i].minEntry.weight=MAX_VAL;
            }
               
        }
        
        Q1[i]=i;
        Qb1[i]=b[i];
        Qb2[i]=0;
        
    }

    if(type ==1) // Unsorted mode
    {
        #ifdef DYN
            #pragma omp parallel for schedule(guided,CHUNK)
        #else
            #pragma omp parallel for schedule(static, CHUNK)
        #endif
        for(int i=0;i<m;i++)
            mark[i]=0;
    }
   
    cout<<"Initialization Done...!!"<<endl;
    
    
    if(verbose)
        t2=omp_get_wtime();
    
    cout << "Start Matching" << Qsize << endl;
    
    while(Qsize>0)
    {
   
        t4=omp_get_wtime();
        //#ifdef DYN
        //    #pragma omp parallel for schedule(guided, CHUNK)
        //#else
        //    #pragma omp parallel for schedule(static, CHUNK)
        //#endif

	    #pragma omp parallel num_threads(numThreads) 
	    {
		    int threadid = omp_get_thread_num();

            #pragma omp for	schedule(guided,CHUNK)		
		    for(int i=0;i<Qsize;i++) 
		    {
		        float min_w,heaviest,weight;
		        int min_idx,bVer,current,next_vertex,sold,eold,skip;
		        int partnerIdx,partner,y,j,done,tid=omp_get_thread_num();
		    
		        #ifdef LBUF
			    int* LQ=TQ[tid];
			    int* lindx=&tindx[tid];
		        #endif
		    
		        current=Q1[i];
		        bVer=Qb1[current];
		        Qb1[current]=0;

		        while(bVer>0){ // For Each vertex we want to find bVer=b partners

			    partner=-1;
			    heaviest=0.0;
			    next_vertex=-1;
			    done=1;

			    sold=start[current];
			    eold=end[current];
			    j=sold;
			    if(type!=1) /// Sorted mode
			    {
			        while(j<end[current]) // Loop over neighbors of the current vertex 
			        {    

				        y = verInd[j].id;         // y is the neighbor of the current vertex
				        weight=verInd[j].weight;  // weight is w(current,y)
				        
                        if(weight<=0)  // If there is no more positive weighted edge
                        {
                            heaviest=-1.0;
                            break;
                        }
                        
                        min_w=S[y].minEntry.weight;
				
				        if(weight <=0 || min_w > weight)
				        {
				            j++;
				            continue; // Next edge
				        }
				
				        while(__sync_lock_test_and_set(&nlocks[y],1)==1);
							   
				        min_w=S[y].minEntry.weight;
				        min_idx=S[y].minEntry.id;

							    
				        if((min_w > weight) || (weight == min_w && current < min_idx))
				        {    
				            __sync_lock_release(&nlocks[y]);
				            j++;
				        }
				        else
				        {
				            heaviest=weight;
				            S[y].AddHeap(heaviest,current);    
								
				            __sync_lock_release(&nlocks[y]); //nlocks[y]--;
				            start[current]=j+1;

				            if(min_idx!=-1 && end[min_idx]!=-1)
				            {    
				                if(verbose)
					                __sync_fetch_and_add(&kickC,1);;
					
                            #ifdef LBUF
					        if(__sync_fetch_and_add(&Qb2[min_idx],1)==0)
					        {   
					            LQ[(*lindx)]=min_idx;
					            (*lindx)++;
					            if((*lindx)==BSIZE)
					            {
						            int start=__sync_fetch_and_add(&Qindx,BSIZE);
						            int count=0;
						            for(int k=start;k<start+BSIZE;k++)
						            {    
						                Q2[k]=LQ[count];
						                count++;
						            }
						            (*lindx)=0;
					            }
					        }
					        #else
					        if(__sync_fetch_and_add(&Qb2[min_idx],1)==0)
					            Q2[__sync_fetch_and_add(&Qindx,1)]=min_idx;
					        #endif
				            } 
				    
				            break;
				        }
			        }   // while(j<end[current])
			    
                    if(heaviest<=0)
			        {
				    
				        start[current]=j;
				        if(end[current]<ver[current+1] && heaviest==0)
				        {
					        end[current]=custom_sort(verInd,start[current],
                                ver[current+1],stepM*b[current],NULL,type);
			    
				            done=0;
					
                        #ifdef LBUF
                        if(__sync_fetch_and_add(&Qb2[current],bVer)==0)
                        {   
                            LQ[(*lindx)]=current;
                            (*lindx)++;
                            if((*lindx)==BSIZE)
                            {
                            int start=__sync_fetch_and_add(&Qindx,BSIZE);
                            int count=0;
                            for(int k=start;k<start+BSIZE;k++)
                            {    
                                Q2[k]=LQ[count];
                                count++;
                            }
                            (*lindx)=0;
                            }
                        }
                        
                        #else
                        if(__sync_fetch_and_add(&Qb2[current],bVer)==0)
                           Q2[__sync_fetch_and_add(&Qindx,1)]=current;
                        #endif
                        
                        }
				        else
				            end[current]=-1;
			        }
			    }
			    else /// Unsorted Mode
			    {
			        j=ver[current];
			    
			        for(;j<ver[current+1] && j!=-1;j++) 
                    { // Loop over neighbors of the current vertex 
				
				        if(mark[j]==1)
				            continue;
				
                        y = verInd[j].id;       // y is the neighbor of the current vertex
                        weight=verInd[j].weight;        // weight is w(current,y)
                        min_w=S[y].minEntry.weight;

                        if((weight<heaviest)|| (weight == heaviest && y < partner)||min_w>weight)
                            continue;

                        if(min_w==weight) 
                        {        
                            skip=0;
                    
                            while(__sync_lock_test_and_set(&nlocks[y],1)==1);
                    
                            min_w=S[y].minEntry.weight;
                            min_idx=S[y].minEntry.id;
                            if((min_w > weight)||(weight == min_w && current < min_idx))
                            skip=1;
                            __sync_lock_release(&nlocks[y]);
                       
                            if(skip==1)
                                continue;
                        }
                       
                        heaviest = weight;  // Store the weight of the heaviest edge found so far
                        partner = y;
                        partnerIdx=j;
                    
                    
                    } // loop over neighbors
                    
                    if (heaviest > 0) // True if there is a new partner
                    {    
                        if(__sync_lock_test_and_set(&nlocks[partner],1)==0) //Locking partner
                        {
                            min_w=S[partner].minEntry.weight;
                            min_idx=S[partner].minEntry.id;

                            if(!S[y].find_id(current) && 
                                    ((heaviest > min_w)||((heaviest==min_w)&&(current>min_idx)))) 
                            {

                                // Need to re check again 
                                // because during time of locking someone else can change the state
                                
                                next_vertex=min_idx;
                                S[partner].AddHeap(heaviest,current);
                            }    
                            else 
                            // Got the lock but someone already changed 
                            // the state so search partner all over again
                                next_vertex=current;
                         
                            mark[partnerIdx]=1;
                            __sync_lock_release(&nlocks[partner]); // Unlocking partner
           
                        }
                        else 
                            // Missed the lock, so someone else may or 
                            // may not changed the state. Hence, RE DO the partner search
                            next_vertex=current;
                  
                        if (next_vertex != -1 && next_vertex!=current)  
                            // True if current vertex just kicked another vertex which is alive
                        {  
                            #ifdef LBUF

                            if(__sync_fetch_and_add(&Qb2[next_vertex],1)==0)
                            {   
                                LQ[(*lindx)]=next_vertex;
                                (*lindx)++;
                                if((*lindx)==BSIZE)
                                {
                                    int start=__sync_fetch_and_add(&Qindx,BSIZE);
                                    int count=0;
                                    for(int k=start;k<start+BSIZE;k++)
                                    {    
                                        Q2[k]=LQ[count];
                                        count++;
                                    }
                                    (*lindx)=0;
                                }
                            }
                            #else
                            if(__sync_fetch_and_add(&Qb2[next_vertex],1)==0)
                                Q2[__sync_fetch_and_add(&Qindx,1)]=next_vertex;
                            #endif
                        }
                    }
                    else
                    end[current]=-1;    
                    //This means the neighbors list is exhausted, 
                    //this vertex will never be considered again (dead).!!
			    }
			
			    if(end[current]==-1)  // if the vertex is dead 
			        break;     // Do not forget to decrease bVer ...!!
			    else
			        if(next_vertex!=current && done==1)
				        bVer--;
                if(type!=1 && done ==0)
                    break;
		    } // while(bVer)
		} // loop over vertices
	}       
 
    #ifdef LBUF
    #pragma omp parallel for
    for(int i=0;i<numThreads;i++)
    { 
        int* LQ=TQ[i];
        int* lindx=&tindx[i];
        if((*lindx)>0)
        {    
            int start=__sync_fetch_and_add(&Qindx,(*lindx));
            int count=0;
            for(int k=start;k<start+(*lindx);k++)
            {        
                Q2[k]=TQ[i][count];
                count++;
            }
            (*lindx)=0;
        }
    }
    #endif
    Qtemp=Q1;
    Q1=Q2;
    Q2=Qtemp;
    Qtemp=Qb1;
    Qb1=Qb2;
    Qb2=Qtemp;
    Qsize=Qindx;
    Qindx=0;
    iter++;
    
    if(verbose)
        cout<<"Iteration: "<<iter-1<<" "<<Qsize<<" "<<omp_get_wtime()-t4<<endl;
    }

    cout<<"Matching Done....!!"<<endl;

    if(verbose)
    {
        t3=omp_get_wtime();

        cout<<"# of Iteration: "<<iter-1<<endl;
        cout<<"Initialization Time: "<<t2-t1<<endl;
        cout<<"Matching Time: "<<t3-t2<<endl;
        cout<<"Total Time: "<<t3-t1<<endl;
    } 

    delete Q1,Q2,Qb1,Qb2,tindx;
    for(int i=0;i<numThreads;i++)
        delete TQ[i];
    delete []TQ;

}// end of bSuitor

