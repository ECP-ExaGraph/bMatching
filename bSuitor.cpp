
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
using namespace std;


int pivotID;
float pivotW;
int ptype;

void Node::print()
{
    int i;
    printf("Printing the heap array\n");
    for(i=0;i<curSize;i++)
    {
        printf("%d %f\n",heap[i].id,heap[i].weight);
    }
}


bool verifyMatching(CSR* g, Node* S, int n)
{
    bool flag=false;
    float weight=0.0;
    int count=0;
    //#pragma omp parallel for schedule(static)
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<S[i].curSize;j++)
        {
            int a=S[i].heap[j].id;
            weight+=S[i].heap[j].weight;
            count++;
            if(a!=-1 && !S[a].find_id(i))
            {    
                cout<<"("<<i<<","<<a<<") :";
                for(int k=0;k<S[i].curSize;k++)
                    cout<<S[i].heap[k].id<<"("<<S[i].heap[k].weight<<")"<<" ";
                cout<<" I ";
                for(int k=0;k<S[a].curSize;k++)
                    cout<<S[a].heap[k].id<<"("<<S[a].heap[k].weight<<")"<<" ";
                cout<<endl;
                
                cout<<"i: ";
                for(int k=0;k<n;k++)
                    if(k!=i && S[k].find_id(i))
                    cout<<k<<" ";
                cout<<endl;
                cout<<"a: ";
                for(int k=0;k<n;k++)
                    if(k!=a && S[k].find_id(a))
                    cout<<k<<" ";
                cout<<endl;
                cout<<"i: ";
                for(int k=g->verPtr[i];k<g->verPtr[i+1];k++)
                    cout<<g->verInd[k].id<<" "<<g->verInd[k].weight<<",";
                cout<<endl;
                cout<<"a: ";
                for(int k=g->verPtr[a];k<g->verPtr[a+1];k++)
                    cout<<g->verInd[k].id<<" "<<g->verInd[k].weight<<",";
                cout<<endl;
                flag=true;
                break;
            }
        }
        if(flag)
            break;
    }

    
    
    if(flag)
        return false;
    else
    {
        if(count==0)
            cout<<"It is an empty matching"<<endl;
        else
            cout<<"Matching Weight: "<<weight/2.0<<endl;
        return true;
    }
}

void printMatching(Node* S, int n)
{

    cout<<"---- Matching Output----"<<endl;
    for(int i=0;i<n;i++)
    {
        cout<<i<<" :";
        for(int j=0;j<S[i].curSize;j++)
            cout<<S[i].heap[j].id<<" ";
        cout<<endl;
    }
    cout<<endl;
}

bool comparator(Edge left, Edge right)
{
    return (left.weight > right.weight || (left.weight==right.weight && left.id > right.id));
}


bool comparatorE(EdgeE left, EdgeE right)
{
    return (left.weight > right.weight || (left.weight==right.weight && left.id > right.id));
}

bool heapComp(Info left, Info right)
{
    return (left.weight<right.weight || (left.weight==right.weight && left.id<right.id));
}


void Node::AddHeap(float wt, int idx )
{
    if(curSize==maxSize)
    {
        if(maxSize>2)
        {
            //heap[0].weight=wt;
            //heap[0].id=idx;

            /// Only heapify one branch of the heap tree
            int small,ri,li,pi=0;
            int done=0;
            
            if(heap[2].weight >heap[1].weight || (heap[2].weight==heap[1].weight && heap[2].id>heap[1].id))
                small=1;
            else
                small=2;
                
            if(wt>heap[small].weight || (wt==heap[small].weight && idx>heap[small].id))
            {
                heap[0].weight=heap[small].weight;
                heap[0].id=heap[small].id;
                heap[small].weight=wt;
                heap[small].id=idx;
            }   
            else
            {
                heap[0].weight=wt;
                heap[0].id=idx;
            }

            pi=small;
            while(!done)
            {
                li=2*pi+1;
                ri=2*pi+2;
                small=pi;

                if(li <maxSize && (heap[li].weight< heap[small].weight || (heap[li].weight == heap[small].weight && heap[li].id < heap[small].id )))
                    small=li;
                if(ri <maxSize && (heap[ri].weight< heap[small].weight || (heap[ri].weight == heap[small].weight && heap[ri].id < heap[small].id)))
                    small=ri;

                if(small != pi)
                {
                    wt=heap[pi].weight;
                    idx=heap[pi].id;

                    heap[pi].weight=heap[small].weight;
                    heap[pi].id=heap[small].id;

                    heap[small].weight=wt;
                    heap[small].id=idx;
                }
                else done=1;

                pi=small;
            }
        }
        else
        {
            if(maxSize != 1 && (wt>heap[1].weight || (wt==heap[1].weight && idx > heap[1].id)))
            {
                heap[0].weight=heap[1].weight;
                heap[0].id=heap[1].id;
                heap[1].weight=wt;
                heap[1].id=idx;
            }
            else
            {
                heap[0].weight=wt;
                heap[0].id=idx;
            }
        }

        minEntry.id=heap[0].id;
        minEntry.weight=heap[0].weight;
    }
    else
    {
        heap[curSize].weight=wt;
        heap[curSize].id=idx;
        curSize++;
        if(curSize==maxSize)
        {    
            sort(heap,heap+curSize,heapComp);
            minEntry.weight=heap[0].weight;
            minEntry.id=heap[0].id;
        }
    }            
}



int pivoting(Edge* verInd, int start, int end, int id, float weight, int* pindx,int step)
{
    int tstart,tend,r,p;
    Edge temp;
    
    tend=end-1;
    tstart=start;
    p=tstart;
    r=tend;
    while ( p < r )
    {
        while ( verInd[p].weight > weight || (verInd[p].weight==weight && verInd[p].id>id))
        {
            //if(p>tend)
                //break;
            p++;
        }
        while ( verInd[r].weight < weight || (verInd[r].weight==weight && verInd[r].id < id) )
        { 
            //if(r<tstart)
                //break;
            r--;
        }
        if(verInd[p].id==verInd[r].id)
            p++;
        else
            if( p < r ) 
            {
                /*temp.id= verInd[p].id;
                temp.weight=verInd[p].weight;

                verInd[p].id  = verInd[r].id;
                verInd[p].weight = verInd[r].weight;

                verInd[r].id = temp.id;
                verInd[r].weight = temp.weight;*/
                p++;
                r--;
            }
    }

    
    /*if(r<start)
        return start;
    else
    {
        if(r==tend)
        {
            sort(verInd+start,verInd+end,comparator);
            return end;
        }
        else
        {
            sort(verInd+start,verInd+r+1,comparator);
            return r+1;
        }
    }*/

    if(r<start)
        (*pindx)=start;
    else
    {
        if(r==tend)
            (*pindx)=end;
        else
            (*pindx)=r+1;
    }
    if(((*pindx)-start)>step)
        step=(*pindx)-start;
    return custom_sort(verInd,start,end,step,NULL,4);
}

int custom_sort(Edge* verInd, int start, int end, int step, int* order, int type)
{
    int part=start+step;
    int tstart, tend,tpart,k,length;
    int id,p,r,tid;
    float weight;
    Edge temp;

    switch(type)
    {
        case 1: break;
        case 2: part=end;
                sort(verInd+start,verInd+end,comparator);
                //csort(verInd,start,end);
                break;
        
        case 3: break;
        case 4: if(part>=end)
                {
                    part=end;
                    sort(verInd+start,verInd+end,comparator);
                }
                else
                {   
                    tend=end-1;
                    tstart=start;
                    k=step+1;
                    while(true)
                    {
                        p=tstart;
                        r=tend;
                        weight = verInd[r].weight;
                        id=verInd[r].id;
                        while ( p < r )
                        {
                            while ( verInd[p].weight > weight || (verInd[p].weight==weight && verInd[p].id>id))
                                p++;
                            while ( verInd[r].weight < weight || (verInd[r].weight==weight && verInd[r].id < id) )
                                r--;
                            
                            if(verInd[p].id==verInd[r].id)
                                p++;
                            else
                                if( p < r ) 
                                {
                                    temp.id= verInd[p].id;
                                    temp.weight=verInd[p].weight;

                                    verInd[p].id  = verInd[r].id;
                                    verInd[p].weight = verInd[r].weight;

                                    verInd[r].id = temp.id;
                                    verInd[r].weight = temp.weight;
                                }
                        }

                        length=r-start+1;
                        if(length==k)
                            break;
                        else
                        {
                            if(length > k)
                                tend=r-1;
                            else
                                tstart=r+1;
                        }
                        if(tstart==tend)
                            break;
                    }

                    //nth_element(verInd+start,verInd+part,verInd+end,comparator);
                    sort(verInd+start,verInd+part,comparator);
                }
                break;

        default: sort(verInd+start,verInd+end,comparator);
    }
    return part;
}

void bSuitor(CSR* G, int* b, Node* S, int algo, bool verbose)
{

    int numThreads, stepM=3,type=4;
    double t1;
    /*********** Memory Allocation *************/
    
    int* nlocks=new int[G->nVer];    //required for all schemes
    int* start=new int[G->nVer];
    int* end=new int[G->nVer];
    char* mark=new char[G->nEdge];    //required for unsorted and part sorted
    
    #pragma omp parallel
    numThreads=omp_get_num_threads();
    
    if(!verbose)
        t1=omp_get_wtime();

    //type 1,2,4 = DU DS DP
    if(algo==0)
        type=1;

    bSuitorBPQD(G,b,nlocks,S,start,end,mark,type,stepM,verbose);
    
    if(!verbose)
    {
        t1=omp_get_wtime()-t1;
        cout<<"Total Matching Time: "<<t1<<endl;
    }
    else
    {
        if(verifyMatching(G,S,G->nVer))
            cout<<"General Matching Verified..!!"<<endl;
        else
            cout<<"We are Screwed..!!"<<endl;
    }

    //printMatching(S,g.nVer);
    
    delete nlocks;
    delete start;
    delete end;
    delete mark;
    
    return;

}
