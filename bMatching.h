
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
 
#ifndef BMATCHING_H
#define BMATCHING_H

#include <omp.h>
#include <cmath>
#include <iostream>
#include "mtxReader.h"
//#include <immintrin.h>
#include <getopt.h>
#include <float.h>
using namespace std;

#define CHUNK 256
#define DYN
#define BSIZE 256
#define LBUF
#define MAX_VAL FLT_MAX/2.0

extern int pivotID;
extern float pivotW;
extern int ptype;

struct Param
{
    int cstart;
    int cend;
    int ostart;
    int oend;
};

struct Info
{
    int id;
    float weight;
};

class Node 
{
    public:
    int maxSize;
    int curSize;
    Info* heap;
    Info minEntry;

    void print();
    int min_id()
    {
        return (curSize>0 && curSize==maxSize)?heap[0].id:-1;
    }

    float min_weight()
    {
        return (curSize>0 && curSize==maxSize)?heap[0].weight:0.0;
    }

    int find_id(int idx)
    {
        for(int i=0;i<curSize;i++)
            if(heap[i].id==idx)
                return 1;
        return 0;
    }

    void Add(float wt, int id);
    void AddHeap(float wt, int id);
    
};

struct Walker
{
    int midx1;
    int midx2;
    float W1;
    float W2;
    float* M1;
    float* M2;
};
//Edge gEdge;

void bSuitorBPQD(CSR* G, int *b, int *nlocks, Node* S, 
        int* start, int* end, char* mark, int type, int stepM, bool verbose);

void bSuitor(CSR* G, int* b, Node* S, int algo, bool verbose);

int custom_sort(Edge* verInd, int start, int end, int part, int* order, int typei);
bool verifyMatching(CSR* g, Node* S, int n);
void printMatching(Node* S, int n);

bool comparator(Edge left, Edge right);
bool comparator(Edge left, Edge right);

struct bmatching_parameters
{
    char* problemname;
    char* bFileName;
    int algorithm;
    int b;
    bool verbose;
    bool bipartite;

    bmatching_parameters();
    void usage();
    bool parse(int argc, char** argv);

};

#endif //BMATCHING_H
