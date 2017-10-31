
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
 
#ifndef MTXREADER_H
#define MTXREADER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cassert>
using namespace std;

struct EdgeE
{
    int head;
    int id;         // Edge tail
    float weight;  // Edge weight
};
struct Edge
{
    int id;         // Edge tail
    float weight;  // Edge weight
};

class CSR
{
    public:
    int nVer;       // number of vertices 
    int nEdge;      // number of edges
    int maxDeg;
    int* verPtr;    // vertex pointer array of size nVer+1
    Edge* verInd;   // Edge array
    bool bipartite; // general graph or bipartite graph
    int rVer;       // The number of vertices on Right for bipartite graph;
    int lVer;       // The number of vertices on left for bipartite graph;
    
    bool readMtxG(char * filename); // reading as a general graph
    bool readMtxB(char * filename); // reading as a bipartite graph
    bool mtxG2csrbin(char* filename, char* outfile); //converting mtx to binary file
    bool mtxB2csrbin(char* filename, char* outfile); //converting mtx to binary file
    bool readCSRbin(char* filename, int opt); // reading binary file of csr format
    bool readCSRbinBipartite(char* filename, int opt); // reading binary file of csr format
    
    CSR():nVer(0),nEdge(0),verPtr(NULL),verInd(NULL),bipartite(false){}
    ~CSR()
    {
        if(verPtr!=NULL)
            delete verPtr;

        if(verInd!=NULL)
            delete verInd;
    }

};

#endif //MTXREADER_H
