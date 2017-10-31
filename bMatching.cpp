
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
#include "mtxReader.h"
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <float.h>

using namespace std;


bmatching_parameters::bmatching_parameters():problemname(NULL),bFileName(NULL),b(5),verbose(false),
    bipartite(false)
{}

void bmatching_parameters::usage()
{
    const char *params =
	"\n\n"
    "Usage: %s -f <problemname> -e <bfilename> -b <bval> -a <algorithm> -v -g\n\n"
	"	-f problemname  : file containing graph. Currently inputs .cbin files\n"
	"	-e bfilename    : Optional input. (currently not implemented)\n"
	"	-b bval         : constant b value if b=0 then randomly generated.\n"
	"	-a algorithm    : Algorithm 0:unsorted 1:partial sorted mode\n"
    "   -t              : bipartite graph \n";
    "   -v              : verbose \n\n";
    fprintf(stderr, params);
}

bool bmatching_parameters::parse(int argc, char** argv)
{
    static struct option long_options[]=
    {
        // These options don't take extra arguments
        {"verbose", no_argument, NULL, 'v'},
        {"help", no_argument, NULL, 'h'},
        {"bipartite", no_argument, NULL, 'g'},
        
        // These do
        {"problem", required_argument, NULL, 'f'},
        {"bVar", required_argument, NULL, 'e'},
        {"b", required_argument, NULL, 'b'},
        {"a", required_argument, NULL, 'a'},
        {NULL, no_argument, NULL, 0}
    };

    static const char *opt_string="vhtf:e:b:a:";
    int opt, longindex;
    opt=getopt_long(argc,argv,opt_string,long_options,&longindex);
    while(opt != -1)
    {
        switch(opt)
        {
            case 'v': verbose=true; 
                      break;
            
            case 't': bipartite=true; 
                      break;

            case 'h': usage(); 
                      return false; 
                      break;

            case 'f': problemname=optarg; 
                      if(problemname==NULL)
                      {
                        cerr<<"Problem file is not speficied"<<endl;
                        return false;  
                      }
                      break;

            case 'e': bFileName=optarg;
                      if(bFileName==NULL)
                      {
                        cerr<<"Problem file is not speficied"<<endl;
                        return false;  
                      }
                      break;
            case 'b': b=atoi(optarg);
                      if(b<0)
                      {
                        cerr<<"the b value can't be negative"<<endl;
                        return false;
                      }
                      break;
            case 'a': algorithm=atoi(optarg);
                      if(algorithm<0)
                      {
                        cerr<<"the a value can't be negative"<<endl;
                        return false;
                      }
                      break;
        }
        opt=getopt_long(argc,argv,opt_string,long_options,&longindex);
    }

    return true;
}
int main(int argc, char** argv)
{
    bmatching_parameters opts;
    if(!opts.parse(argc,argv))
    {
        cout<<"Argument Parsing Done..!!"<<endl;
        return -1;
    }
    double t1;
    int numThreads,best;
    double rt_start = omp_get_wtime();	
    
    /********* Reading the INPUT ******/
    CSR G;
    //G.readCSRbin(opts.problemname,0);
    if(opts.bipartite==false)
        G.readMtxG(opts.problemname);
    else
        G.readMtxB(opts.problemname);

    double rt_end = omp_get_wtime();	
    cout<<"Graph (" << G.nVer << ", " << G.nEdge/2 << 
            ") Reading Done....!! took " << rt_end - rt_start <<endl;
    
    /*********** Memory Allocation *************/
    
    int *b=new int[G.nVer];
    Node* S=new Node[G.nVer];      //Heap data structure
    
    
    #pragma omp parallel
    numThreads=omp_get_num_threads();
    
    
    /************ Assigning bValues **************/
    int avgb=0;
    for(int i=0;i<G.nVer;i++)
    {
        if(opts.b>0)
        {    
            int card=0;
            for(int j=G.verPtr[i];j<G.verPtr[i+1];j++)
                if(G.verInd[j].weight>0)
                    card++;

            if(card>opts.b)
                b[i]=opts.b;
            else b[i]=card;
        }
        else
        {
            //int deg=(int)sqrt(G.verPtr[i+1]-G.verPtr[i]);
            int deg=0;
            for(int j=G.verPtr[i];j<G.verPtr[i+1];j++)
                if(G.verInd[j].weight>0)
                    deg++;
            deg=(int)sqrt(deg);
            if(deg==0)
                b[i]=0;
            else
            {
                if(deg==1)
                    b[i]=1;
                else
                    b[i]=floor(deg);
                    //b[i]=rand()%deg+1;
            }
        }
        avgb+=b[i];
    }
    
    avgb=avgb/G.nVer;

    
    if(opts.b>0)
        cout<<"bValue is constant for each vertex"<<endl;
    else
        cout<<"bValue is randomly generated, b_avg= "<<avgb<<endl;
    

    for(int i=0;i<G.nVer;i++)       
        if(b[i]>0)
            S[i].heap=new Info[b[i]];      //Each heap of size b
        else
            S[i].heap=new Info[1];      //Each heap of size b
        
    
    cout << "Input Processing Done: " << omp_get_wtime() - rt_end <<endl;	

    bSuitor(&G,b,S,opts.algorithm,opts.verbose);
    
	if(opts.verbose)
	{
		ofstream fout;
		fout.open("matching.net",ios::out);
		if(fout.is_open())
		{
			fout<<"*Vertices "<<G.nVer<<endl;
			fout<<"*arcs"<<endl;
			for(int i=0;i<G.nVer;i++)
				for(int j=0;j<S[i].curSize;j++)
					if(i>S[i].heap[j].id)
						fout<<i+1<<" "<<S[i].heap[j].id+1<<" "<<S[i].heap[j].weight<<endl;

		}
		fout.close();
	}
    for(int i=0;i<G.nVer;i++)
        delete S[i].heap;
    
    delete S;
    
    return 0;
}
