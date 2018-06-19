This is a implementation of a half-approximation algorithm, b-Suitor, for computing a b-Matching
of maximum weight in a graph with weights on the edges. b-Matching is a generalization of the 
well-known Matching problem in graphs, where the objective is to choose a subset of M edges in the 
graph such that at most a specified number b(v) of edges in M are incident on each vertex v. Subject
to this restriction we maximize the sum of the weights of the edges in M.

CITATION: 

Khan, A., Pothen, A., Mostofa Ali Patwary, M., Satish, N. R., Sundaram, N., 
Manne, F., & Dubey, P. (2016). Efficient approximation algorithms for weighted b-matching. 
SIAM Journal on Scientific Computing, 38(5), S593-S619.

CONTACT:

Arif Khan,    			ariful.khan@pnnl.gov
Alex Pothen, 			apothen@purdue.edu
Mahantesh Halappanavar,	hala@pnnl.gov


TO COMPILE:

Make

INPUT:

The implementation accepts two inputs: 
i) a graph G in MATRIX MARKET FORMAT (.mtx) and 
ii) EITHER 
    a file that has n lines where n is the number of vertices in G. The i th lines contains
    a single integer representing the b value of i th vertex, i.e., b[i] value. 
    OR
    a constant b for all vertices. if b=0 then the 
    code generates b-value for each vertex v, randomly between 1 and sqrt(d(v)), 
    where d(v) is the degree of v. 

OUTPUT:

For each vertex v, a priority queue holds the final matching. 
The user are requested look into bMatching.h for the priority queue data structure, called NODE. 
An example of the implementation is included in the main function in bMatching.cpp.

HELP:

./bMatching -h 

will show the input argument for the code as follows:

Usage:  -f <problemname> -e <bfilename> -b <bval> -a <algorithm> -v -g

        -f problemname  : file containing graph. Currently inputs .mtx files
        -e bfilename    : Optional input. (currently not implemented)
        -b bval         : constant b value if b=0 then randomly generated.
        -a algorithm    : Algorithm 0:unsorted 1:partial sorted mode (defualt)
        -t              : bipartite graph 
        -v              : verbose

TEST EXAMPLE:

A test problem named Ga10As10H30, a quantum chemistry problem, from Florida Matrix Collection 
has been included with the codes. To experiment with the problem with constant b value = 2,
the user can issue the following command:

./bMatching -f turon_m.mtx -b 2 -v

THE CORRECT OUTPUT:
Graph (189924, 778531) Reading Done....!! took 1.35321
bValue is constant for each vertex
Input Processing Done: 0.010898
Initialization Done...!!
Start Matching189924
Iteration: 1 54135 0.004274
Iteration: 2 19842 0.00208
Iteration: 3 6572 0.001006
Iteration: 4 1703 0.000416
Iteration: 5 533 0.000182
Iteration: 6 93 0.000154
Iteration: 7 41 0.000203
Iteration: 8 13 0.000167
Iteration: 9 4 0.000136
Iteration: 10 1 0.000168
Iteration: 11 1 0.000108
Iteration: 12 0 0.000106
Matching Done....!!
Iterations: 12
Initialization Time: 0.006812
Matching Time: 0.009305
Total Time: 0.016117
Matching Weight: 109378
General Matching Verified..!!

The Iteration lines specifies how many vertices to be processed after that iteration and 
how much time it took for that iteration. The user can also run the same experiment
without the -v command which will generate less information.

Graph (189924, 778531) Reading Done....!! took 1.28027
bValue is constant for each vertex
Input Processing Done: 0.010949
Initialization Done...!!
Start Matching189924
Matching Done....!!
Total Matching Time: 0.015312
