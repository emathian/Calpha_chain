
//  EMILIE MATHIAN
// Compilation:	g++ -o ARP_WARPQ2 ARP_WARPQ2.cpp 
// Execution: ./ARP_WARPQ2 data_q2.txt


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stack>
#include <list>
#include <vector>
#include <queue>
#include <utility>
#include <iostream>

using namespace std;

#define MAX_ATOMS       10000
#define LINE_LENGTH     30
#define Pseudo_valence_angle_min 0
#define Pseudo_valence_angle_max 110
#define NUMATOMS 330

typedef struct {
    double  x, y, z;
} Point;


typedef struct Coord{
    int  i, j;
} Coord;


typedef struct Atom{
    int serial;
    Point   centre;
} Atom;


char* substring(const char* str, size_t begin, size_t len);

double dist_atom2(struct Atom *atom1, struct Atom *atom2);
 /* Returns the Euclidian distance between two atoms. */

float pseudo_valence_angle(struct Atom *A, struct Atom *B, struct Atom *C);
/* Returns the angle according between three atoms.  */


pair<int, vector<int> > bfs(int u,  int adjacency_matrix[][330],  struct Atom *seq ,const int numAtoms);
/*Breadth First Traversal: traversing algorithm.  */

vector<int> longestPathLength( int adjacency_matrix[][330],  struct Atom *seq, int numAtoms);
/* Returns the longest path in a graph.*/


void dfs(int numAtoms, int adjacency_matrix[][330],  struct Atom *seq);
/* Depth-First Search: Tansversal algorithm.*/


void dfsUtil(int u, int numAtoms,int adjacency_matrix[][330],  struct Atom *seq, bool visited[], 
             vector<pair<int, int> > road_used, int parent, int it, vector<int> v_sol);
/* Depth-First Search Recursive loop (called by dfs).*/

/*
 * Declare an array to hold data read from the ATOM records of a PDB file.
 */

Atom    atom[MAX_ATOMS+1];
Coord   coord[MAX_ATOMS*5];


int read_data(const char* filename)
{
    FILE    *stream;
    char    line[LINE_LENGTH];
        char    s_serial[8];
        char    s_x[8];
        char    s_y[8];
        char    s_z[8];

        int     serial;
        double  x;
        double  y;
        double  z;

    int i=0;

        if ( (stream = fopen(filename, "r")) == NULL ) {
                (void) fprintf(stderr, "Unable to open %s\n", filename);
                exit(0);
        }

        while ( fgets(line, LINE_LENGTH, stream) ) {
                
    
                         /* Split the line into its constituent fields.
                         * We are only interested in columns 1-54.
                         */

                        strncpy(s_serial,  &line[0],  7); s_serial[7]  = '\0';
                        strncpy(s_x,       &line[7], 7); s_x[7]       = '\0';
                        strncpy(s_y,       &line[14], 7); s_y[7]       = '\0';
                        strncpy(s_z,       &line[20], 7); s_z[7]       = '\0';

                        /*
                         * Convert the numeric fields to integers or doubles.
                         * The library functions atoi() and atof() are
                         * described in the UNIX manual pages ('man atoi' and
                         * 'man atof').
                         */

                        serial = atoi(s_serial);
                       // resSeq = atoi(s_resSeq);
                        x      = atof(s_x);
                        y      = atof(s_y);
                        z      = atof(s_z);

            /*
             * Copy values to the next element in the atom array.
             */

            if ( ++i > MAX_ATOMS ) {
                (void) fprintf(stderr, "Too many atoms read\n");
                exit(0);
            }
            atom[i].serial = serial;
            atom[i].centre.x = x;
            atom[i].centre.y = y;
            atom[i].centre.z = z;
        }
    //  }
    //}
    return i;
}


void write_pdb_atom(int serial, Point centre)
{
    printf("ATOM  %5d %8.3f%8.3f%8.3f\n",
                                serial,
                                centre.x,
                                centre.y,
                                centre.z);
}


int main(int argc, char ** argv)
{
    int numAtoms;
    int i;
    int j;
    double threshold_min = 3.3;
    double threshold_max = 4.3;
    double d;


    if ( argc<2 ) {
            (void) fprintf(stderr, "usage: atom_array file.pdb\n");
            exit(0);
    }

    numAtoms = read_data(argv[1]);

    // Creation of a list of adjacency
    int c =0;
    for (i=1; i<=numAtoms; ++i) {
        for (j=1; j<=numAtoms; ++j){
            d= dist_atom2(&atom[i], &atom[j]);
            if (atom[i].serial !=  atom[j].serial && d < threshold_max && d > threshold_min){
            
                    coord[c].i = i -1; // We go from 0 to numatom -1 //-1
                    coord[c].j = j -1;// We go from 0 to numatom -1 // -1
                    printf("%f  %d %d \n", d, coord[c].i, coord[c].j );
                    c++;          
            }
        }
    }
    
    
    // Initialisation of an adjacency matrix
    int adjacency_matrix[NUMATOMS][NUMATOMS];
    for (int i = 0; i < numAtoms; ++i)
    {
    	for (int j = 0; j <numAtoms; ++j)
    	{
    		adjacency_matrix[i][j] =0;
    	}
    }


    // Fill the adjacency matrix
    for (int i = 0; i < c; ++i){
    	int in = coord[i].i;
    	int out = coord[i].j;
    	adjacency_matrix[in][out] =1;
    }
  




    // Find the longest chain of alpha-carbon
    // !!! I DIDN'T MANAGE TO RETURN A VECTOR OF U BY DFSUTILS (Illegal instruction: 4)
    /// THAT IS WHY EACH U ADDED TO THE QUEUE IS PRINTED. 
    /// THE RESULTS ARE WRITTEN IN A FILE AND THE LONGEST VECTOR IS SEARCHED ON R. !!!!
    dfs(numAtoms, adjacency_matrix,  &atom[0]);
   
        return 0; 
}


char* substring(const char* str, size_t begin, size_t len) 
{ 
  if (str == 0 || strlen(str) == 0 || strlen(str) < begin || strlen(str) < (begin+len)) 
    return 0; 

  return strndup(str + begin, len); 
} 


double dist_atom2(struct Atom *atom1, struct Atom *atom2) 
{ 
	/* Returns the Euclidian distance between two atoms. */
  double d  = sqrt(pow(atom1->centre.x - atom2->centre.x,2) + pow(atom1->centre.y - atom2->centre.y,2) + pow(atom1->centre.z - atom2->centre.z,2));
  return d;
} 


float pseudo_valence_angle(struct Atom *A, struct Atom *B, struct Atom *C){
	/* Returns the angle according between three atoms.  */

    // AB = B - A
    float AB[3];
    AB[0] = B->centre.x  - A->centre.x;
    AB[1] = B->centre.y  - A->centre.y;
    AB[2] = B->centre.z  - A->centre.z;
    // BC = C -B
    float BC[3];
    BC[0] = C->centre.x  - B->centre.x;
    BC[1] = C->centre.y  - B->centre.y;
    BC[2] = C->centre.z  - B->centre.z;


    // Dot product AB.AC
    float prod_scal = AB[0] * BC[0] + AB[1] * BC[1] + AB[2] * BC[2];

    // Norm AB
    float norm_AB = sqrt(pow(AB[0],2) + pow(AB[1],2) + pow(AB[2],2));
   
    // Norm AC
    float norm_BC = sqrt(pow(BC[0],2) + pow(BC[1],2) + pow(BC[2],2));

    // Angle
    double p =  prod_scal / (norm_AB * norm_BC);
    float theta = acos(p)*(180/M_PI); // Degree


    return theta;
}



pair<int, vector<int> > bfs(int u,  int adjacency_matrix[][330],  struct Atom *seq ,const int numAtoms)
{ 
	/*Breadth First Traversal: traversing algorithm. 
  // Arguments:
    int u : initial point
    adjacency_matrix,
    numAtoms: number of atoms  
  BFS is a tranversing algorithm.
  */


    //  mark all distance with -1 -> Not visited
    int dis[numAtoms]; 
    for (int i = 0; i < numAtoms; ++i)
    {
    	dis[i] = -1;
    }
    // Initialisation of the queue
    queue<int> q; 
    q.push(u); 
  
    //  distance of u from u will be 0 
    dis[u] = 0; 

    // Initialisation of the output vector
  	vector<int> v;
    while (!q.empty()) 
    { 
        int t = q.front();       q.pop(); 
        // Add the new node (atom)
        v.push_back(t);

        //  loop for all adjacent nodes of t 
        for (int i = 0; i < numAtoms; ++i)
        {
        	if(adjacency_matrix[t][i] == 1){
        		int neigh = i;
        		if(dis[i] == -1){ // If the node hasn't been visited 
        			 q.push(neigh); // Add the neighbors of t to the queue
        			 dis[neigh] = dis[t] + 1; 
        		}
        	}
        }
    }

    int maxDis = 0; 
    int nodeIdx; 
  
    //  get farthest node distance and from u
    for (int i = 0; i < numAtoms; i++) 
    { 
        if (dis[i] > maxDis) 
        { 
            maxDis = dis[i]; 
            nodeIdx = i; 
        } 
    } 

    return make_pair(nodeIdx, v); 
} 


vector<int> longestPathLength( int adjacency_matrix[][330],  struct Atom *seq, int numAtoms)
{ 
  /*longestPathLength: 
  // Arguments:
    adjacency_matrix,
    numAtoms: number of atoms  
  In order to find the longest path, the BFS algorithm is calculate for the point 0. 
  Then a second BFS is calculated from the fartest point to 0. This method assures to find the 
  longest path.
  */
  
    // first bfs to find one end point of 
    // longest path 
    pair<int, vector<int> >  t1_pair = bfs(0, adjacency_matrix, &seq[0],numAtoms); 
    int t1 = t1_pair.first;
    //  second bfs to find actual longest path 
    pair<int, vector<int> > t2_pair = bfs(t1, adjacency_matrix, &seq[0], numAtoms);
    return t2_pair.second;
} 



void dfsUtil(int u, int numAtoms, int adjacency_matrix[][330],  struct Atom *seq, bool visited[], 
             vector<pair<int, int> > road_used, int parent, int it, vector<int> v_sol ) 
{ 
	/* Depth-First Search Recursive loop (called by dfs).
	Argument:
	-u current node
	- numAtoms number of atoms
	- adjacency_matrix
	- *seq : list of atoms
	- visited: a node is visited when it has ever been find by dfs.
	- road_used: vector conting each vertex visited.
	- it : iterator
	- v_sol : vector containing each element of the queue. I DIDN'T MANEGE TO RETURN THIS VECTOR*/

    int c = 0; 

    // !!!! THIS STOP CONDITION IS FAR FROM BE PERFECT. IT ALLOWS TO STOP THE RECURSIVE PROCESS WHEN THE ALGORITHM
    // REPEATS IS "BLOCKED" IN A VERTEX !!!
    if (road_used.size() > 4 && road_used[road_used.size() -3].first == road_used[road_used.size() -2].second  
                    && road_used[road_used.size() -3].first == road_used[road_used.size() -1].first 
                    && road_used[road_used.size() -2].second == road_used[road_used.size() -1].first &&
                    road_used[road_used.size() -1].first == road_used[road_used.size() -4].second &&
                    road_used[road_used.size() -2].second == road_used[road_used.size() -4].second &&
                     road_used[road_used.size() -3].first == road_used[road_used.size() -4].second  ){
        for (int i = 0; i < numAtoms; i++) {
            visited[i] = true;
            c++; 
        }
        c--; 
    }
    if (c == numAtoms -1){
    	// Print the solution
        printf("VSOL\n");
        for (int i = 0; i < v_sol.size(); ++i)
        {
            printf("%d \n", v_sol[i] );
        }

        return;
    }
    
    // Mark u as visited
    visited[u] = true; 
  
    // Track the current edge 
    pair<int, int> parent_u = make_pair(parent, u);
    road_used.push_back(parent_u); 

    // Initialisation of the vector solution
    v_sol.push_back(u+1);
   
    // Check for not visited node and proceed with it. 
    for (int i = 0; i < numAtoms; ++i){
        if(adjacency_matrix[u][i] == 1){
            if (!visited[i]){
                 // Calculation of the pseudo valence angle between the current vertex and the following node (i);
                if (road_used.size() > 3 && parent_u.first != parent_u.second
                    && parent_u.second != i
                    && parent_u.second != parent_u.first){
                    
                   float ang  = pseudo_valence_angle(&seq[parent_u.first +1 ], &seq[parent_u.second +1 ], &seq[i + 1 ]);

                   if (ang > Pseudo_valence_angle_min && ang < Pseudo_valence_angle_max){
                    dfsUtil(i, numAtoms, adjacency_matrix, &seq[0], visited, road_used, u, it + 1, v_sol); 
                    
                   }
                }
                else{ // There is less than three elements in road_used
                    dfsUtil(i, numAtoms, adjacency_matrix, &seq[0], visited, road_used, u, it + 1, v_sol); 
                }
               
            }
        }
    }
    // Backtrack 

    for (int i =0; i < road_used.size(); ++i){ 
        pair<int, int> y = make_pair(road_used[i].first, road_used[i].second);
        if (y.second == u) {
       
     
            dfsUtil(y.first, numAtoms, adjacency_matrix, &seq[0] ,visited, 
                   road_used, u, it + 1, v_sol); 
        }
    }
}  

  

void dfs(int numAtoms, int adjacency_matrix[][330],  struct Atom *seq) 
{ 
  /* Depth-First Search Recursive loop (called by dfs).*/

    // Create a array of visited ndoe 
    bool visited[numAtoms]; 
  
    // Vector to track last visited road 
    vector<pair<int, int> > road_used; 
  
    // Initialize all the node with false 
    for (int i = 0; i < numAtoms; i++) 
        visited[i] = false; 
  
    // call the function bfs to find the extremity
    pair<int, vector<int> >  t1_pair = bfs(0, adjacency_matrix, &seq[0], numAtoms); 
    int t1 = t1_pair.first;
    // From this point call dfsUtils
    vector<int> v_sol;
    dfsUtil(t1, numAtoms, adjacency_matrix, &seq[0] ,visited, road_used, -1, 0, v_sol); 

} 


