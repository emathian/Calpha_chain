
//  EMILIE MATHIAN



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
#define Pseudo_valence_angle_min 28
#define Pseudo_valence_angle_max 105
#define NUMATOMS 327

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
 
float pseudo_valence_angle(struct Atom *A, struct Atom *B, struct Atom *C);

pair<int, vector<int> > bfs(int u,  int adjacency_matrix[][327],  struct Atom *seq ,const int numAtoms);

vector<int> longestPathLength( int adjacency_matrix[][327],  struct Atom *seq, int numAtoms);
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
    /*for (i=1; i<=numAtoms; ++i) {
        write_pdb_atom(
            atom[i].serial,
            atom[i].centre);
    }
    */
    for (i=1; i<=numAtoms; ++i) {
        printf("i %d  serial %d x %f y %f z %f \n", i ,atom[i].serial, atom[i].centre.x ,   atom[i].centre.y,  atom[i].centre.z);
    }

 




    int c =0;
    for (i=1; i<=numAtoms; ++i) {
        for (j=1; j<=numAtoms; ++j){
            d= dist_atom2(&atom[i], &atom[j]);
            if (atom[i].serial !=  atom[j].serial && d < threshold_max && d > threshold_min){
            
                    coord[c].i = i -1; // We go from 0 to numatom -1 //-1
                    coord[c].j = j -1;// We go from 0 to numatom -1 // -1
                    //printf("%f  %d %d \n", d, coord[c].i, coord[c].j );
                    c++;          
            }
        }
    }

    // Remove vertices between atom with unexpected pseudo valence angle
    struct Coord cp_coords[MAX_ATOMS*5]; // Copy of seq

    for (int i = 0; i < c; ++i)
    {
       cp_coords[i].i = coord[i].i;
       cp_coords[i].j = coord[i].j;
    }

     for (int i = 0; i < c; ++i){
        printf("(%d %d )\n", cp_coords[i].i , cp_coords[i].j );
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

    for (int i = 0; i < c; ++i){
    	int in = coord[i].i;
    	int out = coord[i].j;
    	adjacency_matrix[in][out] =1;
    }
  




    // Find begining
    //int farest_from0 = bfs(0, adjacency_matrix, numAtoms, 1);

   	vector<int> lp = longestPathLength(adjacency_matrix, &atom[0] ,numAtoms);
    for(int i=0; i< lp.size(); i++){
        printf("res %d \n", lp[i] +1 );
    }


   	
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
  double d  = sqrt(pow(atom1->centre.x - atom2->centre.x,2) + pow(atom1->centre.y - atom2->centre.y,2) + pow(atom1->centre.z - atom2->centre.z,2));
  return d;
} 


float pseudo_valence_angle(struct Atom *A, struct Atom *B, struct Atom *C){

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



pair<int, vector<int> > bfs(int u, int adjacency_matrix[][327], struct Atom *seq ,const int numAtoms) 
{ 
	// return the most distant point from u


    //  mark all distance with -1 
    int dis[numAtoms]; 
    for (int i = 0; i < numAtoms; ++i)
    {
    	dis[i] = -1;
    }

    queue<int> q; 
    q.push(u); 
    int vect_angle[3];
    vector<int> v_angle;
    v_angle.push_back(u);
    //  distance of u from u will be 0 
    printf("u %d\n",u );
    dis[u] = 0; 
  	int length = 0;
	vector<int> v;
    while (!q.empty()) 
    { 
        int t = q.front();       q.pop(); 
        v.push_back(t);
        printf("Top  %d \n", t);
  		
        //  loop for all adjacent nodes of node-t 
        for (int i = 0; i < numAtoms; ++i)
        {
        	if(adjacency_matrix[t][i] == 1){
        		int neigh = i;
                printf("neigh %d\n", neigh );
                v_angle.push_back(neigh);
                if(v_angle.size() >= 3){
                    for (int i= 0 ; i< v_angle.size(); ++i){
                        printf("vangle %d  ", v_angle[i] );
                    }
                    printf("\n");
                    printf("v_angle2 %lu %lu %lu \n", v_angle.size()-3, v_angle.size()-2, v_angle.size()-1);
                   // printf("v_angle2 %d %d %d \n", v_angle[v_angle.size()-2],v_angle[v_angle.size()-1], v_angle[v_angle.size()] );
                    printf("Serial %d %d %d \n", seq[v_angle[v_angle.size()-3]+1].serial , seq[v_angle[v_angle.size()-2]+1].serial, seq[v_angle[v_angle.size()-1]+1].serial);
                    float ang = pseudo_valence_angle(&seq[v_angle[v_angle.size()-3]+1] , &seq[v_angle[v_angle.size()-2]+1], &seq[v_angle[v_angle.size()-1]+1]);
                    printf("ang %f \n", ang );
                }
        		if(dis[i] == -1){
        			 q.push(neigh); 
        			 dis[neigh] = dis[t] + 1; 
        		}
        	}
        }
        printf("LOOP WHILE  \n\n");
    }

       
    int maxDis = 0; 
    int nodeIdx; 
  
    //  get farthest node distance and its index 
    for (int i = 0; i < numAtoms; i++) 
    { 
        if (dis[i] > maxDis) 
        { 
            maxDis = dis[i]; 
            nodeIdx = i; 
        } 
    } 

   printf("\n \n \n");
    return make_pair(nodeIdx, v); 
} 

vector<int> longestPathLength( int adjacency_matrix[][327], struct Atom *seq, int numAtoms) 
{ 

  
    // first bfs to find one end point of 
    // longest path 
    pair<int, vector<int> >  t1_pair = bfs(0, adjacency_matrix, &seq[0], numAtoms); 
    int t1 = t1_pair.first;
  
    //  second bfs to find actual longest path 
    pair<int, vector<int> > t2_pair = bfs(t1, adjacency_matrix, &seq[0], numAtoms);
    return t2_pair.second;
} 

