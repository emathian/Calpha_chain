
//  EMILIE MATHIAN



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stack>
#include <list>
#include <vector>
#include <iostream>

using namespace std;

#define MAX_ATOMS       10000
#define LINE_LENGTH     30



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

int extremity(struct Coord *cp_coords, const int len_list,  int i);
 

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
    double threshold_min = 2.8;
    double threshold_max = 4.8;
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
                    printf("%f  %d %d \n", d ,atom[i].serial, atom[j].serial );
                    coord[c].i = i;
                    coord[c].j = j;
                    c++;          
            }
        }
    }

    struct Coord cp_coords[MAX_ATOMS*5]; // Copy of seq

    for (int i = 0; i < c; ++i)
    {
       cp_coords[i].i = coord[i].i;
       cp_coords[i].j = coord[i].j;
    }


    // Find begining

    int beg_pos = extremity(&cp_coords[0], numAtoms , 0);
  
    stack <vector<struct Coord> > s; 

    // Initialize the stack
    vector<struct Coord>  v1;
   

    v1.push_back(cp_coords[beg_pos]);
    s.push(v1);

    vector<struct Coord>  v_res;
    vector<vector<int > > vv_sol;
    int n =0;
    while (!s.empty ()){
        vector<struct Coord> current_path;
        current_path = s.top();
        s.pop();
        int last_in = current_path.back().i;
        int last_out = current_path.back().j;
        int B = 0;
        for (int ii = 0; ii < c; ++ii)
        {
           
            if (cp_coords[ii].i == last_out && cp_coords[ii].j != last_in)
            {
                struct Coord new_ele;
                new_ele.i = cp_coords[ii].i;
                new_ele.j = cp_coords[ii].j;
                current_path.push_back(new_ele);
                s.push(current_path);
                B =1;
             
            }
        }
            if(B == 0){
    
                int fin = extremity(&cp_coords[0], c , beg_pos+1);
               vector<int> v_sol;
                for(int j = 0; j < current_path.size(); ++j){
                    v_sol.push_back(current_path[j].i);
                    printf("%d\n", current_path[j].i);
                }
                v_sol.push_back( cp_coords[fin].i);
                vv_sol.push_back(v_sol);
                printf("%d \n", cp_coords[fin].i);
                printf("n %d \n", n);
                n++;
            }
        }
      
        
        printf("HERE \n");
        for(int i = 0; i < vv_sol.size(); ++i){
            for (int j = 0; j <vv_sol[i].size(); ++j)
            {
                    printf("%d\n", vv_sol[i][j]);
            }
                
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

int extremity(struct Coord *cp_coords, const int len_list,  int i) 
{ 
    int posL = 0;
    //int i =0;
    int stop =0;
    int occ[len_list];
    while(stop == 0 && i <= len_list){ // c  = Nb atom
        if (i == 0){
            occ[posL] = 1; // Initialisation
        }
        else{
            if (cp_coords[i-1].i != cp_coords[i].i && cp_coords[i+1].i != cp_coords[i].i){
                //res  = cp_coords[i].i;
                stop = 1;
            }
            else{
                if(cp_coords[i].i == cp_coords[i-1].i ){
                    occ[posL]++;
                }
                else{
                    posL++;
                    occ[posL]++;
                }
            }
        }
        i++;
    }
    return i-1;
} 