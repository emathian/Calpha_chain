
//  EMILIE MATHIAN



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

int begining(struct Coord *cp_coords, const int numAtoms);
 

/*
 * Declare an array to hold data read from the ATOM records of a PDB file.
 */

Atom    atom[MAX_ATOMS+1];
Coord   coord[MAX_ATOMS*5];


int read_data(filename)
    char    *filename;
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


void write_pdb_atom(serial, centre)
    int serial;
    Point   centre;
{
    printf("ATOM  %5d %8.3f%8.3f%8.3f\n",
                                serial,
                                centre.x,
                                centre.y,
                                centre.z);
}


int main(argc, argv)
    int argc;
    char    **argv;
{
    int numAtoms;
    int i;
    int j;
    double threshold_min = 2.5;
    double threshold_max = 5;
    double d;

        if ( argc<2 ) {
                (void) fprintf(stderr, "usage: atom_array file.pdb\n");
                exit(0);
        }

    numAtoms = read_data(argv[1]);
    printf("NUM ATOM  %d\n", numAtoms );
    for (i=1; i<=numAtoms; ++i) {
        write_pdb_atom(
            atom[i].serial,
            atom[i].centre);
    }

    /*for (i=1; i<=numAtoms; ++i) {
        printf("i %d  serial %d x %f y %f z %f \n", i ,atom[i].serial, atom[i].centre.x ,   atom[i].centre.y,  atom[i].centre.z);

    }*/

    int c =0;
    for (i=1; i<=numAtoms; ++i) {
        for (j=1; j<=numAtoms; ++j){
            d= dist_atom2(&atom[i], &atom[j]);
            //printf("%f\n", d);
            //printf("%f  %d %d \n", d ,atom[i].serial, atom[j].serial );
            if (atom[i].serial !=  atom[j].serial && d < threshold_max && d > threshold_min){
                    
                    printf("%f  %d %d \n", d ,atom[i].serial, atom[j].serial );
                    coord[c].i = i;
                    coord[c].j = j;
                    c++;
                    //printf(" max %f\n", threshold_max);
            
            }
        }
    }
    

  

    struct Coord cp_coords[MAX_ATOMS*5]; // Copy of seq

    for (int i = 0; i < c; ++i)
    {
       cp_coords[i].i = coord[i].i;
       cp_coords[i].j = coord[i].j;
    }

    for (int i = 0; i < c; ++i)
    {
       printf("Coords: %d / %d \n", cp_coords[i].i, cp_coords[i].j);
    }


    // Find begining
    printf("RES %d\n", begining(&cp_coords[0], numAtoms ));


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

int begining(struct Coord *cp_coords, const int numAtoms) 
{ 
    int posL = 0;
    int i =0;
    int stop =0;
    int occ[numAtoms];
    int res;
    while(stop == 0 && i <= numAtoms){ // c  = Nb atom
        if (i == 0){
            occ[posL] = 1; // Initialisation
        }
        else{
            if (cp_coords[i-1].i != cp_coords[i].i && cp_coords[i+1].i != cp_coords[i].i){
                res  = cp_coords[i].i;
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
    return res;
} 
