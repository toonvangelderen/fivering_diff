#include <stdlib.h>
#include "path.h"



void initIntArray(IntArray *a, size_t initialSize) {
  //https://stackoverflow.com/questions/3536153/c-dynamically-growing-array
  //size_t initialSize is the intial size of the array. for example "5" intergers, do not keep it too small, as the array is doubled if it is too short 
  // a->used is the number of used entries, upon initialization it is zero.
  
  
  if(((a->ints = (int *) malloc(initialSize * sizeof(int) ))==NULL)) {;
       error("could not malloc initialSize * sizeof(int) ");
    }
  
  a->used = 0;
  a->size = initialSize;
  return;

}

void insertIntArray(IntArray *a, int element) {
  // the original array a is elongated if you want to add more elements that is currently has. Just for easiness it is doubled. 
  // NOT handy when the array becomes big!
  // a->used is the number of used entries, because a->ints[a->used++] updates a->used only *after* the array has been accessed.

  // Therefore a->used can go up to a->size 
  if (a->used == a->size) {
    a->size *= 2;
    if( (a->ints = realloc(a->ints, a->size * sizeof(int)))==NULL) {;
       error("could not realloc IntArray->ints ");
    }
    // a->ints = realloc(a->ints, a->size * sizeof(int)); // reallocate the elongated array with the new size.
  }
  a->ints[a->used++] = element;
  return;

}

void freeIntArray(IntArray *a) {
  free(a->ints); //free the memory
  a->ints = NULL;// make pointer zero
  a->used = a->size = 0;   
  return;

}

void removeIthElementIntArray(IntArray *a, int elemi ){
  // elemi stands for the ith element 
  // the array a has a->used  elements currently
  int start =0,i;

  if (a->used < elemi ){
    printf("the max length of the int arrat is %zu, but you want to remove element %d\n",a->used,elemi );
    error("Not possible");
  }

  a->used--;
  for (i=elemi;i<a->used;i++){
      a->ints[i]=a->ints[i+1];
  }
  return;
}

void removeElementXIntArray(IntArray *a, int elemX ){
  // elemX stands for element number X
  // the array a has a->used  elements currently
  int start =0,i;
  int last_int=(int)a->used-1;
  ///because you remove here one on forehand,  a[i]=a[i+1]; will not give problems
  // dprint((int)a->used);
  
  if(a->ints[last_int]==elemX){
      a->used--;
      return;
  }
  else {
    a->used--;
    for (i=0;i<a->used;i++){
      if(a->ints[i]==elemX){
        start=1;
      }
      if (start){
          a->ints[i]=a->ints[i+1];
      }
    }
  }
  
  if (start==0   ){
    // dprint((int)a->used+1);
    dprint(a->ints[last_int]);
    printf("you want to remove elemnt %d \n",elemX);
    printf("the remaining elements are :\n");

    for (i=0;i<a->used;i++){
      printf(" i = %d  ",i);
      dprint(a->ints[i]);
    }
    error("you are trying to remove an element that is not in the list.");
  }

  return;

}

int checkElementXIntArray(IntArray *a, int elemX ){
  // check if elementX is in the array
  // the array a has a->used  elements currently
  int check =0,i;

  for (i=0;i<a->used;i++){
      if(a->ints[i]==elemX){
        check=1;
        break;
      }      
  }

  return check;

}

void printIntArray(IntArray *a ){
  // check if elementX is in the array
  // the array a has a->used  elements currently
  int i;
  for (i=0;i<a->used;i++){
    printf(" item %d has value %d\n", i,a->ints[i]);
     
  }

  return ;

}



/*void memory_allocation(){
  //depending on the type of analysis, method estc you want to create variables 
  //analysis types
  Average *rdf_average, *globalCN_average;

  if (sys.rdfanalysis){
    if((rdf_average = (Average *)malloc(RDFBINS * sizeof(Average))==NULL)) {;
       error("could not malloc rdf_average  RDFBINS* sizeof(Average) ");
    }
    if((globalCN_average = (Average *)malloc(  sizeof(Average))==NULL)) {;
       error("could not malloc globalCN_average   sizeof(Average) ");
    }
  }
  // now these pointers should be passed to somewhere


  int *s_histogram,*r_S_2dhist,*xy_plane_patchangle;
  if(sys.s_distribution){
      if((s_histogram = (int *)malloc(MAXBIN * sizeof(int))==NULL)) {;
         error("could not malloc s_histogram  MAXBIN* sizeof(int) ");
      }

      if((r_S_2dhist = (int *)malloc(MAXBIN * MAXBIN * sizeof(int))==NULL)) {;
         error("could not malloc r_S_2dhist  MAXBIN* MAXBIN *sizeof(int) ");
      }
  }


  
  Statistics bond_length,msd,Sfraction;

  //MC types
  Mcacc  trans,        fintrans,        rot,        finrot;
  Mcacc  trans_cluster,fintrans_cluster,rot_cluster,finrot_cluster ;
  Mcacc  trans_mono,   fintrans_mono,   rot_mono,   finrot_mono;

  //BMD
  if(sys.sim_type==0){ 
    Langevin langevin;
  }


  //Emilio
  Constraint cconst;
  int crossing_check;
  Trajectory segment;


  //advanced sampling: TPS / TIS / replica exchange
  State state[MAXSTATES];
  Replica *replica[MAXREPLICA];
  Path     path;
  Slice *init_state,*fwtrial,*bwtrial,*trial;

  //nearest neighbor variables
  List list;
}
*/

