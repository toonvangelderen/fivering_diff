#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"


/*SETUP / UPDate in constructnn.c*/

// List list;  

// void setup_nnlist();
// void print_cell_list();


void setup_nnlist(){
  /* name a bit misleeding? it creates the cells, not a nn list? */
  int ix,iy,iz,k;

  printf("minimum cellsize = %g\n", MINCELLSIZE ); 
  printf("maximum cells in x =%d\n",NCELLX );
  
  list.cutoff=3.5;
  list.cutoff=4; /* ??*/
  list.cutoff=2; /* rc Simons Potential = +- 1.1 sigma* */

  list.cutoff2= list.cutoff * list.cutoff;
  list.cellx = sys.boxl.x / (0.5*list.cutoff);  /*nr of cells in x,y,z direction (int). determined by the cutoff*/
  list.celly = sys.boxl.y / (0.5*list.cutoff);
  list.cellz = sys.boxl.z / (0.5*list.cutoff);
  
  list.cellsize.x  =sys.boxl.x / list.cellx; /* the size of the cell= boxlength/ (int) nr cells */
  list.cellsize.y  =sys.boxl.y / list.celly;
  list.cellsize.z  =sys.boxl.z / list.cellz;
  list.inv_cellsize.x = 1./ list.cellsize.x; 
  list.inv_cellsize.y = 1./ list.cellsize.y;
  list.inv_cellsize.z = 1./ list.cellsize.z;
  
  dprint( list.cellx);
  dprint( list.celly);
  dprint( list.cellz);
  vprint(sys.boxl);
  vprint( list.cellsize);
  vprint( list.inv_cellsize);


  list.cutoff_buffer2 = 0.5*(list.cutoff-sys.rcutoff); 
  if(list.cutoff_buffer2<=0){
    error("buffer size too small\n");
  }

  gprint(list.cutoff_buffer2);
  list.cutoff_buffer2 *=  list.cutoff_buffer2; 
  gprint(list.cutoff_buffer2);
  // list.bigpart_cutoff2 =2*list.cutoff_buffer2; 

   for(ix=-1;ix<= 1;ix++) 
    for(iy=-1;iy<= 1;iy++) 
      for(iz=-1;iz<= 1;iz++) {
        k = iz + 3*iy + 9*ix + 13;
        list.trans[k].x =ix*sys.boxl.x;
        list.trans[k].y =iy*sys.boxl.x;
        list.trans[k].z =iz*sys.boxl.x;
        vprint(list.trans[k]);
      }

}

void put_parts_in_box(Slice *psl){
  vector dr,r;
  int i; 

  /* r = psl->pts[0].r; 
  dr=nulvec;
  if (r.x < -0.5) dr.x=+1;
  if (r.x >  0.5) dr.x=-1;
  if (r.y < -0.5) dr.y=+1;
  if (r.y >  0.5) dr.y=-1;
  if (r.z < -0.5) dr.z=+1;
  if (r.z >  0.5) dr.z=-1;
  for(i=0;i<sys.nsolute;i++) {
    vector_add(psl->pts[i].r,dr,psl->pts[i].r);
  }*/

  /*  if (sys.colorforce==1) {
 
    for(i=0;i<sys.npart;i++) {
      r = psl->pts[i].r; 
      if (r.x >  0.5) { 
        r.x-=1;
        if (psl->pts[i].type==0) if (ran3() < sys.prob_colorchange) psl->pts[i].type =1;
      }
      if (r.x < -0.5) {
        r.x+=1;
        if (psl->pts[i].type==1) if (ran3() < sys.prob_colorchange) psl->pts[i].type =0;
      }
      if (r.y >  0.5) r.y-=1;
      if (r.y < -0.5) r.y+=1;
      if (r.z >  0.5) r.z-=1;
      if (r.z < -0.5) r.z+=1;
      //    periodic_boundaries(r);
      psl->pts[i].r =r;
      psl->pts[i].dr =nulvec;
    }
    } else {  */

      for(i=0;i<sys.npart;i++) {
        r = psl->pts[i].r; 
        pbc(r,sys.boxl); 
        psl->pts[i].r =r;
        psl->pts[i].dr =nulvec;  
      }

      //}
}

void create_cell_list(Slice *psl){
  /* checks in which cell all particles lie, adds it to the stack and counts the number of particles in the cells*/
  int ix,iy,iz,i,n;
  double half_boxsizex=0.5*sys.boxl.x,half_boxsizey=0.5*sys.boxl.y,half_boxsizez=0.5*sys.boxl.z;

  for(ix=0;ix< list.cellx;ix++) /*set the number particles n in the cells to zero*/
    for(iy=0;iy< list.celly;iy++) 
      for(iz=0;iz< list.cellz;iz++) 
        list.cell[ix][iy][iz].n =0;
  
  for (i=0;i<sys.npart;i++) { /*In which cell does the particle lie?  */
    ix = (int)((psl->pts[i].r.x+half_boxsizex) * list.inv_cellsize.x );
    iy = (int)((psl->pts[i].r.y+half_boxsizey) * list.inv_cellsize.y );
    iz = (int)((psl->pts[i].r.z+half_boxsizez) * list.inv_cellsize.z );

    if ((ix<0) || (ix>=list.cellx)) error("ix not correct\n");
    if ((iy<0) || (iy>=list.celly)) error("iy not correct\n");
    if ((iz<0) || (iz>=list.cellz)) error("iz not correct\n");
        
    n = list.cell[ix][iy][iz].n;
    if (n+1>STACKDEPTH) error("too many particles in cell\n"); /* stackdepth= 30. why 30? */

    /*add the number of particle i to the stack of particlenumbers of cell[ix][iy][iz] and add +1 to the number of particles in cell[ix][iy][iz]*/
    list.cell[ix][iy][iz].stack[n] =i;   
    list.cell[ix][iy][iz].n++;

  }
}

void find_neighbor_using_cells(Slice *psl,int ipart){
  struct image_type {
    int  n;
    int  stack[MAX_IM_NEIGHBORS];
  } image_nn[NSHIFTS],*pil;

  List_item *nnl; /* *nnl=pointer naar de nearest neighbor list*/
  Cell *pjcell;
  vector shift,dr,ri;
  double dr2,dr2a, half_boxsizex=0.5*sys.boxl.x,half_boxsizey=0.5*sys.boxl.y,half_boxsizez=0.5*sys.boxl.z;
  int ix,iy,iz,dix,diy,diz,icx,icy,icz,i,j,k,m,jpart;
  int tmp=0;
  vector pbc,dra,drb;

  for(k=0;k<NSHIFTS;k++) image_nn[k].n=0;

  /*  dprint(ipart);*/
  ri=psl->pts[ipart].r;
  /*determine cell number of particle i*/
  ix = (int)((ri.x+half_boxsizex) * list.inv_cellsize.x ); /* shouldn't i be ipart. and therefore we can write (ri.x+half_boxsize) ?*/
  iy = (int)((ri.y+half_boxsizey) * list.inv_cellsize.y );
  iz = (int)((ri.z+half_boxsizez) * list.inv_cellsize.z );

  /*loop over neighboring cells*/
  for(dix=-GRIDOFFSET;dix<= GRIDOFFSET;dix++) { /* determine the shift of ?? dix from -3 to +3*/
    icx = ix + dix; 
    shift.x=0;
    if (icx <0) { icx+=list.cellx; shift.x =1;} else if (icx >=list.cellx) { icx-=list.cellx; shift.x =-1;}

    for(diy=-GRIDOFFSET;diy<= GRIDOFFSET;diy++) {
      icy = iy + diy;
      shift.y=0;
      if (icy <0) { icy+=list.celly; shift.y =1;} else if (icy >=list.celly) { icy-=list.celly; shift.y =-1;}
      
      for(diz=-GRIDOFFSET;diz<= GRIDOFFSET;diz++) {
        icz = iz + diz;  
        shift.z=0;
        if (icz <0) { icz+=list.cellz; shift.z =1;} else if (icz>=list.cellz) { icz-=list.cellz; shift.z =-1;}
        
        if ((icx<0) || (icx>=list.cellx)) error("icx not correct\n");
        if ((icy<0) || (icy>=list.celly)) error("icy not correct\n");
        if ((icz<0) || (icz>=list.cellz)) error("icz not correct\n");
        
        /*the neighborng cell */
        pjcell = &list.cell[icx][icy][icz];
        k = shift.z + 3*shift.y + 9*shift.x + 13;
        if ((k<0)|| (k>26)) error("k out of range\n");  


        /*loop over the particles that are in the cell. with max = number of particles in the neighboring cell (=pjcell->n)*/
        for(j =0; j<pjcell->n; j++) {
          jpart = pjcell->stack[j];     /* in stack there are particles numbers stored that are in that particular cell*/
          tmp++;
          /*      dprint(jpart);
          vprint(ri);
          vprint(psl->pts[jpart].r);*/

          if (jpart >ipart) {
            vector_minus(ri,psl->pts[jpart].r,dr); /* in currect code ri and dr are in real coordinates*/

            /*      vprint(dr);*/
            pbc.x = (int)(dr.x+half_boxsizex) ;  
            pbc.y = (int)(dr.y+half_boxsizey) ;  
            pbc.z = (int)(dr.z+half_boxsizez) ;  
 
            /*      vprint(pbc);*/
            dr.x += (int)(dr.x+half_boxsizex) ;  
            dr.y += (int)(dr.y+half_boxsizey) ;  
            dr.z += (int)(dr.z+half_boxsizez) ;
            dra =dr;  
            /*vector_times(dr,sys.boxl,dr);  dr is already in real units */

            dr2a = vector_inp(dr,dr);  /* was already in real units*/
            /*      printf("dr2 in pbc = %g\n",dr2a);*/
            vector_minus(ri,psl->pts[jpart].r,dr);  /* was scaled  units*/
            drb =dr;                        /*drb is the distance between ipart and jpart; was scaled units*/  

            /*it used to be scalar_times, but if you use boxl.x!=boxl.y then its vector times
            scalar_times(shift,sys.boxl.x,shift)    shift was from -1 to +1, so turn it into real coordinates*/
            vector_times(shift,sys.boxl,shift)
            vector_add(dr,shift,dr);   

            /*      vprint(shift);*/
            /* vector_times(dr,sys.boxl,dr); dr is already real*/
            dr2 = vector_inp(dr,dr);
            /*      printf("k= %d, dr2=%g\n",k,dr2);*/
            if( fabs(dr2-dr2a)>1e-2) {
              dprint(jpart);
              printf("%d %d %d\n",ix,iy,iz);
              printf("%d %d %d\n",dix,diy,diz);
              printf("%d %d %d\n",icx,icy,icz);
              vprint(ri);
              vprint(psl->pts[jpart].r);
              vprint(dra);
              // drb/=sys.boxl.x; drb was originally in scaled units
              // vprint(drb); 
              vprint(shift);
              vprint(pbc);
              printf("dr2 in pbc = %g\n",dr2a);
              printf("k= %d, dr2=%g\n",k,dr2);
              printf("pcb niet gelijk aan shift");
            }
            if (dr2< list.cutoff2) {
              pil = &image_nn[k];     /* pil=pointer naar the image_nn met (n en stack) van elke NSHIFT (=26) van particle i*/
              if (pil->n > MAX_IM_NEIGHBORS) error("too many neighbors per image\n");
              pil->stack[pil->n] = jpart; /* add particle jpart to the back of the stack*/
              pil->n++;
            }
          }
        }
      }
    }
  }

  /*  dprint(tmp);*/

  /* ga kijken hoe */
  for(k=0;k<NSHIFTS;k++){
    pil = &image_nn[k];
    if(pil->n >0) { /* als er een neighbor in een image zit dan:*/
      nnl = &list.nli[list.nimages];
      list.nimages++;
      if (list.nimages >MAX_IMAGES) error("too many images in list\n");
      nnl->ipart = ipart;
      nnl->image_id =k;
      nnl->first= &list.nlj[list.nneighbors];
      
      if (list.nneighbors + pil->n > MAX_NEIGHBORS)  error("too many neigbors in list\n");
      for (m=0;m<pil->n;m++) {
        list.nlj[list.nneighbors]= pil->stack[m];
        list.nneighbors++;
      }
    }
  }
}

void find_neighbor(Slice *psl,int ipart){
  struct image_type {
    int  n;                          /*number of neighbors of particle ipart in image k */
    int  stack[MAX_IM_NEIGHBORS];    /* what is stacking? is a list of particle numbers that are the neighbors and in image k */
  } image_nn[NSHIFTS],*pil;          /* pil points to the 26 structs that contain info of n and stack*/

  List_item *nnl;                   /*nearest neighbor List_item*/
  Cell *pjcell;
  vector shift,dr,ri;
  double dr2,dr2a,cutoff2, half_boxsizex=0.5*sys.boxl.x,half_boxsizey=0.5*sys.boxl.y,half_boxsizez=0.5*sys.boxl.z;
  int k,m,jpart;
  int tmp=0;
  vector pbc,dra,drb;

  for(k=0;k<NSHIFTS;k++) image_nn[k].n=0; /*reset/put to zero of 26 k images*/

  /*  dprint(ipart);*/
  ri=psl->pts[ipart].r;
  cutoff2=  list.cutoff2;  /* define list.cutoff2. it should be r_v al bit larger than r_c*/

  for(jpart=ipart+1;jpart<sys.npart;jpart++) {
    /* determine in which image particle j lies*/
    vector_minus(ri,psl->pts[jpart].r,dr);
    shift=nulvec;                      
    if (dr.x < -half_boxsizex) shift.x=+1;        /*code uses nonreduced units for xyz. therefore, half_boxsize*/
    if (dr.x >  half_boxsizex) shift.x=-1;
    if (dr.y < -half_boxsizey) shift.y=+1;
    if (dr.y >  half_boxsizey) shift.y=-1;
    if (dr.z < -half_boxsizez) shift.z=+1;
    if (dr.z >  half_boxsizez) shift.z=-1;
    k = shift.z + 3*shift.y + 9*shift.x + 13;
    if ((k<0)|| (k>26)) error("k out of range\n");      

    // scalar_times(shift,sys.boxl.x,shift); turned into vectortimes
    vector_times(shift,sys.boxl,shift)
    vector_add(dr,shift,dr); 
    
    dr2 = vector_inp(dr,dr);
    if (dr2< cutoff2) {
      pil = &image_nn[k];
      if (pil->n > MAX_IM_NEIGHBORS) error("too many neighbors per image\n");
      pil->stack[pil->n] = jpart;
      pil->n++;
    }
  }
        
  /* for each image that contains a neighbor; fill in the (global) nearest neighbor list (nnl) */
  for(k=0;k<NSHIFTS;k++){
    pil = &image_nn[k];
    if(pil->n >0) { /* if there is a particle in image k*/
      nnl = &list.nli[list.nimages];
      list.nimages++;
      if (list.nimages >MAX_IMAGES) error("too many images in list\n");
      nnl->ipart = ipart;
      nnl->image_id =k;
      nnl->first= &list.nlj[list.nneighbors];
      
      if (list.nneighbors + pil->n > MAX_NEIGHBORS)  error("too many neigbors in list\n");
      for (m=0;m<pil->n;m++) {
        list.nlj[list.nneighbors]= pil->stack[m];
        list.nneighbors++;   /* add one to the nneighbors*/
      }
    }
  }
}

void print_cell_list(){
  int ix,iy,iz,i,count=0;
  
  for(ix=0;ix< list.cellx;ix++) 
    for(iy=0;iy< list.celly;iy++) 
      for(iz=0;iz< list.cellz;iz++) {
        printf("cell (%d,%d,%d) bevat %d deeltjes:",ix,iy,iz,list.cell[ix][iy][iz].n);
        for(i=0;i<list.cell[ix][iy][iz].n;i++) printf(" %d",list.cell[ix][iy][iz].stack[i]);
        printf("\n");
        count+= list.cell[ix][iy][iz].n;
      }
  dprint(count);
  if (count!=sys.npart) {
    printf("celllist corrupted %d\n",count);
    printf("exiting");
  }
}

int check_nnlist(Slice *psl){
  vector dr;
  double dr2;
  int i,j,check_count=0;

  for(i=0;i<sys.npart-1;i++) {
    for(j=i+1;j<sys.npart;j++) {
      vector_minus(psl->pts[i].r,psl->pts[j].r,dr); 
      pbc(dr, sys.boxl);
     
      dr2 = vector_inp(dr,dr);  
      if (dr2<list.cutoff2) check_count++;
    }
  }

  return check_count;    
}

void update_nnlist(Slice *psl){
  int i;
  printf("updated nnlist\n");
  put_parts_in_box(psl);
  /*  create_cell_list(psl);*/
  /*  print_cell_list();*/
  list.nimages = list.nneighbors =0;
 
  for (i=0;i<sys.npart;i++) { 
    find_neighbor(psl,i);
    /*find_neighbor_using_cells(psl,i);*/
  }
  list.nli[list.nimages].first= &list.nlj[list.nneighbors];
  // print_nnlist();

  // dprint(list.nneighbors);
 // sys.block_stats.nnlength = list.nneighbors;
  // sys.block_stats.nnup++;
  if (check_nnlist(psl) != list.nneighbors) printf("neighborlist corrupted\n");
  return;
}

void saveupdate_nnlist(Slice *psl){
  int i, nimages, nneighbors;
  vector r;
  Pts *dummy_pts[NPART];

  /*maybe this MC cluster cycle, because I use the nnl made with the copyslice in the MC cluster cycle*/
  Slice *copyslice = malloc(sizeof(Slice));
  memcpy(copyslice, psl, sizeof(Slice));


  

  // put_parts_in_box(psl);
  for(i=0;i<sys.npart;i++) {
      r = psl->pts[i].r; 
      pbc(r,sys.boxl); 
      copyslice->pts[i].r = r;
      copyslice->pts[i].dr = nulvec;  
    }
  
  // list.nimages = list.nneighbors =0;
  nimages = nneighbors =0;
 
  for (i=0;i<sys.npart;i++) { 
    find_neighbor(copyslice,i);
    /*find_neighbor_using_cells(psl,i);*/
  }
  list.nli[list.nimages].first= &list.nlj[list.nneighbors];
  // print_nnlist();

  // dprint(list.nneighbors);
 // sys.block_stats.nnlength = list.nneighbors;
  // sys.block_stats.nnup++;
  if (check_nnlist(psl) != list.nneighbors) printf("neighborlist corrupted\n");
  return;
}


void print_nnlist(){
  int i,k,ipart,*plj,count=0;
  
  printf("i \t shift\t j\n");
  for(i=0;i<list.nimages;i++) {
    ipart = list.nli[i].ipart;
    k = list.nli[i].image_id;
    printf("%d \t %d \t",ipart,k);
    for(plj=list.nli[i].first;plj<list.nli[i+1].first;plj++) {
      printf("%d ",*plj);
      count++;
    }
    printf("\n");
  }
  dprint(count);
}












