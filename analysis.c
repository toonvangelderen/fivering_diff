#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"

// I might want to move this function to somewhere else, maybe in analysis.c?
void linked_structure_icluster(Slice *psl, int icluster, int ipart_ref){
    // if a structure crossed tha pbc, rotating it is straightforward. 
    // with this fucntion you first check if structure crosses pbc ; are there particles 1 isgma from the border? in the 2D plane
    int i,b,current,ipart,boundary_crossing=0,nbonds;
    int particlesleft=cluster.clustersize[icluster];
    vector loc_prev;
    int clusteri_size=cluster.clustersize[icluster];
    vector dr;
    Pts newpos[NPART],p_ref;
    // printf(" copy psl to newpos\n");
    memcpy(newpos, psl->pts, sizeof(psl->pts));

    // first try to center around p_ref, and check if there are particles at the edge of the box
    // ipart_ref = cluster.pic[icluster].stack[ref_random];
    int edge=0;
    // dprint(ipart_ref);

    p_ref = psl->pts[ipart_ref];
    // vprint(p_ref.r);
    //move the structure such that p_ref is in the center, while checking if the new structure is near the boundary after pbc
    for(i=0; i<clusteri_size; i++){
        ipart= cluster.pic[icluster].stack[i];
        
        pbc(newpos[ipart].r,sys.boxl);
        
        if ((0.5*sys.boxl.x-fabs(newpos[ipart].r.x)<1. )| (0.5*sys.boxl.y-fabs(newpos[ipart].r.y)<1.)| (0.5*sys.boxl.z-fabs(newpos[ipart].r.z)<1.) ){
            edge=1;
            // vprint(newpos[ipart].r);
            // vector_add(newpos[ipart].r,p_ref.r,newpos[ipart].r); // so newpos[ipart].r is the chain with p_ref in (0,0,0)
            // printf(" part of structure (particle %d ) near the boundary\n",ipart);
            break; // the structure is near the boundary and may cross it
        }
        

    }
    if (edge==0){
        // printf(" no particles at edge, newpos is correct as it is\n");
        //no particles at edge, newpos is correct as it is
        memcpy( psl->pts, newpos, sizeof(psl->pts));
        return;
    }


    if(edge==1){
        // printf(" edge =1 \n");
        Picl pici=cluster.pic[icluster];
        IntArray particlesleft_list,branchlist; 
        int b, particlesleft=clusteri_size;
        int branch=0,q=0;
        int current,next;
        vector dr;
        

        //initiate the branchlist (empty) and particlesleft_list (filled)
        initIntArray(&branchlist,  (size_t) clusteri_size);
        initIntArray(&particlesleft_list,  (size_t) clusteri_size);
        for(i=0; i<clusteri_size; i++){         
            insertIntArray(&particlesleft_list, pici.stack[i]);
        }
        
        // printf("whole array of particles_left\n");
        // printIntArray(&particlesleft_list);
        // you will need to loop over the bonds, because the boundaries are crossed by the structure
        //start from ipart_ref; 
        current=ipart_ref; 
        newpos[current].r=psl->pts[current].r; // the starting point/particle, check if it is a " branchpoint", you need to go backwards too
        if(psl->pts[current].nbonds>1){ //printf(" particle %d added to branchlist\n",current);
        insertIntArray(&branchlist, current);}
        removeElementXIntArray( &particlesleft_list ,  current); //remove current from particlesleft
        
        // printf("removed current=%d from particles_left\n",current);
        // printIntArray(&particlesleft_list);
        i=0;
        do{
            if (psl->pts[current].nbonds>i){
                next=psl->pts[current].bonds[i];
                if(checkElementXIntArray(&particlesleft_list, next)){
                    // printf(" particle %d (current) and %d (next) are compared\n",current,next);
                    dr=bond_vector(psl, current, next,  1); //check if order of( current, next) is correct (direction of vector)
                    //position next
                    vector_minus(newpos[current].r,dr,newpos[next].r);

                    //add branchpoint if the current particle has more than 2 bonds and is still in the list (to add it only once)
                    if (psl->pts[current].nbonds>2){   
                        // printf(" particle %d added to branchlist\n",current);
                        insertIntArray(&branchlist, current); //save the particle number of the branch
                    }
    
                    // you have given next a new postion, so remove it particle from &particlesleft_list
                    if (checkElementXIntArray(&particlesleft_list,next)){
                        // printf("removed current=%d from particles_left\n",next);
                        removeElementXIntArray( &particlesleft_list ,  next); 
                    }
                    
                    // go to next particle in the cluster, 
                    current=next; 
                    i=0;
                }
                else{
                    i++;
                }
            }
            else{
                //go further the next bond of current branchpoint ; 
                //because current keeps on changing if you walk over the bonds
                if (i==psl->pts[current].nbonds){ // do this if all bonds of current are not in particlesleft_list anymore}
                    if (checkElementXIntArray(&branchlist,current)){
                        removeElementXIntArray( &branchlist ,  current); 
                    }
                    //or if there are no other bonds: go further with new branch points
                    if (branchlist.used > 0){  //if there are branchpoints left
                        current=branchlist.ints[0];//this contains the branchpoints which might have multiple bonds>2
                        // start with the first bond of this nranchpoint partilce, and keep let i run upto a particle that is still in the list
                        i=0;
                        //go to i==nbonds or i is the next neighbor that is still in the particlesleft_list
                        while ( (i<psl->pts[current].nbonds) && (checkElementXIntArray(&particlesleft_list,psl->pts[current].bonds[i])==0)){
                            i++;
                        } 
                        
                    } 
                }
            }
        }while(((int)particlesleft_list.used>0)  || ((int)branchlist.used>0))  ;

        // free the memory!!
        freeIntArray(&branchlist);
        freeIntArray(&particlesleft_list);
    }

    //copy new structure to psl;
    memcpy( psl->pts, newpos, sizeof(psl->pts));
    
    return;
}

void nematic_order_parameter(Slice *psl){
    /*https://journals.aps.org/pra/pdf/10.1103/PhysRevA.31.1776
    https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.98.095701
    2 dimensions
    Si = 1/N_b \sum_j=0^Nb 2|| u_i dot u_j ||^2 -1
    Si local nematic order paramtere of bond i
    N_b number of neighbor particle of bond i within 2 sigma
    u_i is unit vector of bond i  
    */

    /*  you loop over the particle and open list of its bond,
    loop over the bond if the other particle psi->bonds[k]>ipart  //uniquely looop over bonds
    loop over all other particles, jpart=0. and look at the distance cut-off p_cutoff of the particles.
    then look  up the l bond those jpart particles make, check again the cutoff bond_cutoff2. and calculate Si contribution
    track the Si contributions, and if all viscinity bonds are found. Save Si/Nb to pts struct of particle ipart  */

    int ipart, jpart,n,k,l,nbjl,p,k_bo=0,Nb=0;
    Pts *psi,*psj;
    vector rij,uik_bond,ujl_bond,uik_center,ujl_center, rikjl_centers,ri,rj, rk, rl,rjl,uik_bond05,ujl_bond05;
    double r2, r,Si=0., vinp; 
    static double bond_cutoff,bond_cutoff2,p_cutoff,p_cutoff2;
    static int initiate=1;
    
    printf("starting nematic_order_parameter analysis\n");
    if (initiate){
        double max_diameter=0;
        for(p=0;p<sys.nparticle_types; p++){
            if (sys.particletype[p].diameter>max_diameter){
                max_diameter=sys.particletype[p].diameter;
            }
        }
        initiate=0;
        p_cutoff=3.*(max_diameter + sys.rcutoff); //p_cutoff stands for particle cutoff
        p_cutoff2=p_cutoff*p_cutoff;
        bond_cutoff=2.*max_diameter;  
        bond_cutoff2=bond_cutoff*bond_cutoff;

        gprint(p_cutoff);
        gprint(bond_cutoff);

    }
    

    // >>>>>actually loop over clusters, in this way you avoid including bond that are in the same cluster or chain.
    int cid_i,cid_j,pic_i,pic_j;

    for( ipart=0; ipart<sys.npart; ipart++ ){
        psi=&psl->pts[ipart];
        psi->k_bo=0;
    }
    // for( ipart=0; ipart<sys.npart; ipart++ )  {
        //>>>>>loop over the particles in cluster id cid_i (=cluster id i)
    for(cid_i=0;cid_i<psl->nclusters;cid_i++){
        // printf("new cluster %d\n", cid_i); pic_i ((=particles in cluster i))
        for (pic_i=0;pic_i<cluster.clustersize[cid_i]; pic_i++) {
            ipart = cluster.pic[cid_i].stack[pic_i]; 

            psi=&psl->pts[ipart];
            ri=psi->r;
            pbc(ri,sys.boxl);

            /*look up the bonds particle ipart makes*/
            for(k=0;k<psi->nbonds;k++){
                //psi->bonds[k]; this is the particle number of the particle that makes a bond with ipart
                if (psi->bonds[k]<ipart){
                    continue; //avoid double counting; uniquely loop over the bonds 
                }
                    // dprint(psi->bonds[k]);

                rk= psl->pts[psi->bonds[k]].r;
                pbc(rk,sys.boxl);

                vector_minus(ri,rk ,uik_bond); //ri-rj bond vector
                pbc(uik_bond,sys.boxl); // 
                scalar_times(uik_bond,-0.5,uik_bond05);

                uik_bond.z=0; //2D analysis ??
                normvec(uik_bond,uik_bond);

                vector_add(ri,uik_bond05,uik_center); //bond center
                pbc(uik_center,sys.boxl); // again pbc , because the particle i and k could be at the edges of the box
                psi->bond_op[psi->k_bo].r = uik_center; //safe the 3D location of the bond center here, before z-component is set to zero

                uik_center.z=0;


                //<<<<here loop over all  jparticles > ipart+1, because you need to find unique bond combinations; not sure about this yet
                // >>>>loop over the other clusters and their stack
                // for( jpart=0; jpart<sys.npart; jpart++ )  {
                for(cid_j=0;cid_j<psl->nclusters;cid_j++){
                    if (cid_j==cid_i){
                        continue; //skip of course the cluster that you are   already at
                    }
                    //>>>>>loop over the particles in cluster id cid_j
                    for (pic_j=0;pic_j<cluster.clustersize[cid_j]; pic_j++){
                        jpart = cluster.pic[cid_j].stack[pic_j]; 

                        // also not the next bond
                        if ((jpart==psi->bonds[k] )|| (jpart==ipart)  ){
                            continue; //skip the particles that makes bond nj
                        }

                        psj=&psl->pts[jpart];
                        rj=psj->r;
                        pbc(rj,sys.boxl);
                        
                        vector_minus(ri,rj,rij);
                        r2 = vector_inp(rij,rij);
                            
                        // if particle jpart lies in the viscinity, check if it has bonds
                        if(r2<p_cutoff2) { 
                            for(l=0;l<psj->nbonds;l++){
                                // if (psj->bonds[l]<jpart){
                                //     continue;
                                // }
                                // nbjl=psj->bonds[l]; 
                                rl=psl->pts[psj->bonds[l]].r;
                                pbc(rl,sys.boxl);

                                vector_minus(rj, rl,ujl_bond ); //ri-rj
                                pbc(ujl_bond,sys.boxl);
                                scalar_times(ujl_bond,-0.5,ujl_bond05);

                                ujl_bond.z=0; //2D analysis 
                                normvec(ujl_bond,ujl_bond);

                                vector_add(rj,ujl_bond05,ujl_center); //bond center
                                pbc(ujl_center,sys.boxl); // again pbc , because the particle i and k could be at the edges of the box
                                // scalar_divide(ujl_center,2.,ujl_center);
                                // if (rj.x>0.45*sys.boxl.x){
                                //     vprint(rj);
                                //     vprint(rl);
                                //     vprint(ujl_bond);
                                //     vprint(ujl_center);
                                //     printf("\n");
                                // }

                                // if (ujl_center.z>2.){
                                //     printf("      ujl_center.z>2.\n");
                                //     vprint(rj);
                                //     vprint(rl);
                                //     vprint(ujl_bond);
                                //     vprint(ujl_center);
                                //     printf("\n");
                                // }
                                ujl_center.z=0;


                                // are the uik_center and ujl_center in eachothers vicinity?
                                vector_minus(uik_center,ujl_center,rikjl_centers);
                                pbc(rikjl_centers,sys.boxl); // again pbc , because the particle i and k could be at the edges of the box

                                r2=vector_inp(rikjl_centers,rikjl_centers);

                                if(r2<bond_cutoff2) {
                                    //add it to the list of bonds; both rik and rjl 
                                    // add it to the pts of psi 
                                    //Si = 1/N_b \sum_j=0^Nb 2|| u_i dot u_j ||^2 -1
 
                                    vinp=vector_inp(uik_bond,ujl_bond) ;
                                    // gprint(2.*(vector_inp(uik_bond,ujl_bond) )*(vector_inp(uik_bond,ujl_bond) ) -1.);
                                    Si+= (2.*(vinp*vinp) -1.);
                                    Nb+=1;
                                    k_bo=1; //k_bo is?? --> number of bond order properties
                                    // gprint(Si);
                                    // gprint(Nb);
                                    // dprint(k_bo);
                                }
                            }
                        }
                    }
                }
                // safe the order parameter to the bond property struct of particle i;
                if (Nb!=0){
                    psi->bond_op[psi->k_bo].nematic_OP=(double)Si/(double)Nb;
                }
                else{
                    psi->bond_op[psi->k_bo].nematic_OP=0;
                }
                psi->k_bo+=1;
                
                Si=0.;
                Nb=0;
                k_bo=0;
            }
        }
    }

}

double Sfrac(Slice *psl){
    /*the Sfrac (that simon plots in his polymer paper fig 3f)
    N^b_{S>0.4} / N^b 
    in words: the number of bonds with S>0.4 divided by total number of bonds, i.e. fraction of bonds that has S>0.4
    where S is the nematic order parameter in void nematic_order_parameter(Slice *psl)*/
    int  ipart,Nbonds_tot=0,k;
    double Nb_S=0,Sfrac;
    Pts *psi;

    /*loop over all partices, and check whether the bond order parameter is above 0.4*/
    for(ipart=0;ipart<sys.npart;ipart++){
        psi=&psl->pts[ipart];
        Nbonds_tot+=psi->k_bo;
        //loop over the (unique) bonds that have a bond order property
        for(k=0;k<psi->k_bo;k++){
            if (psi->bond_op[k].nematic_OP>0.4){
                Nb_S+=1.;
            }
        }
    }
    if (Nbonds_tot>0){
        Sfrac=(double)Nb_S/(double)Nbonds_tot;
    }
    else{
        Sfrac=0;
    }
    return Sfrac;
}

void timecorrelation_function_P(Slice *psl){
    /*from Ilie2015 eqn 40
    P(t)= 1.5<(vec_u(t)vec_u(0))^2> -0.5
    vec_u is a body fixed vector, taken as the patch vector  and t is time [s]

    for simulation of a single particle, give it 1 patch vector. 
    perform it 100 times and calculate  P(t)= 1.5<(vec_u(t)vec_u(0))^2> -0.5 */

    double P,vecinp0,vecinpt;
    vector vec_ut=psl->pts[0].patchvector[0];
    
    static int initialize =1; 
    static vector vec_u0,vec_prev;

    if (initialize){
        vec_u0=psl->pts[0].patchvector[0];
        vec_prev=psl->pts[0].patchvector[0];
        initialize=0;
    }

    vecinp0=vector_inp(vec_u0,vec_ut);
    vecinpt=vector_inp(vec_prev,vec_ut);
    // gprint(vecinp0);
    // gprint(vecinpt);
    vec_prev=vec_ut;
    /*print to file*/
    FILE *file;

    file = fopen("vector_inp_u0ut.dat","a");
    printf("printing to \"vector_inp_u0ut.dat\"\n ");

    if (file == NULL){
            printf("Error with opening \"vector_inp_u0ut.dat\" \n");
    }
    
    fprintf(file, "%.6lf %.12lf %.12lf  %.12lf %.12lf\n", sys.c_time , vecinpt, vec_ut.x, vec_ut.y,vec_ut.z);
    
    fclose(file);

    return;

}

void MSD(Slice *psl){

    /*root mean squared displacement*/

    static vector old_pos;
    static int initialize=1;
    static double msd_track;
    double vecinp;
    vector dr;

    if (initialize){
        old_pos=psl->pts[0].r;
        initialize=0;
        msd_track=0;
    } 
    vector_minus(psl->pts[0].r,old_pos,dr);
    pbc(dr,sys.boxl)
    dr.z=0; //put it to zero, as experimentalist measure MSD in xy-plane
    
    // vprint(dr);
    running_statistics(&msd,vector_inp(dr,dr));
    msd_track+=(vector_inp(dr,dr));
    old_pos=psl->pts[0].r;

    /*print to file*/
    FILE *file;

    file = fopen("msd.dat","a");
    printf("printing to \"msd.dat\"\n ");

    if (file == NULL){
            printf("Error with opening \"msd.dat\" \n");
    }
    
    fprintf(file, "%.6lf %lf %lf\n", sys.c_time , msd_track,msd.mean);
    
    fclose(file);

    return;
     
}

void xy_plane_patchangle_calc(Slice *psl){
    /*this is to calculate the angle of the patch wrt the xy-plane
    to see if it is uniformly distirbuted
    make a histogram*/
    int i,gamma_box,ipart,x;
    double gamma,rad=PI/180., d_gamma=1./ANGLEBINS;
    vector p1;
    Pts *psi;

    for(i=0;i<psl->nclusters;i++){
         //loop only over monomers and select size==1
        if( cluster.clustersize[i]==1){
            ipart=cluster.pic[i].stack[0];
            psi=&psl->pts[ipart];
            p1=psi->patchvector[0];
            gamma =  fabs(psi->patchvector[0].z);
            if (gamma>1){
                gamma=1;
            }
            gamma_box = floor(gamma/d_gamma);
            
            if ((gamma_box<0 )|| gamma_box>ANGLEBINS ){
                dprint(ipart);
                gprint(gamma);
                dprint(gamma_box);
                vprint(p1);
                gprint(acos(fabs(p1.z)));
                error("gamma_box<0 or gamma_box>ANGLEBINS  in xyplane patchangle");
            }
            //make a histogram of the angles
            xy_plane_patchangle[gamma_box] ++;
        }
    }

    //print to file 
    FILE *file;

    file = fopen("xy_plane_patchangle.dat","w");
    printf("printing to \"xy_plane_patchangle.dat\"\n ");

    if (file == NULL){
            printf("Error with opening \"xy_plane_patchangle.dat\" \n");
    }
    for (x=0; x<ANGLEBINS;x++){
        fprintf(file, "%lf %d\n", (double)(x+0.5)*d_gamma , xy_plane_patchangle[x]);
    }
    
    fclose(file);

    
    return;
    
} 


int bond_check(Slice *psl, int i, int j){
    /* check whether there exists a bond between particle i and j. Based on the treshold value sys.bond_cutoffE
    returns 1 for yes, returns 0 for no*/
    Pts *psi,*psj; 
    int n;  

    psi=&psl->pts[i];
    psj=&psl->pts[j];    

    for (n=0;n<psi->nbonds;n++){
        if (psi->bonds[n]==j){
            return 1;
        }
    }
    for (n=0;n<psj->nbonds;n++){
        if (psj->bonds[n]==i){
            return 1;
        }
    }

    return 0;
}


double bond_distance(Slice *psl, int i, int j, int pbc){
    /* calculates and returns the distance between two particles. specify with pbc if pbc should be on or off*/

    double dr2, r;
    vector dr;

    if(pbc==0 || pbc ==1 ){
        dr = bond_vector(psl, i, j, pbc);
        dr2 =vector_inp(dr,dr);
        r=sqrt(dr2) - sys.particletype[psl->pts[i].ptype].radius-sys.particletype[psl->pts[j].ptype].radius;
    }
    else{
        dprint(pbc);
        error("in bond_distance specify only pbc on (= 1) or off (= 0), no other values allowed");
    }
    return r;
}

void s_distribution_analysis(Slice *psl){
    // loop over the particles and see to which other particls it is bound
    // check the s distance with bond_distance(psl,  i,  j,  1)
    int i,n,j,s_bin,x;
    double s, rangemax=0.02,ds=rangemax/MAXBIN;
    Pts *psi;
    FILE *file;

    for (i=0;i<sys.npart; i++){
        psi = &psl->pts[i];
        for (j=i+1;j<sys.npart; j++){
            s=bond_distance(psl,  i,  j,  1);
            if (s<=rangemax){
                s_bin = floor(s/ds);
                // printf(" for s=%lf it is in bin=%d where ds=%lf defines the bins and s/ds = %lf\n", s,s_bin, ds,s/ds);
                s_histogram[s_bin]+=1;
            }
        }
    }

    file = fopen("s_histogram.dat","w");
    printf("printing to \"s_histogram.dat\"\n ");

    if (file == NULL){
            printf("Error with opening \"s_histogram.dat\" \n");
    }
    for (x=0; x<MAXBIN;x++){
        fprintf(file, "%lf %d\n", (x+0.5)*ds , s_histogram[x]);
    }
    
    fclose(file);
    
    return;
}

void bond_length_average(Slice *psl){
    Pts *psi;
    int ipart, n, jpart;
    double s;
    
    for(ipart=0; ipart<sys.npart; ipart++){
        psi = &psl->pts[ipart];
        for (n=0;n<psi->nbonds;n++){ // you will only walk over the bonds
            jpart = psi->bonds[n];
            if (jpart<ipart){continue;} //skip if particle number jpart < ipart, to avoid double counting
            s=psi->bond_distance[n];
            running_statistics(&bond_length,s);
        }
    }

    return ;
}

vector bond_vector(Slice *psl, int i, int j, int pbc){
    /* calculates and returns the vector rij(=ri-rj) between two particles. specify with pbc if pbc should be on(pbc=1) or off(pbc=0)*/
    vector dr;

    if(pbc==0 || pbc ==1 ){
        vector_minus(psl->pts[i].r,psl->pts[j].r,dr); /* dr of particle j and k */
        if(pbc==1){  pbc(dr, sys.boxl); } 
    }
    else{
        dprint(pbc);
        error("in bond_distance specify only pbc on (= 1) or off (= 0), no other values allowed");
    }
    return dr;  
}




/*void calc_rdf_localCN(Slice *psl) { 
  //, Average *rdf_average ,Average *globalCN_average
    //function to calculate the local coordination number (localCN) and radial distribution function (rdf) of a slice
   
    //loop over all particles, calculate the distances between particle i and the rest. 
    //ocalCN (Coordination Number) =  count number of particles within distance 1.5 sigma per particle
    //globalCN = count number of particles within distance 1.5 sigma on average 
    //RDF(r) = <N(r)>/4PI*r^2*dr*rho 
     
    int ipart, jpart, rdf_i, N;
    double l, r, l2, boxlength;
    MIN(sys.boxl.x,sys.boxl.y,boxlength);
    double rdf_dr=boxlength/RDFBINS, rho=sys.npart/(sys.boxl.x*sys.boxl.y*sys.boxl.z), correction, globalCN=0.0;
    vector lsq, ripart ;
    FILE *globalCNfile;

    double rdf_slice[RDFBINS]={0};

    for (ipart=0; ipart<sys.npart ; ipart++){
        N=0; //set to zero to calculate the localCN per particle
        ripart = psl->pts[ipart].r;
        
        for(jpart=0;jpart<sys.npart; jpart++){
            if(ipart!=jpart){
                // l = bond_distance(psl, ipart, jpart, 1)
                vector_minus(ripart, psl->pts[jpart].r, lsq);         
                pbc(lsq, sys.boxl);                                 //correct for PBC
                l2 = vector_inp(lsq,lsq);                           //length = qsrt(x^2 + y^2 + z^2)
                l= sqrt(l2);
                
                // localCN 
                if(l<1.5*1.017){
                    N+=1;
                }

                // RDF; start by determinening in which bin of the g(r) the distance of particle i and j falls(done by drf_i = floor(l/rdf_dr));
                //then, we count +1. later we will normalize the g(r) with 4pi*r^2*dr*rho and sys.npart (for the mean value of <dn>) 
                rdf_i = floor(l/rdf_dr);
                rdf_slice[rdf_i]+=1;
            }
        }
        psl->localCN[ipart]=N;           // stores the localCN per particle specifically for the slice, the localCB is not yet implemented in e.g. graphics 
        globalCN +=N;
    }

    globalCNfile = fopen("globalCN.dat","a");
    printf("printing to \"globalCN.dat\" :    globalCN  = %1.2f \n",globalCN_average.sum/globalCN_average.n);

    if (globalCNfile == NULL){
            printf("Error with opening \"globalCN.dat\" \n");
    }

    fprintf(globalCNfile, "%d  %.5f\n", sys.ncycle2, globalCN_average.sum/globalCN_average.n);
    fclose(globalCNfile);
    

    //Calculate the average globalCN of the slice; and add it to globalCN_average/
        update_average(globalCN_average, globalCN/sys.npart);

    // here id the normalization of  g(r) done where 
    //    r= rdf_i*rdf_dr; dr = rdf_dr; rho=sys.npart/boxlx^3 -> 4*PI*rdf_i^2*rdf_dr^3*rho
     //   we only need average g(r) of the complete simulation; not of one slice. Therefore, the g(r) is added to rdf_average  
    for(rdf_i=1; rdf_i<RDFBINS/2; rdf_i++){
        r = rdf_i*rdf_dr;
        correction = 4*PI* r*r * rdf_dr* rho*sys.npart; 
        update_average(rdf_average[rdf_i], rdf_slice[rdf_i]/correction);
    }

    return;
}*/

void clustersize_freq_update(Slice *psl){
    /*this function saves the frequency of occurence of the clusterlengths  */

    int i,l;
    int   freq_clusterlength[NPART]={0};
    

    for(i=0;i<psl->nclusters;i++){
        freq_clusterlength[cluster.clustersize[i]-1]++;
    } 

    for(l=0;l<sys.npart;l++){
        /*loops over the lengths which has maximumvalue of sys.npart*/
        // update_average(chain.length_histogram[i], freq_clusterlength[i]);
        
        running_statistics(&cluster.size_histogram.length[l], freq_clusterlength[l]);
        running_statistics(&cluster.size_distribution.length[l], (double)freq_clusterlength[l]/psl->nclusters);

    } 

    return; 
}

void chain_linkedlist(Slice *psl){
    /*To identify the chains, put the particle in a sequence such that you can walk over the chain and calculate e.g. end-to-end distance, FT analysis*/
    /* it works with head, tail and neighbor. You want the find the next in the chain, next to the tail or head, based on the distance but exclude the neighbor you already found*/

    int id,ipart, nclusters=psl->nclusters, n,k, kpart, particlesleft, neighbor, current_tail, added, next_ptsnr,p, current_head;
    Picl pici;

    // printf("Performing chain identification; order the cluster based on position \n");

    // printf("loop over the picl's and make the right sequence of the particles in the doubly linked ChainNode struct\n");
    //loop over the picl's and make the right sequence of the particles in the doubly linked ChainNode struct
    for(id=0; id<nclusters; id++){
        pici=cluster.pic[id];
        chain.head[id]=NULL;
        InsertAtChainTail(pici.stack[0], id);  /* start the chain at the first particle; add accordingly at tail or head;*/ 

        //dummies
        current_tail=pici.stack[0];   
        neighbor=-1;
        particlesleft=cluster.clustersize[id]-1;
    
        do{
            for(p=1;p<cluster.clustersize[id];p++){
                added=0;
                kpart=pici.stack[p];
                if(kpart!=neighbor && kpart!=current_tail){
                    if( bond_check(psl,current_tail,kpart)==1){
                        InsertAtChainTail(kpart,id);
                        neighbor=current_tail;
                        current_tail=kpart;
                        particlesleft-=1;
                        // printf("particle %d was inserted at tail;   neighbor %d; current tail %d; particles left %d\n", kpart, neighbor, current_tail, particlesleft );
                        added=1;
                        break;
                    }
                }
            }   
        }while((added==1 ) && (particlesleft>0));
        /*then add it to the head*/
        if(particlesleft>0){
            // printf("  THERE ARE %d PARTICLES LEFT AFTER TAIL INSERTION\n", particlesleft );
            current_head=chain.head[id]->pts_nr;
            neighbor=chain.head[id]->next->pts_nr;

            do{
                for(p=0;p<cluster.clustersize[id];p++){
                    added=0;
                    kpart=pici.stack[p];
                    if(kpart!=neighbor && kpart!=current_head){
                        // printf("kpart is %d, neighbor is %d, current_head is %d\n", kpart, neighbor, current_head);
                        if( bond_check(psl,current_head,kpart)==1){
                            InsertAtChainHead(kpart,id);
                            neighbor=current_head;
                            current_head=kpart;
                            particlesleft-=1;
                            // printf("particle %d was inserted at head;   neighbor %d; current head %d; particles left %d\n", kpart, neighbor, current_head, particlesleft );
                            added=1;
                            break;
                        }
                    }
                }
            }while((added==1) && (particlesleft>0));
        }
    }    
    // printf("All clusters have been sequenced.\n");
    return;
}

void check_maxbonds(Slice *psl){
    /*checks if any particle exeeds the maximum number of bonds (more than one per site)*/
    int ipart;
    Pts *psi;

    for(ipart=0;ipart<sys.npart;ipart++){
        psi=&psl->pts[ipart];
        if(psi->nbonds>sys.particletype[psi->ptype].nsites){
            dprint(ipart);
            dprint(psi->nbonds);
            dprint(sys.particletype[psi->ptype].nsites);
            error("too many bonds for the particle");
        }
    }
    return;              
}

int cluster_analysis(Slice *psl){
    /* color all particles white, loop over all other particles (j) color them grey; find bonding particles (with dr and SiSj)*/
    /*particles are in a cluster if they have a bond energy <sys.bond_cutoffE*/
    Pts *psj;
    vector dr;
    int label[sys.npart],nwhite,notallgreysgone,ngrey,ncluster,nblack,i,j,k, WHITE=1, BLACK=2, GREY=3, id;
    double r,dr2,s,r_radii,  potential_energy_P_d=0;
    double s_patch;

    //  count=0;
    notallgreysgone=0;
    ngrey=nblack=ncluster = 0;
    for(i=0;i<sys.npart;i++) {
        label[i] = WHITE;
    }

    nwhite = sys.npart;
    for(i=0;i<sys.npart;i++) {
        if (label[i] == WHITE) {
            ncluster++;
            nwhite--;
            ngrey =1;
            label[i]= GREY;
            do {
                notallgreysgone = 0;
                for (j=0;j<sys.npart;j++) {
                    if (label[j]==GREY) {
                        notallgreysgone = 1;
                        ngrey--;
                        nblack++;
                        label[j]=BLACK;
                        psj = &psl->pts[j];
                        psj->cluster_id=ncluster; /* give particle j the cluster_id of particle i*/
                        /*            printf("deeltje %d wordt zwart,nwhite= %d, ngrey= %d nblack= %d, ntot=%d\n",j,nwhite,ngrey,nblack,nwhite+nblack+ngrey);*/
                        for (k=i;k<sys.npart;k++){
                            if ((k!=j) && (label[k] == WHITE)) {
                                vector_minus(psj->r,psl->pts[k].r,dr); /* dr of particle j and k */
                                pbc(dr, sys.boxl); 
                                dr2 =vector_inp(dr,dr);

                                r=sqrt(dr2);
                                r_radii=sys.particletype[psl->pts[k].ptype].radius+sys.particletype[psj->ptype].radius;
                                s=r-r_radii;
                                /* implement here the check if the distance is within the bond distance and the patches form a bond (without NNL)*/
                                if(s<sys.rcutoff && s>0.0) { 
                                    /*this function returns the value of the potential energy to evaluate the existence of a bond*/
                                    potential_energy_P_d= potential_attractive_bond_energy(psl, j, k);

                                    if(potential_energy_P_d<sys.bond_cutoffE){
                                        label[k]=GREY;
                                        nwhite--;
                                        ngrey++;
                                    }
                                }
                            }
                        }
                    }
                }
                if (ngrey==0)   notallgreysgone=0;
            } while (notallgreysgone);
            if (nwhite==0) i= sys.npart;
        }
    }


    return ncluster;
}

void clustersize_identification(Slice *psl){
    /* puts all particles in the correct stack if its clusterid = cluster.pic[id].stack
    and counts how many clusters of a certain size there are in cluster.clustersize. Its start counting at length 1 */
    int id,ipart, nclusters=psl->nclusters, n;
    Picl pici;
    int  clustersize[NPART]={0};

    // printf("loop over all particles en deel ze in in een picl\n" );
    //loop over all particles en deel ze in in een picl 
    for(ipart=0;ipart<sys.npart;ipart++){
        id = psl->pts[ipart].cluster_id-1;
        cluster.pic[id].stack[clustersize[id]]=ipart;
        clustersize[id]++;
    }
    memcpy(cluster.clustersize, clustersize, sizeof(clustersize));

    return;
}

double op_calc(Slice *psl) {
    // Compute order paramater ; emilio's order parameter
    Pts *psi;
    int ipart, jpart,n;
    double s,switch_value, Erep, Ec,totE,Ebond;
    double lambda = 0, lambdatmp;

    // iterate over particles, skipping last one (N-1 patch-patch interactions) - Emilio
    // >> walk over the bonds and and calculate OP , it assumes that before this step the total_energy is calculated - Hannah
    for(ipart=0; ipart<sys.npart; ipart++){
        psi = &psl->pts[ipart];
        for (n=0;n<psi->nbonds;n++){ // you will only walk over the bonds
            jpart = psi->bonds[n];
            if (jpart<ipart){continue;} //skip if particle number jpart < ipart, to avoid double counting
            s=psi->bond_distance[n];

            // if bond is in repusive regime, set s to s_min and calc totE, esle totE=bond energy;
            if(s<sys.s_min){

                Erep=potential_repulsive_energy_sdist(s);
                Ec=potential_attractive_energy_sdist(s);
                Ebond=psi->bond_energy[n];
                switch_value=(Ebond-Erep)/Ec;
                
                totE = sys.Erep_smin + switch_value*sys.Ec_smin;
                
            }
            else{
                totE=psi->bond_energy[n];
            }

            lambdatmp = 1 - (totE/sys.E_smin); // value coinciding with minimum potential
            if(lambdatmp > lambda || lambda == 0){
                lambda = lambdatmp;
            }
        }
    }
    if(lambda>1){
        gprint(lambda);
        error("lambda > 1");
    }
    return lambda;
}

double end_to_end_distance2(Slice *psl, int chainid){
    /*calculate the end-to-end distance of a chain (the chainid given as input variable) */
    ChainNode *tail;
    vector tracker=nulvec, dr ;
    double end_to_end2, end_to_end;
    int n;
    
    tail = chain.head[chainid];

    /*check if the start of the chain has only 1 bond. definition of end of chain */
    if(psl->pts[tail->pts_nr].nbonds!=1){
        // printf("head particle %d has more than one bond?? with:\n",tail->pts_nr);
        // for (n=0;n<psl->pts[tail->pts_nr].nbonds;n++){
        //     printf(" %d      bond %d/%d \n",psl->pts[tail->pts_nr].bonds[n],n,psl->pts[tail->pts_nr].nbonds);
        // }
        // printf("\n\n");
        return 0.;
    }



    while(tail->next != NULL){  
        dr = bond_vector(psl, tail->pts_nr,  tail->next->pts_nr, 1);
        vector_add(tracker, dr, tracker); 
        tail = tail->next; // Go To last Node       
    } 
   
    
    /*check if the end of the chain has only 1 bond*/
    if(psl->pts[tail->pts_nr].nbonds!=1){
        // printf("tail particle %d has more than one bond?? with:\n",tail->pts_nr);
        // for (n=0;n<psl->pts[tail->pts_nr].nbonds;n++){
        //     printf(" %d      bond %d/%d \n",psl->pts[tail->pts_nr].bonds[n],n,psl->pts[tail->pts_nr].nbonds);
        // }
        // printf("\n\n");
        return 0.;
    }


    end_to_end2 = vector_inp(tracker,tracker);

    return end_to_end2;
}



ChainNode* GetNewNode(int x) {
    ChainNode *newNode;
    newNode = (ChainNode*)malloc(sizeof(ChainNode));

    newNode->pts_nr = x;
    newNode->prev = NULL;
    newNode->next = NULL;
    return newNode;
}

void InsertAtChainHead(int x, int chainid) {
    //Inserts a Node at head of doubly linked list
    ChainNode *newNode;
    newNode = GetNewNode(x);
    if(chain.head[chainid] == NULL) {
        chain.head[chainid] = newNode;
        return;
    }
    chain.head[chainid]->prev = newNode;
    newNode->next = chain.head[chainid]; 
    chain.head[chainid] = newNode;
    return;
}

void InsertAtChainTail(int x, int chainid) {
    //Inserts a Node at tail of doubly linked list
    ChainNode *newNode, *temp;
    
    temp = chain.head[chainid];
    newNode = GetNewNode(x);
    if(chain.head[chainid] == NULL) {
        chain.head[chainid] = newNode;
        return;
    }
    while(temp->next != NULL) temp = temp->next; // Go To last Node
    temp->next = newNode;
    newNode->prev = temp;
    return;
}

int particle_in_wall(Slice *psl, int ipart){
    /* returns 1 if particle is in the wall, if not, it returns 0. The wall is defined at 1*/
    if(psl->pts[ipart].r.z<=sys.particletype[psl->pts[ipart].ptype].radius){
        return 1;
    }
    else{
        return 0;
    }
}

void reset_running_statistics(Statistics *stats){
    /* reset the statistics; set values to zero  */

    stats->mean=0.;
    stats->variance2=0.;
    stats->n=0;
    return;
}

double running_variance2(double x_new, double s_old2, double u_old, double u_new, long n_old){
    /*  x_new is the new measurement,
        s_old2 is the old variance^2 and stored in chain.e2e2_sigma2
        u_old is the old mean
        u_new is the new mean
        s_old2 is stored in chain.e2e2_sigma2[chainlength];
        (N-1)S_new^2 =  (N-2)S_old^2 +(x_new-u_new)(x_new-u_old)
        
        return  S_{new}^2 */

    double s_new2, s_new2_N_1;
    long n_new=n_old+1;

    s_new2_N_1 = s_old2*n_old + (x_new-u_new)*(x_new-u_old);
    s_new2 = s_new2_N_1 /(n_new);
    return s_new2;
}

double running_mean(double u_old, double x_new, long n_old){
    /*  x_new is the new measurement,
        u_old is the old mean
        n_old is the old number of measurements;
        n_new = n_old+1
        u_new = u_old + (x_new-u_old)/n_new
        
        return  n_new */
    long n_new = n_old+1;
    double u_new;

    u_new = u_old + (x_new-u_old)/n_new;
    return u_new;
}

void running_statistics(Statistics *oldstats, double x_new){
    /* updates the statistics
    it needs the old statistics (oldstats) and the value of the new measurement
    the new mean value and the new variance is calculated and the statistics are updated.  */

    double u_old =oldstats->mean, s_old2 =oldstats->variance2;
    long n_old = oldstats->n;
    double u_new, s_new2;

    u_new = running_mean(u_old,  x_new,  n_old);
    s_new2 = running_variance2( x_new,  s_old2,  u_old,  u_new,  n_old);

    oldstats->mean=u_new;
    oldstats->variance2=s_new2;
    oldstats->n++;
}

void update_patch_vector_ipart(Slice *psl, int ipart){
    Pts *psi;
    int isite,ptype;
    tensor rotmati;
    
    psi = &psl->pts[ipart];
    rotmati = getrotmatrix(psi->q);
    for( isite=0; isite<sys.particletype[psi->ptype].nsites; isite++) {
        ptype=psi->ptype;
        matrix_x_vector(rotmati,sys.particletype[psi->ptype].site[isite].r,psi->patchvector[isite]); 
    }
    return;
}

void update_patch_vectors(Slice *psl){
    int ipart;
    for(ipart=0;ipart<sys.npart;ipart++){
        update_patch_vector_ipart(psl, ipart);
    }
    return;
}
