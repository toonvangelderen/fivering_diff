// #include <stdlib.h>
// #include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"



void copy_clusterparticles( Slice * ,Slice *,int );
int check_internal_energy(Slice *,Slice *, int );

void propagate_mc(Slice *psl) {
    int i,nnlupdate=0, ipart;
    double movetype,dr2;
    int k,j;
    double r;


    //printf("propagating via mc\n");
    for(i=0; i<sys.npart; i++) {
        movetype =  RandomIntegerRange(0, 2);
        ipart = single_particle_move(psl, movetype);
    }
    return;
}

int single_particle_move(Slice *psl, int r_or_t){
    /*one function for either rotate (0) or translate(1)*/
    int ipart, nbonds_old_tot,n,particle;
    double dr2, Eold, Enew, Edif;
    vector dr, dr_dummy;
    quaternion dq,dqrot;

    /*choose a random particle*/
    ipart=(int)(RandomNumber()*sys.npart);
    
    // old single particle info
    Eold=particle_energy(psl,ipart,0); 
    Pts oldpart=psl->pts[ipart];

    if(psl->energy/psl->energy!=1){
        error("energy nan in before  single particle move ");
    }

    if (r_or_t==0){
        trans.tries++; // might want to add new struct of mc_signle
        /*make a random dr vector*/
        dr=RandomVector(sys.drmax);
        if (sys.gravity>0){
            dr.z/=(sys.gravity*10.);
        }
        /*displace the particle*/
        vector_add(psl->pts[ipart].r, dr, psl->pts[ipart].r);

        /*check if the particle is put into the wall, than always reject. No need to calculate the energy.*/
        if(sys.gravity>0){
            if(oldpart.r.z>1.){
            // if(oldpart[ipart].r.z>1.){
                if(particle_in_wall(psl,ipart)==1){
                    psl->pts[ipart]= oldpart;
                    return -1;
                }
            }
        }
        
        pbc(psl->pts[ipart].r,sys.boxl); 
        
    }
    else{
        rot.tries++;
        /*rotate the quaternion*/
        dq=RandomQuaternionRange(sys.dqmax);
        rotate_quaternion_ipart(  psl,  ipart,  dq );
    }
    /*calculate the new particle energy*/
    Enew=particle_energy(psl,ipart,0);
    Edif = Enew - Eold;

    /*metropolis rule*/
    if(exp(-sys.beta*Edif)<RandomNumber()) {
        psl->pts[ipart]= oldpart;
        return -1;
    }

    // reject too if bond breakage is not allowed
    if (sys.bond_breakage==0){
        if ( oldpart.nbonds != psl->pts[ipart].nbonds){
            /* if the translation caused a bond to break/form,if so do a clusterupdate before performing clustermoves*/
            psl->pts[ipart]=oldpart;
            return -1;
        }
        else{
            for(n=0;n<oldpart.nbonds; n++){
                if(oldpart.bonds[n]!=psl->pts[ipart].bonds[n]){
                    psl->pts[ipart]=oldpart;
                    return -1;
                }
            }
        }
    }

    /*accepted*/
    psl->energy+=Edif;
    if(psl->energy/psl->energy!=1){
        printf("Efid %lf\n", Edif);
        printf("Enew %lf -Eold %lf\n",Enew,Eold );
        error("energy nan single particle move");
    }

    /*add acceptance*/
    if (r_or_t==0){
        trans.acc++;
        particle=-1;
    }
    else{rot.acc++;
        particle = ipart; //return the particle number for updating dr (nearestneighbor)
    }

    /* if the translation/rotation caused a bond to break/form,if so do a clusterupdate before performing clustermoves*/
    if (sys.cluster_MC ){
        if ( oldpart.nbonds != psl->pts[ipart].nbonds){ 
            cluster.update=1;
            return -1;
        }
        else{
            for(n=0;n<oldpart.nbonds; n++){
                if(oldpart.bonds[n]!=psl->pts[ipart].bonds[n]){
                    cluster.update=1;
                }
            }
        }
    }
    return particle;
}

void cluster_propagate_mc(Slice *psl){
    /* do a cluster move instead of single particle move*/

    int n,update=0, i, icluster, ipart;
    double nwhich,dr2;


    //printf("propagating via mc\n");
    for(n=0; n<sys.npart; n++) {
        nwhich = RandomNumber();
        if(nwhich<0.5) {
            icluster = translatepart_cluster_Ekparts(psl);
        }
        else{
            icluster = rotatepart_cluster_Ekparts(psl);
        }
    }

    return;
}

void copy_clusterparticles( Slice *psl_new ,Slice *psl_source,int icluster){
    // copies the particles in icluster from psl_source to psl_new
    int clusteri_size=cluster.clustersize[icluster];
    int ipart;

    /*perform copy on selected particles */
    for(int i=0; i<clusteri_size; i++){
        ipart = cluster.pic[icluster].stack[i];
        psl_new->pts[ipart]=psl_source->pts[ipart];
    } 


    return;
}

int translatepart_cluster_Ekparts(Slice *psl) {

    int icluster, nclusters=psl->nclusters;
    int ipart, i, n;
    int update=0, nbonds_old=0, nbonds_new=0, tot_nbonds;
    double Ebond_new,Ebond_old;
    double Eold=0., Enew=0., Edif;

    
    vector dr;

    /*select a cluster randomly*/
    icluster=(int)(RandomNumber()*nclusters);
    int clusteri_size=cluster.clustersize[icluster];
    int cluster_nbonds = clusteri_size-1;


    Pts p_ref;

    /* differentiate between single and multi particle clusters */
    if(clusteri_size==1){
        dr=RandomVector(sys.drmax_mono);
        trans_mono.tries++;
    }
    else{ // printf("picked a cluster\n");
        dr=RandomVector(sys.drmax_cluster);
        trans_cluster.tries++; 
    }

    if (sys.gravity>0){
        dr.z/=(sys.gravity*50.);
        dr.z=0;
    }

    /*save old energy, nbonds, positions*/
    for(i=0; i<clusteri_size; i++){
        ipart=cluster.pic[icluster].stack[i];
        Eold+=particle_energy(psl, ipart,0) ;
        nbonds_old+=(psl->pts[ipart].nbonds);
    }
   
    
    /* save oldparticles*/
    copy_clusterparticles(psl_old , psl, icluster);
    // memcpy(oldpart, psl->pts, sizeof(psl->pts));

    /*perform translation*/
    for(i=0; i<clusteri_size; i++){
        ipart = cluster.pic[icluster].stack[i];
        vector_add(psl->pts[ipart].r,dr,psl->pts[ipart].r); 

        if(sys.nearest_neighbor!=1){
            //perform pbc if no neighborlist is used, you may just moved the particle ipart out of the box 
            pbc(psl->pts[ipart].r,sys.boxl); 
        }

    } 

    /*check if any particle is put into the wall (z<1.0), than always reject. No need to calculate the energy.*/
    if(sys.gravity>0){
        for(i=0; i<cluster.clustersize[icluster]; i++){
            ipart= cluster.pic[icluster].stack[i];
            if((psl_old->pts[ipart].r.z>=1.) && (particle_in_wall(psl,ipart)==1)){
                // memcpy(psl->pts, oldpart, sizeof(psl->pts));
                copy_clusterparticles( psl, psl_old,icluster);
                return -1; 
            }
        }
    }

    /*calc new  energy, nbonds*/
    /* perform energy caluclation always without using the neighborlist*/
    for(i=0; i<clusteri_size; i++){
        ipart = cluster.pic[icluster].stack[i];
        Enew+=particle_energy(psl, ipart,0) ;
        nbonds_new+=(psl->pts[ipart].nbonds);
    }

    //first check if you created bonds. then already reject due to detailed balance
    if (nbonds_new>nbonds_old){
        // memcpy(psl->pts, oldpart, sizeof(psl->pts));
        copy_clusterparticles( psl, psl_old,icluster);
        // printf("rejected based on  nbonds\n");
        return -1;
    }
    else if (nbonds_new<nbonds_old){
        printf("WARING: a bond has been broken during translation cluster move");
    }

    

    /* check if internal energy has stayed equal*/
    int trans_error= check_internal_energy( psl, psl_old,  icluster);

    if (trans_error){
        error("cluster translation gone wrong. in bond energy\n");
        copy_clusterparticles( psl, psl_old,icluster);
        return -1;
    }

    Edif = (Enew - Eold);

    /*if reject based on energy or nbonds */
    if(exp(-sys.beta*Edif)<RandomNumber() ) {
        // memcpy(psl->pts, oldpart, sizeof(psl->pts));
        copy_clusterparticles( psl, psl_old,icluster);
        return -1;
    }

    psl->energy+=Edif;


    if(clusteri_size==1){
        trans_mono.acc++;
    }
    else{
        trans_cluster.acc++; 
    }

    return icluster;
}

int rotatepart_cluster_Ekparts(Slice *psl) {

    // printf("r");
    int ipart_ref, nclusters=psl->nclusters, icluster, ref_random;
    int  nbonds_old=0, nbonds_new=0, tot_nbonds;
    int ipart,i,n;
    double Eold=0., Enew=0., Edif;
    double Ebond_new, Ebond_old;

   
    
    quaternion dq, dqrot;
    tensor R;
    vector dr, rotvec;

    // printf("\n\ncluster rotate, \n");

    /*select a cluster randomly*/
    icluster=(int)(RandomNumber()*nclusters);
    int clusteri_size=cluster.clustersize[icluster];
    int cluster_nbonds = clusteri_size-1;


    if(clusteri_size==1){
        /*the rotation quaternion and matrix*/
        dq=RandomQuaternionRange(sys.dqmax_mono); 
        R = getrotmatrix(dq); 
        rot_mono.tries++;
        ipart= cluster.pic[icluster].stack[0];

        //old energy
        Eold=particle_energy(psl, ipart,0) ;
        nbonds_old=(psl->pts[ipart].nbonds);
        // memcpy(oldpart, psl->pts, sizeof(psl->pts));
        Pts oldpart =psl->pts[ipart];

        //the rotation
        quat_times(dq,psl->pts[ipart].q,dqrot);
        psl->pts[ipart].q = dqrot;
        update_patch_vector_ipart(psl,ipart); 

        //new energy + nbonds
        Enew=particle_energy(psl, ipart,0) ;
        nbonds_new=(psl->pts[ipart].nbonds);

        Edif = (Enew - Eold);

        /*reject if*/
        if(exp(-sys.beta*Edif)<RandomNumber() || nbonds_old!=nbonds_new) {
            // memcpy(psl->pts, oldpart, sizeof(psl->pts));
            psl->pts[ipart]=oldpart;
            // printf("rejected based on MC or nbonds\n");
            return -1;
        }
        /* acccepted*/
        psl->energy+=Edif;
        if(Edif>=10000){
            printf("exp(-sys.beta*Edif) = %lf\n",exp(-sys.beta*Edif));
            printf("Energy OLD   %lf\n", total_energy(psl));
            error("WARINGN: BIG energydifference when accepted the move?? \n");
        }
        
        rot_mono.acc++;

        return icluster;

    }
    else{
        // dq=RandomQuaternionRange(sys.dqmax_cluster); 
        dq=QuaternionZaxis((2*RandomNumber()-1)*sys.dqmax_cluster); 
        R = getrotmatrix(dq); 
        rot_cluster.tries++; 
    }
    Pts p_ref;

    /*save old energy, nbonds*/
    for(i=0; i<clusteri_size; i++){
        ipart= cluster.pic[icluster].stack[i];
        Eold+=particle_energy(psl, ipart,0) ;
        nbonds_old+=(psl->pts[ipart].nbonds);
    }

    /*save old particles information (bv position, bond energy etc.)*/
    // memcpy(oldpart, psl->pts, sizeof(psl->pts));
    copy_clusterparticles(  psl_old,psl,icluster);

    /*pick randomly the reference particle" ref_random -> [0,0.999> pick particle 0 etc*/
    ref_random= (int)(RandomNumber()*clusteri_size);
    /*pick middle particle*/
    // ref_random=floor(clusteri_size/2);

    ipart_ref = cluster.pic[icluster].stack[ref_random];
    p_ref = psl->pts[ipart_ref];

    // printf("linked_structure_icluster\n");
    linked_structure_icluster(psl,icluster,  ipart_ref);
    //p_ref is not in center yet.
    vector dummy,dummy2;

    // printf("now rotating the structure\n");
    /*rotate the cluster, now you can just use the dr as defined below*/
    for(i=0; i<cluster.clustersize[icluster]; i++){
        ipart= cluster.pic[icluster].stack[i];
        vector_minus(psl->pts[ipart].r,psl->pts[ipart_ref].r,dummy); // dummy points from ipart_ref to ipart
        /* rotate the structure around the reference particle*/
        matrix_x_vector(R, dummy, dummy2);
        vector_add(dummy2,psl->pts[ipart_ref].r,psl->pts[ipart].r); // dummy points from ipart_ref to ipart

        // vprint(dummy);
        // vprint(psl->pts[ipart].r);
        
        if(sys.nearest_neighbor!=1){
            /* perform pbc if no neighborlist is used, you may just moved the particle ipart out of the box */
            pbc(psl->pts[ipart].r,sys.boxl); 
        }   

        /* rotate the quaternion*/
        rotate_quaternion_ipart(  psl,  ipart,  dq );
    }

    // printf("done rotating structure\n");

    /*check if any particle is put into the wall, than always reject. No need to calculate the energy.*/
    if(sys.gravity>0){
        for(i=0; i<cluster.clustersize[icluster]; i++){
            ipart= cluster.pic[icluster].stack[i];
 
            if(psl_old->pts[ipart].r.z>=1.){
                if(particle_in_wall(psl,ipart)==1){
                    printf("particle_in_wall\n");

                    // memcpy(psl->pts, oldpart, sizeof(psl->pts));
                    copy_clusterparticles( psl, psl_old,icluster);
                    return -1;
                }
            }
        }
    }

   

    /*calc new  energy, nbonds*/
    /* perform energy caluclation always without using the neighborlist*/
    for(i=0; i<clusteri_size; i++){
        ipart = cluster.pic[icluster].stack[i];
        Enew+=particle_energy(psl, ipart,0) ;
        nbonds_new+=(psl->pts[ipart].nbonds);
    }

    /*reject based on bond formation*/
    if (nbonds_new>nbonds_old){
        // memcpy(psl->pts, oldpart, sizeof(psl->pts));
        copy_clusterparticles( psl, psl_old,icluster);
        // printf("rejected based on  nbonds\n");
        return -1;
    }
    else if (nbonds_new<nbonds_old){
        printf("WARING: a bond has been broken during rotation cluster move");
    }
    

    /* check if internal energy has stayed equal*/
    int rot_error= check_internal_energy( psl, psl_old,  icluster);

    if (rot_error){
        error("cluster rotation gone wrong. in bond energy\n");
        copy_clusterparticles( psl, psl_old,icluster);
        return -1;
    }
    
    Edif = (Enew - Eold);

    /*reject if*/
    if(exp(-sys.beta*Edif)<RandomNumber() ) {
        // memcpy(psl->pts, oldpart, sizeof(psl->pts));
        copy_clusterparticles( psl, psl_old,icluster);
        return -1;
    }


    /* acccepted*/
    psl->energy+=Edif;
    if(Edif>=10000){
        printf("exp(-sys.beta*Edif) = %lf\n",exp(-sys.beta*Edif));
        printf("Energy OLD   %lf\n", total_energy(psl));
        error("WARINGN: BIG energydifference when accepted the move?? \n");
    }

    // printf("*** ROTATION SUCCESFULL **\n");
    rot_cluster.acc++; 
    
    return icluster;
}
void energy_divergence_check(Slice *psl, char loc[50]){
    // checks for the divergence of the energy due to the use of a running energy in the MC code
    if(fabs(total_energy(psl)-psl->energy)>.001){
        printf("there is a difference in the recalc.  after %s\n", loc);
        printf("recalc total_energy              %lf\n", total_energy(&slice[0]));
        printf("tracked slice[0].energy               %lf\n", slice[0].energy);
        // update_nnlist(&slice[0]);

        if(fabs(total_energy(psl)-psl->energy)>1){
            printf("there is still a difference in the recalc.  after %s\n", loc);
            printf("recalc total_energy              %lf\n", total_energy(&slice[0]));
            printf("tracked slice[0].energy               %lf\n", slice[0].energy);
            error("end calc");

        }

    }
    return;
}


int check_internal_energy(Slice *psl_new,Slice *psl_old, int icluster){
    int clusteri_size=cluster.clustersize[icluster];
    double Ebond_new,Ebond_old;
    int cluster_error=0;

    for(int i=0; i<clusteri_size; i++){  
        int ipart= cluster.pic[icluster].stack[i];

        for(int n=0;n<psl_new->pts[ipart].nbonds;n++){
            if ( psl_old->pts[ipart].nbonds != psl_new->pts[ipart].nbonds){ // "<"bonds are broken? 
                // if (psl_new->pts[ipart].bonds[n] != psl_old->pts[ipart].bonds[n]){ 
                printf("number of bonds has changed from %d (old) to %d (new) for particle %d\n", psl_old->pts[ipart].nbonds,psl_new->pts[ipart].nbonds, ipart);
                // continue;
                for(int m=0;m<psl_new->pts[ipart].nbonds;m++){
            
                    Ebond_new=psl_new->pts[ipart].bond_energy[m];
                    Ebond_old=psl_old->pts[ipart].bond_energy[m];
                    dprint(m);
                    dprint(psl_new->pts[ipart].bonds[m]);
                    dprint(psl_old->pts[ipart].bonds[m]);
                    if(fabs(Ebond_new-Ebond_old) >1e-9){
                        rot_cluster.tries--; 
                        
                        
                        printf("ipart: %d, Ebond_new = %.12lf , Ebond_old=%.12lf,      difference = %.12f\n",ipart,Ebond_new,Ebond_old, fabs(Ebond_new-Ebond_old));
                        
                        cluster_error=1;
                        // printf("rejected based on internal bond_energy\n");
                        // return -1;
                    }

                    printf(" psl the other particle %d has %d bonds \n",psl_new->pts[ipart].bonds[m],psl_new->pts[psl_new->pts[ipart].bonds[m]].nbonds);
                    printf(" psl_old the other particle %d has %d bonds \n",psl_old->pts[ipart].bonds[m],psl_old->pts[psl_old->pts[ipart].bonds[m]].nbonds);
                }
                printf("\n");
                for(int m=0;m<psl_old->pts[ipart].nbonds;m++){
            
                    Ebond_new=psl_new->pts[ipart].bond_energy[m];
                    Ebond_old=psl_old->pts[ipart].bond_energy[m];
                    dprint(m);
                    dprint(psl_new->pts[ipart].bonds[m]);
                    dprint(psl_old->pts[ipart].bonds[m]);
                    if(fabs(Ebond_new-Ebond_old) >1e-9){
                        rot_cluster.tries--; 
                        
                        
                        printf("ipart: %d, Ebond_new = %.12lf , Ebond_old=%.12lf,      difference = %.12f\n",ipart,Ebond_new,Ebond_old, fabs(Ebond_new-Ebond_old));
                        
                        cluster_error=1;
                        // printf("rejected based on internal bond_energy\n");
                        // return -1;
                    }

                    printf(" psl the other particle %d has %d bonds \n",psl_new->pts[ipart].bonds[m],psl_new->pts[psl_new->pts[ipart].bonds[m]].nbonds);
                    printf(" psl_old the other particle %d has %d bonds \n",psl_old->pts[ipart].bonds[m],psl_old->pts[psl_old->pts[ipart].bonds[m]].nbonds);
                }
                printf("\n");
            }
        }
    }
    return cluster_error;
}




void optimizemc() {

    /*single particle moves*/
    static int initiate=1;

    double maxdr;
    MIN(sys.boxl.x,sys.boxl.y,maxdr);


    if (initiate ){
        initiate=0;
        return;
    }
    /*calculate the ratio*/
    trans.ratio=(double)trans.acc/(double)trans.tries;
    rot.ratio=(double)rot.acc/(double)rot.tries;

    /*update the total accepted, tried translate or rotate moves*/
    fintrans.acc+=trans.acc;
    finrot.acc+=rot.acc;
    fintrans.tries+=trans.tries;
    finrot.tries+=rot.tries;
    
    /*set accept. and tries to zero*/
    trans.acc=0;
    trans.tries=0;
    rot.acc=0;
    rot.tries=0;
   
    /*single particle*/
    if(trans.ratio<0.3 && sys.drmax>0.0001) {
        sys.drmax/=1.1;
    }
    else if(trans.ratio>0.7 && sys.drmax<maxdr) {
        if(sys.nearest_neighbor){
            /*  there is a restriction in drmax because if nnl is used, the particles sees other particles upto r=2 and the update treshold is r_treshold=0.4
             the drmax cannot be bigger than r_treshold = sqrt(list.cutoff_buffer2),  because the particle cannot see further than that. 
             Before a particle is moved, there is a check whether the move will cause a crossing of the treshold. If so, the energy w/o nnl is performed to check acceptance.
             After acception, the nnl is updated.  */
            if(sys.drmax<(0.90*sqrt(list.cutoff_buffer2))){
                sys.drmax*=1.1;
            }
            else{
                sys.drmax=(0.90*sqrt(list.cutoff_buffer2));
            }
        }
        else{
            sys.drmax*=1.1;
        }

    }
    if(rot.ratio<0.3 && sys.dqmax>0.0001) {
        sys.dqmax/=1.1;
    }
    else if(rot.ratio>0.7 && sys.dqmax<PI) {
        sys.dqmax*=1.1;
    }
    else if(sys.dqmax>PI){
        sys.dqmax=PI;
    }
    
    /*cluster moves*/
    if(sys.cluster_MC==1){
        /*calculate the ratio translate and rotate cluster */
        trans_cluster.ratio=(double)trans_cluster.acc/(double)trans_cluster.tries;
        rot_cluster.ratio=(double)rot_cluster.acc/(double)rot_cluster.tries;

        /*update the total accepted, tried translate or rotate moves*/
        fintrans_cluster.acc+=trans_cluster.acc;
        finrot_cluster.acc+=rot_cluster.acc;
        fintrans_cluster.tries+=trans_cluster.tries;
        finrot_cluster.tries+=rot_cluster.tries;

        /*update the total accepted, tried translate or rotate moves*/
        fintrans_mono.acc+=trans_mono.acc;
        fintrans_mono.tries+=trans_mono.tries;
        

        /*set accept. and tries to zero*/
        trans_cluster.acc=0;
        trans_cluster.tries=0;
        rot_cluster.acc=0;
        rot_cluster.tries=0;

        if(trans_cluster.ratio<0.3 && sys.drmax_cluster>0.0001) {
            sys.drmax_cluster/=1.1;
        }
        else if(trans_cluster.ratio>0.7 && sys.drmax_cluster<maxdr) {
            sys.drmax_cluster*=1.1;
        }
        if(rot_cluster.ratio<0.1 && sys.dqmax_cluster>0.00001) {
            sys.dqmax_cluster/=1.1;
        }
        else if(rot_cluster.ratio>0.7 && sys.dqmax_cluster<179) {
            sys.dqmax_cluster*=1.1;
        }
        else if(sys.dqmax_cluster>180.){
            sys.dqmax_cluster=179.;//PI;
        }

        /*calculate the ratio translate and rotate single particle in clustermove (mono) */
        rot_mono.ratio=(double)rot_mono.acc/(double)rot_mono.tries;
        trans_mono.ratio=(double)trans_mono.acc/(double)trans_mono.tries;

        /*update the total accepted, tried translate or rotate moves*/
        fintrans_mono.acc+=trans_mono.acc;
        fintrans_mono.tries+=trans_mono.tries;
        finrot_mono.acc+=rot_mono.acc;
        finrot_mono.tries+=rot_mono.tries;

        /*set accept. and tries to zero*/
        trans_mono.acc=0;
        trans_mono.tries=0;
        rot_mono.acc=0;
        rot_mono.tries=0;

        if(trans_mono.ratio<0.3 && sys.drmax_mono>0.0001) {
            sys.drmax_mono/=1.1;
        }
        else if(trans_mono.ratio>0.7 && sys.drmax_mono<maxdr) {
            sys.drmax_mono*=1.1;
        }
        if(rot_mono.ratio<0.3 && sys.dqmax_mono>0.0001) {
            sys.dqmax_mono/=1.1;
        }
        else if(rot_mono.ratio>0.7 && sys.dqmax_mono<PI) {
            sys.dqmax_mono*=1.1;
        }
        else if(sys.dqmax_mono>PI){
            sys.dqmax_mono=PI;
        }


    }


    return;
}

void printstatusmc(Slice *psl) {

    double en,l;
    static int initiate=1;


    if (initiate ){
        initiate=0;
        return;
    }
   
    printf("Energy system %lf\n",slice[0].energy);

    trans.ratio=(double)trans.acc/(double)trans.tries;
    rot.ratio=(double)rot.acc/(double)rot.tries;
    printf("Single particle moves:\n");
    printf("Accepted %10ld translations from %10ld tries, ratio %lf\n",trans.acc,trans.tries,trans.ratio);
    printf("Accepted %10ld rotations    from %10ld tries, ratio %lf\n",rot.acc,rot.tries,rot.ratio);
    printf("drmax = %lf\n",sys.drmax);
    printf("dqmax = %lf\n",sys.dqmax);
    printf("\n");

    if(sys.cluster_MC==1){
        trans_cluster.ratio=(double)trans_cluster.acc/(double)trans_cluster.tries;
        rot_cluster.ratio=(double)rot_cluster.acc/(double)rot_cluster.tries;

        printf("Cluster moves:\n");
        printf("Accepted %10ld translations from %10ld tries, ratio %lf\n",trans_cluster.acc,trans_cluster.tries,trans_cluster.ratio);
        printf("Accepted %10ld rotations    from %10ld tries, ratio %lf\n",rot_cluster.acc,rot_cluster.tries,rot_cluster.ratio);
        printf("drmax_cluster = %lf\n",sys.drmax_cluster);
        printf("dqmax_cluster = %lf\n",sys.dqmax_cluster);
        printf("\n");

        trans_mono.ratio=(double)trans_mono.acc/(double)trans_mono.tries;
        rot_mono.ratio=(double)rot_mono.acc/(double)rot_mono.tries;

        printf("Mono (cluster) moves:\n");
        printf("Accepted %10ld translations from %10ld tries, ratio %lf\n",trans_mono.acc,trans_mono.tries,trans_mono.ratio);
        printf("Accepted %10ld rotations    from %10ld tries, ratio %lf\n",rot_mono.acc,rot_mono.tries,rot_mono.ratio);
        printf("drmax_mono = %lf\n",sys.drmax_mono);
        printf("dqmax_mono = %lf\n",sys.dqmax_mono);
        printf("\n");

    }

    return;
}





