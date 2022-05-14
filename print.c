#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "path.h"

void printstatusbmd();
void r_S_hist_2D(Slice *);
// here all the functino for terminate block

void terminate_block() {

    int count_npart=0, monomer=0,emptycluster=0;
    int n, id,nbonds_old,nbonds_new;
    double e2e2;
        
    static int dipatch=1;
    for (n=0;n<sys.nparticle_types;n++){
        if (sys.particletype[n].nsites>2){
            dipatch=0;
            break;
        }
    }

    nbonds_old=slice[0].nbonds;
    slice[0].energy=total_energy(&slice[0]);
    printf("the total energy is %lf\n", slice[0].energy);
    nbonds_new=slice[0].nbonds;
    
    printenergy();
    conf_output(&slice[0]);

    // identify clusters and count bonds
    slice[0].nclusters = cluster_analysis(&slice[0]);
    
    
    for(id=0;id<sys.npart;id++){
        /* add up all clustersizes; it should be equal to npart*/
        count_npart +=cluster.clustersize[id];
        if(cluster.clustersize[id]==0){
            emptycluster ++;
        }
        if(cluster.clustersize[id]==1){
            monomer ++;
        }
    }
    printf("\n *****\n");
    printf("sys.npart     = %d          = count_npart = %d \n", sys.npart,  count_npart);
    printf("emptyclusters = %d            = nbonds      = %d       in %d clusters with %d monomers\n",emptycluster , slice[0].nbonds, slice[0].nclusters, monomer);
    printf("nbonds + nclusters = %d +  %d  = %d   =  sys.npart     = %d\n", slice[0].nbonds, slice[0].nclusters , slice[0].nclusters+slice[0].nbonds, sys.npart);
    printf(" *****\n");

    // printf("start terminating block\n");
    // if(sys.rdfanalysis==1){
    //     calc_rdf_localCN(&slice[0] );
    // }
    if(sys.s_distribution==1){
        // r_S_hist_2D(&slice[0]);
        // s_distribution_analysis(&slice[0]);
        print_energy_s_theta1_theta2(&slice[0]);

    }
    
    if(sys.sim_type==0){         //bmd
        printstatusbmd(); // print the energy

        if ((sys.bond_breakage==0) & (nbonds_old != nbonds_new) & sys.c_time > 0 & (sys.s_distribution!=1)){ 
            dprint(sys.bond_breakage);
            dprint(nbonds_old);
            dprint(nbonds_new);
            gprint(sys.c_time);
            dprint(sys.s_distribution);
            error("bond has broken and it's not allowed (sys.bond_breakage==0), stop the simulation");
        }
       
    }
    else if(sys.sim_type==2) {   //mc   
        printstatusmc(&slice[0]);     
        optimizemc();
    }

    if(sys.xy_print==1){
        print_xy_positions(&slice[0]);
    }

    if(cluster.analysis==1){
        if (sys.adjacency){
            print_adjacency_matrix_coordinates(&slice[0]);
        }
        
        clustersize_identification(&slice[0]);
        clustersize_freq_update(&slice[0]);
        print_statistics_file(&cluster.size_histogram);
        print_statistics_file(&cluster.size_distribution);

        bond_length_average(&slice[0]); 
        printf("bond length average =  %.12lf \n", bond_length.mean);
        printf("bond length var^2   =  %.12lf \n", bond_length.variance2); 

        // make  radius of gyration instead of end-to-end
        //rewrite this part still!!
        // do chain analysis if there are dipatch particles only
        
        // xy_plane_patchangle_calc(&slice[0]);
        if (sys.npart==1){
            timecorrelation_function_P(&slice[0]);
            MSD(&slice[0]);
        }
        if (dipatch==1){
            printf("Performing analysis on chain lengths  \n");

            /*identify the chains, which particles belong to the chains, what is the frequency of occurance of each length, and print*/
            // clustersize_identification(&slice[0]);            
            chain_linkedlist(&slice[0]); /*we need linked list for calculating e2e2*/
        
            /* save end_to_end_distance in a average variable
            loop over the chainid (nclusters) and calculate the end-to-end-distance
            acumulate the average e2e2 for each chainlength*/

            for(id=0;id<slice[0].nclusters;id++){
                // dprint(cluster.clustersize[id]);
                if(cluster.clustersize[id]>1){
                    // print_chain_order(&slice[0],id);
                    /*this function return 0 if the chain is a ring*/
                    e2e2 = end_to_end_distance2(&slice[0], id);
                    // printf("e2e2 is %lf\n", e2e2);
                    if(e2e2>1e-15){
                        /*if not a ring*/
                        running_statistics(&chain.e2e2.length[cluster.clustersize[id]-1], e2e2);
                    }     
                         
                }            
            }

            nematic_order_parameter(&slice[0]);
            double Sfraci=Sfrac(&slice[0]);
            running_statistics(&Sfraction ,Sfraci);
            print_statistics_N1(Sfraction,getName(Sfraction));

            print_statistics_file(&chain.e2e2);
        }
    }

    
    return;
}

void printstatusmc_II(Mcacc *type){

    type->ratio=(double)type->acc/(double)type->tries;  
    printf("Accepted %21ld translations from %21ld tries, ratio %lf\n",type->acc,type->tries,type->ratio);
    return;
}

void finalstat() {

    int i,istate,jstate,irep,id;
    // printf("starting with finalstat\n");
    printf("\nFinished Simulation\n");


    if(sys.sim_type==0) {
        printstatusbmd(); // print the energy
    }
    else if(sys.sim_type==2) {
        printf("\n*** MC stats **\n");
        fintrans.ratio=(double)fintrans.acc/(double)fintrans.tries;
        finrot.ratio=(double)finrot.acc/(double)finrot.tries;
        
        printf("Single particle moves:\n");
        // printstatusmc(&fintrans);
        // printstatusmc(&finrot);

        printf("Accepted %21ld translations from %21ld tries, ratio %lf\n",fintrans.acc,fintrans.tries,fintrans.ratio);
        printf("Accepted %21ld rotations    from %21ld tries, ratio %lf\n",finrot.acc,finrot.tries,finrot.ratio);

        if(sys.cluster_MC==1){
            fintrans_cluster.ratio=(double)fintrans_cluster.acc/(double)fintrans_cluster.tries;
            finrot_cluster.ratio=(double)finrot_cluster.acc/(double)finrot_cluster.tries;
            printf("Cluster moves:\n");
            printf("Accepted %21ld translations from %21ld tries, ratio %lf\n",fintrans_cluster.acc,fintrans_cluster.tries,fintrans_cluster.ratio);
            printf("Accepted %21ld rotations    from %21ld tries, ratio %lf\n",finrot_cluster.acc,finrot_cluster.tries,finrot_cluster.ratio);

            fintrans_mono.ratio=(double)fintrans_mono.acc/(double)fintrans_mono.tries;
            finrot_mono.ratio=(double)finrot_mono.acc/(double)finrot_mono.tries;
            printf("Mono (cluster) moves:\n");
            printf("Accepted %21ld translations from %21ld tries, ratio %lf\n",fintrans_mono.acc,fintrans_mono.tries,fintrans_mono.ratio);
            printf("Accepted %21ld rotations    from %21ld tries, ratio %lf\n",finrot_mono.acc,finrot_mono.tries,finrot_mono.ratio);
        }
    }
    else if(sys.sim_type==4){
        int p;
        FILE *fp;

        if((fp=fopen("ffs_output/freq.dat","w"))==NULL)
            error("Warning: can not open files\n");


        fprintf(fp, "path count\n");
        for(p=0; p<sys.ncycle1; p++){
            fprintf(fp,"%d %d\n",p,state[0].freq[p]);
        }
        fclose(fp);
    }
    else{
        error("use only simtype 2 for MC");
    }
    

    //  print RDF to file if rdfanalysis is turned on
    // if(sys.rdfanalysis==1){
    //     printrdf();
    // }


    
    return;
}

void print_adjacency_matrix_coordinates(Slice *psl){
    /* the adjacency matrix is a NxN matrix with zero's and one's
    0 if no bond between particle i and j, and a 1 if there is a bond.
    If the total energy is calculated, the bonds are tracked.  see potential_energy(Slice *psl)
    We will first try to print the whole complete adjacency matrix.
    And evaluate later if the files are not getting too big. THere is a lot of redundant data
    We might want to only print which particles have a bond, and later construct the adjacency matrix in python
    */
    int ipart, jpart,n;
    FILE *fp;
    Pts *psi;

    if((fp=fopen("adjacency_matrix_coordinates.dat","a"))==NULL) {;
        printf("Warning: can not open adjacency_matrix.dat\n");
    }
    else {
        /*we loop over the particles. 
        in psi->bonds[psi->nbonds]=jpart the particle numbers are printed which have a bond with psi
        start with n=0, there are no bond found with ipart yet*/
        for(ipart=0; ipart<sys.npart; ipart++) {
            psi = &psl->pts[ipart];
            for(n=0; n<psi->nbonds; n++){
                jpart=psi->bonds[n];
                if(jpart>ipart){
                    fprintf(fp,"%d,%d\n", ipart,jpart); 
                }
            }
        }
    }
    fprintf(fp,"\n");    
    fclose(fp);
    return;
}


void printstatusbmd() {
    double energy,en,op,l;
    char chr;

    printf("printstatusbmd\n");

    printf("cumulative time  %lf s\n", sys.c_time);

    return;
}



void conf_output(Slice *psl) {

    int ipart, isite;
    FILE *fp;


    if((fp=fopen("conf.out","w"))==NULL) {;
        printf("Warning: can not open conf.out\n");
    }
    else {
        fprintf(fp,"%d %.6f  %.6f  %.6f \n", sys.npart, sys.boxl.x,sys.boxl.y,sys.boxl.z);
        for(ipart=0; ipart<sys.npart; ipart++) {
            fprintf(fp,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", psl->pts[ipart].r.x, psl->pts[ipart].r.y, psl->pts[ipart].r.z,
                                 psl->pts[ipart].q.q0, psl->pts[ipart].q.q1, psl->pts[ipart].q.q2, psl->pts[ipart].q.q3);
        } }
    fclose(fp);
    return;
}


void conf_input(Slice *psl) {

    int ipart, isite,npart,nsites;
    vector boxl;
    FILE *fp;

    if((fp=fopen("conf.inp","r"))==NULL) {;
        printf("Warning: can not open conf.inp\n");
    }
    else {
        fscanf(fp,"%d %lf %lf %lf\n", &npart, &boxl.x, &boxl.y, &boxl.z);
        for(ipart=0; ipart<npart; ipart++) {
            fscanf(fp,"%lf %lf %lf ", &psl->pts[ipart].r.x, &psl->pts[ipart].r.y, &psl->pts[ipart].r.z);
            fscanf(fp,"%lf %lf %lf %lf\n", &psl->pts[ipart].q.q0, &psl->pts[ipart].q.q1, &psl->pts[ipart].q.q2, &psl->pts[ipart].q.q3);
        }
    }
    fclose(fp);

    if(npart!=sys.npart) {
        error("Error: number of particles in system not same as in conf.inp\n");
        
    }
    if((boxl.x!=sys.boxl.x) | (boxl.y!=sys.boxl.y) | (boxl.z!=sys.boxl.z)) {
        vprint(boxl);
        vprint(sys.boxl);
        error("Error: boxlength in system not same as in conf.inp\n");
    }

    return;
}




void print_chain_order(Slice *psl, int chainid){
    ChainNode *tail;   
    int ibond;

    tail = chain.head[chainid];
    while(tail->next != NULL){ 
        printf("%d, ", tail->pts_nr );      
        tail = tail->next; // Go To last Node
    } 
    printf("\n");

    tail = chain.head[chainid];
    while(tail->next != NULL){ 
        printf("%d, ", psl->pts[tail->pts_nr].nbonds );      
        tail = tail->next; // Go To last Node
    } printf("\n");

    tail = chain.head[chainid];
    while(tail->next != NULL){ 
        for(ibond=0;ibond<psl->pts[tail->pts_nr].nbonds; ibond++){
            printf("%d, ", psl->pts[tail->pts_nr].bonds[ibond] );    
        }
        printf("\n\n");  
        tail = tail->next; // Go To last Node
    } printf("\n");

    return;
}

void printenergy(){
    /*Printing the energy to a file. Each printed energy is after ncycle2 steps */
    FILE *energyfile;
    char energyfilename[100];
    char* a = "energy";
    char* extension = ".dat";
    double cov;

    cov=sys.npart/(sys.boxl.x*sys.boxl.y);

    snprintf( energyfilename, sizeof( energyfilename ), "%s_dT%2.2lf_cov%.2lf_%s", a, sys.dT, cov, extension );
    energyfile = fopen(energyfilename,"a");


    if (energyfile == NULL){
            printf("Error with opening \"energy.dat\" file");
    }
    else{    
        if(sys.sim_type==0){
            fprintf(energyfile, "%8.12lf %8.12lf\n", sys.c_time, slice[0].energy);  }
        else{
                fprintf(energyfile, "%8.12lf\n", slice[0].energy);   }  
    }
    fclose(energyfile);

    return;

}

void r_S_hist_2D(Slice *psl){
    /* only bonds are analyzed*/
    Pts *psi,*psj;
    int ipart, n, jpart,r_bin,S_bin,x,y;
    double dr=sys.rcutoff/MAXBIN, dS=1./MAXBIN,S=dS*0.5,R=dr*0.5;
    double s,Erep,Ebond,Ec,switch_value,Derj_pref;
    FILE *file;
    // printf("\n\n*********inside r_S_hist_2D***************\n");
    // total_energy(psl);

    for(ipart=0; ipart<sys.npart; ipart++){
        psi = &psl->pts[ipart];
        // printf("particle %d has %d bonds with: ", ipart,psi->nbonds);
        // for (n=0;n<psi->nbonds;n++){
        //     printf("%d (n=%d), ", psi->bonds[n],n);
        // } 
        // printf("\n");  
        for (n=0;n<psi->nbonds;n++){ // you will only walk over the bonds
            jpart = psi->bonds[n];
            if (jpart<ipart){continue;} //skip if particle number jpart < ipart, to avoid double counting
            // printf("\n");
            s=psi->bond_distance[n];
            Derj_pref=sys.particletype[psi->ptype].Derj_scaling[psl->pts[jpart].ptype];
            Erep=potential_repulsive_energy_sdist(s)*Derj_pref;
            Ec=potential_attractive_energy_sdist(s)*Derj_pref;
            Ebond=psi->bond_energy[n];
            switch_value=(Ebond-Erep)/Ec;


            if (switch_value>1){
                error("S is bigger than 1? in r_S_hist_2D");
            }
            r_bin=floor(s/dr);
            S_bin=floor(switch_value/dS);

            // gprint(s);
            // gprint(dr);
            // dprint(r_bin);
            // gprint(switch_value);
            // gprint(dS);
            // dprint(S_bin);
            r_S_2dhist[r_bin][S_bin]++;
            // printf("\n");
        }
    }
    file = fopen("r_S_2dhistogram.dat","w");
    printf("printing to \"r_S_2dhistogram.dat\"\n ");

    if (file == NULL){
            printf("Error with opening \"r_S_2dhistogram.dat\" \n");
    }


    for (x=0; x<MAXBIN;x++){
        S=dS*0.5;
        for (y=0; y<MAXBIN;y++){
            fprintf(file, "%lf,%lf,%d\n",  R,S,r_S_2dhist[x][y]);
            S+=dS;
        }
        R+=dr;
    }
    
    fclose(file);

    return;

}





void print_energy_s_theta1_theta2(Slice *psl){
    /*Printing the energy to a file. Each printed energy is after ncycle2 steps */
    FILE *energyfile;

    double   s_treshold=0.5;
    Pts *psi, *psj;
    int i,j, isite, jsite;
    vector rij,p1,p2,u1;
    double rad_to_degree=180./PI;
    double s_bb,s_bp1,s_bp2,s_pp,cosphi_i,cosphi_j,cosphi_ij;
    double r,r2,s,r_radii,Derj_pref;
    double Erep_pp,Erep_bb;
    double Eatr,S,Erep,Ebond;
    double treshold_anglei,treshold_anglej;



    if((energyfile=fopen("ctime_s_theta1_theta2_energy.dat","a"))==NULL) {;
            printf("Error with opening \"print_energy_s_theta1_theta2.dat\" file");
    }
    else{  
        for (i=0;i<sys.npart; i++){
            for (j=i+1;j<sys.npart; j++){
                psi = &psl->pts[i];
                psj = &psl->pts[j];
                vector_minus(psi->r,psj->r,rij);
                pbc(rij,sys.boxl); /*do perform pbc when using no nnl*/

                r2 = vector_inp(rij,rij);
                r=sqrt(r2);
                scalar_divide(rij,r,u1);
                r_radii=sys.particletype[psi->ptype].radius+sys.particletype[psj->ptype].radius;
                s=r-r_radii;
                if(s<s_treshold && s>0.0) { 
                    treshold_anglei=sys.particletype[psj->ptype].cosphi;
                    treshold_anglej=sys.particletype[psi->ptype].cosphi;

                     /*  calculate here the s, s_pp,s_bp and the angles that are use for the switch function*/
                    orientation_parameters(psl,i,j, &s_bb,&s_bp1,&s_bp2,&s_pp, &cosphi_i, &cosphi_j, &cosphi_ij);
                    if (fabs(s-s_bb)>1e-10){
                        gprint(s);
                        gprint(s_bb);
                        error("something went wrong in the calculation of the distance between big particles");
                    }

                    //the derjaguin scaling of this combintaion of particle types
                    Derj_pref=sys.particletype[psi->ptype].Derj_scaling[psj->ptype];

                    // patch attraction
                    Eatr=potential_attractive_energy_sdist(s_pp)*Derj_pref;
                    if ((cosphi_i<treshold_anglei) || (cosphi_j<treshold_anglej)){
                        S=0;
                    }
                    else{
                        S=S_value(psl,cosphi_i,cosphi_j, nulvec,nulvec, nulvec, i, j); // I can use nulvec here, because I take S90 anyway
                    }
                    Ebond=Eatr*S;

                    //repulsions; I only include patch-patch  and bulk-bulk repulsion for now
                    
                    Erep_pp=potential_repulsive_energy_sdist(s_pp)*Derj_pref; //patch-patch repulsion
                    Erep_bb=potential_repulsive_energy_sdist(s_bb); //bulk-bulk repulsion
                    
                    Erep=Erep_pp+Erep_bb;
                            
                    fprintf(energyfile, "%8.6lf %8.8lf %8.8lf %8.8lf %8.8lf %8.8lf %8.8lf %8.8lf %8.8lf %8.8lf\n", sys.c_time,  s_bb,s_pp, acos(cosphi_i)*rad_to_degree, acos(cosphi_j)*rad_to_degree, Erep+Eatr*S, S, Erep_pp,Erep_bb, Eatr );     
                        
                    
                }
            }
        }
    }
    fclose(energyfile);

    return;

}


void print_pos_sites(Slice *psl){
    //looop over the particles and print its position and its patches
    int ipart,p, s;
    Pts *psi;

    for(p=0;p<sys.nparticle_types;p++){
        printf("%d particles have patch: \n",  sys.particletype[p].nparticles);
        for(s=0; s<sys.particletype[p].nsites;s++){
            vprint(sys.particletype[p].site[s].r);
        }
    }
    return;
}


// analysis
// void printrdf() {
//     /* this function prints the RDF, which is stored in rdf_averag, to the file "rdf.dat". The rdf per slice is calculated in calc_rdf_localCN() */
//     int i;
//     FILE *rdffile;
//     double maxl;
//     MIN(sys.boxl.x,sys.boxl.y,maxl);

//     rdffile = fopen("rdf.dat","w");
//     printf("\nPrinting the RDF to rdf.dat\n");

//     if (rdffile == NULL){
//             printf("Error with opening \"rdf.dat\" file");
//     }
//     else{
//         for(i=0; i<RDFBINS/2; i++){
//             fprintf(rdffile, "%2.4lf  %.5lf\n", i*maxl/RDFBINS, rdf_average[i].sum/rdf_average[i].n);
//         }
//     }
//     fclose(rdffile);
//     return;
// }

void print_statistics_N1(Statistics stats_name, char filename1[100]){
    FILE *file;
    char *pt,line[NPART], filename[100];
    int i=0, dummy;
    char* extension = ".dat";
  
    snprintf(filename, sizeof(filename), "%s%s", filename1,extension);
    
    
    if ((file = fopen(filename,"w"))==NULL){
        printf("%s \n", filename);
        error("input: can't be opened \n");
    }
    else{
        printf("printing statistics to %s\n", filename );
        fprintf(file, "%8.12lf %8.12lf %ld\n",  stats_name.mean, stats_name.variance2,stats_name.n); 
            
    }
    fclose(file);

    return;
}

void print_statistics_file(StatsLength *stats_name){
    FILE *file;
    char *pt,line[NPART], filename[100];
    int i=0, dummy;

    memcpy(filename,stats_name->filename,sizeof(filename));  
    
    if ((file = fopen(filename,"w"))==NULL){
        printf("%s \n", filename);
        error("input: can't be opened \n");
    }
    else{
        printf("printing statistics to %s\n", filename );
        for(i=0;i<sys.npart;i++){ /*the length histogram starts with length=1; the minimum length of a chain. therefore i+1*/ 
            fprintf(file, "%d %8.12lf %8.12lf %ld\n", i+1, stats_name->length[i].mean, stats_name->length[i].variance2,stats_name->length[i].n); 
            // printf("%d %8.12lf %8.12lf %d\n", i, stats_name->length[i].mean, stats_name->length[i].variance2,stats_name->length[i].n);
        }
    }
    fclose(file);

    return;
}

void print_xy_positions(Slice *psl){
    // x   y   fx  fy  particle    time
    int ipart;
    static int linecount_xypos=0;
    static double cycle_mc=0;
    FILE *posfile;
    char filename[100];
    char* a = "xypos";
    char* extension = ".csv";

    vector v1,r_ref,dv,dv_norm,dr,pi_shifted,dr_ipart;
    vector yaxis=nulvec,xaxis=nulvec,drcheck;
    vector v[sys.npart],v_new[sys.npart];
    double costheta,rotangle;
    quaternion q;
    tensor R;
    vector rotvec;

    yaxis.y=1.;
    xaxis.x=1.;


    snprintf( filename, sizeof( filename ), "%s%s", a,  extension );
    
    if( access( filename, F_OK ) == -1 ) {
        posfile = fopen(filename,"w");
        if (posfile == NULL){
            printf("Error with creating \"%s\" file",filename);
        }
        else{
            printf("creating the file %s",filename);
            fprintf(posfile,",x,y,fx,fy,particle,time\n");
        }
        fclose(posfile);
    }
    
    posfile = fopen(filename,"a");
    if (posfile == NULL){
        printf("Error with opening \"%s\" file",filename);
    }
    else{
        printf("\nPrinting the xy positions to %s\n",filename);
        /* rotate the chain before printing, such that the chain is oriented along the x-axis. This is better for the analysis tool of Simon*/
        /*  translate all particles such that v0=(0,0,0) 
            take the angle between  the x-value of particle 0 along x-axis; the y-value should be around 0.
            and the vector between particle 0 and N-1
            use this angle to do the rotation via the quaternion in z-axis QuaternionZaxis(degree) */
       
        // calculate ther rotation matrix R
        // /choose particle 0 as reference
        // shift the reference frame. walk through the chain
        r_ref=psl->pts[0].r;
        v[0]=nulvec;
        // vprint(v[0]);
        for(ipart=1; ipart<sys.npart; ipart++) {
            //the vector difference between particle i and j
            vector_minus(psl->pts[ipart].r,psl->pts[ipart-1].r,dr);
            pbc(dr,sys.boxl);
            vector_add(v[ipart-1],dr,v[ipart]);
            // vprint(v[ipart]);
        }
        
        // calculate the rotation matrix
        dv.x=v[sys.npart-1].x;
        dv.y=v[sys.npart-1].y;
        dv.z=0.;

        scalar_divide(dv,sqrt(vector_inp(dv,dv)),dv_norm); 

        /*build in if angle>180(== costheta<0, the rotation becomes is 360-angle)*/
        costheta=vector_inp(dv_norm,xaxis);
        if (dv_norm.y>0.){
            //rotations are done tegen de klok in
            rotangle=360.-acos(costheta)*180./PI;

        }
        else{
            rotangle=acos(costheta)*180./PI;
        }
        // gprint(rotangle);
        q=QuaternionZaxis(rotangle);
        R = getrotmatrix(q);
        
        //perform the rotarion
        for(ipart=0; ipart<sys.npart; ipart++) {
            matrix_x_vector(R, v[ipart], v_new[ipart]);


            //print to file        
            if (sys.sim_type==0){
                fprintf(posfile,"%d,%.5lf,%.5lf,0.0,0.0,%d,%.5lf\n", linecount_xypos,v_new[ipart].x, v_new[ipart].y,ipart ,cycle_mc);
                linecount_xypos++;
            }   
            else if(sys.sim_type==2){
                fprintf(posfile,"%d,%.5lf,%.5lf,0.0,0.0,%d,%.5lf\n", linecount_xypos ,v_new[ipart].x, v_new[ipart].y,ipart ,sys.c_time);
                linecount_xypos++;
            }         
        }
    } 
    
    cycle_mc+=0.1;
    fclose(posfile);
    return;

}


