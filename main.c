#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <dirent.h> 
// #include <omp.h> // for paralelization


// clean these global variables such that they are not global but dynamic according to the type of simulation you do
// create these variables in heap

//truly global 
System sys;
vector nulvec={0};
SaccentInt s_int; //switch function variables
Slice *slice,*psl_old;
Cluster cluster;
Chain chain;


//analysis types
// Analysis_ptrs analysis_ptrs;
//Average rdf_average[RDFBINS], globalCN_average;
int crossing_check,s_histogram[MAXBIN], r_S_2dhist[MAXBIN][MAXBIN],xy_plane_patchangle[ANGLEBINS];
Statistics bond_length,msd,Sfraction;

//MC types
Mcacc  trans,        fintrans,        rot,        finrot;
Mcacc  trans_cluster,fintrans_cluster,rot_cluster,finrot_cluster ;
Mcacc  trans_mono,   fintrans_mono,   rot_mono,   finrot_mono;

//BMD
Langevin langevin;
Slice *init_state,*fwtrial,*bwtrial,*trial;


//Emilio
Constraint cconst;
int crossing_check;
Trajectory segment;


//advanced sampling: TPS / TIS / replica exchange
State state[MAXSTATES];
Replica *replica[MAXREPLICA];
Path     path;

//nearest neighbor variables
List list;








int main(int argc, char **argv) {

    int icycle,jcycle;
    double energy;
    FILE *emptyfile, *rdffile;
    time_t t0 = time(0);
    
    setup_simulation();
    
    if(sys.empty_files==1){
        emptyfiles();
    }
    
    plotpotential();
    
    // if(sys.graphics==1) {
    //     Init_Graphics(argc,argv);
    // }

    for(icycle=0; icycle<sys.ncycle1; icycle++) {
        time_t block2_0 = time(0);
        for(jcycle=0; jcycle<sys.ncycle2; jcycle++) {
            if(sys.sim_type==0) {
                bmdcycle();
            }
            else if(sys.sim_type==2) {
                mccycle();
            }
            else if(sys.sim_type==4){
                //what is crossing check? what does it do
                crossing_check = 1;
                ffscycle(icycle, jcycle);
                // if(crossing_check==0)
                    // break;
            }
            else {
                error("choose only MC = sys.sim_type==2 or bmd = sys.sim_type==0");
                // associationcycle();
            }
        }
        printf("\nBLOCK %d\n",icycle);
        if(sys.sim_type != 4){
            terminate_block();
        }
        time_t block2_1 = time(0);
        double datetime_diff2 = difftime(block2_1, block2_0);
        printf("This block took %lf  [sec]\n", datetime_diff2);
    }
    // if(sys.planalysis==1){
    //     samplepl(1);
    // }
    dprint(icycle);
    dprint(jcycle);

    finalstat();
    time_t t1 = time(0);
    double datetime_diff = difftime(t1, t0);/* in seconds*/
    printf("total time [min] %lf\n", datetime_diff/60.); 
    return 0;
}



void bmdcycle() {
    /*sys.sim_type=0*/
    int istep, update=0, ipart;
    double dr2;

    for(istep=0;istep<10;istep++){
        propagate_bd(&slice[0]);
    }
    // forcecheck_nn(&slice[0]);
    /* is this needed?
    if(sys.npart==1) {
        //angle distribution and translational/rotational diffusion constant calculation
        sampleangles(0); 
        samplediffusion(0);
    }

    if(sys.npart==2) {
        samplefreeenergy(0);
        sampleenergy(0);
        samplemdrate(0);
    }*/

    return;
}

void mccycle() {
    /*sys.sim_type=2*/

    int i,ipart,istep, nclusters,nbonds_before ;
    double which, dr2,r; 
    int nnlupdate=0,j,k;
    vector dr;
 
  
    /* first identify the clusters; the number of clusters will not change in a clustermove, only the positions of them*/
    /* we first make the clusters with cluster_analysis(psl), followed by identifying the lengths (with clustersize_identification), and linked-list 
    such that the clusters are available via head[chainid] and they are sequenced */

    // printf("nbonds_after %d \n", nbonds_after);
    // printf("start MC\n");
    if(sys.cluster_MC==1){

        Slice *copyslice = malloc(sizeof(Slice));

        for(istep=0;istep<10;istep++){
            /* choose single or cluster movement*/
            which = RandomNumber();

            /* cluster moves, where nbonds and nclusters will not change*/
            /* the number of cluster should stay constant during clustermoves ivm detailed balance*/
            if(which<0.05){

                /*save old positions to evaluate after clusterpropagation if an update of the nnl is neccessary*/
                if(sys.nearest_neighbor==1){
                    memcpy(copyslice, &slice[0], sizeof(Slice));
                }
                nbonds_before= slice[0].nbonds;

                /* expensive function, do this function only when the number of bonds has changed. */
                if (cluster.update==1){
                    // printf("****performing cluster analysis*****\n");
                    slice[0].nclusters = cluster_analysis(&slice[0]); 
                    clustersize_identification(&slice[0]);
                    // chain_linkedlist(&slice[0]); /* linked listst are neccessary for check the bond lengths in the rotation cluster move*/
                    cluster.update=0;
                }
            
                /* *** the cluster move*** */
                cluster_propagate_mc(&slice[0]);

                // if(sys.nearest_neighbor==1){
                //     for(ipart=0; ipart<sys.npart; ipart++){
                //         /* add the dr difference between new-old, dr is not fixed */
                //         vector_minus(slice[0].pts[ipart].r, copyslice->pts[ipart].r, dr);
                //         vector_add(slice[0].pts[ipart].dr,dr,slice[0].pts[ipart].dr);
                //         dr2 = vector_inp(slice[0].pts[ipart].dr , slice[0].pts[ipart].dr);

                //         if ((dr2 > list.cutoff_buffer2)) {  nnlupdate=1;   }
                //     }
                // }
                // printf("a cluster move is performed, the energy is         %lf\n", slice[0].energy);
                energy_divergence_check(&slice[0],"cluster move");

                /*check nclusters is constant, it is constant if nbonds is constant*/
                if(nbonds_before != slice[0].nbonds){
                    int totbonds=0;
                    for(ipart=0;ipart<sys.npart;ipart++){
                        totbonds +=slice[0].pts[ipart].nbonds;
                    }
                    printf("tot bonds results from sum over all slice[0].pts[ipart].nbonds= %d \n", totbonds );
                    printf("nbonds before : %d  after : %d \n",nbonds_before, slice[0].nbonds);
                    error("ERROR:: nclusters not constant during clustermoves");
                }

                /*checks if any particle exeeds the maximum number of bonds (more than one per site)*/
                // printf("after cluster moves\n");
                check_maxbonds(&slice[0]);
            }
            
            else{
                /* nearest neighbor list cannot be used in cluster monte carlo
                if(nnlupdate==1){
                    update_nnlist(&slice[0]);
                    // printf("Updated the neighborlist\n");
                    nnlupdate=0;
                } */

                propagate_mc(&slice[0]);   
                // printf("a single particle move is performed, the energy is %lf\n", slice[0].energy);
                energy_divergence_check(&slice[0],"single particle moves");
                /*checks if any particle exeeds the maximum number of bonds (more than one per site)*/
                check_maxbonds(&slice[0]);
 
            }
            
        }
        free(copyslice);
    }
    else{
        for(istep=0;istep<1;istep++){
            propagate_mc(&slice[0]);

            /*checks if any particle exeeds the maximum number of bonds (more than one per site)*/
            // printf("after single particle moves\n");
            check_maxbonds(&slice[0]);
            // update_nnlist(&slice[0]);
            conf_output(&slice[0]);
            energy_divergence_check(&slice[0],"single particle moves");
                
        }   
    }

    
  
    return;
}

void ffscycle(int init_config, int interface) {
    /*sys.sim_type=4*/
    int trial_counter, success=0,timestep, k;
    double lambda0, lambdai, lambdaj,lambda_prev;

    printf("\nStarting FFS cycle:\n");
    lambda0 = state[0].lambda[0];
    // lambda_prev = state[0].lambda[interface-1];
    lambdai = state[0].lambda[interface];
    lambdaj = state[0].lambda[interface+1];

    // reached state B
    if(interface==state[0].nrep-1){
        double rb_weight=1, ratio;
        for(trial_counter=1; trial_counter<state[0].nrep; trial_counter++){
            rb_weight *= (double) state[trial_counter].nsucc/sys.K;
        }

        ratio = rb_weight / state[0].weight;
        printf("        Wn / Wo : %lf / %lf =  %lf\n", rb_weight, state[0].weight, ratio);

        if(ratio >= 1){
            state[0].weight = rb_weight;
            save_ffscycle(init_config);
            state[0].freq[init_config] += 1;
            state[0].current_path = init_config;
            printf("    -> new path accepted! Writing to FFS_paths.dat...");
        }
        else{
            double r = RandomNumber();
            if(r < ratio){
                state[0].weight = rb_weight;
                save_ffscycle(init_config);
                state[0].freq[init_config] += 1;
                state[0].current_path = init_config;
                printf("    -> new path accepted! Writing to FFS_paths.dat...");
            }
            else{
                state[0].freq[state[0].current_path] += 1;
                printf("    -> new path rejected!");
            }
        }
    }
    // run K trials from lambdai
    else{
        printf("Running trials from interface %d (lambdai : %lf)...\n", interface+1, lambdai);
        if (interface>2){
            lambda_prev = state[0].lambda[interface-2];
        }else{
            lambda_prev = state[0].lambda[0];
        }
        printf("Until crossing lambda_prev : %lf or lambdaj : %lf\n", lambda_prev, lambdaj);
        
        /*reads the coordinates+quaternions of the last frame/segment of interface i
        and puts them in psl*/
        read_coordinates_of_interface(init_state, interface); /*used to be called ffs_init by emilio*/
        Segment *segment_list[sys.K];

        /* paralelized part of code with omp*/
        #pragma omp parallel for private(trial_counter, trial)
        for(trial_counter=0; trial_counter<sys.K; trial_counter++){
            int ipart, time;
            double lambda_now, lambda_prev;
            trial = (Slice *)calloc(1,sizeof(Slice));
            Segment *start = (Segment *)calloc(1,sizeof(Segment));
            Segment *current;
            printf("running trial %d...", trial_counter);

            // set up initial configuration
            for(ipart=0;ipart<sys.npart;ipart++){
                trial->pts[ipart] = init_state->pts[ipart];
                start->pts[ipart] = init_state->pts[ipart];
            }

            time=0;
            lambda_now = op_calc(trial); // initial OP
            //lambda_now = lambdai;
            start->lambda=lambda_now;
            current = start;
            // run until crossing either lambda0 or lambdaj
            do{
                current->next = (Segment *)calloc(1,sizeof(Segment));
                current = current->next;
                lambda_prev = lambda_now;

                propagate_bd(trial);
                //lambda_now += 0.05*RandomGaussianNumber();
                lambda_now = op_calc(trial);
                time++;
                current->t=langevin.ninter*langevin.timestep*time;
                current->lambda=lambda_now;
                //printf("lambda_now : %lf (thread %d)\n", lambda_now, omp_get_thread_num());

                for(ipart=0;ipart<sys.npart;ipart++){
                    current->pts[ipart] = trial->pts[ipart];
                }

            } while(lambda_now > lambda0 && lambda_now < lambdaj);

            // if trial crossed next interface, store path
            if(lambda_now > lambdaj){
                success++;
                segment_list[success-1] = start;
                // printf("    -> trial %d succesfully crossed interface %d!\n", trial_counter, interface+2);
            }
        }

        int traj;
        Segment *current, *next;
        if(success==0){
            printf("    -> no successful trials\n");
            crossing_check = 0;

            for(traj=0; traj<success; traj++){ //free allocated memory
                current = segment_list[traj];
                while(current->next != NULL){
                    next = current->next;
                    free(current);
                    current = next;
                }
            }
        }
        else if(success>0){
            state[interface+1].nsucc=success;
            printf("    -> %d successful trials\n", success);
            k = (int)floor(RandomNumber()*success);
            state[interface+1].path = segment_list[k];
            printf("    -> Chose trial %d out of %d initial configurations\n", k, state[interface+1].nsucc);

            for(traj=0; traj<success; traj++){ //free allocated memory
                if(traj==k)
                    continue;

                current = segment_list[traj];
                while(current->next != NULL){
                    next = current->next;
                    free(current);
                    current = next;
                }
            }

        }
    }
}

// void mainloop_for_graphics() {

//     static int icycle,jcycle,initial,wait=0;

//     if(initial==0) {
//         initial=1;
//         icycle=jcycle=1;
//         if(sys.bond_op){  wait=1;}
//     }
//     if(sys.bond_op){    
//         clustersize_identification(&slice[0]);
//         clustersize_freq_update(&slice[0]);
//         nematic_order_parameter(&slice[0]);
//     }
//     if(sys.sim_type==0) {
//         bmdcycle();
//     }
//     if(sys.sim_type==2) {
//         mccycle();
//     }
//     // printf("mainloop_for_graphics\n" );
    

//     if(jcycle%sys.ncycle2==0) {
//         printf("\nBLOCK %d\n",icycle);
//         terminate_block();
        
//         if((sys.graphics==1 )&& (sys.snapshot==1)){
//             if (wait>1){
//                 saveScreenshotToFile( 700,800);
//                 exit(1);
//             }
//             wait+=1;
//         }
//         jcycle=0;
//         if(icycle%sys.ncycle1==0) {
//             printf("Finished Simulation\n");
//             finalstat();
//             exit(1);
//         }
//         icycle++;
//     }

//     jcycle++;

//     return;
// }

void emptyfiles(){
    struct dirent *de;
    int status;

    DIR *dr = opendir(".");
    if (dr == NULL)  // opendir returns NULL if couldn't open directory 
    { 
        printf("Could not open current directory. No files are deleted. continue" ); 
        return ; 
    } 
  
    // Refer http://pubs.opengroup.org/onlinepubs/7990989775/xsh/readdir.html 
    // for readdir() 
    while ((de = readdir(dr)) != NULL){
        // printf("%s\n", de->d_name); /* prints all de files in the current directory*/
        if (strstr(de->d_name, ".dat") != NULL){
            status = remove(de->d_name);

            if (status == 0){
                printf("deleted successfully:     %s\n", de->d_name);
            }
            else{
                printf("Unable to delete the file\n");
            }
        }
    }
  
    closedir(dr);     
    return ;
}

void error(char *msg){
  printf("error: %s\n",msg);
  exit(0);
}























