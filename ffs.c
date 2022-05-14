
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"

void conf_input_ffs(Slice *);

void read_coordinates_of_interface(Slice *psl, int interface){
    int ipart;
    double lambda_now;
    
    if(interface==0){
        conf_input_ffs(psl);
        read_particletypes(psl);
        update_patch_vectors(psl);
    }
    else{
        Segment *current;
        current = state[interface].path;

        while(current->next != NULL){
            current = current->next;
        }
        
        memcpy(psl->pts, current->pts, sizeof(Pts[NPART]));
    }
    printf("initial energy = %lf\n",total_energy(psl) );
    lambda_now = op_calc(psl);
    printf("initial lambda : %lf\n", lambda_now);
}

// saves FFS path to output file
void save_ffscycle(int init_config){
    int interface, ipart;
    vector posi;
    quaternion quati;
    Pts *psl;
    FILE *pp,*cp;

    // setup correct path file
    char str[25] = "ffs_output/path";
    char crossp[27] = "ffs_output/crossp";
    char num[4];
    sprintf(num, "%d", init_config);
    strcat(str, num);
    strcat(crossp, num);
    strcat(str, ".dat");
    strcat(crossp, ".dat");

    if((pp=fopen(str,"w"))==NULL || (cp=fopen(crossp, "w"))==NULL)
        error("Warning: can not open ffs files\n");

    // store path[init_config] to file
    fprintf(pp, "PATH NUMBER %d | | | | | | |\n", init_config);

    //why is interface_max - sys.ncycle2?
    for(interface=1; interface<sys.ncycle2; interface++){
        Segment *current, *next;
        current = state[interface].path;
        fprintf(cp, "%d %lf %lf\n", interface, state[0].lambda[interface-1], ((double) state[interface-1].nsucc) / ((double) sys.K));
        while(current->next != NULL){
            for(ipart=0; ipart<sys.npart; ipart++){
                posi=current->pts[ipart].r;
                quati=current->pts[ipart].q;
                fprintf(pp,"%d %lf %lf %lf %lf %lf %lf %lf %.10f %.10f\n",ipart, posi.x,posi.y,posi.z,quati.q0,quati.q1,quati.q2,quati.q3,current->t, current->lambda);
            }
            next = current->next;
            free(current);
            current = next;
        }
    }
    fprintf(cp, "%d %lf %lf\n", interface, state[0].lambda[sys.ncycle2-1], ((double) state[sys.ncycle2-1].nsucc) / sys.K);
    fclose(pp);
    fclose(cp);
}

void conf_input_ffs(Slice *psl){
    int ipart, isite,npart,nsites;
    vector boxl;
    FILE *fp;

    static int ffs_setup;
    char str[50]="conf";
    char num[4];
  
    sprintf(num, "%d", ffs_setup);
    strcat(str, num);
    strcat(str, ".inp");
    ffs_setup++;

    if((fp=fopen(str,"r"))==NULL)  {;
        printf("Warning: can not open %s\n",str);
    }
    else {
        fscanf(fp,"%d %lf %lf %lf\n", &npart, &boxl.x,  &boxl.y,  &boxl.z);
        for(ipart=0; ipart<npart; ipart++) {
            fscanf(fp,"%lf %lf %lf ", &psl->pts[ipart].r.x, &psl->pts[ipart].r.y, &psl->pts[ipart].r.z);
            fscanf(fp,"%lf %lf %lf %lf\n", &psl->pts[ipart].q.q0, &psl->pts[ipart].q.q1, &psl->pts[ipart].q.q2, &psl->pts[ipart].q.q3);
        }
    }
    fclose(fp);

    if(npart!=sys.npart) {
        dprint(sys.npart);
        dprint(npart);
        error("Error: number of particles in system not same as in conf.inp\n");
        
    }
    if((boxl.x!=sys.boxl.x) | (boxl.y!=sys.boxl.y) | (boxl.z!=sys.boxl.z)) {
        vprint(boxl);
        vprint(sys.boxl);
        error("Error: boxlength in system not same as in conf.inp\n");
    }
    return;

}
