#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"

void convert_bulk_to_patch_variables(vector *,vector *,vector *, int , int , double ,double *,double *,double *);
void calc_angles( vector *,vector *,vector *, double *,double *,double *);


double total_energy(Slice *psl) { 
    /*  this function calculates the total energy (= gravitational energy + bond energies) of the system.
    A loop over all particles is done.
    */
    Pts *psi,*psj;
    vector rij,rnorm;
    int ipart,jpart,tot_bonds=0;
    double Eparticle,Etot=0;
    
  
    
    if (sys.nearest_neighbor){
        // perform energy calculation with nearest neighbor list
        // printf("with nn \n");
        Etot=total_energy_neighborlist(psl);
    }
    else{
        for( ipart=0; ipart<sys.npart; ipart++ )  {
            //put all number of bonds to zero
            psl->pts[ipart].nbonds =0;
        }
        //calculate patchy Simons potential interaction between all particles
        for( ipart=0; ipart<sys.npart; ipart++ )  {
           
            jpart=ipart+1;
            Eparticle=particle_energy(psl, ipart, jpart);  // particle energy 
            Etot+=Eparticle;
            tot_bonds+=psl->pts[ipart].nbonds;
        }
        // printf("particle energies done\n");
        psl->nbonds=tot_bonds/2;
    }

    // gprint(Etot);
    // error("quit");

    return Etot;
}

double harmonic_oscillator(double k, double x){
    return 1./2*k*x*x;
}

double wall_interaction_ipart(Slice *psl, int ipart){
    /* if the wall has an interaction with the patches
    there is a harmonic oscillator that governs the interaction strength with its minimum below 0*/
    Pts *psi;
    psi=&psl->pts[ipart];
    int isite;
    vector z_axis=nulvec,x_axis=nulvec,y_axis=nulvec,p1;
    double Ewall_tmp=0,Ewall,z,angle_max=sys.particletype[psi->ptype].cosdelta;
    double zcut=sys.particletype[psi->ptype].zcut+sys.particletype[psi->ptype].radius*1.1,angle,cosangle,S,E_height,pi180=PI/180.;
    z_axis.z=1;
    x_axis.x=1;
    y_axis.y=1;

   
    // only if the particle is in the vicinity of the wall: 
    if (psi->r.z<=zcut & psi->r.z>0){
        // find thenpatch that is pointing downwards. take the innerproduct with patch vector if 
        for(isite=0;isite<sys.particletype[psi->ptype].nsites;isite++){
            p1=psi->patchvector[isite];
            if (p1.z>0){
                continue;
            }
            cosangle=-vector_inp(p1,z_axis); 
            
            if (cosangle<angle_max){
                continue;
            }
            // printf("wall angle ipart %d = %lf\n", ipart,angle);

            S= Saccent( cosangle, psi->ptype);
            
            // E_height=gravitational_energy_ipart(psl, ipart); // gravitational energy is zero at the wall
            z=fabs(p1.z-zcut);
            E_height=-1.*harmonic_oscillator(sys.wall_int,z);  
            // E_height=0;
            Ewall=E_height*S;

            // gprint(Ewall);
            if (Ewall<Ewall_tmp){
                Ewall_tmp=Ewall;

                printf("*** ipart %d with isite %d***\n",ipart,isite);
                gprint(Ewall);
                gprint(S);
                angle=acos(cosangle);
                gprint(angle/pi180);
                printf("\n");


                // // vprint(p1);
                // // vprint(z_axis);
                // // gprint(sqrt(vector_inp(p1,p1)));

                // // printf("angle x axis\n");
                // // cosangle=vector_inp(p1,x_axis); 
                // // gprint(cosangle);
                // // gprint(angle=acos(cosangle)/pi180);
                // // printf("angle y axis\n");
                // // cosangle=vector_inp(p1,y_axis); 
                // // gprint(cosangle);
                // // gprint(angle=acos(cosangle)/pi180);


            }
        }
    }
    
    // gprint(Ewall_tmp);
    return Ewall_tmp ;//-Ewall_max;
}

double particle_energy(Slice *psl, int ipart, int jpart) { 
    /*  this function calculates the particle energy of particle ipart.
    this includes gravity if on.
    It starts counting from particle jpart.
    Use jpart=0 in MC code, while in total_energy you use jpart=ipart+1;

    For MC you need to count bonds, 
    for the calculation of a single particle you need to known if nr bonds changed (fot the cluster update)
    Only info about that specific particle is needed

    If total energy is calculated, all the particles should contain info about bonden (so if p1-p2 bond, than both should know that)
    */

    Pts *psi, *psk;
    vector rik;
    
    double Ebond,Erep, Ebond_tot=0, Etot=0, Erep_tot=0, Egrav_tot=0, Ewall=0,s_patch;
    double s_bb,s_bp1,s_bp2,s_pp,cosphi_i,cosphi_j,cosphi_ij;
    double r,r2,s,r_radii,Derj_pref;
    double Eatr,S;
    int kpart;

    psi = &psl->pts[ipart];
     // printf("the pointer of particle %d %p\n", ipart, &psi);
    // printf("partilce energy %d\n",ipart);
    if (jpart==0){
        //if you start from MC (jpart==0) then set nbonds of ipart to zero first
        //if you start from total_energy it is already set to zero in that function
        psl->pts[ipart].nbonds =0;
    }   
    // if (ipart>1000){printf("bonds=0\n");}
    /* gravity*/
    if(sys.gravity>0){
     //   for( ipart=0; ipart<sys.npart; ipart++ )  {
        Egrav_tot+=gravitational_energy_ipart(psl, ipart);
        if (sys.wall_int!=0){
            Ewall+=wall_interaction_ipart(psl,ipart);
        }
    }
    // printf("partilce energy gravity %f\n",Egrav_tot);
    /*kpart is de teller, bereken de energy van ipart, start counting vanaf jpart*/
    for( kpart=jpart; kpart<sys.npart; kpart++ ) {
  
        if(ipart==kpart){
            continue;
        }
        psk = &psl->pts[kpart];

        vector_minus(psi->r,psk->r,rik);
        pbc(rik,sys.boxl); /*do perform pbc when using no nnl*/


        r2 = vector_inp(rik,rik);
        r=sqrt(r2);
        r_radii=sys.particletype[psi->ptype].radius+sys.particletype[psk->ptype].radius;
        s=r-r_radii;
        // if (ipart>1000){gprint(s);}
        if(s<sys.rcutoff && s>0.0) { 
            

            /*  calculate here the s, s_pp,s_bp and the angles that are use for the switch function*/
            orientation_parameters(psl,ipart,kpart, &s_bb,&s_bp1,&s_bp2,&s_pp, &cosphi_i, &cosphi_j, &cosphi_ij);
            if (fabs(s-s_bb)>1e-10){
                gprint(s);
                gprint(s_bb);
                error("something went wrong in the calculation of the distance between big particles");
            }
            
            //the derjaguin scaling of this combintaion of particle types
            Derj_pref=sys.particletype[psi->ptype].Derj_scaling[psk->ptype];

            // patch attraction
            Eatr=potential_attractive_energy_sdist(s_pp)*Derj_pref;
            double cosphi_i_treshold=sys.particletype[psi->ptype].cosphi;
            double cosphi_j_treshold=sys.particletype[psk->ptype].cosphi;
            if ((cosphi_i<cosphi_i_treshold) || (cosphi_j<cosphi_j_treshold)){
                S=0;
            }
            else{
                S=1; // S_value(psl,cosphi_i,cosphi_j, nulvec,nulvec, nulvec, ipart, jpart); // I can use nulvec here, because I take S90 anyway
            }
            Ebond=Eatr*S;

            //repulsions; I only include patch-patch  and bulk-bulk repulsion for now
            double Erep_pp,Erep_bb;
            Erep_pp=potential_repulsive_energy_sdist(s_pp)*Derj_pref; //patch-patch repulsion
            Erep_bb=potential_repulsive_energy_sdist(s_bb); //bulk-bulk repulsion
            
            Erep=Erep_pp+Erep_bb;
            Erep_tot += Erep; 

            // printf("The distance s=%lf , s_patch=%.5lf, Ebond=%.5lf, Erep=%.5lf,\n", s, s_patch,Ebond,Erep);

            if(Ebond<sys.bond_cutoffE){
                psi->bonds[psi->nbonds]=kpart;
                psi->bond_energy[psi->nbonds]= Erep+Ebond;
                psi->bond_distance[psi->nbonds]=s;
                psi->nbonds++; 
                // printf("%p  of nbonds %d\n", &psi->nbonds,psi->nbonds);
                // if (ipart>1000){printf(", and now has total %d bonds\n", psi->nbonds);}
                if( (jpart>0 )&& (kpart>ipart)){
                    /*in case of calculating the total energy, also add the bond to the other particle*/
                    psk->bonds[psk->nbonds]=ipart;
                    // printf("       and ipart %d is added to kpart %d who had already %d bonds \n",ipart,kpart,psk->nbonds);
                    psk->bond_energy[psk->nbonds]= Erep+Ebond;
                    psk->bond_distance[psk->nbonds]= s;
                    psk->nbonds++;
                }
            }
            //add it to the potential energy of the patch
            Ebond_tot += Ebond ;
        }
        else if(s<=0.0){
            // printf("overlapping particles\n");
            Ebond_tot +=1e6;
        }
    }


    Etot = Ebond_tot + Erep_tot + Egrav_tot + Ewall;
    return Etot;
}


void convert_bulk_to_patch_variables(vector *u1,vector *p1,vector *p2, int ptype1, int ptype2, double r_bulk,double *r_patch_ptr,double *cosphi_i,double *cosphi_j){
    double r_patch2,r_patch,d2,d1,theta_ij,theta_i,theta_j,rad_to_degree=180./PI;
    vector pi,pj, pi_norm=*p1,pj_norm=*p2,r_pp,r=*u1;
    scalar_times(r,r_bulk,r);

    d1=sys.particletype[ptype1].d;
    d2=sys.particletype[ptype2].d;

    calc_angles( u1,p1, p2,&theta_i,&theta_j,&theta_ij);
   
    scalar_times(pi_norm,d1,pi);
    scalar_times(pj_norm,d2,pj);
    vector_minus(pj,r,pj);

    //vector of r_patch and its length
    vector_minus(pi,pj,r_pp) //r_pp points form j->i
    r_patch2= vector_inp(r_pp,r_pp);
    r_patch=sqrt(r_patch2);
    *r_patch_ptr=r_patch;

    //angle phi of particle i  & j;
    *cosphi_i=-vector_inp(pi_norm,r_pp)/r_patch;
    *cosphi_j=vector_inp(pj_norm,r_pp)/r_patch;

    // printf("the angles of the bulk are theta_i=%.5lf , theta_j=%.5lf, theta_ij=%.5lf, and r_patch=%.5lf \n",theta_i,theta_j,theta_ij,r_patch);

    return;
}

void calc_angles( vector *u1,vector *p1,vector *p2, double *theta_i,double *theta_j,double *theta_ij){
    /* calclulates the angles between a interparticle vector u1 and two patch vectors p1 and p2*/
    double costheta_i,costheta_j,costheta_ij,rad_to_degree=180./PI;
    vector e1,e2,e3,proj_u1_p1,proj_u1_p2;
    vector u11=*u1,p11=*p1,p22=*p2,u2,u3;

    costheta_i=-vector_inp(u11,p11); 
    costheta_j=vector_inp(u11,p22); 
    // vprint(u11);
    // vprint(p11);
    // vprint(p22);

    if (costheta_i>=1 || costheta_j>=1){
        *theta_i=(costheta_i<=1.)? acos(costheta_i)*180./PI:0;
        *theta_j=(costheta_j<=1.)? acos(costheta_j)*180./PI:0;
        *theta_ij=0.;
    } 
    else{  


    // projection p1 onto u1, and create (unit vector) e2 the perpendicular projection of p1 
    scalar_times(u11,-costheta_i,proj_u1_p1); 
    vector_minus(p11,proj_u1_p1,u2);
    scalar_divide(u2,sqrt(vector_inp(u2,u2)),e1);
    // projection p2 onto u11, and create (unit vector) e3 the perpendicular projection of p2 
    scalar_times(u11,costheta_j,proj_u1_p2); 
    vector_minus(p22,proj_u1_p2,u3);
    scalar_divide(u3,sqrt(vector_inp(u3,u3)),e2);

    /*calc angle3={0,180}, you need it in angles, because the interpolation needs angles*/
    costheta_ij=vector_inp(e1,e2);

    *theta_i=acos(costheta_i);
    *theta_j=acos(costheta_j);
    *theta_ij=acos(costheta_ij);
    }

    return;
}

   
vector find_directing_patchvector(Slice *psl, int ipart,vector u1){
    int isite; 
    double tracker_i=-1,cositheta;
    vector p1_select=nulvec;

    

    for( isite=0; isite<sys.particletype[psl->pts[ipart].ptype].nsites; isite++ ) {
        vector p1=psl->pts[ipart].patchvector[isite];
        cositheta=vector_inp(u1,p1); 
        
        if(cositheta<tracker_i) {
            continue;
        }
        tracker_i=cositheta;
        p1_select=p1;
    }

    return p1_select;
}

void orientation_parameters(Slice *psl, int ipart, int jpart, double *s_bb,double *s_bp1,double *s_bp2,double *s_pp, double *cosphi_i, double *cosphi_j, double *cosphi_ij){
    /* this function calculates the bulk0bulk distance s, the bulk-patch distance s_bp, the patch-pathc distance s_pp,
    the theta and phi angles */

    Pts *psi,*psj;
    
    double r,r2;
    double r_patch,s_patch;
    vector u1,min_u1,rij, p1,p2;

    psi=&psl->pts[ipart];
    psj=&psl->pts[jpart];

    vector_minus(psi->r,psj->r,rij);
    pbc(rij,sys.boxl);
    r2 = vector_inp(rij,rij);
    r=sqrt(r2);
    scalar_divide(rij,r,u1);  

    //find smallest combination of the patch vector angles
    scalar_times(u1,-1,min_u1); //make here -u1, because you need to flip the u1 vector, to find p1 (and patch angle) correctly.
    p1= find_directing_patchvector(psl, ipart,  min_u1);
    p2= find_directing_patchvector(psl, jpart,  u1);
    convert_bulk_to_patch_variables(&u1,&p1,&p2,  psi->ptype,  psj->ptype, r,&r_patch,cosphi_i,cosphi_j);


    *s_pp=r_patch-sys.particletype[psi->ptype].small_radius-sys.particletype[psj->ptype].small_radius;     
    *s_bb=(r-sys.particletype[psi->ptype].radius-sys.particletype[psj->ptype].radius);
    // *s_pp=sdist_patch;
    *s_bp1=10.; //I put it to 10 , because I don't know if im going to use it. 
    *s_bp2=10.; //I put it to 10 , because I don't know if im going to use it. 
    *cosphi_ij=0; // I put it here to 90 degrees, I dont think I will need it in the future as we are using S90 in swichfunction


    return;
}

double potential_attractive_bond_energy(Slice *psl, int ipart, int jpart){
    /*this fucntion calculates Vc*S and returns this energy
    you already know that s is small enough < rccutoff
    */
   
    double s_bb,s_bp1,s_bp2,s_pp,cosphi_i,cosphi_j,cosphi_ij;
    double Ebond,S,Eatr;
    double cosphi_i_treshold=sys.particletype[psl->pts[ipart].ptype].cosphi,cosphi_j_treshold=sys.particletype[psl->pts[jpart].ptype].cosphi;

    //energy and S based on s_patches distance, and phi's
    // step 1: convert s_bulk and theta's to s_patches and phi's

    /*  calculate here the s, s_pp,s_bp and the angles that are use for the switch function*/
    orientation_parameters(psl,ipart,jpart, &s_bb,&s_bp1,&s_bp2,&s_pp, &cosphi_i, &cosphi_j, &cosphi_ij);
    
    int ptype1=psl->pts[ipart].ptype,ptype2=psl->pts[jpart].ptype;
    // patch attraction
    Eatr=potential_attractive_energy_sdist(s_pp)*sys.particletype[ptype1].Derj_scaling[ptype2];;
    
    if ((cosphi_i<cosphi_i_treshold) || (cosphi_j<cosphi_i_treshold)){
        S=0;
    }
    else{
        S=1; // S_value(psl,cosphi_i,cosphi_j, nulvec,nulvec, nulvec, ipart, jpart); // I can use nulvec here, because I take S90 anyway
    }
    Ebond=Eatr*S;


    if(Ebond!=0.){
        if (Ebond/Ebond!=1 ){
            gprint(S);
            gprint(cosphi_i);
            gprint(cosphi_j);
            gprint(Eatr);
            gprint(s_bb);

            error("Ebond is nan");
        }
    }
    
    return Ebond;
}

void energycheck_nn(Slice *psl, double en_wo, double en){
    int i;

    // printf("end\n\n");
    if(en-en_wo>=1e-6){
        printf("WARNING energy %lf not equal to energy with neighborlist %lf ", en_wo, en );
        printf("with a difference of     %lf\n", en_wo-en);

        // en=potential_energy(psl);
        if(fabs(en-en_wo)>=0.9  ){
            error("The difference is larger than 1.0 kT; quitting using the nearestneighbor list\n");
        }
    }
    return;
}

double gravitational_energy_ipart(Slice *psl, int ipart){
    double gravitational_energy=0.0, zinv,zinv2,zinv6,zinv12,psiz;
    double fg,zcut,b_zc,radius;
   
    Pts *psi;
    psi = &psl->pts[ipart];

    
    zcut=sys.particletype[psi->ptype].zcut;
    fg=sys.particletype[psi->ptype].fg;
    b_zc=sys.particletype[psi->ptype].b_zc;
    
    radius=sys.particletype[psi->ptype].radius;
    
    psiz =psi->r.z-radius;


    if(psiz<0){
        psiz+=sys.boxl.z;
    }
    if(psiz>=zcut){ 
        gravitational_energy = fg*psiz - b_zc; /*V_g = Fg*z - b*/
    }
    else{
        zinv = 1./psiz;
        zinv2 = zinv*zinv;
        zinv6 = zinv2*zinv2*zinv2;
        zinv12 = zinv6*zinv6;
        gravitational_energy = 4.*sys.epsilongravLJ*(zinv12-zinv6+1./4.);
    }
    // gprint(zinv12);
    // gprint(zinv6);
        

    return gravitational_energy;


}

double total_energy_neighborlist(Slice *psl){
     
    /* the potential energy using the neighbor list
      the neighbor list is a structure that works with 26 images */
 
    Pts *psi,*psj;
    vector rij,ri;
    double r,r2,s,r_radii,s_patch;
    double Erep,Erep_tot=0,Ebond,Ebond_tot=0 ,Egrav_tot=0,Etot;

    int i,k,a,j,count,*b,tot_bonds=0,ipart;

    // printf("*** with negihborlist***\n");

    /* gravity*/
    
    for( ipart=0; ipart<sys.npart; ipart++ )  {
        psl->pts[ipart].nbonds =0;
        if(sys.gravity>0){
            Egrav_tot+=gravitational_energy_ipart(psl, ipart);
        }
    }

    //calculate patchy Simons potential interaction between the neighboring particles 
    /* use neighborlist from here on*/
    count=0;

    for(a = 0;a<list.nimages;a++) {
        i = list.nli[a].ipart;
        k = list.nli[a].image_id;

        if ((k<0)|| (k>=NSHIFTS)) error("force:k out of range\n");   
        psi = &psl->pts[i];

        vector_add(psi->r,list.trans[k],ri) ;
        
        for(b=list.nli[a].first; b<list.nli[a+1].first; b++){
            j = *b;                     /* b is a pointer to the next particle number in the image */
            psj =&psl->pts[j];  
            vector_minus(psi->r,psj->r,rij);
            // pbc(rij,sys.boxl); /*do perform pbc when using no nnl*/
            
            r2 = vector_inp(rij,rij);
            r=sqrt(r2);
            r_radii=sys.particletype[psi->ptype].radius+sys.particletype[psj->ptype].radius;
            s=r-r_radii;
            // gprint(s);

            if(s<sys.rcutoff && s>0.0) { 
                //calculate repulsive energy(=isotropic,i.e. not dep on patch angles)
                //I need to think about how I want to make the repulstion work.
                // Erep=potential_repulsive_energy_sdist(s);

                // Erep_tot += Erep; 
                /*this calculates: Vc(r)*S */
                Ebond=potential_attractive_bond_energy(psl,i,j);
                Erep=potential_repulsive_energy_sdist(s_patch);
                Erep_tot += Erep;

                // printf("The distance s_patch=%.5lf, Ebond=%.5lf, Erep=%.5lf,\n", s_patch,Ebond,Erep);

                if((Ebond<sys.bond_cutoffE)  ) { //&& (i>j)
                    psi->bonds[psi->nbonds]=j;
                    psi->bond_energy[psi->nbonds]= Erep+Ebond;
                    psi->bond_distance[psi->nbonds]=s;
                    psi->nbonds++; 

                    psj->bonds[psj->nbonds]=i;
                    psj->bond_energy[psj->nbonds]= Erep+Ebond;
                    psj->bond_distance[psj->nbonds]= s;
                    psj->nbonds++;

                    tot_bonds++;
                }
            
                //add it to the potential energy of the patch
                Ebond_tot += Ebond ;
            }
            else if(s<0.0){
                // printfr_{min}(psj->r);

                Ebond_tot +=1e6;
            }
            count++;
        }
    }
    if (count!= list.nneighbors) error("nlist corrupted in energy"); 
    // printf("tot bonds are %d \n\n",tot_bonds);
    psl->nbonds=tot_bonds;
    Etot = Ebond_tot + Erep_tot + Egrav_tot;

    return Etot;
}
