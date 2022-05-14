#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"

void wall_interaction_force_ipart(Slice *, int );

void calculate_total_forces(Slice *psl) {
    // error("force calculation not adjusted to the new potential yet. ");
    Pts *psi;
    int ipart,jpart;

    //initialize all forces to zero and update patch vectors
    for( ipart=0; ipart<sys.npart; ipart++ )  {
        psi = &psl->pts[ipart];
        psi->f = nulvec;
        psi->t = nulvec;
        //update of patch vectors is aparte function
        update_patch_vector_ipart(psl, ipart);
    }

    //calculate and apply gravitational forces
    if(sys.gravity>0){
        for( ipart=0; ipart<sys.npart; ipart++ ) {
            gravitational_force_ipart(psl,ipart);

            if (sys.wall_int!=0){
                wall_interaction_force_ipart(psl,  ipart);
            }
        }
    }

    //calculate forces based on force definition is simonsparameters.c
    if (sys.nearest_neighbor==1){
        calculate_forces_neighborlist(psl);
        
    }
    else{
        for( ipart=0; ipart<sys.npart; ipart++ ) {
            for( jpart=ipart+1; jpart<sys.npart; jpart++ ) {
                //just give here nulvec as vector, it rij will be calculated if sys.nearest_neighbor!=1
                single_bond_force(psl,ipart,jpart, nulvec);
            }
        }
    }

    
    return ;
}

double F_harmonic_oscillator(double k, double x){
    return -1*k*x;
}

void wall_interaction_force_ipart(Slice *psl, int ipart){
    /* if the wall has an interaction with the patches
    there is a harmonic oscillator that governs the interaction strength with its minimum below 0*/
    Pts *psi;
    psi=&psl->pts[ipart];
    int isite;
    vector z_axis=nulvec,x_axis=nulvec,y_axis=nulvec,p1;
    double Fwall_tmp=0,Fwall,z,angle_max=sys.particletype[psi->ptype].cosdelta;
    double zcut=sys.particletype[psi->ptype].zcut+sys.particletype[psi->ptype].radius*1.1,angle,cosangle,S,F_height;
    z_axis.z=1;
    x_axis.x=1;
    y_axis.y=1;

    vector rcrosspi,piperpr;
    double F_harm, E_harm, fmag, Umag,fmagP,dSdcostheta,fmagrinv;

   
    // only if the particle is in the vicinity of the wall: 
    if (psi->r.z<=zcut & psi->r.z>0){
        // find thenpatch that is pointing downwards. take the innerproduct with patch vector if 
        for(isite=0;isite<sys.particletype[psi->ptype].nsites;isite++){
            p1=psi->patchvector[isite];
            if (p1.z>0){
                continue;
            }
            cosangle=-vector_inp(p1,z_axis); 
            // angle=acos(cosangle);
            if (cosangle<angle_max){
                continue;
            }
            
            // F_height=gravitational_force_ipart(psl, ipart); // gravitational energy is zero at the wall
            z=fabs(p1.z-zcut);

            Umag=fabs(harmonic_oscillator(sys.wall_int,z));
            fmag=fabs(F_harmonic_oscillator(sys.wall_int,z));
            
            S= Saccent( cosangle, psi->ptype);

            fmagP = S*fmag;
            
            scalar_mintimes(z_axis,fmagP,psi->f);

            /*the torque = p x rnorm , so here rcrosspi = -torque (due to direction of rnorm)*/
            vector_cross(z_axis,p1,rcrosspi);

            /*the direction of the force that bring patch back to center*/
            /*rnorm points form j to i*/
            /*pi is projected onto pi*/
            vector_cross(z_axis,rcrosspi,piperpr);
            
            //"sliding" force, the one perperdicular to rnorm
            //CALCULATE DEL U / DEL COSTHETA_i derivative to i
            
            dSdcostheta = dSaccent_dcosangle(cosangle,psi->ptype);
            fmag = fabs((Umag)*dSdcostheta);
            fmagrinv=fmag/z;

            scalar_plustimes(piperpr,fmagrinv,psi->f);

            /*rcrosspi is -torque and fmag is -F, therefore plustimes */
            scalar_plustimes(rcrosspi,fmag,psi->t);

        }
    }
    
    
    return  ;
}

void calculate_forces_neighborlist(Slice *psl){
    /*uses the neighbor list to loop over the particles and calculates the single bond force
    the gravitational force and patch vectors are already calculated at this point*/

    Pts *psi,*psj;
    vector rij,ri;  
    int ipart,i,k,a,j,*b;

    list.count=0;
    for(a = 0;a<list.nimages;a++) {
        i = list.nli[a].ipart;
        k = list.nli[a].image_id;
        if ((k<0)|| (k>=NSHIFTS)) error("force:k out of range\n");   
        psi = &psl->pts[i];
        vector_add(psi->r,list.trans[k],ri) ;

        for(b=list.nli[a].first; b<list.nli[a+1].first; b++){
            j = *b;
            psj =&psl->pts[j];  
            vector_minus(ri,psj->r,rij); 
            // vector_times(rij,sys.boxl,rij); //??
            // pbc(rij,sys.boxl); // is this necessary? Im not sure?

            // give the rij to the single_bond_force 
            single_bond_force( psl,  i,  j,rij);

            list.count++;
        }
    }
    if (list.count!= list.nneighbors) error("nlist corrupted in force"); 

    return;
}

void single_bond_force(Slice *psl, int ipart, int jpart, vector rij){
    /*the forces acting on the particle: brownian motion and bond forces (no gravitational force)
    the vector rij is needed if you use nearestneighbors
    paper: 10.1080/00268970601075238 Expressions for forces and torques in molecular simulations using rigid bodies, allen & germano */

    Pts *psi, *psj;
    int isite,jsite;

    vector rnorm;
    vector rcrosspi,rcrosspj,piperpr,pjperpr;

    double r2, r, r_radii,s, rinv;
    double fmag, Umag,S,fmagP,fmagrinv,dSdcostheta;
    double cositheta,cosjtheta,cosdelta_i,cosdelta_j;

    psi = &psl->pts[ipart];
    psj = &psl->pts[jpart];
    
    if (sys.nearest_neighbor!=1){
        vector_minus(psi->r,psj->r,rij); //rij points j-->i
        pbc(rij,sys.boxl);
    }

    r2 = vector_inp(rij,rij);
    r=sqrt(r2);
    rinv=1./r;
    r_radii=sys.particletype[psi->ptype].radius+sys.particletype[psj->ptype].radius;
    s=r-r_radii;


    if(s<0.0){
        printf("s = %lf <0.0\n", s);
        printf("r2 = %lf \n", r2);
        printf("r = %lf \n", r);
        vprint(rij);
        dprint(ipart);
        vprint(psi->r);
        dprint(jpart);
        vprint(psj->r);
        error("Overlapping particles!!!");
    } 

    // printf("the calculation the force between particles %d and %d starts \n",ipart,jpart );
    // printf("they have forces: \n");
    // vprint(psi->f);
    // vprint(psj->f);

    if(s<sys.rcutoff && s>0.0){ 
        fmag=fabs(bond_repulsive_force(s));/*F=-dVrep/ds >= 0*/

        scalar_divide(rij,r,rnorm);
        // vprint(rnorm);

        // q = sys.Arep*exp_app(-(s)/sys.dl);  // enumerator
        // g = 0.58 + 3.2*(s);           // denominator
        // fmag = q/(sys.dl*g) + (q*3.2)/(g*g);    // total repulsive force

        /*add the repulsive force (=a positive number) positively to i and negatively to j (as rij points toward from j to i)*/
        /*>>>>>SHOULDNT IS BE rnorm?? instead of rij<<<<<*/
        scalar_plustimes(rnorm,fmag,psi->f);
        scalar_mintimes(rnorm,fmag,psj->f);
        /*repulsive force done*/

        //calculate angular part of patch force for each site
        for( isite=0; isite<sys.particletype[psi->ptype].nsites; isite++ ) {
            vector p1=psi->patchvector[isite];
            cositheta=-vector_inp(rnorm,p1);
            cosdelta_i=sys.particletype[psi->ptype].cosdelta;

            if(cositheta<cosdelta_i  ) {
                continue;
            }

            for( jsite=0; jsite<sys.particletype[psj->ptype].nsites; jsite++ ) {
                vector p2=psj->patchvector[jsite];
                cosjtheta=vector_inp(rnorm,p2);
                cosdelta_j=sys.particletype[psj->ptype].cosdelta;

                if(cosjtheta<cosdelta_j) {
                    continue;
                }

                /*if you are here. you found the patches that are connected*/
                
                // /* attractive energy and force of the patches wihtout angular correction */
                Umag = fabs(potential_attractive_energy_sdist(s)) ;  // Umag is a negative number; do absolute value
                fmag = fabs(bond_attractive_force(s)); // Fc <=0 , thus fabs() 


                // Gets correct switching function. 0<=S<=1
                S = S_value(psl,  cositheta,  cosjtheta, p1,  p2, rnorm ,  ipart,  jpart);
                // the magnitude of the force on the patches. take absolute to get magnitude
                fmagP = S*fmag;

                /*the effective attractive force(fmagP). The particle attract, so mintimes on i and plustimes on j */
                // scalar_mintimes(rnorm,fmagP,psi->f);
                // scalar_plustimes(rnorm,fmagP,psj->f);

                //the direction along rnorm
                scalar_mintimes(rnorm,fmagP,psi->f);
                scalar_plustimes(rnorm,fmagP,psj->f);

                // printf("patch radial interaction \n");
                // vprint(psi->f);
                // vprint(psj->f);
                // vector cross is right hand rule, wijsv. = a , middel v. = b, duim: a x b, 
                /*the torque = p x rnorm , so here rcrosspi = -torque (due to direction of rnorm)*/
                vector_cross(rnorm,p1,rcrosspi);
                vector_cross(rnorm,p2,rcrosspj);

                /*the direction of the force that bring patch back to center,ie. aligning with rnorm*/
                /*rnorm points form j to i*/
                /*pi is projected onto pi*/
                vector_cross(rnorm,rcrosspi,piperpr);
                vector_cross(rnorm,rcrosspj,pjperpr);

                //"sliding" force, the one perperdicular to rnorm see eqn 4 of the Allen paper
                //CALCULATE DEL U / DEL COSTHETA_i derivative to i
                dSdcostheta = dS_dcostheta(psl,cositheta, ipart,cosjtheta,jpart);
                fmag = fabs((Umag)*dSdcostheta);
                fmagrinv=fmag*rinv;

                scalar_plustimes(piperpr,fmagrinv,psi->f);
                scalar_mintimes(piperpr,fmagrinv,psj->f);

                //there is also a costheta_ij force actually 


                /*rcrosspi is -torque and fmag is -F, therefore plustimes */
                scalar_plustimes(rcrosspi,fmag,psi->t);
                

                // DEL U / DEL COSTHETA_j derivative to j
                dSdcostheta  = dS_dcostheta(psl,cosjtheta,jpart,cositheta, ipart);
                fmag = fabs((Umag)*dSdcostheta);
                fmagrinv=fmag*rinv;
                scalar_mintimes(pjperpr,fmagrinv,psi->f);
                scalar_plustimes(pjperpr,fmagrinv,psj->f);

                //torque on other particle
                scalar_mintimes(rcrosspj,fmag,psj->t);

            }
        }
    }
    if(s<0.0 && sys.lincs==0){
        /*im not sure if this if statement should result an error or print statement. look into it*/
        printf("s = %lf <0.0\n", s);
    }
    

    return;
}

void lincs_solve(Slice *psl, double sol[], vector B[sys.npart], double A[sys.npart][cconst.cmax], double rhs[2][sys.npart]){
    /*written in such a way that it only works on a consecutive chain */
    int w, i, j,n,tol, iters;
    Pts *psi,*psj;
    double sol_check, sol_err;
    vector ir, jr;
    w = 2;
    do{
        sol_err = 0.000000001;
        for(i=0; i<sys.npart; i++){
            psi= &psl->pts[i];
            rhs[w-1][i] = 0;
            for(n=0;n<psi->nbonds;n++){
                j = psi->bonds[n];
                if (j<i){
                    continue;
                }
                if(psi->nbonds==1){
                    rhs[w-1][i] += A[i][0]*rhs[2-w][j];
                }
                else if(psi->nbonds==2){
                    rhs[w-1][i] += A[i][0]*rhs[2-w][j];
                    rhs[w-1][i] += A[i][1]*rhs[2-w][j];
                }
                sol_check = sol[i];
                sol[i] += rhs[w-1][i];
                if( fabs(sol_check - sol[i]) > sol_err ){
                    sol_err = fabs(sol_check - sol[i]);
                }
            }
        }
        w = 3 - w;
    } while(sol_err > 0.000000001);

    for(i=0; i<sys.npart; i++){
        psi= &psl->pts[i];
        for(n=0;n<psi->nbonds;n++){
            j = psi->bonds[n];
            if (j<i){
                continue;
            }
            scalar_times(B[i], cconst.cm*sol[i], ir);

            vector_minus(psl->pts[i].r, ir, psl->pts[i].r);
            if(sys.nearest_neighbor){
                vector_minus(psi->dr, ir, psi->dr);
            }

            vector_add(psl->pts[j].r, ir, psl->pts[j].r);
            if(sys.nearest_neighbor){
                vector_add(psl->pts[j].dr, ir, psl->pts[j].dr);
            }
            
        }
    }
    return;
}

void lincs(Slice *psl_old, Slice *psl){

    int i, j,q,n;
    Pts *psi, *psj;
    double r, p;
    double A[sys.npart][cconst.cmax], sol[sys.npart], rhs[2][sys.npart];
    vector ir, jr, rij;
    vector B[sys.npart];

    // set B matrix : gradient of constraints
    for(i=0; i<sys.npart; i++){
        psi= &psl->pts[i];
        for(n=0;n<psi->nbonds;n++){
            j = psi->bonds[n];
            if (j<i){
                continue;
            }
            ir = psi->r;
            jr = psl->pts[j].r;

            vector_minus(ir, jr, rij);
            r = sqrt(vector_inp(rij,rij));
            scalar_divide(rij, r, B[i]);
        }      
    }

    // set A matrix : normalized constraint coupling
    for(i=0; i<sys.npart; i++){
        q=0;
        A[i][q]=0;
        psi= &psl->pts[i];
        for(n=0;n<psi->nbonds;n++){
            j = psi->bonds[n];
            if (j<i){
                continue;
            }

            if(psi->nbonds==1){
                A[i][q] = cconst.cm2*vector_inp(B[i], B[j]);
            }
            else if(psi->nbonds==2){
                A[i][q] = cconst.cm2*vector_inp(B[i], B[j]);
            }
            q++;

            ir = psi->r;
            jr = psl->pts[j].r;
            vector_minus(ir, jr, rij);
            r = sqrt(vector_inp(rij,rij));

            rhs[0][i] = cconst.cm*(vector_inp(B[i], rij) - cconst.d);
            sol[i] = rhs[0][i];
        }
    }

    // SOLVE
    lincs_solve(psl, sol, B, A, rhs);

    // correction for rotational lengthening
     for(i=0; i<sys.npart; i++){
        q=0;
        psi= &psl->pts[i];
        for(n=0;n<psi->nbonds;n++){
            j = psi->bonds[n];
            if (j<i){
                continue;
            }
            ir = psi->r;
            jr = psl->pts[j].r;
            vector_minus(ir, jr, rij);

            p = sqrt(2.*cconst.d2 - vector_inp(rij,rij));
            rhs[0][i] = cconst.cm*(cconst.d - p);
            sol[i] = rhs[0][i];
        }
    }

    //SOLVE
    lincs_solve(psl, sol, B, A, rhs);

    return;
}

void gravitational_force_ipart(Slice *psl, int ipart){
    double zinv,zinv2,zinv6,zinv7,zinv12,zinv13,psiz,fgrav;
    double fg,zcut,b_zc,radius;
    vector z=nulvec;

    Pts *psi;
    psi = &psl->pts[ipart];

    zcut=sys.particletype[psi->ptype].zcut;
    fg=sys.particletype[psi->ptype].fg;
    b_zc=sys.particletype[psi->ptype].b_zc;
    radius=sys.particletype[psi->ptype].radius;
    psiz =psi->r.z-radius;
   
    /*the z unitvector; used for the giving gravity a direction*/
    z.z=1.; 
    if(sys.gravity>0){

        if(psiz<0){
            psiz+=sys.boxl.z;
        }
        if(psiz>=zcut){ 
            fgrav = fg; /*V'_g(z) = -Fg  minus goes later*/
        }
        else{
            zinv = 1./psiz;
            zinv2 = zinv*zinv;
            zinv6 = zinv2*zinv2*zinv2;
            zinv7 = zinv6*zinv;
            zinv12 = zinv6*zinv6;
            zinv13 = zinv12*zinv;
            fgrav = 4.*sys.epsilongravLJ*(-12.*zinv13+6.*zinv7); /* LJ12-6' = 4*epsilonLJ*(-12/r13+6/r7)*/
        }
        
    }

    scalar_mintimes(z,fgrav,psl->pts[ipart].f);
    return ;
}



void forcecheck_nn(Slice *psl){
    /* testing the force difference with and without using neighborlist*/
    vector forcedif=nulvec, tracker=nulvec;
    int i, ipart, jpart;
    vector af[sys.npart], bf[sys.npart];
    vector ar[sys.npart], br[sys.npart];
    double df2;

    printf("checking the force of nn is good \n");
    for(i=0;i<sys.npart;i++){
        bf[i]=nulvec;
        af[i]=nulvec;
        br[i]=nulvec;
        ar[i]=nulvec;

        psl->pts[i].f=nulvec;
    }

    //calculate forces based on force definition is simonsparameters.c
    for( ipart=0; ipart<sys.npart; ipart++ ) {
        for( jpart=ipart+1; jpart<sys.npart; jpart++ ) {
            //just give here nulvec as vector, it rij will be calculated if sys.nearest_neighbor!=1
            single_bond_force(psl,ipart,jpart, nulvec);
        }
    }
    
    for(i=0;i<sys.npart;i++){
        bf[i]=psl->pts[i].f;
        br[i]=psl->pts[i].r;
        psl->pts[i].f=nulvec;
    }

    calculate_forces_neighborlist(psl);
    for(i=0;i<sys.npart;i++){
        af[i]=psl->pts[i].f;
        ar[i]=psl->pts[i].r;
    }

    for(i=0;i<sys.npart;i++){
        vector_minus(af[i], bf[i], forcedif );
        vector_add(forcedif, tracker, tracker)
    }
    // vprint(tracker);
    df2 = vector_inp(tracker,tracker);
    if(fabs(df2)>0.9){
        printf("The df2 is %.12lf \n", vector_inp(tracker,tracker));
        error("force of neighborlist not good");
    }
    else{
        printf("Force is good, difference df2 is %lf\n", df2);
    }
    return;
}

void propagate_bd(Slice *psl)  {

    int ipart,inter, update=0;
    Pts *psi;
    
    Particletype *ptypei;
    vector theta,f,t,u,evi;
    quaternion qu1,qu2,qu3,qprime,qua;
    tensor rotmat;
    quattensor Bmat;
    double dt,lambdaq,fnorm, dr2,ctime;
    double Eparticle, s_dist;

    Slice *copyslice = malloc(sizeof(Slice));
    

    for( inter=0; inter<langevin.ninter; inter++) {
        if(sys.lincs==1){ // save old constrained positions
            // memcpy(copyslice, psl, sizeof(Slice *));
            memcpy(copyslice, psl, sizeof(Slice **));
            // save_old_positions(psl_old, psl);
        }
        /* calculate the forces on all particles*/
        calculate_total_forces(psl); 

        for( ipart=0; ipart<sys.npart; ipart++) {
            psi = &psl->pts[ipart];
            ptypei=&sys.particletype[psi->ptype];
            // printf("start propagate particle %d\n", ipart );
            // printf("its position is :\n");
            // vprint(psi->r);
            // printf("its force is :\n");
            // vprint(psi->f);

            //translational part due to force
            //DT*F*dt*Beta = muT*F*dT see eqn 9 first term of Ilie2015
            /* scalar_times(psi->f,langevin.dtD,f);
             the above line this was written, it uses langevin.dtD=sqrt(2.0*langevin.timestep*sys.temp) instead of
            langevin.dtBeta=langevin.timestep*sys.beta . see eqn 9 the first term on rhs of Ilie2015*/
            
            scalar_times(psi->f,langevin.timestep,f);
            scalar_times(f,sys.mobilityT,f);
            vector_add(psi->r,f,psi->r);
            // vprint(psi->r);
            if(sys.nearest_neighbor){ vector_add(psi->dr,f,psi->dr); }

            //random translation
            theta=RandomBrownianVector(langevin.dtD);
            scalar_times(theta,sys.sqrtmobilityT,theta);
            vector_add(psi->r,theta,psi->r);
            if(sys.nearest_neighbor){ vector_add(psi->dr,f,psi->dr); }
            
            
            //get matrices for q(t)
            //rotmat=A (eqn 12 Ilie2015), Bmat= matrix for quaternion rotation: Baalpha
            rotmat = getrotmatrix(psi->q); 
            Bmat = getquatmatrix(psi->q);

            if(ptypei->activity==1){
                evi=nulvec;
                // deterministic self-propulsion for translation
                matrix_x_vector(rotmat,ptypei->active_force.r,evi);
                scalar_times(evi,sys.mobilityT*ptypei->active_force.F0*langevin.dtBeta,evi);
                vector_add(psi->r,evi,psi->r);
                if(sys.nearest_neighbor){ vector_add(psi->dr,f,psi->dr); }
            }

            if(sys.nearest_neighbor!=1){
                pbc(psi->r,sys.boxl);
            }

            //rotational part due to force. eqn 13 first term is qu1
            //Baalpha * (muR) * rotmatA * torque * delt = qu1
            scalar_times(psi->t,langevin.dtBeta,t);
            //do not forget, multiply with the inverse rotation matrix to convert to body-fixed torque
            matrixT_x_vector(rotmat,t,u);
            /*unit of mobility_R, rad/s? anwser: I think it is in units of sin(a/2) where a is the rotation angle.

            first, vector(torque) x scalar -> vector(3x1) , matrix(3x3) x vector(3x1) -> vector(u), 
            vector u (3x1)x scalar mobilityR (1x1)-> vector(3x1)  , matrix Bmat (4x3) x vector(3x1) -> quaterion (4x1)
            it leads to a quaternion q.0 = cos(a/2) q.1 =sin(a/2)i, q.2 = sin(a/2)j, q.3=sin(a/2)k where (i,j,k) is the rotating axis.
            as you start with a vector --> (i,j,k) and you multiply it with
             */
            scalar_times(u,sys.mobilityR,u);
            /* ok here you multiply u (torque x timestep x drag) */
            quatmatrix_x_vec(Bmat,u,qu1);

            
            //Baalpha * (muR) * theta = qu2
            /*make a random rotating axis first, it is not unit? */
            theta=RandomBrownianVector(langevin.dtD);
            scalar_times(theta,sys.sqrtmobilityR,theta);
            //Bmat is the matrix representation of the quaternion of the particle at time t
            quatmatrix_x_vec(Bmat,theta,qu2);

            //qprime = qu1+qu2+q(t) for the first time
            quat_add(qu1,qu2,qprime);
            quat_add(qprime,psi->q,qprime);

            //find lambdaq, I guess it is th smallest...check which one keeps q(t+delt)^2=1
            //woohoo it works for the smaller lambdaq, maybe also simply for the bigger...
            lambdaq=langrange_multiplier_quat(qprime, psi->q);
            if(lambdaq>1e20) {
                //langrange multiplier did not work, simply renormalize quaternion
                scdivide_quat(qprime,sqrt(quat_inp(qprime,qprime)),psi->q);
            }
            else {
                //lambdaq*q(t) = qu3
                sctimes_quat(psi->q,lambdaq,qu3);
                //q(t+delt) = qprime+qu3
                quat_add(qprime,qu3,psi->q);
            }

            if(sys.nearest_neighbor){ 
                dr2 = vector_inp(psi->dr , psi->dr);
                if ((dr2 > list.cutoff_buffer2)) {
                    update=1;
                }
            }
        }

        if(sys.lincs==1){ // LINCS algorithm
            // printf("performing LINCS\n");
            lincs(copyslice, psl);
        }

        if (sys.nearest_neighbor==1 && update) {
            // printf("Updating the neighborlist\n");
            update_nnlist(psl);
            update =0;
        }
    }
    // cumulative time 
    sys.c_time+=(langevin.timestep*langevin.ninter);

    free(copyslice);
    return;
}

