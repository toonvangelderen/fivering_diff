#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"


double find_minimum_of_potential(void);
double find_trunc_of_Saccent_1(int);
double find_zcut(double );
void gravitational_parameters(void);

void read_input(void);
vector check_read_unitvec(vector );
void init_simtype(Slice *);
void print_input(void);
void setup_delta(void);
void derivative_check(Slice *);
void print_particle_properties(void);
void check_input_with_MAXDEFS(void);
void distances_bulkpatch_particles(void);
void check_boxlengths(void);

double Derjaguin_prefactor(double , double ); 


double Derjaguin_prefactor(double radius1, double radius2){
    /*take it relative wrt Rp=1.0micron= 0.3125 sigma*/
    double sigma=3.7,Derjaguin_prefactor,Rp,fitted_Derj;
    
    // printf("calculating Derjaguin_prefactor\n");

    if((radius1<1e-2) ||(radius2<1e-2)  ){
        gprint(radius1);
        gprint(radius2);
        error("radius of patch particle too small");
    }
     
     Rp=0.5/sigma ;
     fitted_Derj=(Rp*Rp)/(Rp+Rp);
     Derjaguin_prefactor=(radius1*radius2)/(radius1+radius2);


    return Derjaguin_prefactor/fitted_Derj;

}

void derivative_check(Slice *psl){
    /* this is to compare (and check) the angular part and (radial) of derivative the potential
    The derivative is numerically approximated at various delta_cos and the difference with the exact solution.

    I check the derivative at max_iter_q=100 points (0.25 degrees from eachother)
    so at theta=0,1,2,3,4,..., 25 degrees
    I calculate the energy at theta, theta+delta_cos etc where delta_theta= E-6, E-5, E-4, E-3, E-2, E-1
    Calculate the force numerically by evaluating the energy difference F = -dE/dr 
    note thet dr = 0.5PI * delta_theta / 180  

    There are 2 particles,
    rotate one particle, keep the other fixed*/


    double delta_theta=0.1,cosphi;
    double dV_NUM, dV_exact,max_angle=25.0, dangle=0.25;
    int max_iter_p = 10, max_iter_q=(int)max_angle/dangle;
    
    int p,p0=0,q, old_nsites,gravity_psl,npart_psl, ipart, pt;
    vector old_site, rnorm,p1,p2,rij;
    quaternion dq,dqrot;
    double theta_check,y_prev,dcostheta, Fnum, Fexact,E0,E1,theta=0,r,cos1theta,cos2theta;
    double Eatr, S1,S0;

    Pts *psi, *psj;
    Particletype old_pts_types[PTYPES];

    
    //create a temporary slice for the 2 particles to perform the force check
    Slice *copyslice = malloc(sizeof(Slice));
    memcpy(copyslice, psl, sizeof(Slice));

    // save old settings
    gravity_psl=sys.gravity;
    npart_psl=sys.npart;

    //turn gravity off and set npart to 2
    sys.gravity=0;
    sys.npart=2;

    //loop over the particle types
    for (pt=0; pt<sys.nparticle_types;pt++){
        // printf("/nCalculating forces numerically for patchtype %d:\n", pt);
        // printf("Fnum   =- (E1-E0)/dcosphi \n");
        // printf("Fexact =- dE/dcosphi      \n");
        // // printf("            Theta         |    delta Theta     |     F numerically  |       F exact      |  Fexact-Fnum   |     E0 - E1     |          E0        |        S0       |          E1        |        S1       |\n");
        // printf("   Theta_0  | delta Theta |    check Theta_1     |     F numerically  |       E 1      |       F exact      |  Fexact-Fnum   \n");
        // // printf("__________________________________________________________________________________________________________________________________________________________________________________________________________\n");
        // printf("________________________________________________________________________________________________________________________________\n");

        //save nsites, and first site of particle type
        old_nsites=sys.particletype[pt].nsites;
        old_site=sys.particletype[pt].site[0].r;
        

        //put the nsites to 1 and site as (0 1 0)
        sys.particletype[pt].nsites=2;
        sys.particletype[pt].site[0].r=nulvec;
        sys.particletype[pt].site[0].r.y=1;

        y_prev=-0.48*sys.boxl.y;
       // setup the 2 particles, which both have 1 patch
        for( ipart=0; ipart<sys.npart; ipart++) {
            psi=&copyslice->pts[ipart]; 
            psi->ptype=pt;

            psi->r=nulvec;
            psi->r.y = y_prev + (sys.particletype[psi->ptype].radius+sys.s_min) ;                    
            y_prev = psi->r.y+sys.particletype[psi->ptype].radius;
            
            //this makes the patches in opposite directions
            if(ipart==0) {
                psi->q.q0 = 1.0;
                psi->q.q1 = 0.0;
                psi->q.q2 = 0.0;
                psi->q.q3 = 0.0;
            }
            if(ipart==1) {
                psi->q.q0 = 0.0;
                psi->q.q1 = 1.0;
                psi->q.q2 = 0.0;
                psi->q.q3 = 0.0;
            }
        }


        vector_minus(copyslice->pts[0].r,copyslice->pts[1].r,rij);
        r = sqrt(vector_inp(rij,rij));
        scalar_divide(rij,r,rnorm);

        // first loop changes the cos / composition
        theta=0;
        Eatr = potential_attractive_energy_sdist(sys.s_min);
        for(q=0; q<max_iter_q;q++){
            
            // Calculate the energy at delta_theta=0
            update_patch_vectors(copyslice);

            p1=copyslice->pts[0].patchvector[0];
            p2=copyslice->pts[1].patchvector[0];

            cos1theta=-vector_inp(rnorm,p1);
            cos2theta=vector_inp(rnorm,p2);

            cosphi=cos(theta/180.*PI);


            S0=S_value(copyslice,cos1theta,cos2theta, p1,p2, rnorm, 0, 1);
            E0=S0*Eatr;

            // Fexact=dSaccent_dcosangle(costheta, pt);
            Fexact=dS_dcostheta(copyslice,cos1theta, 0,cos2theta,1);
            // Fexact = dSaccent_dcosangle(cos1theta,pt)*Saccent(cos2theta,pt);

           
            delta_theta=0.1;
            //second loop changes delta_theta= 0.001, 0.01, 0.1 etc
            //I put p0 and max_iter_p to 6 and 7, because I only want to check at delta theta=1E-6
            p0=5;
            max_iter_p=p0+1;

            for(p=0;p<max_iter_p; p++){
                //so you start with a conformation at a specific angle 0, 0.25, etc
                //change the angle according to delta_theta , calculate the energy and rotate it back for the next delta_theta
                if (p<p0){ 
                    delta_theta/=10.;
                    continue;
                }
                //rotate particle 0 along x-axis 
                dq=QuaternionXaxis(delta_theta);
                quat_times(dq,copyslice->pts[0].q,dqrot);
                copyslice->pts[0].q = dqrot; 
                update_patch_vectors(copyslice);

                p1=copyslice->pts[0].patchvector[0];
                cos1theta=-vector_inp(rnorm,p1);
                theta_check = acos(cos1theta)/PI*180;

                S1=S_value(copyslice,cos1theta,cos2theta, p1,p2, rnorm, 0, 1);
                E1=S1*Eatr;
                
                dcostheta= cos((theta)/180*PI) - cos((theta+delta_theta)/180*PI) ;
                Fnum = - (S1-S0)/dcostheta;
                
                if (Fexact- Fnum>1e-2){
                    printf("   Theta_0  | delta Theta |    check Theta_1     |     F numerically  |       E 1      |       F exact      |  Fexact-Fnum   \n");
                    printf("    %4.2f       1E-%2d       %20.15f %20.15f %20.15f %20.15f %20.15f \n", theta, p+1 , theta_check, Fnum, E1 ,Fexact,Fexact- Fnum ); //, E0-E1, E0, S0, E1, S1);

                    error("Fexact- Fnum>1e-3 in the angular force");
                }
                // printf("%20.15f %20.15f %20.15f %20.15f %16.15f %16.15f %20.15f %16.15f %20.15f %16.15f\n", theta, dcostheta , Fnum, Fexact,Fexact- Fnum , E0-E1, E0, S0, E1, S1);
                // printf("    %4.2f       1E-%2d       %20.15f %20.15f %20.15f %20.15f %20.15f \n", theta, p+1 , theta_check, Fnum, E1 ,Fexact,Fexact- Fnum ); //, E0-E1, E0, S0, E1, S1);

                //rotate back particle 0 along x-axis 
                dq=QuaternionXaxis(-1.*delta_theta);
                quat_times(dq,copyslice->pts[0].q,dqrot);
                copyslice->pts[0].q = dqrot; 
                update_patch_vectors(copyslice);

                delta_theta/=10.;
            }
            // printf("________________________________________________________________________________________________________\n");
            dq=QuaternionXaxis(dangle);
            quat_times(dq,copyslice->pts[0].q,dqrot);
            copyslice->pts[0].q = dqrot; 
            update_patch_vectors(copyslice);
            theta=theta+dangle;
            // go from old to new position
        }
        sys.particletype[pt].site[0].r=old_site;
        sys.particletype[pt].nsites=old_nsites;
    }

    sys.gravity=gravity_psl;
    sys.npart=npart_psl;
    /* I checked if original particles are still there, (they are in psl, not in copyslice)*/
    //free the memory of the copyslice
    free(copyslice);

    return;
}

double find_minimum_of_potential(void){
    double s, Utot, Urep, Uattr, Utot_tracked=0, s_min=1.1, Uattr_tracked=0, Urep_tracked=0;
    int r;

    for (r=1;r<100000;r++){
        s=r/100000. *0.10;

        Urep = potential_repulsive_energy_sdist(s);
        Uattr = potential_attractive_energy_sdist(s);
        Utot = Urep+Uattr;
        if (sys.S_fixed>0){
            Utot*=sys.S_fixed;
        }

        if (Utot<Utot_tracked){
            s_min=s;
            Utot_tracked=Utot;
            Uattr_tracked=Uattr;
            Urep_tracked=Urep;

        }
    }
    sys.Erep_smin=Urep_tracked;
    sys.Ec_smin=Uattr_tracked;
    sys.E_smin=Utot_tracked;
    printf(" the colloid-colloid distance with minimum energy sys.E_smin=%lf (sys.Ec_smin=%lf, sys.Erep_smin=%lf) is %lf \n", sys.E_smin, sys.Ec_smin, sys.Erep_smin, s_min );
    return s_min;
}

double find_trunc_of_Saccent_1( int ptype){
    /* find the angle at which S' is almost zero,i.e. < treshold*/
    double theta, cosphi,S_new,treshold=5e-4,max_angle=90.0;
    int i,imax=100000000;

    for(i=0;i<imax; i++){
        theta=(double)i/imax*max_angle;
        cosphi=cos(theta/180.*PI);        
        S_new=Saccent(cosphi,ptype);

        if (S_new<treshold){
            gprint(S_new);
            return cosphi;
        }
    }
    gprint(treshold);
    error("cosphi not found.");
    
    return -1;
}

double find_zcut(double fg){
    /*fg_LJ = 4.*sys.epsilongravLJ*(-12.*zcutinv12*zcutinv+6.*zcutinv6*zcutinv);
    zcut between 1 and 1.2*/

    int z;
    double z_search,zmax=1e8;
    double fg_LJ;
    double zcutinv, zcutinv2,zcutinv6,zcutinv12;

    for(z=0;z<(int)zmax;z++){
        z_search = (double)z/zmax*0.2+1.;
        zcutinv = 1./z_search;
        zcutinv2 = zcutinv*zcutinv;
        zcutinv6 = zcutinv2*zcutinv2*zcutinv2;
        zcutinv12 = zcutinv6*zcutinv6;
        fg_LJ = 4.*sys.epsilongravLJ*(-12.*zcutinv12*zcutinv+6.*zcutinv6*zcutinv);
        if(fabs(fg - fg_LJ)<=1e-4){
            gprint(fg);
            gprint(fg_LJ);
            gprint(z_search);
            return z_search;
        }
    }
    error("zcut not found");
    return 0;
}

void gravitational_parameters(void){
    double g,delta_rho_kg_m3,fg_LJ, boltzmann, T, ktsigma;
    double zcutinv, zcutinv2,zcutinv6,zcutinv12;
    Particletype *ptypen;
    double fg,zcut,b_zc,r_micron;
    int n;
    
    
    // constants
    g=9.80665; /*gravitation acceleration of eath [m/s^2]*/
    boltzmann=1.38064852e-5; /*joule. the e-5 comes from e-23 *(original value of Boltzmann) times (e-6)^3=e-18 (because in sys.fg we use sys.r_micron^3) */
    T=33.8+273.15;/*K, this an is approximate for all measurements*/
    sys.sigma=3.2e-6;/*sigma is reducing measure for distance 3.2 mum= the colloidal diameter*/

    //some particle indep paramters
    ktsigma = (boltzmann * T)/sys.sigma;
    printf("boltzmann %.12lf T %.12lf sigma %.12lf\n", boltzmann,T,sys.sigma);
    printf("kt/sigma %.12lf \n",ktsigma);
    
    //some particle dep paramters
    for (n=0;n<sys.nparticle_types;n++){
        ptypen=&sys.particletype[n];
        if((ptypen->nsites==2 )|| (ptypen->nsites==1 && ptypen->activity==1)){
            ptypen->delta_rho_kg_m3=60.568460120322576 ; /* the density diffence between the lutide-water solution and the colloid [kg/m^3]*/
        }
        else if((ptypen->nsites==4 )|| (ptypen->nsites==3)){
            ptypen->delta_rho_kg_m3=63.5699387858425;
        }
        else{
            dprint(ptypen->nsites);
            error("npatch not corresponding to 2 or 4. Indicate delta_rho_kg_m3 in the code");
        }
    }

    //calculate the gravity per particle, 
    for (n=0;n<sys.nparticle_types;n++){
        ptypen=&sys.particletype[n];
        //the radius of the particle in micron
        r_micron=ptypen->radius*sys.sigma*1e6;
        

        fg = 4./3.*PI* (r_micron*r_micron*r_micron)*ptypen->delta_rho_kg_m3*g/(ktsigma)*sys.gravity; /* V_g = fg*z [kT/sigma*] at 307K*/
        // printf("fg = %lf [kT/sigma]\n", fg);
        zcut=find_zcut(fg);
        zcutinv = 1./(zcut);
        zcutinv2 = zcutinv*zcutinv;
        zcutinv6 = zcutinv2*zcutinv2*zcutinv2;
        zcutinv12 = zcutinv6*zcutinv6;

        fg_LJ = 4.*sys.epsilongravLJ*(-12.*zcutinv12*zcutinv+6.*zcutinv6*zcutinv);
        if(fabs(fg - fg_LJ)>=0.01){
            printf("the difference of Fg(line)-fgLJ = %lf -%lf  = %lf", fg , fg_LJ, fg - fg_LJ);
            error("gravity force not equal at zcut");
        }
        /* b_zc should include correction of radius. The wall is at 0, the center of the colloid is than at z=radius*/
        b_zc = fg*(zcut) - 4.*sys.epsilongravLJ*(zcutinv12-zcutinv6+1./4) ;
        
        ptypen->fg=fg;
        ptypen->zcut=zcut;
        ptypen->b_zc=b_zc;

    }

    return;
}

void check_boxlengths(void){

    if (sys.boxl.y<1e-5 & sys.boxl.x>1e-5){
        //its a cubic box
        sys.boxl.y = sys.boxl.z = sys.boxl.x;
    }
    else if(sys.boxl.y>1e-5 & sys.boxl.x>1e-5){
        sys.boxl.z = sys.boxl.x;
    }
    else if(sys.boxl.x<1e-5){
        error("box length is not read correctly");
    }
    return; 
}

void setup_derj_scaling_params(void){
    /*calculate the scaling parameter that scale with the radii of the two spheres
    the potential parameter (simon's) are fitted for spheres with radius of 1 micron
    each combination of two particletype have a  Derjaguin prefactor. save these prefactors in the particle type structure
    you can access them via  sys.particletype[i].Derj_scaling[j]*/
    int itype,jtype;
    double Derj_pref;

    for (itype=0;itype<sys.nparticle_types;itype++){
        for (jtype=itype;jtype<sys.nparticle_types;jtype++){
            Derj_pref=Derjaguin_prefactor(sys.particletype[itype].small_radius, sys.particletype[jtype].small_radius);
            if (Derj_pref<1e-1){
                gprint(Derj_pref);
                error("stop, Derjaguin_prefactor very low?");
            }
            sys.particletype[itype].Derj_scaling[jtype]=Derj_pref;
            sys.particletype[jtype].Derj_scaling[itype]=Derj_pref;
        }
    }
    return;
}

void init_model(Slice *psl) {
    /*startis with defining the potential based on rwet and sc,
    then determines the minimum of the potential and de bond energy threshold
    */
    double bond_treshold_fraction=0.001;

    check_boxlengths();
    
    sys.temp = 1.0/sys.beta;
    cluster.update=1; //why is this? 

    if (sys.nearest_neighbor==1 & sys.sim_type==2 & sys.cluster_MC==1){
        error("nearest_neighbor, cluster MC and MC do not go together!");
    }
    
    /*setting up the parameters of the potentials and Derjaguin prefactors*/
    setup_Simons_potential_parameters();


    /*take minimum bond length for bond average*/
    printf("\nLooking for the minimum of the potential...");
    sys.s_min=find_minimum_of_potential();
    printf("found s_min = %lf.\n", sys.s_min);
    if (sys.s_min>1.0){
        // printf("bond cutoff = sys.bond_cutoffE = - (sys.A/sys.xi) * exp(-(sys.s_min-1.)*(sys.s_min-1.)/sys.xi) * 0.50 ; \n");
        error("define sys.s_min such that the minimum energy and the bond definition ( == 1 procent of minimum ) can be calculated ");
    }
    /*bond is defined if energy between particles is less that bond_treshold_fraction*minimum depth */
    sys.bond_cutoffE =(potential_repulsive_energy_sdist(sys.s_min)+potential_attractive_energy_sdist(sys.s_min))* bond_treshold_fraction;
    if (sys.S_fixed>0){
        sys.bond_cutoffE *= sys.S_fixed;
    }

    //first read in the particle, specifies particle types etc
    printf("Setting up the positions and sites for the particles...");
    setup_positions_sites(psl);
    printf("done\n");
    
    //patchwidth variables
    printf("Setting up the patch widths definitions for the particles...");
    setup_delta();
    distances_bulkpatch_particles();
    printf("done\n");
    
    // derjaguin scaling is related to the radius of curvature of the particle (types)
    setup_derj_scaling_params();

    printf("printing potentials ... ");
    plotpotential();
    printf("done\n");


    /*Calculate the gravitational force fg and correction b_zc*/
    if(sys.gravity>0){
        printf("\nSetting GRAVITY PARAMETERS...");
        gravitational_parameters();
        printf("done\n");
    }
    
    /*rcut^2*/
    if (sys.rcutoff>0.5){
        gprint(sys.rcutoff);
        error("sys.rcutoff too big make it smaller than 0.5");
    }
    sys.rcutoffsq=sys.rcutoff*sys.rcutoff;
   
    /*does sim_type specific definitions*/
    init_simtype(psl);


    return;
}

void setup_delta(void){
    // walk over the particle types and identify what kind of switch function is has, and what the delta cutoff is (i.e. where is S' 0)
    int ptype;
    Particletype *part_type;
    double mini_cosphi=1.;

    for(ptype=0;ptype<sys.nparticle_types; ptype++){
        part_type=&sys.particletype[ptype];
        if(sys.S_fixed>0){
            part_type->cosphi=cos(part_type->delta_degree/180.*PI);
        }
        // if delta if >0 then use old switch as S'
        if (s_int.s_accent==0 ){
            printf("OLD SWITCH\n");
            part_type->s_accent=0;
            part_type->cosphi=cos(part_type->delta_degree/180.*PI); // theta cutoff
            part_type->oneover_cosphi = 1.0/(1.0-part_type->cosphi);
        }
        if ( s_int.s_accent==2 ){
            printf("LINEAR SWITCH\n");
            part_type->s_accent=2;
            part_type->cosphi=cos(part_type->delta_degree/180.*PI); // theta cutoff
        }
        //if delta<0 (<1e-5) then use integrated switch function
        else if (s_int.s_accent==1){
            printf("NEW SWITCH\n");
            part_type->s_accent=1;
            part_type->cosphi=find_trunc_of_Saccent_1(ptype);
            // part_type->delta_degree is needed to calculate the distance between bulk and patch particle
        }

        if (part_type->cosphi<mini_cosphi){
            mini_cosphi = part_type->cosphi;
        }
    }

    return;
}

void init_simtype(Slice *psl){
    int  id;
    /*if bmd or ffs*/
    if (sys.sim_type==0 || sys.sim_type==4){
        //do here the force check, you need s_min
        printf("\nChecking the derivatives/forces...");
        derivative_check(psl);
        // gprint(sys.npart); works 
        printf("done\n");
        /*langevin/brownian md parameters*/
        // CHANGE NAME OM dtD and dTBETA
        printf("\nDefining the langevin timestep and mobility parameters... ");
        langevin.dtD=sqrt(2.0*langevin.timestep*sys.temp);
        langevin.dtBeta=langevin.timestep*sys.beta;
        // code works for istropic rotational and translational diffusion
        sys.sqrtmobilityT = sqrt(sys.mobilityT); //eqn 5,  
        sys.sqrtmobilityR = sqrt(sys.mobilityR); // eqn 4 Ilie2015
        printf("done\n");
    }
    else if(sys.sim_type==2) {
        //Monte Carlo
        printf("\nSetting up the MC parameters... ");
        trans.acc=0;
        fintrans.acc=0;
        rot.acc=0;
        finrot.acc=0;
        trans.tries=0;
        fintrans.tries=0;
        rot.tries=0;
        finrot.tries=0;
        sys.drmax=0.05;
        sys.dqmax=0.5;
        if(sys.cluster_MC==1){
            sys.drmax_cluster=0.5;
            sys.dqmax_cluster=20.;
            sys.drmax_mono=5.0;
            sys.dqmax_mono=0.5;
        }  
        printf("done\n");  
    }
     // Rosenbluth Forward Flux Sampling (RBFFS)
    if(sys.sim_type==4){
        int n_int;

        FILE *fplam;
        printf("Reading in interfaces...\n");

        if((fplam = fopen("interfaces.inp", "r"))==NULL)
            error("interfaces.inp not found!\n");

        fscanf(fplam, "%d", &state[0].nrep);
        if(state[0].nrep != sys.ncycle2)
            error("number of interfaces defined at top of file does not match sys.ncycle2");

        for(n_int=0; n_int<state[0].nrep; n_int++){
            fscanf(fplam, "%lf", &state[0].lambda[n_int]);
            printf("interface %d : %lf\n", n_int, state[0].lambda[n_int]);
        }

        if(n_int != sys.ncycle2)
            error("number of interfaces does not match number of values");

        fclose(fplam);

        printf("\nFFS setup: %d initial crossings, %dtrial runs, %d interfaces\n", sys.ncycle1, sys.K, sys.ncycle2);

        state[0].weight=0.000000000000000000000001; // set first RB weight
        init_state = (Slice *)calloc(1,sizeof(Slice));
    }

    if(sys.lincs==1){ // Set up constraint constants
        printf("setting the paramters for lincs\n");
        if(sys.sim_type!=0 && sys.sim_type!=4)
            error("sim_type must be 0 (bmd) or 4 (FFS)!");

        // psl_old=(Slice *)calloc(1,sizeof(Slice)); //??why this 
        cconst.cmax=2;
        cconst.d = sys.s_min+1.; // fixed distance
        cconst.d2 = cconst.d*cconst.d;
        cconst.K = sys.npart-1;
        cconst.cm = 1. / sqrt(2.); // inverse mass term
        cconst.cm2 = cconst.cm*cconst.cm;


    }

    if(sys.nearest_neighbor){
        
    }
    if(cluster.analysis==1){
        // chain. = (Statistics *)calloc(NPART,sizeof(Statistics));
        printf("\nCode performs cluster analysis \n");
        printf("the total energy is %lf\n", total_energy(&slice[0]));
        printf("now cluster analysis\n");
        slice[0].nclusters = cluster_analysis(&slice[0]); 
        clustersize_identification(&slice[0]);
        chain_linkedlist(&slice[0]);
        printf("Initial # bonds is %d \n", slice[0].nbonds);
        bond_length_average(&slice[0]);
          
        printf("bond length average =  %.12lf \n", bond_length.mean);
        printf("bond length var^2   =  %.12lf \n", bond_length.variance2);
        
        strcpy(cluster.size_histogram.filename,"clustersize_histogram.dat");
        strcpy(cluster.size_distribution.filename,"clustersize_distribution.dat");

        strcpy(chain.e2e2.filename,"e2e2.dat");
        
        if(sys.empty_files==0){
            read_statistics_file(&cluster.size_histogram);
            read_statistics_file(&cluster.size_distribution);
            read_statistics_file(&chain.e2e2);
        }   
    }

    return;
}

void plotpotential(void){
    FILE *fp;
    int i,j,n, t;
    double r,s,z, V=0 , V_LJ, rinv, rinv3, rinv6,rinv12,rinv24, V_g, Vrep=0, Vattr=0;
    double Eattr=0, Erep=0, E=0;
    double theta, cosphi,S_old, S_new,S,dS,cosphi_n; 
    
    if ((fp = fopen("pot.dat","w"))==NULL){
        printf("output: can't open pot.dat\n");
        return;
    }
    else {
        for(i=1;i<20000; i++){
            s=i/100000.;
            
            Vrep = potential_repulsive_energy_sdist(s);
            Vattr = potential_attractive_energy_sdist(s);
            if (sys.S_fixed>0){
                Vattr *=sys.S_fixed; 
            } 
            
            V= Vrep + Vattr;
            
            fprintf(fp, "%30.15lf %30.15lf %30.15lf %30.15lf\n", s, V, Vrep, Vattr);//, E, Erep, Eattr);
        }
    fclose(fp);
    }

    if ((fp = fopen("2nd_der.dat","w"))==NULL){
        printf("output: can't open 2nd_der.dat\n");
        return;
    }
    else {
        for(i=1;i<20000; i++){
            s=i/100000.;
            
            Vrep = second_der_Vrep(s);
            Vattr = second_der_Vc(s);
            if (sys.S_fixed>0){
                Vattr *=sys.S_fixed; 
            }           
            V= Vrep + Vattr;
            
            fprintf(fp, "%30.15lf %30.15lf %30.15lf %30.15lf\n", s, V, Vrep, Vattr);//, E, Erep, Eattr);
        }
    fclose(fp);
    }

    if ((fp = fopen("force.dat","w"))==NULL){
        printf("output: can't open force.dat\n");
        return;
    }
    else {
        for(i=1;i<20000; i++){
            s=i/100000.;
            
            Vrep = bond_repulsive_force(s);
            Vattr = bond_attractive_force(s);
            if (sys.S_fixed>0){
                Vattr *=sys.S_fixed; 
            }     
            V= Vrep + Vattr;
            
            fprintf(fp, "%30.15lf %30.15lf %30.15lf %30.15lf\n", s, V, Vrep, Vattr);//, E, Erep, Eattr);
        }
    fclose(fp);
    }


    if ((fp = fopen("switchfunction.dat","w"))==NULL){
        printf("output: can't open switchfunction.dat\n");
        return;
    }
    else {
        
        for(i=0;i<10000; i++){
            theta=i/10000.*50.0;
            cosphi=cos(theta/180.*PI);
            fprintf(fp, "%30.15lf ",theta);
            for(n=0;n<sys.nparticle_types;n++){
                cosphi_n=sys.particletype[n].cosphi;

                if(cosphi<cosphi_n) {
                    S=0;
                    dS=0;
                }
                else{
                    S=Saccent(cosphi,n); 
                    dS=dSaccent_dcosangle(cosphi, n);
                }   
                

                fprintf(fp, "%30.15lf %30.15lf",S,dS);
            }
            
            fprintf(fp, "\n");//, E, Erep, Eattr);
        }
    fclose(fp);
    }


    if(sys.gravity>0){
        if ((fp = fopen("gravity.dat","w"))==NULL){
            printf("output: can't open gravity.dat\n");
            return;
        }
        else {
            
            for(i=0;i<10000; i++){
                r=i/10000.*3+1;
                
                fprintf(fp, "%30.15lf ", r);

                for(n=0;n<sys.nparticle_types;n++){
                    z=r-sys.particletype[n].radius;

                    if(z>=sys.particletype[n].zcut){ 
                        V_g = sys.particletype[n].fg*z - sys.particletype[n].b_zc ; /*V_g = Fg*z - b*/
                    }
                    else{
                        rinv = 1/z;
                        rinv3= rinv*rinv*rinv;
                        rinv6=rinv3*rinv3;
                        rinv12= rinv6* rinv6;
                        rinv24 = rinv12*rinv12;
                        V_g = 4*sys.epsilongravLJ*(rinv12-rinv6+1./4); /* LJ12-6 = 4*b(1/r12-1/r6)*/
                    }

                    fprintf(fp, "%30.15lf ", V_g);
                }
                fprintf(fp, "\n");
            }
            
        fclose(fp);
        }
    }


    return;
}

#define MAXLINE 100
void read_input(void) {

    FILE *fp;
    char *pt,line[MAXLINE];


    if((fp = fopen("path.inp","r"))==NULL) {
        printf("ERROR: could not read input file\n");
        exit(1);
    }

    printf("\nReading input from path.inp\n");

    while(fgets(line,MAXLINE, fp) != NULL) {
        pt = strtok(line," ");
        if( strcmp(pt,"sim_type")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.sim_type);
        } else if( strcmp(pt,"bond_breakage")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.bond_breakage);
        } else if( strcmp(pt,"start_type")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.start_type);
        } else if( strcmp(pt,"particle_setup")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.particle_setup);
        } else if(strcmp(pt,"ncycle1")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.ncycle1);
        } else if( strcmp(pt,"ncycle2")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.ncycle2);
        } else if( strcmp(pt,"graphics")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.graphics);
        } else if( strcmp(pt,"snapshot")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.snapshot);
        } else if( strcmp(pt,"bond_op")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.bond_op);
        } else if( strcmp(pt,"cluster_MC")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.cluster_MC);
        } else if( strcmp(pt,"nearest_neighbor")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.nearest_neighbor);
        } else if( strcmp(pt,"mobilityT")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.mobilityT);
        } else if( strcmp(pt,"mobilityR")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.mobilityR);
        } else if( strcmp(pt,"ninter")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&langevin.ninter);
        } else if( strcmp(pt,"timestep")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&langevin.timestep);
        }  else if( strcmp(pt,"K")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.K);
        } else if( strcmp(pt,"lincs")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.lincs);
        } else if( strcmp(pt,"rdfanalysis")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.rdfanalysis);
        } else if( strcmp(pt,"cluster_analysis")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&cluster.analysis);
        } else if( strcmp(pt,"adjacency")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.adjacency);
        } else if( strcmp(pt,"s_histogram")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.s_distribution);
        } else if( strcmp(pt,"empty_files")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.empty_files);
        } else if( strcmp(pt,"xy_print")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.xy_print);
        } else if( strcmp(pt,"gravity")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.gravity);
        }  else if( strcmp(pt,"S_fixed")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.S_fixed);
        } else if( strcmp(pt,"epsilongravLJ")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.epsilongravLJ);
        } else if( strcmp(pt,"switch_method")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.switch_method);
        } else if( strcmp(pt,"switch_expc")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&s_int.switch_expc);
        } else if( strcmp(pt,"switch_expd")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&s_int.switch_expd);
        } else if( strcmp(pt,"switch_expe")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&s_int.switch_expe);
        } else if( strcmp(pt,"switch_expf")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&s_int.switch_expf);
        } else if( strcmp(pt,"switch_expg")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&s_int.switch_expg);
        } else if( strcmp(pt,"switch_exph")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&s_int.switch_exph);
        } else if( strcmp(pt,"switch_expi")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&s_int.switch_expi);
        } else if( strcmp(pt,"switch_variable")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&s_int.switch_variable);
        } else if( strcmp(pt,"s_accent")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&s_int.s_accent);
        } else if( strcmp(pt,"dT")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.dT);
        } else if( strcmp(pt,"r_wetting")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.r_wetting);
        } else if( strcmp(pt,"surface_charge")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.surface_charge);
        } else if( strcmp(pt,"wall_int")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.wall_int);
        } else if( strcmp(pt,"npart")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.npart);
        } else if( strcmp(pt,"rcutoff")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.rcutoff);
        } else if( strcmp(pt,"beta")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.beta);
        } else if( strcmp(pt,"boxl")==0) { //boxlx is not read as boxlx, because of old inputfiles using boxl
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.boxl.x);
        } else if( strcmp(pt,"boxly")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.boxl.y);
        } else if( strcmp(pt,"boxlz")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.boxl.z);
        } else if( strcmp(pt,"nchains")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.nchains);
        } else if( strcmp(pt,"chaingap")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.chaingap);
        } 
                // } else if( strcmp(pt,"timestep")==0) {
        //     pt = strtok(NULL," ");
        //     sscanf(pt,"%lf",&langevin.timestep);
        // } else if( strcmp(pt,"ninter")==0) {
        //     pt = strtok(NULL," ");
        //     sscanf(pt,"%d",&langevin.ninter);
        // }  else if( strcmp(pt,"nshoot")==0) {
        //     pt = strtok(NULL," ");
        //     sscanf(pt,"%d",&path.nshoot);
        // } else if( strcmp(pt,"nrepswap")==0) {
        //     pt = strtok(NULL," ");
        //     sscanf(pt,"%d",&path.nrepswap);
        // } else if( strcmp(pt,"nstateswap")==0) {
        //     pt = strtok(NULL," ");
        //     sscanf(pt,"%d",&path.nswapstates);
        // } else if( strcmp(pt,"nreverse")==0) {
        //     pt = strtok(NULL," ");
        //     sscanf(pt,"%d",&path.nreverse);
        // } else if( strcmp(pt,"stateswapbias")==0) {
        //     pt = strtok(NULL," ");
        //     sscanf(pt,"%d",&path.stateswapbias);
        // } else if( strcmp(pt,"fixedbias")==0) {
        //     pt = strtok(NULL," ");
        //     sscanf(pt,"%d",&path.fixedbias);
        // }
        else if( strcmp(pt,"\n")==0) {
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"SIMULATION\n")==0) {
            printf("Reading simulation parameters\n");
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"GRAPHICS\n")==0) {
            printf("Reading graphics parameters\n");
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"MC\n")==0) {
            printf("Reading MC parameters\n");
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"BMD\n")==0) {
            printf("Reading BMD parameters\n");
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"HARMONIC_OSCILLATOR\n")==0) {
            printf("Reading HARMONIC parameters\n");
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"FFS\n")==0) {
            printf("Reading FFS parameters\n");
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"ANALYSIS\n")==0) {
            printf("Reading Analysis parameters\n");
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"GRAVITY\n")==0) {
            printf("Reading gravity parameters\n");
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"SWITCH_FUNCTION\n")==0) {
            printf("Reading switch function parameters\n");
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"Simons_Potential\n")==0) {
            printf("Reading potential parameters\n");
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"SYSTEM\n")==0) {
            printf("Reading system parameters\n");
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"K\n")==0) {
            printf("Reading K parameters\n");
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"TIS\n")==0) {
            printf("Reading TIS parameters\n");
        } else if( strcmp(pt,"\n")==0) {
            // pt = strtok(NULL," ");
        } else {
            printf("Keyword unknown: %s\n",pt);
        }
    }
    
    check_input_with_MAXDEFS();

    
    fclose(fp);

    printf("Done reading path.inp\n");

    return;
}

void check_input_with_MAXDEFS(){
    /*in path.h there are max definitions given. 
    Check after reading the input file if any of these max values are crossed*/

    if (sys.npart>NPART){
        dprint(sys.npart);
        dprint(NPART);
        error("sys.npart > NPART");
    }

    return;

}

void distances_bulkpatch_particles(void){
    /*the d distance, is distance between bulk and patch particle */
    double dp_c, dp_p;
    double theta_c , theta_p ;
    double R_c , R_p ;
    double degree_to_rad=PI/180.;
    double AB,AC,BC,S_AS_B;
    int n;
    Particletype *part_type;
    printf("Caluclating the distance between bulk and patch particle\n");
    for(n=0;n<sys.nparticle_types;n++){
        part_type=&sys.particletype[n];

        theta_c=part_type->delta_degree*degree_to_rad;
        gprint(part_type->delta_degree);

        R_c=part_type->radius;
        R_p=part_type->small_radius;
        
        /*first check of the input patchwidths and diameters correspond to a correct system*/
        dp_c= sin(0.5*theta_c)*R_c;
        // dp_p= 2.*sin(0.5*theta_p)*R_p;
        gprint(dp_c);
        gprint(theta_c);
        gprint(R_p);

        /*calcutate the distance from from center big colloid to center small patch colloid*/
        AC=cos(0.5*theta_c)*R_c;
        BC=sqrt(R_p*R_p-dp_c*dp_c);
        gprint(AC);
        gprint(BC);
        // gprint(AB);
        AB = AC-BC;
        S_AS_B= AB+ R_p -R_c ;
        S_AS_B=0;
        gprint(AB);
        gprint(S_AS_B);
        if((AB<0) || (AB>R_c)|| (S_AS_B<0)){
            error("AB<0 or AB>Rc not possible. it is a length from center big colloid to center small patch colloid");
        }
        
        
        part_type->d=(AB-S_AS_B); // [\sigma]
        gprint(part_type->d);
        double d_cbsp= part_type->d+R_p;
        gprint(d_cbsp);
        gprint(sys.s_min);
        gprint(potential_repulsive_energy_sdist(sys.s_min));
        part_type->d_cbsp=d_cbsp; // [\sigma]
        part_type->cosdelta=cos(theta_c);
        gprint(part_type->d_cbsp);
        gprint(theta_c);

        gprint(part_type->cosphi);
        gprint(part_type->cosdelta);
    }
    
    return;
    
}

void read_particletypes(Slice *psl){
    /*read sites and diameter from particle.inp (used to be sites.inp)*/
    /* the file is structured as follows:
    NPARTICLE_TYPES
    NSITES NPARTICLES DIAMETER DELTA ACTIVITY
    SITE_1 
    SITE_N
    ACTIVITY_VECTOR
    etc.
    */
    FILE *file;
    int isite, ptype,a,p,ipart=0;
    int amount, totalparticles=0;
    vector new;
    Particletype *part_type;
    double activity_p;
    sys.dipatch_only=1;

    
    if(sys.particle_setup==1) {
        printf("Reading particles setup directly from particles.inp\n");
        
        if ((file = fopen("particles.inp","r"))==NULL){
            error("input: can't open particles.inp \n");
        }
        else {
            /*first line of sites.inp contains the number of particle types*/
            fscanf(file,"%d\n",&sys.nparticle_types);
            if (sys.nparticle_types>PTYPES){
                dprint(sys.nparticle_types);
                dprint(PTYPES);
                error("there is a boundary of maximum number of particle types");
            }
            /*such that you know how many to read in*/
            for(ptype=0;ptype<sys.nparticle_types; ptype++){
                /*the first line of the particle type contains:
                number of sites, number of particles, diameter in sigma*/
                part_type=&sys.particletype[ptype];
                fscanf(file,"\n%d %d %lf %lf %lf %lf\n",&part_type->nsites, &part_type->nparticles, &part_type->diameter,&part_type->small_diameter, &part_type->delta_degree, &activity_p);
                part_type->radius=part_type->diameter/2.;
                part_type->small_radius=part_type->small_diameter/2.;

                if (part_type->nsites!=2){
                    sys.dipatch_only=0;
                }

                /*track if all particles have an assigned patch/diameter definition*/
                amount =  part_type->nparticles;
                totalparticles += part_type->nparticles;

                /*read the sites, and make sure these vectors are unit vectors*/
                for(isite=0; isite<part_type->nsites; isite++) {
                    fscanf(file,"%lf %lf %lf\n",
                            &new.x,
                            &new.y,
                            &new.z);
                    part_type->site[isite].r= check_read_unitvec(new);
                }
                if (activity_p>0){
                    part_type->activity=1;
                    part_type->active_force.F0=activity_p;
                    fscanf(file,"%lf %lf %lf\n",
                            &new.x,
                            &new.y,
                            &new.z);
                    part_type->active_force.r=check_read_unitvec(new); 
                }
                do{ /*give the individual particles in the slice a particle type*/
                    psl->pts[ipart].ptype = ptype;
                    // printf("psl->pts[%d].ptype= %d\n",ipart,psl->pts[ipart].ptype);

                    amount--;
                    ipart++;
                }
                while(amount>0 && ipart<sys.npart);
            }
            /* check if all particles have an assigned structure (patches and diameter)*/
            if(totalparticles!=sys.npart){
                error("totalparticles  !=sys.npart  ");
            }
            fclose(file);
        }
    }
    else{
        error("sys.particle_setup can only be 1");
    }
    return;
}

vector check_read_unitvec(vector source){
    /*used when reading in vectors from particle.inp, check if the vector is a finite size (>0).
    it normalizes the vector before saving it*/
    /*if the source vector is super small e.g. (0,0,0), it is probably an error*/
    vector source_new;
    if ( fabs(source.x) <1e-5 && fabs(source.y) <1e-5 && fabs(source.z)<1e-5){
        gprint(source.x);
        gprint(source.y);
        gprint(source.z);
        error("patch site not well defined.");
    }
    /*normalize and return*/
    normvec(source,source_new);

    return source_new;
}

void read_statistics_file(StatsLength *stats_name){
    FILE *file;
    char *pt,line[NPART], filename[100];
    int i=0,j=0, dummy;

    memcpy(filename,stats_name->filename,sizeof(filename));
    // printf("filname %s\n", filename);
    if ((file = fopen(filename,"r"))==NULL){
        printf("%s \n", filename);
        error("input: can't be opened \n");
    }
    else{
        printf("reading %s ...", filename);
    }
   
    while(fgets(line,1000, file) != NULL) {
        sscanf(line,"%d %lf %lf %ld", &dummy, &stats_name->length[i].mean, &stats_name->length[i].variance2,&stats_name->length[i].n);
        i++;
    }
    printf("done\n");
   
    return;
}

void setup_positions_sites(Slice *psl) {

    double r2,s,x=0.,y=0.,z=0.;
    int ipart, jpart, overlap, wall,n;
    vector dr;
    Pts *psi, *psj;
    double radius_i;

    read_particletypes(psl);
    // print_particle_properties();
    // int  fract=floor(sys.npart/100);
    // dprint(fract);

    /*place them random; do not put the particle at z=[0,1] if gravity is on*/
    if(sys.start_type==0){ 
        printf("Randomly placing particles with random orientation\n");
        for( ipart=0; ipart<sys.npart-1; ipart++) {

            psi=&psl->pts[ipart];
            radius_i=sys.particletype[psi->ptype].radius;
            psi->q=RandomQuaternion();
            
            //give it a random position
            do {
                overlap=0;
                if (sys.gravity>0){
                    do{
                        psi->r=RandomVector3D(sys.boxl);
                    }while((psi->r.z<1.12+radius_i ) || (psi->r.z>1.12+2.7*radius_i));
                }
                else{psi->r=RandomVector3D(sys.boxl);}

                // check for overlap
                for( jpart=0; jpart<ipart; jpart++) {

                    psj=&psl->pts[jpart];
                    vector_minus(psj->r,psi->r,dr);
                    pbc(dr,sys.boxl);
                    r2=vector_inp(dr,dr);

                    s=sqrt(r2)-radius_i-sys.particletype[psj->ptype].radius;
                    if(s<sys.s_min*0.9) {
                        // printf("overlap\n");
                        overlap=1;

                    }


                }

            } while(overlap==1);

        }
    }

    /*read from input*/
    if(sys.start_type==1) {
            printf("Reading from conf.inp\n");
            conf_input(&slice[0]);
    }

    /* place them in a chain*/
    if(sys.start_type==2){
            printf("Placing %d particles in a chain, with s_min=%.5lf\n", sys.npart,sys.s_min);
            double y_previousparticle=-0.48*sys.boxl.y;
            if((double)sys.boxl.x/sys.npart<=sys.s_min){
                error("ERROR: Too many particles to put in a 1D chain. Adjust boxl or npart!");
            }
            else{
                for( ipart=0; ipart<sys.npart; ipart++) {
                    psi=&psl->pts[ipart]; 
                    radius_i=sys.particletype[psi->ptype].d_cbsp;
                    // gprint(radius_i);
                    psi->r=nulvec;
                    psi->r.z =1.15+radius_i;
                    psi->r.y = y_previousparticle + (radius_i+sys.s_min) ;                    
                    psi->q.q0=1.;
                    psi->q.q1=0.;
                    psi->q.q2=0.;
                    psi->q.q3=0.;
                    y_previousparticle = psi->r.y+radius_i;
                }
            }  
    }

    if(sys.start_type==4){
        int maxchain_length=floor(sys.npart/sys.nchains)+1,mod_ipart;

        printf("Placing %d particles in %d chains\n", sys.npart,sys.nchains);
        if((double)(sys.boxl.y<=maxchain_length*(sys.s_min+1) )|| (double)(sys.boxl.x)<=sys.chaingap*(sys.nchains)){
            printf("ERROR: Too many particles to put in a 1D chain %lf<%lf or %lf<%lf . Adjust boxl or npart! \n", (double)sys.boxl.y,maxchain_length*(sys.s_min+1),sys.boxl.x,sys.chaingap*(sys.nchains+1));
            exit(1);
        }
        else{
            // printf("%d = %d/ %d\n", sys.npart/sys.nchains,sys.npart,sys.nchains);
            // gprint(sys.s_min);
            n=0;
            // printf("\nnew chain n=%d\n", n);

            for( ipart=0; ipart<sys.npart; ipart++) {
                mod_ipart=ipart%(maxchain_length+1);
                psi=&psl->pts[ipart]; 
                psi->r.x = -0.48*sys.boxl.x + sys.chaingap*n;
                psi->r.z = 1.13+sys.particletype[psi->ptype].d_cbsp;
                psi->r.y = -0.48*sys.boxl.y + (sys.s_min+2.*(sys.particletype[psi->ptype].d_cbsp))*mod_ipart ;   
                // printf("psi->r.y = -0.48*sys.boxl.y + sys.s_min+2*(sys.particletype[psi->ptype].radius)*mod_ipart ;\n" );                
                // printf("psi->r.y = -0.48*s%.5lf + %.5lf+2*(%.5lf)*%d ;\n", sys.boxl.y,sys.s_min,sys.particletype[psi->ptype].radius,mod_ipart);                
                psi->q.q0=1.;
                psi->q.q1=0.;
                psi->q.q2=0.;
                psi->q.q3=0.;
                // vprint(psi->r);
                if(ipart%maxchain_length==0 & ipart>0){
                    n++;
                    // printf("new chain n=%d as ipart modulus sys,nchains =  %d\n",n,ipart%sys.nchains );

                }
                // }
            }  
        }  
    }
    
    // printf("Setting the positions for all the particles is done. now printing\n");
    // for( ipart=0; ipart<sys.npart; ipart++) {
    //     psi=&psl->pts[ipart]; 
    //     // vprint(psi->r);
    // }
    /*setup al the patch vectors*/
    update_patch_vectors(psl);

    if(sys.nearest_neighbor==1){
        /* create the neighbor list here*/
        printf("Setting up the neighbor list\n");
        setup_nnlist();
        printf("making the neighbor\n");
        update_nnlist(psl);
    }
    return;
}

void setup_simulation() {

 
    printf("Setting up the system\n\n");

    printf("initializing the random numbers\n");
    InitializeRandomNumberGenerator(Random_urandom());
    printf("Random_urandom() %d \n", Random_urandom());
    printf("the time(0l) is %ld\n\n", time(0l));

    printf("Allocating memory for the slice\n");
    slice = (Slice *)calloc(1,sizeof(Slice));
    psl_old=(Slice *)calloc(1,sizeof(Slice)); //used in cluster MC
    // memory_allocation();

    printf("\nReading the input\n");
    read_input();
    printf("\nInit the model\n");
    init_model(&slice[0]);
    printf("\nInput read from path.inp:\n");
    print_input();


    //calculates and prints energy, bonding etc
    printf("\n   calculating energy and bonding of initial configuration:\n");
    terminate_block();



    printf("Setting up the system is done\n");

    return;
}

void print_input(void) {
    
    
    printf("Simulation\n");
    if(sys.sim_type==0){
        printf("********* BMD *********\n");
        printf("mobilityT       %.15lf\n", sys.mobilityT);
        printf("mobilityR       %.15lf\n", sys.mobilityR);
        printf("timestep        %.1lf [ns] / %.1lf [ps]\n", langevin.timestep*1e9, langevin.timestep*1e12);
        printf("ninter           %d\n", langevin.ninter);
        printf("nearest_neighbor %d\n", sys.nearest_neighbor);
    }
    else if(sys.sim_type==1){
        printf("********* TIS + BD *********\n");
        printf("mobilityT       %.15lf\n", sys.mobilityT);
        printf("mobilityR       %.15lf\n", sys.mobilityR);
        printf("timestep        %.1lf [ns] / %.1lf [ps]\n", langevin.timestep*1e9, langevin.timestep*1e12);
        printf("ninter           %d\n", langevin.ninter);
        printf("nearest_neighbor %d\n", sys.nearest_neighbor);
    }
    else if(sys.sim_type==2){
        printf("********* MC *********\n");
        printf("cluster_MC %d\n", sys.cluster_MC);
        printf("bond_breakage %d\n", sys.bond_breakage);

    }
    else if(sys.sim_type==4){
        printf("********* FFS *********\n");
        printf("mobilityT       %.15lf\n", sys.mobilityT);
        printf("mobilityR       %.15lf\n", sys.mobilityR);
        printf("timestep        %.1lf [ns] / %.1lf [ps]\n", langevin.timestep*1e9, langevin.timestep*1e12);
        printf("ninter           %d\n", langevin.ninter);
        printf("nearest_neighbor %d\n", sys.nearest_neighbor);
    }

    printf("checks           %d\n", sys.checks);
    // printf("graphics         %d\n", sys.graphics);
    printf("ncycle1          %d\n", sys.ncycle1);
    printf("ncycle2          %d\n", sys.ncycle2);


    printf("\nSystem\n");
    printf("npart             %d\n", sys.npart);
    printf("#particle types   %d\n", sys.nparticle_types);
    printf("beta              %lf\n", sys.beta);
    printf("boxl.x            %lf\n", sys.boxl.x);
    printf("boxl.y            %lf\n", sys.boxl.y);
    printf("boxl.z            %lf\n", sys.boxl.z);

    if(sys.start_type==2 ||sys.start_type==4 ){
        printf("nchains           %d\n", sys.nchains);
        printf("chaingap x        %lf\n", sys.chaingap);
    }

    printf("\nPotential Simon\n");
    printf("rwet              %.3lf\n", sys.r_wetting);
    printf("surface_charge    %.3lf\n", sys.surface_charge);
    printf("dl                %.12lf\n", sys.dl);
    printf("dT                %lf\n", sys.dT);
    printf("xi(dT)            %.18lf\n", sys.xi);
    printf("A(dT)             %.18lf\n", sys.A);
    printf("Arep              %.18lf\n", sys.Arep);
    printf("rcutoff           %.3lf\n", sys.rcutoff);
    printf("gravity           %lf\n", sys.gravity);
    if (sys.wall_int!=0){
        printf("wall_interaction           %lf   if >0 attraction, <0 repulsion; scales equal to gravity\n", sys.wall_int);
    }
   
    printf("\nParticle Properties:\n");
    print_particle_properties();

    printf("lincs          %d\n", sys.lincs);
    if(sys.lincs==1){
        printf("\nConstraints\n");
        printf("cmax              %d \n",cconst.cmax);
        printf("d2                %lf \n",cconst.d2);
        printf("cm                %lf \n",cconst.cm);
        printf("cm2               %lf \n",cconst.cm2);
        printf("K                 %d \n",cconst.K);
        printf("distance          %lf \n", sys.s_min+1.);
    }
    

   
    printf("\nAnalysis\n");
    printf("RDF                %d\n", sys.rdfanalysis);
    printf("xy_print                %d\n", sys.xy_print);
    
    printf("s_histogram        %d\n", sys.s_distribution);
    printf("cluster_analysis   %d\n", cluster.analysis);
    printf("adjacency printing   %d\n",sys.adjacency); 
    printf("bond average       %lf\n", sys.bond_avg );
    printf("bond def. E        %lf\n", sys.bond_cutoffE);

    printf("\n");

    return;
}

void print_particle_properties(){
    int s,n;
    Particletype *part_type;

     for(n=0;n<sys.nparticle_types;n++){
        
        part_type=&sys.particletype[n];
        printf("\n* particle_type %d \n",n);
        printf("%d particle(s) have %d patche(s) with pw=%.2lf degrees: \n", part_type->nparticles, part_type->nsites, part_type->delta_degree);
        printf("and patch vectors :\n" );

        for(s=0; s<part_type->nsites;s++){
            printf("   ( %10.5lf  %10.5lf  %10.5lf ) \n", part_type->site[s].r.x, part_type->site[s].r.y, part_type->site[s].r.z);
            // vprint(part_type->site[s].r);
        }
        if(part_type->activity==1){
            printf("particle has active force with size %lf:\n",part_type->active_force.F0);
            printf("   ( %10.5lf  %10.5lf  %10.5lf ) \n", part_type->active_force.r.x, part_type->active_force.r.y, part_type->active_force.r.z);
        }

        if (sys.S_fixed>0){
            printf("S=Si*Sj is fixed and chosen as %.2lf\n", sys.S_fixed );
        }
        else{
            if (part_type->s_accent==0){
                printf("S old:  ");
                printf("Saccent = 0.5*(1.0-cos(PI*[cosphi-cosphi_particle]/[1-cosphi_particle])) \n");
                printf("        = 0.5*(1.0-cos(PI*[cosphi-%10.5lf]/%10.5lf)) \n",part_type->cosphi, part_type->oneover_cosphi);
            }
            if (part_type->s_accent==1){
                printf("S new:   ");
                printf("Saccent = exp[-cx^2 -dx^3 -ex^4 -fx^5 -gx^6 -hx^7 -ix^8] \n");
                printf("  c                 %.18lf\n", s_int.switch_expc);
                printf("  d                 %.18lf\n", s_int.switch_expd);
                printf("  e                 %.18lf\n", s_int.switch_expe);
                printf("  f                 %.18lf\n", s_int.switch_expf);
                printf("  g                 %.18lf\n", s_int.switch_expg);
                printf("  h                 %.18lf\n", s_int.switch_exph);
                printf("  i                 %.18lf\n", s_int.switch_expi);
                if (s_int.switch_variable==0){
                    printf("x in degrees\n");
                }
                else if (s_int.switch_variable==1){
                    printf("x in cosangle\n");
                }
            }
            if (part_type->s_accent==2){
                printf("S linear:   ");
                printf("Saccent = 1-x/theta_p \n");
                printf("Saccent = 1-x/%5.3lf \n",part_type->delta_degree);
            }
            printf("the threshold of calculating S:\n");
            printf("  cosphi          %lf  \n", part_type->cosphi);
            printf("  cosdelta          %lf  \n", part_type->cosdelta);

            printf("  Derjaguin_prefactor     %lf  \n", part_type->Derj_scaling[n]);
            printf("  delta             %lf  \n", part_type->delta_degree);
            printf("  bulk radius            %lf [sigma]\n", part_type->radius);
            printf("  patch radius           %lf [sigma]\n", part_type->small_radius);
            printf("  patch -bulk distance   %lf [sigma]\n", part_type->d);

            if(sys.gravity>0){
                printf("and gravity parameters:\n");
                printf("  gravity factor    %lf \n", sys.gravity);
                printf("  zcut              %lf [sigma]\n", part_type->zcut);
               
                printf("  Fg                %lf [[kT/sigma]]\n", part_type->fg);
                printf("  b_zc              %lf [sigma]\n", part_type->b_zc);
            }
        }
    }
    printf("switch_method      %d (0 is S0, 1 is S90, 2 is lincomb, 3 is linear)\n\n", sys.switch_method);
    
   return;
}






