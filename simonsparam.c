#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"

/*SIMONS POTENTIAL PARAMETERS*/
//Vrep

// void yukawa_parameters(double radius){
//     double lambda_B, z, dr_center;
//     double eps=25.391640195240345, kb=1.3806e-23,e=1.6022e-19/* e in coulomb*/ ,eps0=8.8542e-12 /*farad/m*/, T=273.15+33.8 /*K*/;
    
//     sys.z=4.*PI*radius*sys.surface_charge*1e18;
//     sys.lambda_B=e*e/(4.*PI*eps0*eps*kb*T);
//     sys.A_yuk=sys.z*sys.z*sys.lambda_B*(1/(1+radius/sys.dl))*(1/(1+radius/sys.dl));
    
//     return;

// }

// double yukawa_rep(double s){
//     /*s is given in units of sigma =3.2 mum*/
//     double  Delta,dr_center;
    

//     dr_center = (s+1.)*3.2;
    
//     return sys.A_yuk*exp(-s/sys.dl)/(dr_center);
// }


double potential_repulsive_energy_sdist(double s){
    /* this is for R2*/
    double erep;


    erep= sys.Arep *exp(-s/sys.dl);
    return erep;
}
//Frep

double bond_repulsive_force(double s){
    // the repulsive force of the bond 
    double Vrep;

    Vrep=potential_repulsive_energy_sdist( s);

    return Vrep/sys.dl;
}
//d^2Vrep/ds^2
double second_der_Vrep(double s){
    /*2nd derivative of Vrep ->> d^2Vrep/ds^2*/
    double Vrep,snd_der;

    Vrep=potential_repulsive_energy_sdist( s);
    snd_der = Vrep/(sys.dl*sys.dl);

    return snd_der;
}


//Vc
double potential_attractive_energy_sdist(double s){
    /* this is for R2*/
    double Uatr;

    Uatr = -(1/sys.xi)*exp(-sys.A*(s/sys.xi)-sys.B*((s*s)/(sys.xi*sys.xi)));

    return Uatr;
}

//Fc
double bond_attractive_force(double s){
    // the attractive force of the bond 
    double Vc;

    Vc=potential_attractive_energy_sdist( s);

    return 2.*s/(sys.xi*sys.xi)*Vc;
}

//d^2Vc/ds^2
double second_der_Vc(double s){
    /*2nd derivative of Vc ->> d^2Vc/ds^2*/
    double Fc, Vc,snd_der, s2_xi2;

    
    s2_xi2= (s*s)/(sys.xi*sys.xi);

    Fc=bond_attractive_force(s);
    Vc=potential_attractive_energy_sdist( s);
    // snd_der = -2.*Vc/(sys.xi*sys.xi*sys.xi) + (2.*s)/sys.xi*(1.*Fc);
    snd_der = Vc * (4.*s2_xi2-2./(sys.xi*sys.xi));

    return snd_der;
}

void setup_Simons_potential_parameters(){
    /*determine here A and xi. You choose :
    choose dt only as:  0.12 upto 0.22
    choose r_wetting only as:  0.40 upto 0.56
    choose surface_charge only as:  -0.10 upto -0.38

    these parameters are based on a system with:
    salt concentration:             1mM MgSO4   = 4.0 
    volume fractie lutidine =       25%vol      
    DP-A particles met diameter=    3.2 micron 
    patch curvature R=              1.0 micron      
        */
    double dt,rwet,sc;
    double A,xi,Arep,Apart1,Apart2,Bpart1,Bpart2,xipart1,xipart2;

    printf("    SIMONS POTENTIAL PARAMETERS:\nthese parameters are based on a system with:\n");
    printf("salt concentration:             1mM MgSO4   = 4.0 \n");
    printf("volume fractie lutidine =       25\%%vol      \n");
    printf("DP-A particles met diameter=    3.2 micron \n");
    printf("patch curvature R=              1.0 micron  \n");

    printf("these parameters are fitted between:\n");
    printf("rwet in [0.30,0.60]\n");
    printf("sc   in [-0.50,-0.01] \n");
    printf("dt   in [0.10,0.22]\n");

    dt=sys.dT;
    rwet=sys.r_wetting;
    sc=sys.surface_charge;

    if(dt<0.02 || dt>0.22){
        printf("sys.dT= %lf  \n", dt);
        error(" choose only dt between 0.10 and 0.22  ");
    }
    if(rwet<0.30 || rwet>0.60){
        printf("sys.r_wetting= %lf  \n", rwet);
        error(" choose only rwet between 0.30 and 0.60");
    }
    if(sc>-0.01 || sc<-0.50){
        gprint(sys.surface_charge);
        error("choose surface_charge between -0.50 and -0.01\n");
    }
    
    sys.dl= 0.000750174488696; // 0.0008673892526; // 0.002775645608174734; //  // = 2.775645608174734e-03 micron/ 3.2 micron

    /*from files: 
    fit_A_poly_4x4dt_patchrad1.00mum
    fit_xi_poly_4x4dt_patchrad1.00mum.txt
    fit_Arep_poly_3x0dt_patchrad1.00mum.txt*/

    /* sys.sigma_mum=1.0; //[mum]
    Apart1=0.056782502496567314  +rwet*0.260215676538175344  +rwet*rwet*-2.087133590729501886  +rwet*rwet*rwet*14.759587560733907097  ;    
    Apart2=0.011502266054040660  +dt*5.278593996230584118  +dt*dt*-29.955521677043037698  +dt*dt*dt*51.198280765408576087  ;
    sys.A=Apart1*Apart2;//[dimensionless]

    xipart1=0.004056363352597567  +rwet*-0.079426720163945003  +rwet*rwet*0.031859894515664684  +rwet*rwet*rwet*-0.011776312522848912  ;   
    xipart2=-0.486451801907159032  +dt*1.964137077233496065  +dt*dt*-3.757782035401259879  +dt*dt*dt*0.711773304860235934 ;
    sys.xi=xipart2*xipart1; //   [dimensionless]
    
    sys.Arep=0.000084830341705028  +sc*0.002376164488813220  +sc*sc*1295462.677729546325281262    ;*/

    /*in directory update_II_200826_1.0mum*/
    // Apart1=-0.015061862283229362  +rwet*2.179971050601182458  +rwet*rwet*-9.919260032900107049  +rwet*rwet*rwet*51.479504724036686980  ;
    // Apart2=-0.226173154804393844  +dt*7.246060901785060793  +dt*dt*-60.600308320350059432  +dt*dt*dt*221.602339052154405863  +dt*dt*dt*dt*-304.832948001714555630  ;

    // sys.A=Apart1*Apart2;// 

    // xipart1=0.219318405905894714  +rwet*-6.415286942580221030  +rwet*rwet*0.684058978296809062  +rwet*rwet*rwet*0.685914887127846429  ;
    // xipart2=-0.007218874637425805  +dt*0.071211357421138557  +dt*dt*-0.501629147203333958  +dt*dt*dt*1.870137054939489119  +dt*dt*dt*dt*-2.783413465699470013  ;

    // sys.xi=xipart2*xipart1;
    
    // sys.Arep=-0.000060954325068766  +sc*-0.000689805950692273  +sc*sc*1295484.309107528068125248  ;

    /*in directory update_II_200826_1.0mum*/
    Apart1=78.19914394228547  +rwet*-140.0594240946654  +rwet*rwet*-243.59516521844142  +rwet*rwet*rwet*229.37356574989568 +rwet*rwet*rwet*rwet*1369.7626413407813 +rwet*rwet*rwet*rwet*rwet*-1676.854890665089;
    Apart2=1.8963340474698267  +dt*-92.2592429895963  +dt*dt*1569.7394089464121  +dt*dt*dt*5563.972698891612  +dt*dt*dt*dt*-488619.8292050912 +dt*dt*dt*dt*dt*4268201.831152699;

    // sys.A=Apart1*Apart2;
    sys.A=2.33836617;

    Bpart1=48.83809400200591  +rwet*-287.0119049995045  +rwet*rwet*568.9030139284323  +rwet*rwet*rwet*-379.4007580217976  ;
    Bpart2=-174.91845934764606  +dt*22438.195763523367  +dt*dt*-867329.5791137232  +dt*dt*dt*15215741.64998858  +dt*dt*dt*dt*-95902211.79087652  ;

    // sys.B=Bpart1*Bpart2;
    sys.B=6.59350900;

    xipart1=1.5426239734853422  +rwet*-3.747131745035845  +rwet*rwet*-4.799611306922257  +rwet*rwet*rwet*20.603788487585664 +rwet*rwet*rwet*rwet*-12.568351277865146 +rwet*rwet*rwet*rwet*rwet*-4.03778471787294;
    xipart2=0.5004919117123887  +dt*-10.832915022158021  +dt*dt*434.1108941227658  +dt*dt*dt*-15641.039533469375  +dt*dt*dt*dt*295784.0017711761 +dt*dt*dt*dt*dt*-1993035.9130637303;

    // sys.xi=xipart2*xipart1;
    sys.xi=0.03777731;
    
    // sys.Arep=0.0001852178  +sc*-0.0043217593  +sc*sc*643757.4642065160  ;
    sys.Arep=5214.366733;    
    // sys.Arep=-0.000060954325068766  +sc*-0.000689805950692273  +sc*sc*1295484.309107528068125248  ;

    // sys.xi*=(3.2) ;//to make it into units of mum
    // sys.A*=(3.2) ;// to make it into units of mum
    // sys.dl*=3.2;


    gprint(sys.A);
    gprint(sys.xi);
    gprint(sys.Arep);
    return;


}
