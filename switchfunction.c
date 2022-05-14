#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"


double S_value(Slice *psl, double cositheta, double cosjtheta, vector p1, vector p2, vector u1 , int ipart, int jpart){
    /* the definition of S and S' define the effective switch funciton S;
    S =max(S'i,S'j) = S0        option 0 of sys.switch_method
    S =S'i*S'j  = S90           option 1
    S = lincomb SO, S90, S180   option 2 
    S' is old switch            option 0 of particletype[ptype].s_accent;
    S' is integrated switch     option 1*/

    double Si, Sj, S;
    
    if (sys.S_fixed>0){
        S=sys.S_fixed;
        return S;
    }

    // first calc Si' and S'j
    Si = Saccent( cositheta,  psl->pts[ipart].ptype);
    Sj = Saccent( cosjtheta,  psl->pts[jpart].ptype);

    // printf("The values of Si=%.5lf and Sj=%.5lf \n", Si,Sj);
    // then create S
    if (sys.switch_method==0){
        S=S0(Si,Sj);
    }
    else if(sys.switch_method==1){
        S=S90(Si,Sj);
    }
    else if(sys.switch_method==2){ //lin comb including S180
        S=switch_method_2(Si,Sj,cositheta,cosjtheta,u1,p1,p2);
    }
    else{
        error("choose switch_method = 0 (S0), 1(S90), 2(lincomb), 3(S fixed)");
    }

    if((S<-1e-1) || (S-1.>1e-3 && S>1.) || ((S/S!=1) && S>0) ){
        printf("WARNING: S is either S<-0.10 or S>1.001 or S=Nan\n");
        printf(" Si        = %20.18lf\n", Si);
        printf(" Sj        = %20.18lf\n", Sj);
        printf(" S         = %20.18lf\n", S);
        printf(" cositheta = %20.18lf\n", cositheta);
        printf(" cosjtheta = %20.18lf\n", cosjtheta);

        printf(" S0(Si,Sj) = %20.18lf\n", S0(Si,Sj));
        printf(" S90(Si,Sj) = %20.18lf\n", S90(Si,Sj));
        printf(" S180(Si,Sj) = %20.18lf\n", S180(Si,Sj));

    }
    return S;
}

double Saccent(double costheta, int ptype ){
    double Saccent;

    if (sys.S_fixed>0){
        Saccent=sys.S_fixed;
        return Saccent;
    }

    if( sys.particletype[ptype].s_accent==0){
        // printf("CALC old_S\n"); the "old" one from the naterue paper
        Saccent=old_S(costheta,ptype);
    }
    else if (sys.particletype[ptype].s_accent==1){
        // printf("CALC new_S\n"); //the integrated one 
        Saccent=new_S(costheta);     
    }
    else if (sys.particletype[ptype].s_accent==2){
        // printf("CALC new_S\n"); the linear one (used in wertheim)
        Saccent=S_lin(costheta,ptype);     
    }
    else{
        error("particletype[ptype].s_accent should be 0 (for Sold) and 1(for Snew) and 2 for linear S");
    }
    return Saccent;

}

double S_lin(double costheta, int ptype){
    double theta_p=sys.particletype[ptype].delta_degree,  Slin, angle;
    
    angle=(costheta<1.)? acos(costheta)*180./PI:0;
    Slin = 1.-angle/theta_p;
    return Slin;
}


double old_S(double costheta, int ptype){
    double Sold, cosdelta_particle = sys.particletype[ptype].cosdelta,oneover_cosdelta = sys.particletype[ptype].oneover_cosdelta;

    Sold=0.5*(1.0-cos(PI*(costheta-cosdelta_particle)*oneover_cosdelta));
    return Sold;
}

double new_S(double cosangle){
    /* switch_expbcd is a switch function as integrated with the patch integration scheme,
    and fitted to a curve exp(cx^2 + dx^3 + ex^4 + fx^5 + gx^6) where x =1-cos(theta) and theta the angle
    the first parameter in the exponent is c, because, a=0 (such that f(0)=1) and b=0, because f'(0)=0
    this is to make the switching function continuous in the fist derivative */
    double angle1, angle2, angle4;
    double cx2, dx3, ex4, fx5,gx6,hx7,ix8;
    double a, S;

    if(s_int.switch_variable==0){
        //use degree
        if (cosangle<=-1){
            angle1=0;
        }
        angle1=(cosangle<1.)? acos(cosangle)*180./PI:0;
    }
    else if (s_int.switch_variable==1){
        //use cosangel
        angle1=1.-cosangle;
        
    }
    angle2=angle1*angle1;
    angle4=angle2*angle2;

    cx2=angle2* s_int.switch_expc;
    dx3=angle2*angle1* s_int.switch_expd;
    ex4=angle4* s_int.switch_expe;
    fx5=angle4*angle1* s_int.switch_expf;
    gx6=angle4*angle2* s_int.switch_expg;
    hx7=angle4*angle2*angle1* s_int.switch_exph;
    ix8=angle4*angle4* s_int.switch_expi;
    a=0.-cx2-dx3-ex4-fx5-gx6-hx7-ix8;
    S=exp(a);
    if (S>0 & S/S!=1){
        gprint(a);
        gprint(cx2);
        gprint(dx3);
        gprint(ex4);
        gprint(fx5);
        gprint(gx6);
        gprint(hx7);
        gprint(ix8);
        gprint(cosangle);
        gprint(angle1);
        dprint(s_int.switch_variable);
        gprint(S);
    }

    return S;
}

double S0(double Si,double Sj){
    //kleinste waarde voor S bepaalt
    return min(Si,Sj);
}

double S90(double Si,double Sj){
    return Si*Sj;
}

double S180(double Si,double Sj){
    double S;
   
    S=S90(Si,Sj) +0.7*(S90(Si,Sj)-S0(Si,Sj));
    return S;
}

double switch_method_2(double Si, double Sj, double cositheta, double cosjtheta, vector u1, vector p1, vector p2){
    /*the linear combination method*/
    double angle1, angle2, angle3,S,costheta3;
    double lambda1, lambda2;
    vector u2,proj_u1_p1,u3,proj_u1_p2,e1,e2;

    // include theta3 for a linear combination of cube S0, sphere S90, and the combi S180.
        /*theta3 [0,180], subtract projection of rij of p1 and p2 and take vector_inp*/

    if(cositheta==1.0 || cosjtheta==1.0){
        // if one of the angles is zero, the rotation around angle3 is invariant. Thus, angle3 is arbitraty 
        S=max(Si,Sj);
    }
    else{
        // error("do not use switch method3")
        /*do gram schmidt orthonormalization to find the projections*/
        /*u1=rnorm*/

        // projection p1 onto u1, and create (unit vector) e2 the perpendicular projection of p1 
        scalar_times(u1,-cositheta,proj_u1_p1); 
        vector_minus(p1,proj_u1_p1,u2);
        scalar_divide(u2,sqrt(vector_inp(u2,u2)),e1);
        // projection p2 onto u1, and create (unit vector) e3 the perpendicular projection of p2 
        scalar_times(u1,cosjtheta,proj_u1_p2); 
        vector_minus(p2,proj_u1_p2,u3);
        scalar_divide(u3,sqrt(vector_inp(u3,u3)),e2);

        /*calc angle3={0,180}, you need it in angles, because the interpolation needs angles*/
        costheta3=vector_inp(e1,e2);

        /* angle3=nan if costheta3>1 or costheta3<-1 (can sometimes happen)*/
        if((costheta3>=1) || (costheta3<=-1)){
                if(costheta3>=1.) angle3=0.;
                if(costheta3<=-1.) angle3=180.;
        }
        else{   angle3=acos(costheta3)*180./PI;
        }

        
        if(angle3>=0. && angle3<=90.0){
            lambda1=angle3/90.;
            S=lambda1*S90(Si,Sj)+(1.-lambda1)*S0(Si,Sj);

        }
        else if(angle3>90.0 && angle3<=180.0){
        // else if(costheta3<0 && costheta3>=-1){
            lambda2=(angle3-90.)/90.;
            S=(1.-lambda2)*S90(Si,Sj)+lambda2*S180(Si,Sj);
        }
        else{
            printf("****esle: \n");
            printf("cositheta = %.18lf \n", cositheta);
            printf("cosjtheta = %.18lf \n", cosjtheta);
            printf("costheta3 = %.18lf \n", costheta3);
            gprint(angle3);

            error("angle3 is bigger than 180??");
        }
    }
    return S;
}

double dS_dcostheta(Slice *psl, double cositheta, int ipart, double cosjtheta, int jpart ){
    /*write the dS/dcostheta there for various switch function 
    take the derivative to cositheta

    S(i,j) or S(i) something else
    dS/dcosi = S(j)* dS(i)/dcosi*/

    double dS_dcostheta;
    int ptypei, ptypej;
    ptypei=psl->pts[ipart].ptype;
    ptypej=psl->pts[jpart].ptype;


    if (sys.S_fixed>0){
        dS_dcostheta=0;
        return dS_dcostheta;
    }
    
    if (sys.switch_method==0){
        dS_dcostheta=dS0_dcostheta(cositheta,ptypei,cosjtheta,ptypej);  
    }
    else if(sys.switch_method==1){
        dS_dcostheta=dS90_dcostheta(cositheta,ptypei,cosjtheta,ptypej);
    }
    else if(sys.switch_method==2){
        error("im not sure how the derivative of this one goes. Im will probably not use it, make it if  you want to.");
    }
    else{
        error("choose switch_method = 0,1,2");
    }


    return dS_dcostheta;
}

double dSaccent_dcosangle(double cosangle, int ptype){
    /* the derivative wrt to cosangle of Saccent 
    S'=exp(a(x))
    dS'/dcostheta   = dS'/da * da/dx * dx/dcostheta = 
                    = S' * b * dx/dcostheta 
    with a=0-cx2-dx3-ex4-fx5-gx6-hx7-ix8 and x dependent on choice of switch_varialble*/

    double angle1,angle2,angle4;
    double b,dSaccent_dcosangle,dx_dcostheta;
    double cx1_2, dx2_3, ex3_4, fx4_5, gx5_6,hx6_7,ix7_8;

    // make distinction between particletype
    if (sys.particletype[ptype].s_accent==0){
        dSaccent_dcosangle=dSold_dcostheta(cosangle,ptype);
    }
    else if (sys.particletype[ptype].s_accent==1){
        dSaccent_dcosangle= dSnew_dcostheta(cosangle, ptype);
    }
    else if (sys.particletype[ptype].s_accent==2){
        error("if s_accent=2 the switchfunction is linear and discontinuous");
    }

    

    return dSaccent_dcosangle;
}

double dSnew_dcostheta(double cosangle, int ptype){
    /* the derivative wrt to cosangle of Saccent 
    S'=exp(a(x))
    dS'/dcostheta   = dS'/da * da/dx * dx/dcostheta = 
                    = S' * b * dx/dcostheta 
    with a=0-cx2-dx3-ex4-fx5-gx6-hx7-ix8 and x dependent on choice of switch_varialble*/

    double angle1,angle2,angle4;
    double b,dSaccent_dcosangle,dx_dcostheta;
    double cx1_2, dx2_3, ex3_4, fx4_5, gx5_6,hx6_7,ix7_8;

    if(s_int.switch_variable==0){
        if (cosangle<=-1){
            angle1=180.;
        }
        if (cosangle>=1.){           
            angle1=0.;
            dx_dcostheta=0;
        }
        else{                        
            angle1=(cosangle<1.)? acos(cosangle)*180./PI:0;
            dx_dcostheta = -1./sin(angle1*PI/180.)*(180./PI);
        }

        /*dx/dcostheta= dx/dtheta * dtheta/dcostheta = dx/dtheta * (dcostheta/dtheta)^(-1)
        x=theta */
    }
    else if (s_int.switch_variable==1){
        angle1=1.-cosangle;
        //x=1-costheta so dx/dcostheta = -1
        dx_dcostheta=-1.;
        
    }

    angle2=angle1*angle1;
    angle4=angle2*angle2;

    /* b=da/dx  */
    cx1_2=2.* s_int.switch_expc*angle1;
    dx2_3=3.* s_int.switch_expd*angle2;
    ex3_4=4.* s_int.switch_expe*angle2*angle1;
    fx4_5=5.* s_int.switch_expf*angle4;
    gx5_6=6.* s_int.switch_expg*angle4*angle1;
    hx6_7=7.* s_int.switch_exph*angle4*angle2;
    ix7_8=8.* s_int.switch_expi*angle4*angle2*angle1;

    b = 0.-cx1_2-dx2_3-ex3_4-fx4_5-gx5_6-hx6_7-ix7_8;
    
    dSaccent_dcosangle = Saccent(cosangle, ptype) * b * dx_dcostheta;

    // gprint(Saccent(cosangle, ptype));
    // gprint(b);
    // gprint(dx_dcostheta);
    return dSaccent_dcosangle;
}

double dSold_dcostheta(double costheta , int ptype){
    double dSold_dcostheta;
    double cosdelta_particle = sys.particletype[ptype].cosdelta,oneover_cosdelta = sys.particletype[ptype].oneover_cosdelta;

    dSold_dcostheta=HALFPI*sin(PI*(costheta-cosdelta_particle)*oneover_cosdelta)*oneover_cosdelta;
    // gprint(dSold_dcostheta);
    // gprint(cosdelta_particle);
    // gprint(oneover_cosdelta);
    return dSold_dcostheta;
}

double dS0_dcostheta(double cosangle1, int ptype1, double cosangle2, int ptype2){
    /*derivative of S0 wrt cosangle1 */
    double dS0_dcostheta;

    /*check which angle is biggest (i.e. cosangle smallest) and use that one as theta*/
    if(cosangle2<cosangle1){
        dS0_dcostheta=0;
    }
    else{
        dS0_dcostheta=dSaccent_dcosangle(cosangle1, ptype1)*Saccent(cosangle2,ptype2);
    }
    
    return dS0_dcostheta;
}

double dS180_dcostheta(double cosangle1, int ptype1, double cosangle2, int ptype2){
    double dS180_dcostheta,dS90=dS90_dcostheta(cosangle1,ptype1,cosangle2,ptype2);
    double dS0 = dS0_dcostheta(cosangle1,ptype1,cosangle2,ptype2);
   
    dS180_dcostheta=dS90 +0.7*(dS90-dS0);
    return dS180_dcostheta;
}

double dS90_dcostheta(double cosangle1, int ptype1, double cosangle2, int ptype2){
    /*derivative of S90 wrt cosangle1
    with S90 = S'(cosangle1) * S'(cosangle2)  */
    double dS90_dcostheta;

    /*dS90_dcostheta = dS'/dcostheta1 *S'(costheta2)*/
    dS90_dcostheta = dSaccent_dcosangle(cosangle1,ptype1)*Saccent(cosangle2,ptype2);
    
    return dS90_dcostheta;
}
