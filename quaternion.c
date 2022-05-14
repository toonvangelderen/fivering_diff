#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"

// vector and matrix operations
void rotate_quaternion_ipart( Slice *psl, int ipart, quaternion dq ){
    // function that rotates the quaternion and updates the patch vectors of particle ipart 

    quaternion dqrot;
    
    // very strange; quat_times needs to go via dqrot, and cannot be directly via psl->pts[ipart].q..
    quat_times(dq,psl->pts[ipart].q,dqrot);
    psl->pts[ipart].q = dqrot;
    update_patch_vector_ipart(psl,ipart); 
    return;
}

double langrange_multiplier_quat(quaternion qprime, quaternion q) {

    double lambdaq,a,b,c,det,root1,root2;

    //solving for equation 15 Ilie Briels den Otter
    //lambdaq^2 + 2*lambdaq*qprime*q + qprime*qprime == 1

    a = 1.0;
    b = 2.0*quat_inp(q,qprime);
    c = quat_inp(qprime,qprime)-1.0;

    det = b*b - 4.0*a*c;
    if(det<0.0) {
        //printf("Warning: determinant lower than 0, can not find roots\nWhat to do?\n");
        //printf("b value %lf c value %lf\n",b,c);
        return BIGNUM;
    }
    else {
        //hmmm which root to choose?
        root1=(-b + sqrt(det))/(2.0*a);
        root2=(b + sqrt(det))/(2.0*a);
    }

    //go for minimum for now...makes more sense
    //if this does not work try just renormalizing it...
    
    if(fabs(root1)<fabs(root2)) {
        lambdaq=root1;
    }
    else {
        lambdaq=root2;
    }

    //printf("lambdaq %lf\n",lambdaq);
    return lambdaq;
}

quaternion quatVecToVec(vector vec1, vector vec2) {

    quaternion qrot;
    vector rotvec,xvec,yvec;
    double angle,s,invs,l;
    //unit cos(angle)
    angle=vector_inp(vec1,vec2);
    if(angle>=1) {
        qrot.q0=1;
        qrot.q1=0;
        qrot.q2=0;
        qrot.q3=0;
    }

    else if(angle < -0.999999) {
        //xvec.x=1;
        //xvec.y=0;
        //xvec.z=0;
        //vector_cross(xvec,vec2,rotvec);   
        //l=sqrt(vector_inp(rotvec,rotvec));
        //if(l==0) {
        //  yvec.x=0;
        //  yvec.y=1;
        //  yvec.z=0;
        //  vector_cross(yvec,vec2,rotvec);
        //  l=sqrt(vector_inp(rotvec,rotvec));
        //}
        //scalar_divide(rotvec,l,rotvec);
        //qrot.r=PI;
        //qrot.u=rotvec;
        //scdivide_quat(qrot,sqrt(quat_inp(qrot,qrot)),qrot);
        qrot.q0=0;
        qrot.q1=0;
        qrot.q2=1;
        qrot.q3=0;
    }

    else {
        s=sqrt(2.*(1.+angle));  
        invs=1./s;
        vector_cross(vec1,vec2,rotvec);
        qrot.q0=0.5*s;
        scalar_times(rotvec,invs,rotvec);
        qrot.q1=rotvec.x;
        qrot.q2=rotvec.y;
        qrot.q3=rotvec.z;
        scdivide_quat(qrot,sqrt(quat_inp(qrot,qrot)),qrot);
    }
                

    return qrot;
}

tensor getrotmatrix(quaternion q) {

    tensor mat;
    double q0,q1,q2,q3,q0sq,q1sq,q2sq,q3sq;

    q0=q.q0;
    q1=q.q1;
    q2=q.q2;
    q3=q.q3;

    q0sq=q0*q0;
    q1sq=q1*q1;
    q2sq=q2*q2;
    q3sq=q3*q3;

    mat.x.x = q0sq + q1sq - q2sq - q3sq ;  mat.x.y = 2.0*(q1*q2 - q0*q3)       ;  mat.x.z = 2.0*(q1*q3 + q0*q2)       ;
    mat.y.x = 2.0*(q1*q2 + q0*q3)       ;  mat.y.y = q0sq - q1sq + q2sq - q3sq ;  mat.y.z = 2.0*(q2*q3 - q0*q1)       ;
    mat.z.x = 2.0*(q1*q3 - q0*q2)       ;  mat.z.y = 2.0*(q2*q3 + q0*q1)       ;  mat.z.z = q0sq - q1sq - q2sq + q3sq ; 

    return mat;
}

quattensor getquatmatrix(quaternion q) {
    /* eqn 14 of Illie2015, but without the first column with is q0 q1 q2 q3 so only the right part */
    quattensor qmat;
    double q0,q1,q2,q3;

    q0=0.5*q.q0;
    q1=0.5*q.q1;
    q2=0.5*q.q2;
    q3=0.5*q.q3;

    qmat.w.x = -q1 ;  qmat.w.y = -q2 ;  qmat.w.z = -q3 ;
    qmat.x.x =  q0 ;  qmat.x.y = -q3 ;  qmat.x.z =  q2 ;
    qmat.y.x =  q3 ;  qmat.y.y =  q0 ;  qmat.y.z = -q1 ;
    qmat.z.x = -q2 ;  qmat.z.y =  q1 ;  qmat.z.z =  q0 ;

    return qmat;
}

quaternion RandomQuaternion(void) {

    double sigma1, sigma2, s, theta1, theta2;
    quaternion q;

    s = RandomNumber();
    sigma1 = sqrt(1.0-s);
    sigma2 = sqrt(s);
    theta1 = 2.0*PI*RandomNumber();
    theta2 = 2.0*PI*RandomNumber();

    q.q0 = sigma2 * cos(theta2);
    q.q1 = sigma1 * sin(theta1);
    q.q2 = sigma1 * cos(theta1);
    q.q3 = sigma2 * sin(theta2);

    return q;
}


quaternion RandomQuaternionRange(double dqmax) {
    /*see  https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation for explination*/
    quaternion q;
    double alpha,sina;
    vector rotaxis;

    alpha = dqmax*RandomNumber();

    rotaxis = RandomUnitVector();

    sina = sin(alpha/2.0);

    q.q0 = cos(alpha/2.0);
    q.q1 = sina*rotaxis.x;
    q.q2 = sina*rotaxis.y;
    q.q3 = sina*rotaxis.z;

    return q;
}

quaternion QuaternionXaxis( double degree) {

    quaternion q;
    double alpha,sina;

    alpha = degree/180.*PI;

    sina = sin(alpha/2.0);

    q.q0 = cos(alpha/2.0);
    q.q1 = sina;
    q.q2 = 0.;
    q.q3 = 0.;

    return q;
}

quaternion QuaternionYaxis( double degree) {

    quaternion q;
    double alpha,sina;

    alpha = degree/180.*PI;

    sina = sin(alpha/2.0);
    
    q.q0 = cos(alpha/2.0);
    q.q1 = 0.;
    q.q2 = sina;
    q.q3 = 0.;

    return q;
}

quaternion QuaternionZaxis( double degree) {

    quaternion q;
    double alpha,sina;

    alpha = degree/180.*PI;

    sina = sin(alpha/2.0);

    q.q0 = cos(alpha/2.0);
    q.q1 = 0.;
    q.q2 = 0.;
    q.q3 = sina;

    return q;
}