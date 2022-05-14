#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>


/*------------------INLINE SUBSTITUTIONS-------------------------------------*/

//Binary definitions
#define TRUE 1
#define FALSE 0
#define CANCELLED -99

//not so necessary defitions
#define P_INDEX(ptr)  (ptr - pts)
#define F_COMP(a,b) (fabs( a-b ) < 1e-10)
#define MIN(a,b,tmp)    (  (tmp=a) < b  ? tmp : b )
#define MAX(a,b,tmp)    (  (tmp=a) > b  ? tmp : b )
#define MESSAGE(a) printf("Message:"#a"\n")
#define SIGN(a)   ( 2*(a>0) -1) 

#define max(a,b) \
    ({ __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })
#define min(a,b) \
    ({ __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })

//Math shorthands + special numbers
#define PI           3.141592653589793
#define HALFPI       1.5707963267948966
#define SQRT2        1.4142136
#define BIGNUM       1e99
#define EPS          1e-10
#define NVT          1
#define cubic(x) ((x)*(x)*(x))
#define square(x) ((x)*(x))

//Print shorthands
#define dprint(expr) printf("  "#expr " = %d\n",expr)
#define gprint(expr) printf("  "#expr " = %20.15g\n",expr)
#define vprint(expr) printf("  "#expr " = ( %16.14g %16.14g %16.14g ) \n",expr.x, expr.y, expr.z)
#define getName(var)  #var

//Vector definitions
#define vector_inp(a, b) (a.x * b.x + a.y * b.y + a.z * b.z )
#define vector_cross(a,b,h)   h.x = a.y * b.z - a.z * b.y; h.y = a.z * b.x - a.x * b.z; h.z = a.x * b.y - a.y * b.x;
#define vector_times(a,b,h)   h.x = a.x * b.x; h.y = a.y * b.y; h.z = a.z * b.z;
#define vector_divide(a,b,h)   h.x = a.x / b.x; h.y = a.y / b.y; h.z = a.z / b.z;
#define vector_plustimes(a,b,h)   h.x += a.x * b.x; h.y += a.y * b.y; h.z += a.z * b.z;
#define vector_mintimes(a,b,h)   h.x -= a.x * b.x; h.y -= a.y * b.y; h.z -= a.z * b.z;
#define scalar_times(a,b,h)  h.x = a.x *b; h.y =a.y * b; h.z = a.z * b;  
#define scalar_divide(a,b,h)  h.x = a.x /b; h.y =a.y / b; h.z = a.z / b;  
#define scalar_plustimes(a,b,h)  h.x += a.x *b; h.y +=a.y * b; h.z += a.z * b;  
#define scalar_mintimes(a,b,h)  h.x -= a.x *b; h.y -=a.y * b; h.z -= a.z * b;  
#define vector_add(a,b,h)  {h.x = a.x +b.x; h.y =a.y + b.y; h.z = a.z + b.z;  }
#define vector_minus(a,b,h) { h.x = a.x - b.x; h.y =a.y - b.y; h.z = a.z - b.z;  }
#define matrix_x_vector(m,a,h) {\
    h.x = m.x.x*a.x + m.x.y*a.y + m.x.z*a.z;\
    h.y = m.y.x*a.x + m.y.y*a.y + m.y.z*a.z;\
    h.z = m.z.x*a.x + m.z.y*a.y + m.z.z*a.z;}
#define matrixT_x_vector(m,a,h) {\
    h.x = m.x.x*a.x + m.y.x*a.y + m.z.x*a.z;\
    h.y = m.x.y*a.x + m.y.y*a.y + m.z.y*a.z;\
    h.z = m.x.z*a.x + m.y.z*a.y + m.z.z*a.z;}
#define normvec(a,b) {b.x=a.x/sqrt(vector_inp(a,a)); b.y=a.y/sqrt(vector_inp(a,a)); b.z=a.z/sqrt(vector_inp(a,a));}

//Quaternion definitions
#define quat_add(a,b,h) { h.q0=a.q0+b.q0; h.q1=a.q1+b.q1; h.q2=a.q2+b.q2; h.q3=a.q3+b.q3; }
#define quat_minus(a,b,h) { h.q0=a.q0-b.q0; h.q1=a.q1-b.q1; h.q2=a.q2-b.q2; h.q3=a.q3-b.q3;}
#define quat_inp(a,b) (a.q0*b.q0 + a.q1*b.q1 + a.q2*b.q2 + a.q3*b.q3)
#define sctimes_quat(a,b,h) { h.q0=a.q0*b; h.q1=a.q1*b; h.q2=a.q2*b; h.q3=a.q3*b;}
#define scdivide_quat(a,b,h) { h.q0=a.q0/b; h.q2=a.q2/b; h.q3=a.q3/b; h.q3=a.q3/b;}
#define quat_times(a,b,h) {\
    h.q0 = a.q0*b.q0 - a.q1*b.q1 - a.q2*b.q2 - a.q3*b.q3;\
    h.q1 = a.q0*b.q1 + a.q1*b.q0 + a.q2*b.q3 - a.q3*b.q2;\
    h.q2 = a.q0*b.q2 - a.q1*b.q3 + a.q2*b.q0 + a.q3*b.q1;\
    h.q3 = a.q0*b.q3 + a.q1*b.q2 - a.q2*b.q1 + a.q3*b.q0;}
#define quat_inverse(a,h) {h.q0 = a.q0; h.q1 = -a.q1; h.q2 = -a.q2; h.q3 = -a.q3;}
#define quatmatrix_x_vec(m,a,h) {\
    h.q0 = m.w.x*a.x + m.w.y*a.y + m.w.z*a.z;\
    h.q1 = m.x.x*a.x + m.x.y*a.y + m.x.z*a.z;\
    h.q2 = m.y.x*a.x + m.y.y*a.y + m.y.z*a.z;\
    h.q3 = m.z.x*a.x + m.z.y*a.y + m.z.z*a.z;}

//Update definitions
#define update_average(a,b) {a.now =b; a.sum += a.now; a.sumsq += a.now*a.now; a.n++;}
#define update_blockaver(a) {a.sum /= a.n; a.sumsq = sqrt(a.sumsq/a.n - a.sum*a.sum);}
#define update_finaver(a,b) {a.sum += b.sum; a.sumsq += b.sum*b.sum;  b.sum = 0; b.sumsq = 0; b.now=0; b.n=0;}
#define update_blockmc(a,b) {a.acc += b.acc; a.tries += b.tries;  b.acc = 0; b.tries = 0;}


#define SWAP(a,b,c) {c=a;a=b;b=c;}

//Periodic Boundary conditions
#define pbc(a,b) {\
    if(a.x>(0.5*b.x)) a.x-=b.x;\
    if(a.x<(-0.5*b.x)) a.x+=b.x;\
    if(a.y>(0.5*b.y)) a.y-=b.y;\
    if(a.y<(-0.5*b.y)) a.y+=b.y;\
    if(a.z>(0.5*b.z)) a.z-=b.z;\
    if(a.z<(-0.5*b.z)) a.z+=b.z;}

    //a.x -= b.x*rint(a.x/b.x);\
    //a.y -= b.y*rint(a.y/b.y);\
    //a.z -= b.z*rint(a.z/b.z);}

#define NPART 1000
#define NSITES 8
#define PTYPES 5  //Site type (patchy particle type)
#define MAXSLICES 100000
#define MAXFRAMES 300000
#define MAXSTATES 3
#define MAXREPLICA 10 //was 30 
#define MAXSETS 0
#define NACC 4
#define NSTAT 2
#define MAXBIN 100
#define RDFBINS 100
#define NINTERVALS 250
#define ANGLEBINS 100
#define MAXN0 200
#define MAXK 1000
#define MAXBONDS NSITES
#define MAXTOTBONDS NSITES*NPART

#define MAXPARTS            NPART
#define MAX_IMAGES          (MAXPARTS *8) /* is MAXPARTS NPART?*/
#define MAX_NEIGHBORS       (MAXPARTS *300)
#define MAX_IM_NEIGHBORS    300
#define STACKDEPTH          30
#define NSHIFTS             27

#define GRIDOFFSET          3
#define LJ_RZERO            1.5
#define MAXBOXSIZE          20
#define MINCELLSIZE         (LJ_RZERO/GRIDOFFSET)
// #define NCELLX              ((int) (MAXBOXSIZE / MINCELLSIZE)) gives errors??
#define NCELLX              40
#define NCELLY              (NCELLX)
#define NCELLZ              (NCELLX)


/*---------------------------------------------------------------------------*/
/*------------------STRUCTURE DEFINITIONS------------------------------------*/

typedef struct IntArraytype {
  // a list type of intergers
  int *ints;
  size_t used;
  size_t size;

} IntArray ;

typedef struct vector_type {

    double        x,y,z;

} vector;


typedef struct quaternion_type {

    double        q0,q1,q2,q3;

} quaternion;


typedef struct tensor_type {

    vector        x,y,z;

} tensor;


typedef struct quattensor_type {

    vector        w,x,y,z;

} quattensor;

typedef struct bondproperty_type{
  /* used for showing nematic order parameter in the graphics for Simons paper
  the bonds are listed with a maximum length of MAXTOTBONDS = NPART*MAXBONDS
  the list contains:
   a location (x,y,z)
   an orientation (q0,q1,q2,q3) in quaterion
   and the nematic order parameter S
  */

  vector       r;

  quaternion   q;

  double       nematic_OP;

} BondProperty;

typedef struct particle_type {

    vector        r,
                  f,
                  t,
                  dr;

    vector        patchvector[NSITES]; //the current ortientation of the patches

    int           cluster_id,
                  ptype, // ptype stands for particle type, with sys.particletype[ptype] you see the characteristics 
                  bonds[MAXBONDS],
                  nbonds,
                  k_bo; //number of bond order properties

    double        bond_energy[MAXBONDS],
                  bond_distance[MAXBONDS];

    quaternion    q;

    BondProperty  bond_op[MAXBONDS];

} Pts;


typedef struct site_type {
    
    vector        r;

    int           switch_type;

} Site;

typedef struct activeforce_type {
  /*The active force with strength F0 and direction r [unit vector]*/
  double        F0;

   vector       r;

} ActiveForce;

typedef struct saccentint_type {
  /*parameters for the integrated Saccent*/
  int             switch_variable,  //0 if the paramters in angles and 1 if in cosangle (only type 1)
                  s_accent ;        //if 0, use old settings where only partilces.inp is dep to chose Sold or Snew_S.
                                    // if 2,use Slin  

  double          switch_expc,          // swich function coefficients (only type 1)
                  switch_expd,
                  switch_expe,
                  switch_expf,
                  switch_expg,
                  switch_exph,
                  switch_expi;

} SaccentInt;

typedef struct particletype_type {
  /*there are diffent kind of particle type
  e.g. dipatch particle with diameter 3.2 miron = 1 sigma
  or tetrapatch particle with diameter 3.7 micron
  in particle type the number of patches plus diameter (-> gravity) is specified*/

  int           nparticles,
                nsites,
                activity, //0 no 1 yes (then read in active force)
                s_accent; // def of S' -> 0 (old switch; broad) 1 (integrated switch, dT dep)

  Site          site[NSITES];

  double        diameter,/*of bulk [sigma]*/
                radius, /*of bulk[sigma]*/
                small_radius, /*of patch partilce[sigma]*/
                small_diameter, /*of patch partilce[sigma]*/
                delta_degree,     /*patch width cutoff in degrees (i.e. where S' is 0)*/
                cosdelta,         /*patch width cutoff in cosdelta*/
                cosphi,         /*patch width cutoff in cosdelta*/
                oneover_cosdelta,  /*used if s_accent==0 */
                oneover_cosphi,  /*used if s_accent==0 */
                delta_rho_kg_m3, /* the density diffence between the lutide-water solution and the colloid [kg/m^3]*/
                fg,           /*gravitational force*/
                zcut,         /*cut off value of z in sigma* = 2mum*/
                b_zc,         /*V_g = LJ(z) or fg*z - b_zc */
                d_cbsp,           // this is the distance (d_) from center bulk particle (cb) to surface of patch particle (sp)
                d,           // this is the distance between center of bulk and patch particle
                Derj_scaling[PTYPES]; // each combination of particle type has a derj scaling factor. Read as: the scaling with particletype 0 is Derj_scaling[0]=alpha

  //each particle(type) can only have one active force, it you want different active forces, specify different particles types
  ActiveForce   active_force; 
  
} Particletype;

typedef struct bondinfo_type {
  /* the struct that saves the bond infromation, */

  double  s,  //surface-surface distance
          anglei, // angle 1 
          anglej,
          Erep, //repulsive energy
          Ec, // attractive casimir 
          S;   // swithcing function

} BondInfo;

// typedef struct bondbreak_type {
//   /* the struct that saves the coordinates, bond information and timeframe*/
//   vector       r[NPART];       // contains coordinates about location, patchvectors etc.
//   quaternion   q[NPART];

//   BondInfo  bondinfo[NPART];  // contains coordinates about bond energies etc.

//   double    ctime;            // timestamp of the timeframe
//   int       nbonds;

// } BondBreakage;

typedef struct statistics_type {
  /* struct for statistics: running mean and variance
  mean = u_n = u_{n-1} + (x_n-u_{n-1})/n
  variance2 = sigma2_n = (sigma2_{n-1} + (x_n-u_n)(x_n-n_{n-1}))/(n-1)*/

    double        mean,
                  variance2;
    long           n;

} Statistics;

typedef struct statslength_type {
  /* struct for statistics: running mean and variance
  mean = u_n = u_{n-1} + (x_n-u_{n-1})/n
  variance2 = sigma2_n = (sigma2_{n-1} + (x_n-u_n)(x_n-n_{n-1}))/(n-1)*/

    char          filename[100];

    Statistics    length[NPART];

} StatsLength;

typedef struct mcacc_type {

    char          name[100];

    long int            acc,
                        tries;

    double        ratio;

} Mcacc;

typedef struct aver_type {

    char          name[100];

    double        now,
                  n,
                  sum,
                  sumsq;

} Average;


typedef struct bonds_mat {

    double        mat[NPART][NPART][NSITES][NSITES];

} Bonds;


typedef struct hist1d_type {
            
    double        *bin,
                  dbin,
                  offset,
                  weight;
    int           maxbin;
            
    char          name[100];

} Hist1d;


typedef struct hist2d_type {

    double        **bin,
                  dbin1,
                  dbin2,
                  offset1,
                  offset2,
                  weight;
    int           maxbin1,
                  maxbin2;
            
    char          name[100];

} Hist2d;



typedef struct dop_type {
    
    int           **flag_dop,
                  *flag_norm,
                  npaths,
                  maxbin1,
                  maxbin2,
                  type;

    double        **dop, 
                  ***vf,
                  **escape,
                  *norm,
                  **dopen,
                  *dopennorm,
                  weight,
                  offset1,
                  offset2,
                  dbin1,
                  dbin2;

    char          name[40];

    FILE          *pathAB;

} Dop;

typedef struct slice_type {

    Pts           pts[NPART];

    double        energy;
    
                  // rdf[RDFBINS],
                  // localCN[NPART];

    int           nclusters,
                  nbonds;

} Slice;

typedef struct segment_type {

    Pts           pts[NPART];

    double        t,          /*time*/
                  lambda;     /*order parameter*/

    struct segment_type*       next;

} Segment;


typedef struct trajectory_type {

    double      lambda[MAXBIN];

    Slice*      path[MAXBIN];

} Trajectory;

typedef struct replica_type {


    double        lambda, 
                  weight,
                  logcrossprob,
                  lnweight,
                  dos;
  
    int           pathlen,
                  index, 
                  numpathstot,
                  numpathsnextint,
                  swapindex,
                  ntotal,
                  navlen,
                  type;  

    long int      avlen,
                  avlensq;

    //Dop           dop[MAXSTATES];

    Hist2d        fe;

    vector        string;


    //FILE          *pathfilep[2];

} Replica;

typedef struct state_type {

    Replica       srep[MAXREPLICA];

    double        weight,
                  lnweight,
                  dos,
                  crosshist[MAXREPLICA][MAXBIN],
                  pathtypenumbers[MAXREPLICA][MAXSTATES][MAXREPLICA],
                  problambdaAA[MAXREPLICA][MAXBIN],
                  problambdaAB[MAXREPLICA][MAXBIN],
                  flux0,
                  flux1,
                  rate[MAXSTATES],
                  lambda[MAXREPLICA],
                  logcrossprob,
                  min,
                  mindist,
                  scalefactor,
                  volume_op;

    int           nrep, 
                  type_mat[MAXREPLICA][2],
                  mstis_mat[MAXSTATES],
                  freq[MAXN0],
                  current_path,
                  maxpaths,
                  ntotal,
                  nflux0,
                  nflux1,
                  fcount,
                  nsucc;

    Slice         target;

    Segment       *path;

} State;


typedef struct langevin_type {
  
    int           ninter; //number of langevin steps per propagate_bd

    double        timestep,
                  dtD,
                  dtBeta,
                  friction;
} Langevin;

typedef struct lookup_table {

    int         intervals;

    double      ktable[NINTERVALS + 1],
                stable[NINTERVALS + 1];

    Statistics  rng_check;

} Lookuptable;

typedef struct chainnode_type { 
  /* here the particles are ordered according to bonding.
  with next and prec you go to the next or previous particle in the chain
  works only for chains,i.e. npatches =2.*/

    int           pts_nr; 
    
    struct chainnode_type   *next, 
                            *prev; // Pointer to previous node in DLL 

} ChainNode;




// typedef struct stats_type {

//     Average       aver[NSTAT];

//     Mcacc         mcacc[NACC];

// } Stats;

typedef struct chain_type {

    // Average       length_histogram[NPART];
                  // e2e2_avg[NPART]; /*make both into an Statistics*/

    StatsLength   e2e2;

    struct chainnode_type      *head[NPART];


} Chain;

typedef struct particles_in_cluster_list_type{
  /* stack,  i.e. a list of particles that are in this cluster, these are not ordered according to the bonds*/
    int           stack[NPART];

} Picl;

typedef struct cluster_type{
  /* the cluster type contain information about all clusters, 
  there are maximally NPART clusters (if all are monomers) 
  in pic (particles in clusterlist) you use the clusterid -> pic[clusterid] 
  to find the stack */

    int           analysis,
                  clustersize[NPART], 
                  update; 

    struct particles_in_cluster_list_type     pic[NPART];

    StatsLength   size_histogram,
                  size_distribution;

} Cluster;

typedef struct path_type {

    int           nslices,
                  ninter,
                  type,
                  nstates,
                  nshoot,
                  nreverse,
                  nrepswap,
                  current_initial_state,
                  current_final_state,
                  nswapstates,
                  initial_state,
                  final_state,
                  current_replica,
                  current_set,
                  stateswapbias,
                  nreplica,
                  fixedbias,
                  maxlength;
    
    double        scalefactor, 
                  current_gsen,
                  enbond,
                  energy;

    FILE          *filepath[MAXSTATES][MAXSTATES],
                  *filepatherr[MAXSTATES],
                  *fileswap[MAXSTATES],
                  *filestateswap,
                  *fprc;

    Dop           dop[3];

    // Stats         block_stats[MAXSTATES][MAXREPLICA],
    //               final_stats[MAXSTATES][MAXREPLICA];       


} Path;

typedef struct analysisptr_type {


  int           adjacency,
                *s_distribution,
                *rdfanalysis,
                *bond_op;


} Analysis_ptrs;

typedef struct system_type {

    int           adjacency,
                  s_distribution,
                  rdfanalysis,
                  dipatch_only,
                  bond_breakage,
                  particle_setup, 
                  cluster_MC,
                  checks,
                  // cc, changed this to lincs
                  empty_files,
                  graphics,
                  K, //K is number trialruns from from lambdai (ffs)
                  lincs,
                  nchains,
                  ncycle1,
                  ncycle2,
                  nearest_neighbor,
                  npart,
                  npatch,
                  bond_op,
                  nparticle_types,
                  start_type,
                  sim_type,
                  snapshot,
                  switch_method,     // definition of full S.  0 :S=max(thetai, thetaj)) or 1: S=Si*Sj, 2 S=(lin comb of S0, S90, 180)
                  warmup,
                  xy_print;
                  

    double        energy,
                  beta,
                  bond_avg,
                  temp,        
                  bond_cutoffE,
                  rcutoff,
                  rcutoffsq,
                  mobilityT,
                  sqrtmobilityT,
                  mobilityR,
                  sqrtmobilityR,
                  oneover_cosdelta,
                  oneover_cosdeltarev,
                  // bondbreakage_treshold, // hard coded to 0.5; this is the trheshold to which the sim is stopped
                  lambda,
                  fSE,
                  gravity,
                  drmax,         //for MC 
                  drmax_cluster,
                  drmax_mono,
                  dqmax,
                  dqmax_cluster,
                  dqmax_mono,
                  c_time,   // cumulatice time, its not really a system variable but ok
                  r_wetting,
                  s_min,
                  surface_charge,
                  dT,
                  A,
                  B,
                  xi,
                  Arep,
                  chaingap,
                  dl,
                  sigma, /*[mum]*/
                  epsilongravLJ,
                  S_fixed,
                  Erep_smin,
                  Ec_smin,
                  E_smin,
                  wall_int;
                  
    vector        boxl;
                  
    Particletype  particletype[PTYPES];

    FILE          *filep;

} System;

typedef struct constraint_handler {
    int         cmax,
                K;

    double      d,
                d2,
                cm,
                cm2;
} Constraint;

typedef struct cell2_type {
 
  int                n                 ,
                     stack[STACKDEPTH] ;
} Cell;

typedef struct nlist_item_type {

  int                ipart             ,
                     image_id          ,
                     *first           ;
} List_item;

typedef struct list_type {

  int                cellx             ,
                     celly             ,
                     cellz             ;  

  vector             cellsize          ,
                     inv_cellsize      ,
                     trans[NSHIFTS]    ;

  List_item          nli[MAX_IMAGES]    ;

  int                nlj[MAX_NEIGHBORS] ,
                     nneighbors         ,
                     nimages           ,
                     count             ;

  double             cutoff            ,
                     cutoff2           ,
                     cutoff_buffer2    ,
                     bigpart_cutoff2   ;

  Cell               cell[NCELLX][NCELLY][NCELLZ];

} List;



                    
/*---------------------------------------------------------------------------*/
/*------------------GLOBALLY DEFINED VARIABLES-------------------------------*/


// extern Average rdf_average[RDFBINS], globalCN_average;
extern Chain chain;
extern Cluster cluster;
extern Constraint cconst;
extern int crossing_check,s_histogram[MAXBIN],r_S_2dhist[MAXBIN][MAXBIN], xy_plane_patchangle[ANGLEBINS];
extern Mcacc  trans,        fintrans,        rot,        finrot;
extern Mcacc  trans_cluster,fintrans_cluster,rot_cluster,finrot_cluster ;
extern Mcacc  trans_mono,   fintrans_mono,   rot_mono,   finrot_mono;
extern Langevin langevin;
extern List list;
extern Path     path;
extern Replica *replica[MAXREPLICA];
extern Slice *slice,*psl_old,*init_state,*fwtrial,*bwtrial,*trial;
extern Statistics bond_length,msd,Sfraction;
extern State state[MAXSTATES];
// extern Stats nulstat;
extern System sys;
extern Trajectory segment;
extern vector nulvec;
extern SaccentInt s_int;

// extern BondBreakage bondbreakage[MAXFRAMES];

/* PUT THESE IN main.c TOO!*/


//truly global 
// extern System sys;
// extern vector nulvec={0};
// extern SaccentInt s_int; //switch function variables
// extern Slice *slice,*psl_old;
// extern Cluster cluster;
// extern Chain chain;
// extern Analysis_ptrs analysis_ptrs;







/*---------------------------------------------------------------------------*/
/*------------------GLOBALLY DEFINED FUNCTIONS-------------------------------*/

/*analysis.c*/
extern void nematic_order_parameter(Slice *psl);
extern double Sfrac(Slice *);
extern void timecorrelation_function_P(Slice *);
extern void MSD(Slice *);
extern void xy_plane_patchangle_calc(Slice *);
extern int bond_check(Slice *, int , int );
extern double bond_distance(Slice *, int , int , int);
extern void bond_length_average(Slice *psl);
extern void s_distribution_analysis(Slice *);
extern vector bond_vector(Slice *, int , int , int);
extern void calc_chainangles(Slice *);
extern void calc_tangentangles(Slice *, int);
extern void calc_rdf_localCN(Slice *);
extern void clustersize_freq_update(Slice *);
extern void chain_linkedlist(Slice *);
extern void check_maxbonds(Slice *);
extern int cluster_analysis(Slice *);
extern void clustersize_identification(Slice *psl);
extern double end_to_end_distance2(Slice *, int);
extern ChainNode* GetNewNode(int);
extern void InsertAtChainHead(int, int);
extern void InsertAtChainTail(int, int );
extern double running_variance2(double, double, double, double, long);
extern double running_mean(double, double, long);
extern void reset_running_statistics(Statistics *);
extern void running_statistics(Statistics *, double );
extern int particle_in_wall(Slice *psl, int ipart);
extern double op_calc(Slice *);
extern void update_patch_vector_ipart(Slice *, int );
extern void update_patch_vectors(Slice *);
extern void print_xy_positions(Slice *);
extern void linked_structure_icluster(Slice *, int , int );

/*constructnn.c*/
extern void setup_nnlist();
extern void put_parts_in_box(Slice *);
extern void create_cell_list(Slice *);
extern void find_neighbor_using_cells(Slice *,int);
extern void find_neighbor(Slice *,int );
extern int check_nnlist(Slice *);
extern void update_nnlist(Slice *);
extern void print_nnlist();

/*draw.c*/
extern void reset_center(Slice *);
extern void Init_Graphics();
extern void readslice(Slice *);
extern int readpath();

/*energy.c*/
extern double total_energy( Slice * );
extern double particle_energy( Slice * , int,int);
extern void energycheck_nn(Slice *, double ,double );
extern double potential_attractive_bond_energy(Slice *, int , int );
extern double gravitational_energy_ipart(Slice *, int );
extern double total_energy_neighborlist(Slice *);
extern double harmonic_oscillator(double , double );
extern void orientation_parameters(Slice *, int , int , double *,double *,double *, double *,double *, double *, double *);
extern vector find_directing_patchvector(Slice *, int ,vector );

/*ffs.c*/
extern void save_ffscycle(int );
extern void read_coordinates_of_interface(Slice *, int );
/*force.c*/
extern void gravitational_force_ipart(Slice *, int );
extern void single_bond_force(Slice *, int , int,vector);
extern void propagate_bd(Slice *);
extern void calculate_total_forces(Slice *);
extern void calculate_forces_neighborlist(Slice *);
extern void save_old_positions(Slice *, Slice *);
// extern void lincs_solve(Slice *, double , vector , double , double );
extern void lincs(Slice *, Slice *);
extern void gravitational_force_ipart(Slice *, int );
extern void forcecheck_nn(Slice *);

/* init.c*/
extern void setup_simulation( void );
extern void init_model(Slice *);
extern void setup_positions_sites( Slice * );
extern void read_input(void);
extern void read_particletypes(Slice *);
extern void read_statistics_file(StatsLength *);
extern void plotpotential(void);

/* main.c*/
extern void terminate_block();
extern void finalstat(); 
extern void bmdcycle();
extern void mccycle();
// extern void tiscycle();
extern void ffscycle(int, int);
extern void mainloop_for_graphics();
extern void emptyfiles();
extern void error(char *);

/*mc.c*/
extern void propagate_mc( Slice * );
extern void cluster_propagate_mc(Slice * );
extern void energy_divergence_check(Slice *, char [50]);
extern void optimizemc();
extern void printstatusmc(Slice *);
extern int single_particle_move(Slice *, int );
extern int single_particle_move_neighborlist(Slice *, int );
extern int rotatepart_cluster_Ekparts(Slice *);
extern int translatepart_cluster_Ekparts(Slice *);
/*function below are actually redundant; they are not used anymore*/
extern int translatepart(Slice *);
extern int translatepart_Ekpart(Slice *);
extern void rotatepart(Slice *);
extern void rotatepart_Ekpart(Slice *);
extern int rotatepart_cluster(Slice *);
extern int translatepart_cluster(Slice *);



/*switchfunction.c*/
extern double S_value(Slice *, double , double ,vector , vector , vector , int, int);
extern double Saccent(double , int  );
extern double old_S(double , int);
extern double new_S(double );
extern double S_lin(double , int );
extern double S0(double ,double );
extern double S90(double ,double );
extern double S180(double ,double );
extern double switch_method_2(double , double , double , double , vector , vector , vector );
extern double dS_dcostheta(Slice *, double , int , double , int  );
extern double dSaccent_dcosangle(double, int );
extern double dSold_dcostheta(double, int  );
extern double dSnew_dcostheta(double , int );
extern double dS0_dcostheta(double, int  ,double, int  );
extern double dS180_dcostheta(double, int  ,double, int  );
extern double dS90_dcostheta(double, int  ,double, int  );



/*print.c*/
extern void print_input(void);
extern void printrdf();
extern void conf_output(Slice *);
extern void conf_input(Slice *);
extern void print_pos_sites(Slice *);
extern void printenergy();
extern void print_statistics_file(StatsLength *);
extern void print_statistics_N1(Statistics , char [100]);
extern void print_chain_order(Slice *, int );
extern void print_adjacency_matrix_coordinates(Slice *);
extern void saveScreenshotToFile( int,int);
extern void print_energy_s_theta1_theta2(Slice *);




/*random.c*/
extern double RandomNumber(void);
extern void InitializeRandomNumberGenerator(double);
extern double RandomGaussianNumber();
extern vector RandomBrownianVector(double);
extern vector RandomVector(double);
extern double RandomVelocity(double);
extern vector RandomUnitVector(void);
extern int Random_urandom(void);
extern vector RandomVector3D(vector ) ;
extern double RandomNumberRange(double , double );
extern int RandomIntegerRange(int  , int );


// simonsparam.c
extern void setup_Simons_potential_parameters();
extern double potential_repulsive_energy_sdist(double );
extern double bond_repulsive_force(double );
extern double second_der_Vrep(double );
extern double potential_attractive_energy_sdist(double );
extern double bond_attractive_force(double );
extern double second_der_Vc(double );


/*quaternion.c*/
extern void rotate_quaternion_ipart( Slice *, int , quaternion  );
extern double langrange_multiplier_quat(quaternion, quaternion);
extern quaternion quatVecToVec(vector, vector);
extern tensor getrotmatrix(quaternion);
extern quattensor getquatmatrix(quaternion);
extern quaternion RandomQuaternion();
extern quaternion RandomQuaternionRange(double);
extern quaternion QuaternionXaxis( double );
extern quaternion QuaternionYaxis( double );
extern quaternion QuaternionZaxis( double );

/*storage.c*/

extern void initIntArray(IntArray *, size_t );
extern void insertIntArray(IntArray *, int ); 
extern void freeIntArray(IntArray *);
extern void removeIthElementIntArray(IntArray *, int  );
extern void removeElementXIntArray(IntArray *, int  );
extern int checkElementXIntArray(IntArray *, int );
extern void memory_allocation(void);
extern void printIntArray(IntArray * );


