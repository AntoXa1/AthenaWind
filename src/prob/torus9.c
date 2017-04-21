
#include "copyright.h"
/*============================================================================*/
/*! \file hkdisk.c
 *  \brief Problem generator for Hawley Krolik disk (Specifically, GT4)
 */
/*============================================================================*/

#include <math.h>
#include <stdio.h>

#include <unistd.h>

#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

//****************************************

#define HAWL
#define RADIAL_RAY_TEST



//#define SOLOV
//#if !defined(HAWL)
//#define SOLOV
//#endif


//****************************************
/* make all MACHINE=macosxmpi */







static void inX1BoundCond(GridS *pGrid);
static void diode_outflow_ix3(GridS *pGrid);
static void diode_outflow_ox3(GridS *pGrid);

void plot(MeshS *pM, char name[16]);
//void plot(MeshS *pM, char name[16]);

static void calcProblemParameters();
static void printProblemParameters();

#ifdef XRAYS
//void optDepthFunctions(MeshS *pM);

void optDepthFunctions(GridS *pG);


Real updateEnergyFromXrayHeatCool(const Real dens, const Real Press, const Real dt, const Real xi_in, int whatToDo);
void xRayHeatCool(const Real dens, const Real Press, const Real xi_in, Real*, Real*,  const Real dt);
Real rtsafe_energy_eq(Real, Real, Real, Real, int*);
#endif

// Input parameters
static Real q, r0, r_in, rho0, e0, dcut, beta, seed, R_ib;
static Real nc0, F2Fedd, fx =0.5, fuv =0.5;
static Real rRadSrc[2]={0.,0.}; //coordinates of the radiation source
static Real Ctor, Kbar, nAd, q1, xm, rgToR0,Tx,Rg,MBH_in_gram,Tgmin,
  Dsc,Nsc,Rsc,Usc,Psc, Esc,Time0,Tvir0,Lamd0,Evir0, Ledd_e, inRadOfMesh,
  Lx, Luv, a1;



static double
  CL = 2.997925E10,
  GRV = 6.67384E-8,
  QE=4.80325E-10,
  MSUN = 1.989E33,
  M2Msun = 10E7,
  GAM53 = 5./3.,
  KPE = 0.4,
  PC = 3.085678E18,
  ARAD = 7.56464E-15,
  RGAS = 8.31E7,
  MP = 1.672661E-24,
  M_U = 1.660531E-24,
  RSUN = 6.95E10,
  SGMB,MSOLYR,
  YR = 365.*24.*3600.,
  M_MW = 1., //mean mol. weight
  tiny = 1.E-10;
//	PI = 3.14159265359;

//----------------------------


#ifdef SOLOV
static Real
  Kbar_Sol,
  Br0_Sol = 1.,
  Br1_Sol = 6.,

  a1_Sol = 1,
  a2_Sol = 1,
  b1_Sol = 3.,

  w0_Sol = 2.,
  w1_Sol = 0.;
#endif
/* ============================================================================= */
/*                       ray tracing functions                                      */
/* ============================================================================= */
void testRayTracings(GridS *pG);
void traceGridCell(GridS *pG, Real *res, int ijk_cur[3], Real xpos[3], Real* x1x2x3, const Real
		   cartnorm[3], const Real cylnorm[3], short*);
void coordCylToCart(Real*, const Real*,const Real, const Real );
void cartVectorToCylVector(Real*, const Real*, const Real, const Real);

void vectBToA(Real*, const Real*, const Real*);

void normVectBToA(Real* v, const Real* , const Real*, Real* );

Real absValueVector(const Real*);
  
/* ============================================================================= */
/*                      boundary conditions                                      */
/* ============================================================================= */

#ifdef MPI_PARALLEL
/* MPI send and receive buffers */
static double **send_buf = NULL, **recv_buf = NULL;
static MPI_Request *recv_rq, *send_rq;
#endif /* MPI_PARALLEL */

//  \brief Generic void function of Grid. 
/* boundary condition function pointers. local to this function  */
static VGFun_t apply_ix1 = NULL, apply_ox1 = NULL;  
static VGFun_t apply_ix2 = NULL, apply_ox2 = NULL;
static VGFun_t apply_ix3 = NULL, apply_ox3 = NULL;

void bvals_tau_init(MeshS *pM);
void bvals_tau(DomainS *pD);

#ifdef MPI_PARALLEL
static void pack_tau_ix1(GridS *pG);
static void pack_tau_ox1(GridS *pG);
static void pack_tau_ix2(GridS *pG);
static void pack_tau_ox2(GridS *pG);
static void pack_tau_ix3(GridS *pG);
static void pack_tau_ox3(GridS *pG);

static void unpack_tau_ix1(GridS *pG);
static void unpack_tau_ox1(GridS *pG);
static void unpack_tau_ix2(GridS *pG);
static void unpack_tau_ox2(GridS *pG);
static void unpack_tau_ix3(GridS *pG);
static void unpack_tau_ox3(GridS *pG);
#endif /* MPI_PARALLEL */

static void boundCondOptDepthLike_ix1(GridS *pG);
static void boundCondOptDepthLike_ox1(GridS *pG);
static void boundCondOptDepthLike_ix2(GridS *pG);
static void boundCondOptDepthLike_ox2(GridS *pG);
static void boundCondOptDepthLike_ix3(GridS *pG);
static void boundCondOptDepthLike_ox3(GridS *pG);
 






//extern x1GravAcc;

/*! \fn static Real grav_pot(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential */


void getIndexOfCellByZ(GridS *pG, Real z, int  *k){
  *k= (int)( (z - pG->MinX[2])/pG->dx3 -0.5 ) + pG->ks;
  /* printf(" %d, \n",*k); getchar(); */
  }
//	*px1 = pG->MinX[0] + ((Real)(i - pG->is) + 0.5)*pG->dx1;
//	*px2 = pG->MinX[1] + ((Real)(j - pG->js) + 0.5)*pG->dx2;
//	*px3 = pG->MinX[2] + ((Real)(k - pG->ks) + 0.5)*pG->dx3;


void apause(){
  printf("pause.. press any key ss... ");
  getchar();
}


#ifdef XRAYS


void optDepthFunctions(GridS *pG){
  // it is assumed that the source is located on the axis of symmetry

  //	GridS *pG=pM->Domain[0][0].Grid;

  Real r, t, z, dl, sToX, sToZ,
    tau=0,dtau=0,tau_ghost,xi=0,x1,x2,x3, ro,rad,
    nnorm, norm[2], sfnorm[2],
    rCur[2],rNextX[2],rNextZ[2],colDens,
    den, xBndMax,zBndMin,zBndMax,x_is,z_is,res, sToX1, sToZ1, rsph;

  int i,k,is,ie,js,je,ks,ke, il, iu, jl,ju,kl,ku,ip,jp,kp,knew,
    my_id=0;
  
  float sgnf=1.;
  
  is = pG->is;
  ie = pG->ie;
  js = pG->js;
  je = pG->je;
  ks = pG->ks;
  ke = pG->ke;
  
  il = is - nghost*(ie > is);
  jl = js - nghost*(je > js);
  kl = ks - nghost*(ke > ks);
  
  iu = ie + nghost*(ie > is);
  ju = je + nghost*(je > js);
  ku = ke + nghost*(ke > ks);

#ifdef MPI_PARALLEL
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

  printf("optdepth,  my_id = %f \n",(float)my_id);

  //infinite loop for parallel debugging
  int  mpi1=1;
  // while( mpi1==1 );
 #endif 
  
    
  for (kp=ks; kp<=ke; kp++) {   // z
    for (jp=js; jp<=je; jp++) { // phi

      tau = pG->tau_e[ kp ][ jp ][ 1 ];
     
      
      pG->tau_e[kp][jp] [ is ] = tau;
      
      // tau = 0;      
      for (ip = is; ip<=ie; ip++) { //R

	//	pG->tau_e[kp][jp][ip] =  1.; //pG->tau_e[kp][jp][1];
	
	// Printf(" %d %d %d %d %d %d \n",il, jl, kl, iu,  ju,  ku );
	// getchar();
	//	tau=0.;
	//pG->tau_e[kp][jp][ip]=0.;	
       
	cc_pos(pG,ip,jp,kp,&x1,&x2,&x3);

	nnorm = sqrt(pow(x1-rRadSrc[0], 2)+pow( x3-rRadSrc[1], 2));
	norm[0]=(x1-rRadSrc[0])/nnorm;
	norm[1]=(x3-rRadSrc[1])/nnorm;
	
	colDens =0.;
	x_is = pG->MinX[0] + 0.5*pG->dx1;
	
  sfnorm[0] =  copysign(1., norm[0] )*fmax(fabs(norm[0]),tiny);
	z_is = norm[1]/sfnorm[0]*(x_is - rRadSrc[0]) + rRadSrc[1];
      
	/* z_is = norm[1]/ fmax(norm[0],tiny)*(x_is - rRadSrc[0]) + rRadSrc[1]; */

	
	rCur[0]=x_is;
	rCur[1]=z_is;
	xBndMax = pG->MaxX[0];
      	zBndMin = pG->MinX[2];
      	zBndMax = pG->MaxX[2];
	
	rsph =  sqrt(pow(rCur[0] - rRadSrc[0], 2) + pow(rCur[1] - rRadSrc[1] ,2));

	//	tau = 0.;

	k = (int) ((z_is - pG->MinX[2])/pG->dx3 + pG->ks);

	i = is;
	while (i < ie)	{
	  sfnorm[0] =  copysign(1., norm[0] )*fmax(fabs(norm[0]),tiny);    
	  sfnorm[1] =  copysign(1., norm[1] )*fmax(fabs(norm[1]),tiny);	 	  
	    
	  // crossing face, xf,i+1
	  rNextX[0] = pG->MinX[0] + (i+1 - pG->is)*pG->dx1;

	  rNextX[1]= norm[1]/sfnorm[0]*(rNextX[0]-rRadSrc[0])+rRadSrc[1];

	  if( rNextX[0]>pG->MaxX[0] ||  rNextX[1]>pG->MaxX[2] ||
	      rNextX[1]<pG->MinX[2] ){
	    break;
	  }

	  sToX = sqrt(pow(rNextX[0] - rCur[0], 2) + pow(rNextX[1] - rCur[1] ,2));
		       
	  sToX1 = sqrt(pow(rNextX[0] - rRadSrc[0], 2) + pow(rNextX[1] - rRadSrc[1] ,2));

	  knew = k + ( copysign(1., norm[1]) );
			
	  rNextZ[1] = pG->MinX[2] + (  (Real) (knew - pG->ks)  )*pG->dx3;

	  if (norm[1] != 0){
	    rNextZ[0] =(norm[0]/sfnorm[1]) * ( rNextZ[1]-rRadSrc[1] )+rRadSrc[0];
	  } else{
	    rNextZ[0] =norm[0]/sfnorm[1]*(rNextZ[1]-rRadSrc[1])+rRadSrc[0];
	  }

	  sToZ = sqrt(pow(rNextZ[0] - rCur[0], 2) + pow(rNextZ[1] - rCur[1] ,2));
	  sToZ1 = sqrt(pow(rNextZ[0] - rRadSrc[0], 2) + pow(rNextZ[1] - rRadSrc[1] ,2));

	  dl=0.;
	  if (sToZ>sToX){
	    // Lz>Lx, k is the same
	    dl = sToX;
	    i+=1;
	    rCur[0]=rNextX[0];
	    rCur[1]=rNextX[1];
	  } else{
	    //Lx>=Lz
	    dl = sToZ; //up or down
	    k=knew;
	    rCur[0]=rNextZ[0];
	    rCur[1]=rNextZ[1];
	  }
	  den = (pG->U[k][jp][i].d < 1.05* rho0) ? tiny : pG->U[k][jp][i].d;

	  den=1.;

	  
	  dtau = dl * den;

	  tau  += dtau;
	  
	  // printf("%f %d, %d , %d , %d  \n", tau, i ,ip,jp,kp );
	  
	  if ( k >= ke || k<ks || i==ip ){
	    //	    printf("%f \n", pG->tau_e[kp][jp][ip] );
	    break;
	  }



	  
	} /* i-loop:  -- end of ray tracing -- */
	
	

	//pG->tau_e[kp][jp][ip] =pG->tau_e[kp][jp][0]*(float)my_id; 
	//printf("%f \n", pG->tau_e[kp][jp][ip] );	
	//	Pg->tau_e[kp][jp][ip]=10;//(float)my_id;

	
	//test, using the spherical radius as a benchmark
	/*   	uncommnet below */
 	pG->tau_e[kp][jp][ip] = Rsc*KPE*Dsc*tau;
	rsph =  sqrt(pow(x1 - rRadSrc[0], 2) + pow(x3 - rRadSrc[1] ,2));

	pG->tau_e[kp][jp][ip] =  tau;

	pG->tau_e[kp][jp][ip] =  tau/rsph;//pG->tau_e[ kp ][ jp ][ 1 ];//1.;

	//	printf("%f %d %d %d \n", pG->tau_e[ kp ][ jp ][ 1 ], kp, jp, ip);  
	
	
	if ((pG->U[k][jp][i].d > tiny) || (rsph > tiny)){

	  xi = Lx / (Nsc* fmax(pG->U[ kp ][ jp ][ ip ].d , tiny))/ pow( Rsc*rsph, 2 );

	}
	else{
	  printf("quite possibly something is negative in opticalProps ");
	  apause();
	}
	pG->xi[kp][jp][ip]  = xi;

	pG->xi[kp][jp][ip]  *=  exp(- (pG->tau_e[kp][jp][ip]) );    

      } //ip = x
      
      	
    }	// jp = phi
  }  //kp = z

}


Real updateEnergyFromXrayHeatCool(const Real dens, const Real Press, const Real dt, const Real xi_in, int whatToDo)
{
  Real Pnew,dP, t_th,Hx, Cx,  dPmax,res,dHxdT, Tg, rts;
  int status;
  dPmax = 0.1;

  if(isnan(Press)){
    printf("Press is nan %e %e \n", dens, Press);
    return (0.);
  }

  xRayHeatCool(dens, Press, xi_in, &Hx, &dHxdT, dt);

  t_th=Press/Gamma_1/fmax( fabs(Hx), tiny);
  Pnew = Press + dt* Hx/Gamma_1;

  //	dP = Pnew-Press;
  //	 if( fabs(dP/Press ) > dPmax){
  //		 dP =fabs(Press) * copysign(dPmax,  dP );
  //	 }
  //	 Hx = dP/dt/Gamma_1;

  //		if (Hx>0) printf("%e %e %e %e %e\n", t_th, dt, Hx, Press, xi_in);
  //	printf("%e \n", Gamma_1 ); //apause();



  //	if ( t_th < dt){
  Pnew = rtsafe_energy_eq(Press, dens, xi_in, dt, &status);
  //if (status==1) {printf("stop 2: %e %e res= %e t= %e H= %e stat= %d \n", Press, Pnew, res, t_th, Hx, status); /* apause(); */}

  dP = Pnew-Press;

  //		 if (dP >0) printf("%f \n", dP);

  if( fabs(dP/Press ) > dPmax){
    dP =fabs(Press) * copysign(dPmax,  dP );
    //		  printf("%e %e %e %e %f \n",max_delta_e, eold,  dt*Hx, de/eold, Tx);
  }
  Hx = dP/dt/Gamma_1;
  //	 }


  if(Pnew <= 0. ) Hx = 0.;

  //		 printf("negative pressure detected");


  //	 Hx = (abs(Hx)>0.0001 ) ?  copysign(0.0001, Hx) : Hx;
  //	 printf( " %e %e \n", Press, Pnew );
  //	 Tg = fabs( Press * Esc * M_MW) / (dens*Dsc*RGAS );
  //	 if(xi_in>100.) printf("%e %e %e %e %e %e \n", Hx, dens, Press,Pnew, Tg, xi_in);

  //	 Hx = 0.;
  return(-Hx); /* returns cooling NOT  heating rate */

}

void xRayHeatCool(const Real dens, const Real Press, const Real xi_in,
		  Real *Hx, Real *dHxdT, const Real dt)
{
  Real max_delta_e = 0.05, Gx,Gc,Lbr,lbr1,lbr2,dGcdT,dGxdT,dLbrdT,dHxde,Tg,xi,dP,
    res=0.,
    epsP=0.5,
    delt=1.,

    c1 = 8.9e-36,
    c2 = 1.5e-21,
    al0 = 3.3e-27,
    al1 = 1.7e-18,
    al2 = 1.3e5,
    bta = -1.;
  Real const Tgmin = 10.;

  Tg = fabs( Press * Esc * M_MW) / (dens*Dsc*RGAS );



  //	 return(Press/Gamma_1);
  //  egas = Evir0 * eij
  //  Tg = fabs( (gam -1.)*m_mw *egas/ (nd* nc0*m_u) /Rgas )

  xi = xi_in;

  //	xi = fmax(xi_in, 1.);

  //	 if(xi_in < 1.) {*Hx=0.; return;}

  Gc =  0.;
  Gx =  0.;
  Lbr = 0.;
  Gx= c2*pow(xi,0.25) /pow(Tg,0.5)*(1. - Tg/Tx);
  Gc = c1* xi*(Tx - 4.*Tg);

  lbr1= al0*sqrt(Tg);
  lbr2= ( al1*exp(-al2/Tg)*( pow(xi,bta) ) /sqrt(Tg)+1.e-24)*delt;



  if (Tg <= 1.e3 ){
    lbr2 =  ( al1*exp(-al2/Tgmin)*( pow(xi,bta) ) /sqrt(Tgmin)+1.e-24)*delt;
  }

  //	 lbr2 = 0.;

  Lbr=lbr1 + lbr2; //(erg cm^3/s)


  dP = dt* fabs(Lbr)/Gamma_1;
  if( fabs ( dP / Press ) > epsP) Lbr = fabs(Press) *epsP /dt/Gamma_1;

  //	 Lbr=0.;
  //	 if(xi_in < 0.01) Lbr = 0.;

  *Hx = Gc +Gx -Lbr;


  dGcdT = - 4.*c1*xi;
  dGxdT = - 0.5*c2/Tx * pow(xi,0.25) *(Tx/Tg+1.)*pow (Tg, -0.5);
  dLbrdT = 0.5* al0 *pow(Tg,-0.5) + delt*al1*pow(xi,bta) *exp(-al2/Tg)*
    ( al2 -  0.5* Tg  )  * pow(Tg,-2.5);
  *dHxdT = dGcdT + dGxdT - dLbrdT;
  dHxde = *dHxdT*Gamma_1*M_MW*Evir0/Dsc/RGAS/Lamd0 * dens * pow( nc0, 2);

  *Hx = *Hx/Lamd0 * pow(dens*nc0, 2);

  //	 if(xi_in>100.)  printf("%e \n", *Hx );
  //	 printf(" %e %e %e %e %e %e %e\n",Tg, Tx, xi, *Hx, Lbr, Gc, Gx); //apause();
  //	 printf("%e %e %e %e \n", dens, Press, Tg, xi);

  if(isnan(*Hx)){
    printf("Hx or/and Cx is nan %e %e %e %e \n", dens, Press, Tg, xi);
    *Hx = 0.;		//apause();
    return;
  }

  //	 *Hx = 0.;

  //	 *Hx = (Tg < Tx) ? *Hx : 0.;
  //	 printf("%e %e %e %e %e \n", Hx, Tg, xi, dens, Press);
}
#endif  /* XRAYS */


void plot(MeshS *pM, char name[16]){
  int i,j,k,is,ie,js,je,ks,ke,nx1,nx2,nx3,nr,nz, il, jl, kl, iu, ju, ku;
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  
  GridS *pG=pM->Domain[0][0].Grid;
  
  printf("WORKING DIR: %s \n" ,  getcwd(cwd, sizeof(cwd) ) );

#ifdef MPI_PARALLEL
  FILE *f = fopen("./athena_plot.tmp", "w");
#else
  FILE *f = fopen("./athena_plot.tmp", "w");
#endif
  
  if (f == NULL){

    printf("Error opening file!\n");
    exit(1);
  }
  is = pG->is; ie = pG->ie;
  js = pG->js; je = pG->je;
  ks = pG->ks; ke = pG->ke;
  
  il = is - nghost*(ie > is);
  jl = js - nghost*(je > js);
  kl = ks - nghost*(ke > ks);
  
  iu = ie + nghost*(ie > is);
  ju = je + nghost*(je > js);
  ku = ke + nghost*(ke > ks);

  nr= (ie-is)+1;
  nz=(ke-ks)+1;
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;

  //	printf("%d %d %d %d %d %d %d\n", nx1,nx2,nx3,nghost, (ie-is)+1 ,(je-js)+1,ie); apause();

  /* for(k=4;k<=10;k++) printf( "%f \n", pG->tau_e[k][10][20] ) ; */
  /* getchar(); */
 
  j=js;
  fprintf(f, "%d  %d\n", nr, nz );

  for (k=ks; k<=ke; k++) {
    /* for (j=js; j<=je; j++){	  */
    for (i=is; i<=ie; i++){
      
      if (strcmp(name, "d") == 0){
	fprintf(f, "%f\n", pG->U[k][j][i].d );
      }
      if (strcmp(name, "E") == 0){
	fprintf(f, "%f\n", pG->U[k][j][i].E );
      }
#ifdef XRAYS
      if (strcmp(name, "tau") == 0){

	/* printf("here -------------------- %f, %d, %d, %d \n", */
	/*        pG->tau_e[k][j][i],  k,j,i); */

	
	fprintf(f, "%f\n", pG->tau_e[k][j][i] );
      }
      if (strcmp(name, "xi") == 0){
	//				 fprintf(f, "%f\n", log10(pG->xi[k][j][i]) );
	fprintf(f, "%f\n", (pG->xi[k][j][i]) );
      }
#endif
    }
    fprintf(f,"\n");
  }
  fclose(f);
  system("./plot_from_athena.py");
}



static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  Real rad;


  rad = sqrt( SQR(x1) + SQR(x3) );

  return -1.0/(rad-rgToR0);

  // return 0.0;
}

/*! \fn static Real grav_acc(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational acceleration */

static Real grav_acc(const Real x1, const Real x2, const Real x3) {

  Real rad,res;

  rad = sqrt( SQR(x1) + SQR(x3) );

  res = (1.0/(SQR(rad-rgToR0)))*(x1/rad);

  return (res);

  //  return (1.0/(SQR(rad-2.0)))*(x1/rad); old

}

//Private functions (x1,x2,x3) = (R,p,z)

/*! \fn Real density(Real x1, Real x2, Real x3)
 *  \brief Calculates the density at x1, x2, x3*/

#define VP(R) pow(R/xm, 1.-q)

Real density(Real x1, Real x2, Real x3) {
  Real rad, tmp, d, res;
  rad = sqrt( SQR(x1) + SQR(x3));


#ifdef SOLOV
  tmp = (Ctor + xm/(rad-rgToR0)  - a1*pow(x1,-q1)) /(nAd+1)/Kbar;
  res= Br0_Sol + pow(x1,2)/2. + 1/sqrt(pow(x1,2) + pow(x3,2)) -
    (Br1_Sol + (w0_Sol*pow(x1,2))/2.)*(a2_Sol*pow(-1 + pow(x1,2),2) + (b1_Sol + a1_Sol*(-1 + pow(x1,2)))*pow(x3,2));
  res /= ((nAd+1)*Kbar);
  d = pow(res, nAd)*( res >=0.01 );
#endif


#ifdef HAWL
  tmp = (Ctor + (xm*1.0/( rad-rgToR0 )) - pow(xm, 2.0*q )*pow( x1, -q1 )/q1)/(nAd+1)/Kbar;
  d = pow(tmp, nAd)*( x1>=r_in );
#endif


  d = MAX(d, rho0);

  //  printf("in density()  %f %f %f %f %f %f\n", r_in, nAd, x1, rad, tmp, d); pause();
  //	  printf("%f %f %f %f %f %f\n", r_in, nAd, x1, rad, tmp, d); pause();
  //if (res>0.){
  //      printf("Kbar = %f %f\n", Kbar, Kbar_Sol);
  //  	  printf("%f %f %f %f\n", x1, tmp, res, d);
  // }
  //  d = pow(temp/Kbar,n)*(x1>=r_in);
  //  d = MAX(d,rho0);
  //  d = x1 < r_in ? rho0 : d;

  return d;

}

/*! \fn Real Volume(Grid *pG, int i, int j, int k)
 *  \brief Calculates the volume of cell (i,j,k) */
Real Volume(GridS *pG, int i, int j, int k) {
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  // Volume_ijk = R_i * dR * dPhi * dZ
  return x1*pG->dx1*pG->dx2*pG->dx3;
}

VOutFun_t get_usr_out_fun(const char *name)
{
  return NULL;
}

/*----------------------------------------------------------------------------*/


/* problem:   */
//#define VP(R) (  (sqrt(r0)/(r0-2)) * pow(R/r0,-q+1)  )

//#define VP(R) pow(xm,q)*pow(R, 1.-q)


#define VP(R) pow(R/xm, 1.-q)



void problem(DomainS *pDomain)
{
  GridS *pG=(pDomain->Grid);

#ifdef XRAYS
  CoolingFunc = updateEnergyFromXrayHeatCool;
#endif


  int i,j,k;
  int is,ie,js,je,ks,ke,nx1,nx2,nx3;
  int il,iu,jl,ju,kl,ku;
  Real rad, IntE, KinE, MagE, rbnd, dz;
  Real x1,x2,x3, x1i, x2i, x3i,rhoa;
  Real Pgas, Pb, TotPgas, TotPb;
  Real scl, ***Ap, ***Bri, ***Bzi;
  Real divB=0.0, maxdivB=0.0;



  //  ------------------------------
 //  Soloviev solution

#ifdef SOLOVо
  Real Psi_Sol;
#endif

  //  ------------------------------


  int my_id;

  is = pG->is; ie = pG->ie;
  js = pG->js; je = pG->je;
  ks = pG->ks; ke = pG->ke;
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;
  il = is - nghost*(ie > is);
  jl = js - nghost*(je > js);
  kl = ks - nghost*(ke > ks);
  iu = ie + nghost*(ie > is);
  ju = je + nghost*(je > js);
  ku = ke + nghost*(ke > ks);


  if ((Ap = (Real***)calloc_3d_array(nx3+1, nx2+1, nx1+1, sizeof(Real))) == NULL) {
    ath_error("[HK-Disk]: Error allocating memory for vector pot\n");
  }

  if ((Bri = (Real***)calloc_3d_array(nx3+1, nx2+1, nx1+1, sizeof(Real))) == NULL) {
    ath_error("[HK-Disk]: Error allocating memory for Bri\n");
  }

  if ((Bzi = (Real***)calloc_3d_array(nx3+1, nx2+1, nx1+1, sizeof(Real))) == NULL) {
    ath_error("[HK-Disk]: Error allocating memory for Bzi\n");
  }

  /* Read initial conditions */
  F2Fedd = par_getd("problem","F2Fedd");
  fx = par_getd("problem","fx");
  fuv = par_getd("problem","fuv");
  nc0 = par_getd("problem","nc0");
  q      = par_getd("problem","q");
  r0     = par_getd("problem","r0");
  r_in   = par_getd("problem","r_in");
  rho0   = par_getd("problem","rho0");
  e0     = par_getd("problem","e0");
  seed   = par_getd("problem","seed");

#ifdef MHD
  dcut = par_getd("problem","dcut");
  beta = par_getd("problem","beta");
#endif


  //printf("%f \n", F2Fedd); getchar();


  Real Tg = fabs( Gamma_1*M_MW *e0/ (rho0*Dsc) /RGAS );

  //	printf(" %e \n", Tg); apause();

  calcProblemParameters();

  printProblemParameters();

  //    Ctor = pow(r_in,-q1)/q1 - 1.0/(r_in-rgToR0);
  //    Kbar = (Ctor + 1.0/(1.-rgToR0) - 1./q1) /(nAd+1);
  //    Ctor = pow(xm, 2.0*q) * pow(r_in,-q1)/q1 - xm*1.0/(r_in-rgToR0);
  //    Kbar = (Ctor + xm*1.0/(xm-rgToR0) - pow(xm, 2.0) /q1) /(nAd+1);


  //   	    cc_pos(pG, is+5,js,ks,&R_ib,&x2,&x3);
  //   	    cc_pos(pG, is,js,ks,&x1,&x2,&x3)	;
  //        printf(" %e  %e \n ", R_ib, x1); apause();
  //    	   Kbar_Sol = (C3_Sol + (6. + 3.*b1_Sol + a2_Sol*w0_Sol)/6.);
  //    	   Kbar_Sol /= (nAd+1);



  //    printf(" %e %e %e %e %e %e %e\n", q, Rsc, Dsc, Usc, Time0,Tvir0,Lamd0); apause();

  cc_pos(pG,is,js,ks,&x1,&x2,&x3);
  inRadOfMesh=x1;

  /* int  mpi1=1; */
  /* while( mpi1==1 ); */
 

  
  // Initialize data structures, set density, energy, velocity and Ap
  // Loop over all cells (including ghosts)
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {

	//    	    printf("%d\n", ju);
	//    	  	printf("%f \n", rho0); apause();

	cc_pos(pG,i,j,k,&x1,&x2,&x3);
        rad = sqrt(SQR(x1) + SQR(x3));

        x1i = x1 - 0.5*pG->dx1;
        x2i = x2 - 0.5*pG->dx2;
        x3i = x3 - 0.5*pG->dx3;

        pG->U[k][j][i].d = rho0;
        pG->U[k][j][i].M1 = 0.0;


	pG->U[k][j][i].M2 = rho0 * sqrt(rad)/(rad-rgToR0);

	//        pG->U[k][j][i].M2 = VP(x1)*rho0;
	//        pG->U[k][j][i].M2 = 0.0;

        pG->U[k][j][i].M3 = 0.0;
        pG->U[k][j][i].E   = e0;

#ifdef MHD
        Ap[k][j][i] = 0.0;
        pG->B1i[k][j][i] = 0.0;
        pG->B2i[k][j][i] = 0.0;
        pG->B3i[k][j][i] = 0.0;
        pG->U[k][j][i].B1c = 0.0;
        pG->U[k][j][i].B2c = 0.0;
        pG->U[k][j][i].B3c = 0.0;
#endif

        // Set up torus

	//        rbnd =1./xm/( pow(xm, 2.0*q) *pow(x1, -q1)/q1 - Ctor ) + rgToR0;

        rbnd =xm/( a1*pow(x1, -q1) - Ctor ) + rgToR0;

        if (rbnd<=0) { printf("rbnd<=0: %e \n", rbnd); apause();}

        r_in = 0.5;
	//        if ( (x1 >= r_in) && (rad <= rbnd) ) { //Checks to see if cell is in torus

	if (x1 >= r_in) { //Dummy

	  rhoa = density(x1, x2, x3);

	  if (rhoa>0.1){
	    //        		printf("%f \n", rhoa);
	  }

	  IntE = pow(rhoa, GAM53)*Kbar/Gamma_1;

	  // Add pressure fluctuations to seed instability
	  IntE = IntE*(1.0 - seed*sin(x2));

	  if ((IntE >= e0) && (rhoa >= rho0)) {
	    //If the values are above cutoff, set up the cell

            pG->U[k][j][i].d = rhoa;
            pG->U[k][j][i].M2 = VP(x1)*pG->U[k][j][i].d;

#ifdef SOLOV
	    Psi_Sol= a2_Sol*pow(-1 + pow(x1,2),2) + (b1_Sol + a1*(-1 + pow(x1,2)))*pow(x3,2);
	    pG->U[k][j][i].M2 = sqrt( ( 1.-w0_Sol* Psi_Sol  > 0 ) ) *pG->U[k][j][i].d;
#endif


            pG->U[k][j][i].E = IntE;

            //Note, at this point, E holds only the internal energy.  This must
            //be fixed later
          }

        }
      }
    }
  }

  // Calculate density and set up Ap if appropriate
  for (k=kl; k<=ku+1; k++) {
    for (j=jl; j<=ju+1; j++) {
      for (i=il; i<=iu+1; i++) {

	cc_pos(pG,i,j,k,&x1,&x2,&x3);

	rad = sqrt(SQR(x1) + SQR(x3));

	x1i = x1 - 0.5*pG->dx1;
        x2i = x2 - 0.5*pG->dx2;
        x3i = x3 - 0.5*pG->dx3;

        rhoa = density(x1i,x2i,x3i);

	//printf("%f \n", dcut); apause();

#ifdef MHD
        if (rhoa >= dcut) {
          // Ap = (density-cutoff)^2

          Ap[k][j][i] = SQR(rhoa-dcut);

        }
	//        else {
	//        	  Ap[k][j][i] = x1i;
	//        }

#endif //MHD

      }
    }
  }


#ifdef MHD
  // Set up interface magnetic fields by using Ap
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {

        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        x1i = x1 - 0.5*pG->dx1;
        x2i = x2 - 0.5*pG->dx2;
        x3i = x3 - 0.5*pG->dx3;

        // Br = -dAp/dz
        pG->B1i[k][j][i] = -(Ap[k+1][j][i] - Ap[k][j][i])/pG->dx3;

        // Bz = (1/R)*d(R*Ap)/dr
        pG->B3i[k][j][i] = (Ap[k][j][i+1]*(x1i+pG->dx1) - Ap[k][j][i]*x1i)/(x1*pG->dx1);

        //Bt = d(Ap)/dz - d(Ap)/dr
	//        pG->B1i[k][j][i] = 0.;
	//        pG->B3i[k][j][i] = 0.;

        dz =  copysign(pG->dx3, x3);

	//        pG->B2i[k][j][i] = fabs( ( Ap[k+1][j][i] - Ap[k][j][i])/dz  - (Ap[k][j][i+1] - Ap[k][j][i])/pG->dx1);

	//        if ( pG->B2i[k][j][i] != 0. ){
	//        		printf("%f %f\n", Ap[k+1][j][i],  Ap[k][j][i] ); // pG->B2i[k][j][i]);
	//        }


	//        // non-zero B in the empty space
	//	    cc_pos(pG,i,j,k,&x1,&x2,&x3);
	//	    rad = sqrt(SQR(x1) + SQR(x3));
	//	    x1i = x1 - 0.5*pG->dx1;
	//	    x2i = x2 - 0.5*pG->dx2;
	//    		x3i = x3 - 0.5*pG->dx3;
	//    		rhoa = density(x1i,x2i,x3i);
	//        if (rhoa <  dcut) {
	//        	pG->B3i[k][j][i] = (Ap[k][j][i+1]*(x1i+pG->dx1) - Ap[k][j][i]*x1i)/(x1*pG->dx1);
	//        pG->B3i[k][j][i] *= dcut;
	//        }

#ifdef SOLOV
	a1_Sol = 1.;
	a2_Sol = 1.;
	b1_Sol = 3.;
	pG->B1i[k][j][i] = (-2*(b1_Sol + a1_Sol*(-1 + pow(x1i,2)))*x3i)/x1i;   //Bx
	pG->B3i[k][j][i] =4*a2_Sol*(-1 + pow(x1i,2)) + 2*a1_Sol*pow(x3i,2); //Bz
#endif

      }
    }
  }
  //apause();
  // Calculate cell centered fields
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {

        cc_pos(pG,i,j,k,&x1,&x2,&x3);

        if (i==iu)
          pG->U[k][j][i].B1c = pG->B1i[k][j][i];
        else

          pG->U[k][j][i].B1c = 0.5*((x1-0.5*pG->dx1)*pG->B1i[k][j][i]+(x1+0.5*pG->dx1)*pG->B1i[k][j][i+1])/x1;


        if (j==ju)
          pG->U[k][j][i].B2c = pG->B2i[k][j][i];
        else
          pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j+1][i] + pG->B2i[k][j][i]);

        if (k==ku)
          pG->U[k][j][i].B3c = pG->B3i[k][j][i];
        else
          pG->U[k][j][i].B3c = 0.5*(pG->B3i[k+1][j][i] + pG->B3i[k][j][i]);

	//        printf("B1c, B2c, B3c %E %E %E \n", pG->U[k][j][i].B1c, pG->U[k][j][i].B2c, pG->U[k][j][i].B3c);

      }
    }
  }
#endif //MHD


#ifdef MHD
  // Calculate scaling factor to satisfy beta, specifically Pgas and Pb per tile
  // Don't loop over ghosts
  Pgas = 0.0;
  Pb   = 0.0;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {

	Pgas += (Gamma-1)*pG->U[k][j][i].E*Volume(pG,i,j,k);

	Pb += 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c)
		   + SQR(pG->U[k][j][i].B3c))*Volume(pG,i,j,k);
      }
    }
  }
#endif //MHD

#ifdef MPI_PARALLEL
  if(MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD, &my_id))
    ath_error("[main]: Error on calling MPI_Comm_rank in torus9.c\n");

  MPI_Reduce(&Pgas, &TotPgas, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Pb, &TotPb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  if (my_id == 0) {
    printf("Total gas pressure = %f\n", TotPgas);
    printf("Total magnetic pressure = %f\n", TotPb);
  }
  
  MPI_Bcast(&TotPgas, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&TotPb, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#else

  TotPgas = Pgas;
  TotPb = Pb;

#endif //PARALLEL


#ifdef MHD
  //calculate and use scaling factor so that the correct beta is ensured
  scl = sqrt(TotPgas/(TotPb*beta));

  printf("Using magnetic scaling factor %f\n", scl);
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
	
        pG->U[k][j][i].B1c *= scl;
        pG->U[k][j][i].B3c *= scl;
        pG->B1i[k][j][i]   *= scl;
        pG->B3i[k][j][i]   *= scl;

        cc_pos(pG,i,j,k,&x1,&x2,&x3);
	x1i = x1 - 0.5*pG->dx1;
	x2i = x2 - 0.5*pG->dx2;
	x3i = x3 - 0.5*pG->dx3;
	rhoa = density(x1i,x2i,x3i);

	//        if (rhoa <  dcut) {
	//        	printf("%f %f \n",
	//		pG->U[k][j][i].B1c,
	//		pG->U[k][j][i].B3c );
	//        }
      }
    }
  }
#endif


  // Fix E to hold the total energy (Kin + Int + Mag) (Loop over ghosts)
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {

	KinE = 0.5*( SQR(pG->U[k][j][i].M1)+ SQR(pG->U[k][j][i].M2)
	      + SQR(pG->U[k][j][i].M3))/(pG->U[k][j][i].d);
        MagE = 0.0;

#ifdef MHD
        MagE = 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c)
		    + SQR(pG->U[k][j][i].B3c));
#endif
        pG->U[k][j][i].E += KinE + MagE;
      }
    }
  }
  /* Enroll the gravitational function and radial BC */


  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;
  
  /* setting pointers to user BC */
  bvals_mhd_fun(pDomain, left_x1,  inX1BoundCond );
  bvals_mhd_fun(pDomain, left_x3,  diode_outflow_ix3 );
  bvals_mhd_fun(pDomain, right_x3,  diode_outflow_ox3);



#ifdef XRAYS
  /* optical depth is calculated for a given domain, boundary conditions are applied in */
  
  //  optDepthFunctions(pG);
  //bvals_tau(pDomain);
  
  /* plot(pDomain, "tau"); */
  /* plot(pDomain, "xi"); */

#endif


    


  //  for (k=kl; k<=ku; k++) {
  //    for (j=jl; j<=ju; j++) {
  //      for (i=il; i<=iu; i++) {
  //
  //		pG->U[k][j][i].d =0.1;
  //		pG->U[k][j][i].E  = e0;
  //
  //		pG->U[k][j][i].M1 = 0.;
  //
  //		cc_pos(pG,i,j,k,&x1,&x2,&x3);
  //
  //		pG->U[k][j][i].M2 = 1./x1;
  //
  ////				pG->U[k][j][i].M3 = 0.;
  //
  //      }
  //
  //    }
  //  }
  //  x1GravAcc = grav_acc;

  //  set_bvals_fun(left_x1,  disk_ir_bc);
  //  set_bvals_fun(right_x1,  disk_or_bc);
  //  set_bvals_fun(left_x3,  diode_outflow_ix3);
  //  set_bvals_fun(right_x3,  diode_outflow_ox3);


  //  plot(pG, "d");

  return;
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 * current() - computes x3-component of current
 * Bp2()     - computes magnetic pressure (Bx2 + By2)
 *----------------------------------------------------------------------------*/

//void problem_write_restart(GridS *pG, DomainS *pD, FILE *fp)

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
//void problem_read_restart(GridS *pG, DomainS *pD, FILE *fp)
{
  GridS *pG=pM->Domain[0][0].Grid;

  //  q = par_getd("problem","q");
  //  r0 = par_getd("problem","r0");
  //  r_in = par_getd("problem","r_in");
  //  rho0 = par_getd("problem","rho0");
  //  e0 = par_getd("problem","e0");
  //	seed = par_getd("problem","seed");

  //#ifdef MHD
  //  dcut = par_getd("problem","dcut");
  //  beta = par_getd("problem","beta");
  //#endif

  /* Enroll the gravitational function and radial BC */
  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;

  calcProblemParameters();
  printProblemParameters();

  bvals_mhd_fun(pG, left_x1,  inX1BoundCond );
  bvals_mhd_fun(pG, left_x3,  diode_outflow_ix3 );
  bvals_mhd_fun(pG, right_x3,  diode_outflow_ox3);

  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

//Gasfun_t get_usr_expr(const char *expr)
//{
//  return NULL;
//}

//void Userwork_in_loop(Grid *pG, Domain *pDomain)  !A.D.

void Userwork_in_loop (MeshS *pM)
{

  int i,j,k,is,ie,js,je,ks,ke, nx1, nx2, nx3, il, iu, jl,ju,kl,ku;
  int prote, protd, my_id;
  int nl, nd;
  Real IntE, KinE, MagE=0.0, x1, x2, x3, DivB,Pgas,rad, dns,as_lim;

  Real di, v1, v2, v3, qsq,p, asq,b1,b2,b3, bsq,rho_max;

  static Real TotMass=0.0;

  //A.D.
  GridS *pG=pM->Domain[0][0].Grid;

  DomainS pD = pM->Domain[0][0];

  is = pG->is; ie = pG->ie;
  js = pG->js; je = pG->je;
  ks = pG->ks; ke = pG->ke;
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;
  il = is - nghost*(ie > is);
  jl = js - nghost*(je > js);
  kl = ks - nghost*(ke > ks);
  iu = ie + nghost*(ie > is);
  ju = je + nghost*(je > js);
  ku = ke + nghost*(ke > ks);

  // Verify divB
  protd = 0;
  prote = 0;

#ifdef MHD
  DivB = compute_div_b(pG);
#endif


#ifdef XRAYS



/* bvals_tau(&pD); */
 
 testRayTracings(pG);
  
/* optDepthFunctions(pG); */
 
/* //for(k=4;k<=10;k++) printf("%f \n", pG->tau_e[k][10][10]) ; */
/* bvals_tau(&pD); */

#endif

//plot(pM, "d");

 
//MPI_Comm_rank(MPI_COMM_WORLD, &my_id); 
//if (my_id != 2)

   plot(pM, "tau");

//plot(pM, "xi");
  //plot(pM, "E");


  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {

	cc_pos(pG,i,j,k,&x1,&x2,&x3);

	//        rad = sqrt(SQR(x1) + SQR(x3));
	//        if (rad< R_ib){
	//             pG->U[k][j][i].d = rho0;
	//             pG->U[k][j][i].E = e0;
	//             pG->U[k][j][i].M1= 0.;
	//             pG->U[k][j][i].M2= 0.;
	//             pG->U[k][j][i].M3= 0.;
	//        }


        KinE = 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2)
		    + SQR(pG->U[k][j][i].M3))/(pG->U[k][j][i].d);
	
#ifdef MHD

        MagE = 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c)
		    + SQR(pG->U[k][j][i].B3c));
#endif

        IntE = pG->U[k][j][i].E - KinE - MagE;
       
        if (isnan(pG->U[k][j][i].d) || isnan(pG->U[k][j][i].E)) {

	  printf("At pos isnan: (%d,%d,%d) (%f,%f,%f), Den1 = %f, E1 = %f\n", i, j, k, x1,x2,x3,
		 pG->U[k][j][i].d, pG->U[k][j][i].E);
	  printf("KinE1, MagE1, IntE (%f,%f,%f), \n", KinE,
		 MagE, IntE);
	  apause();

	  /* printf("%f %f %f \n", , );  */

          pG->U[k][j][i].d = rho0;
        } //end isnan check

	//        printf("%f %f\n", IntE, x1);
	//        printf("%f %f %f\n", pG->U[k][j][i].d, pG->U[k][j][i].M2, x1);

	//printf("%f %f \n", rho0, dcut); getchar();

        di = 1.0/(pG->U[k][j][i].d);
        v1 = pG->U[k][j][i].M1*di;
        v2 = pG->U[k][j][i].M2*di;
        v3 = pG->U[k][j][i].M3*di;
        qsq = v1*v1 + v2*v2 + v3*v3;
        b1 = pG->U[k][j][i].B1c
          + fabs((double)(pG->B1i[k][j][i] - pG->U[k][j][i].B1c));
        b2 = pG->U[k][j][i].B2c
          + fabs((double)(pG->B2i[k][j][i] - pG->U[k][j][i].B2c));
        b3 = pG->U[k][j][i].B3c
          + fabs((double)(pG->B3i[k][j][i] - pG->U[k][j][i].B3c));
        bsq = b1*b1 + b2*b2 + b3*b3;
        p = MAX(Gamma_1*(pG->U[k][j][i].E - 0.5*pG->U[k][j][i].d*qsq
			 - 0.5*bsq), TINY_NUMBER);
        asq = Gamma*p*di;

        rho_max = 100.;
        as_lim = 10.;

        //check for very low density
        pG->U[k][j][i].d = fmax(rho0,  pG->U[k][j][i].d );

        if( (sqrt(asq)>as_lim) || (pG->U[k][j][i].d > rho_max) ) {
	  dns = fmin(rho_max,  pG->U[k][j][i].d);

	  //check for very high density
	  pG->U[k][j][i].d = dns;

	  KinE = 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2)
		      + SQR(pG->U[k][j][i].M3  ))/dns;
     
	  IntE = e0;  // Set this as well to keep c_s reasonable

	  pG->U[k][j][i].E = IntE + KinE + MagE;

	  protd++;
	  prote++;
        }
        else if (IntE < e0) {
          pG->U[k][j][i].E = e0 + KinE + MagE;
          prote++;
        }


        if (sqrt(asq)>as_lim){

	  printf("\n At pos 'asq' (%d,%d,%d) (%f, %f, %f), Den1 = %e, E1 = %e\n", i, j, k, x1,x2,x3,
		 pG->U[k][j][i].d, pG->U[k][j][i].E);
	  printf("asq ,bsq, P, ro : %e %e %e  %e\n", asq, bsq, p, pG->U[k][j][i].d);
	  printf("%e %e", rho0, dcut);
        }



	//        else if (IntE > 1.) {
	//
	//          	pG->U[k][j][i].E = e0 + KinE + MagE;
	//
	//          	printf("torus9 d ,  IntE, KinE, MagE: %f %f %f %f\n", pG->U[k][j][i].d ,  IntE, KinE, MagE);
	//          	printf("At pos (%d,%d,%d) (%f,%f,%f), Den = %f, E = %f\n", i, j, k, x1,x2,x3,pG->U[k][j][i].d, pG->U[k][j][i].E);
	//          	prote++;
	//
	//        } //end  ro and E check


	//        IntE = pG->U[k][j][i].E - KinE - MagE;
	//        Pgas  = (Gamma-1)*IntE;
	//        xRayHeatCool(pG->U[k][j][i].d, Pgas, pG-> dt, pG->xi[k][j][i] );


	//      IntE = pG->U[k][j][i].E - KinE - MagE;
      }
    }
  }

#ifdef MPI_PARALLEL
 if (my_id == 0) {
    printf("\tDivergence @ Orbit %2.3f = %e\n",(pG->time)/61.6, DivB);
    if ((protd+prote) > 0) {
      printf("\tProtection enforced (D=%d,E=%d), Cumulative Mass = %2.5f\n", protd, prote,TotMass);
    }
  }
#else
 printf("\tDivergence @ Orbit %2.3f = %e\n",(pG->time)/61.6, DivB);
  if ((protd+prote) > 0) {
    printf("\tProtection enforced (D=%d,E=%d), Cumulative Mass = %2.5f\n", protd, prote,TotMass);

  }
#endif //PARALLEL

}  //end of Userwork_in_loop



void Userwork_after_loop(MeshS *pM)
{

  //	free_3d_array(pG->xi);

  
  //	free some memory

}

#define MAXIT 100
#define MAXIT_INTERV_EXPAND  5
//xRayHeatCool(const Real dens, const Real Press, const Real xi_in,
//									Real Hx, Real dHxdT)

#ifdef XRAYS
Real rtsafe_energy_eq(Real Press, Real dens, Real xi, Real dt, int* status)
{
  //returns new Pressure
  int j;
  Real df,dx,dxold,f,fh,fl,rat;
  Real temp,xh,xl,rts, x1, x2, x0,dHdT,Hx1,Hx2,
    xacc=1e-4;
  int i;

  x0 = Press;
  x1= Press;
  x2=x1;
  *status = 1;

  for (i =1; i<=MAXIT_INTERV_EXPAND; i++) {	/* loop 1 */
    rat  = 2.;
    x1 /= rat;
    x2 *= rat;
    
    xRayHeatCool(dens, x1, xi, &Hx1, &dHdT, dt);
    xRayHeatCool(dens, x2, xi, &Hx2, &dHdT, dt);

    fl = x1 - x0 -  dt* Hx1/Gamma_1;
    fh = x2 - x0 - dt* Hx2/Gamma_1;

    if (fl == 0.0) return x1;
    if (fh == 0.0) return x2;
    if (fl < 0.0) {
      xl=x1;
      xh=x2;
    } else {
      xh=x1;
      xl=x2;
    }

    if(fl*fh<0.) { /* in the interval  */
      //			printf("cond. met \n");
      *status=1;
      //			printf("%e %e \n", t_th, dt);
      //			apause();
      //			break;

      rts=0.5*(x1+x2);
      dxold=fabs(x2-x1);
      dx=dxold;

      xRayHeatCool(dens, rts, xi, &f, &df, dt);
      f =  rts - x0 -  dt* f/Gamma_1;
      df = 1. - df/Gamma_1;

      for (j=1;j<=MAXIT;j++) { /* loop 2 */
	if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)
	    || (fabs(2.0*f) > fabs(dxold*df))) {
	  dxold=dx;
	  dx=0.5*(xh-xl);
	  rts=xl+dx;
	  if (xl == rts) return rts;
	} else {
	  dxold=dx;
	  dx=f/df;
	  temp=rts;
	  rts -= dx;
	  if (temp == rts) return rts;
	}
	if (fabs(dx) < xacc) return rts;

	xRayHeatCool(dens, rts, xi, &f, &df, dt);
	f = rts - x0 - dt* f/Gamma_1;

	if (f < 0.0) {
	  xl=rts;
	} else {
	  xh=rts;
	}
      } /* loop 2 */

      /* Cap the number of iterations but don't fail */
      printf("%f \n", rts);
      apause();
    }

    else{  /* not in the interval  */

      if ( i==MAXIT_INTERV_EXPAND ){
	status = 0;
	//				printf("cond. has not been met \n");
      }

    } /* interval check */

  } /* for loop 1 */

  return rts;

  //
  //	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
  //		ath_error("interval cannot be expanded further in torus9:rtsafe_energy_eq\n");

}
#endif


#undef MAXIT
#undef MAXIT_INTERV_EXPAND


static void inX1BoundCond(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
  Real x1,x2,x3,rad,rsf,lsf;


#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif


  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {


#ifdef show_debug_messages
	if( (pGrid->U[k][j][i].d < rho0) ||  (pGrid->U[k][j][i].d > 3.) ) {
	  printf("negative density at BC: %e \n", pGrid->U[k][j][i].d );
	}
#endif

	cc_pos(pGrid,is,j,k,&x1,&x2,&x3);

	//    	  printf("%f %f %f\n", pGrid->U[k][j][i].M1, pGrid->U[k][j][i].M2, pGrid->U[k][j][i].M3);
	//   	  printf("%d \n", nghost ); apause();

	rad = sqrt(SQR(x1) + SQR(x3));

	//       printf("%f %f\n", x1, x3);

	// printf("%d %d\n", nghost, is); getchar();

	pGrid->U[k][j][is-i] = pGrid->U[k][j][is];
	pGrid->U[k][j][is-i].M1 = MIN(pGrid->U[k][j][is-i].M1,0.0);

	//    	  	  if (pGrid->U[k][j][is].M1 > 0.) {
	//      	    		pGrid->U[k][j][is-i].M1 = 0.;
	//    	  	  }

	pGrid->U[k][j][is-i].M2 =  pGrid->U[k][j][is].d * sqrt(rad) / (rad-rgToR0);
	pGrid->U[k][j][is-i].M2 =  pGrid->U[k][j][is].d * sqrt(  1 / (rad-rgToR0 ) );




	//    	  pGrid->U[k][j][is-i].d =  fmax(rho0, pGrid->U[k][j][is].d) ;
	//    	      pGrid->U[k][j][is-i].d =  rho0;
	//    	      pGrid->U[k][j][is-i].E =  e0;

	//       printf("%f \n", pGrid->U[k][j][is-i].E);
      }
    }
  }

#ifdef MHD
  /* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost-1; i++) {

	pGrid->B1i[k][j][is-i] = pGrid->B1i[k][j][is];

	rsf = (x1+0.5*pGrid->dx1);
	lsf = (x1-0.5*pGrid->dx1);

	pGrid->B1i[k][j][is-i] = pGrid->B1i[k][j][is]*rsf/lsf;


      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j][is];

        pGrid->B2i[k][j][is-i] = 0;

      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {

        pGrid->B3i[k][j][is-i] = pGrid->B3i[k][j][is];

	//        pGrid->B3i[k][j][is-i] = 0;

      }
    }
  }
#endif /* MHD */

  return;
}


static void diode_outflow_ix3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[ks-k][j][i] = pGrid->U[ks][j][i];
        pGrid->U[ks-k][j][i].M3 = MIN(pGrid->U[ks-k][j][i].M3,0.0);
      }
    }
  }



#ifdef MHD
  /* B1i is not set at i=is-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[ks-k][j][i] = pGrid->B1i[ks][j][i];
      }
    }
  }

  /* B2i is not set at j=js-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[ks-k][j][i] = pGrid->B2i[ks][j][i];
      }
    }
  }

  /* B3i is not set at k=ks-nghost */
  for (k=1; k<=nghost-1; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[ks-k][j][i] = pGrid->B3i[ks][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_ox3(GridS *pGrid)
 *  \brief OUTFLOW boundary conditions, Outer x3 boundary (bc_ox3=2) */

static void diode_outflow_ox3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->U[ke+k][j][i] = pGrid->U[ke][j][i];

        pGrid->U[ke+k][j][i].M3 = MAX(pGrid->U[ke+k][j][i].M3,0.0);

      }
    }
  }

  //  pGrid->U[k][j][is-i].M1 = MIN(pGrid->U[k][j][is-i].M1,0.0);

#ifdef MHD
  /* B1i is not set at i=is-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pGrid->B1i[ke+k][j][i] = pGrid->B1i[ke][j][i];
      }
    }
  }

  /* B2i is not set at j=js-nghost */
  for (k=1; k<=nghost; k++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B2i[ke+k][j][i] = pGrid->B2i[ke][j][i];
      }
    }
  }

  /* k=ke+1 is not a boundary condition for the interface field B3i */
  for (k=2; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[ke+k][j][i] = pGrid->B3i[ke][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}




static void calcProblemParameters(){
  MSOLYR = 1.989e33/YR;
  SGMB = ARAD*CL/4;

  MBH_in_gram = M2Msun*MSUN;
  Rg = 2*GRV*MBH_in_gram/pow(CL,2);
  r0*=Rg;

  Nsc=nc0;
  Tx=1.0e8; //Compton temperature
  Tgmin = 2.3;

  Dsc = nc0*M_MW*	M_U;
  xm =1.;
  r_in *= xm;

  Rsc = xm* r0;
  rgToR0 = Rg/Rsc;
  //   	rgToR0 = 0.;

  //   	printf("Rsc = %e \n", Rsc); apause();

  Usc = sqrt(GRV*MBH_in_gram/Rsc);
  Time0=sqrt( pow(Rsc,3) /GRV/MBH_in_gram);
  Evir0=Dsc*GRV*MBH_in_gram/Rsc;  //(erg/cm^3)
  Psc = Evir0;
  Esc = Psc;

  Tvir0=GRV*MBH_in_gram/RGAS/Rsc*M_MW;
  Lamd0=Evir0/Time0;  // (erg/cm^3/s)
  Ledd_e = 4.*PI *CL*GRV*MBH_in_gram/KPE;
  Lx = fx * F2Fedd*Ledd_e;
  Luv= fuv *F2Fedd*Ledd_e;

  //  torus parameters

  nAd = 1.0/Gamma_1;
  q1 = 2.0*(q-1.0);

#ifdef HAWL

  a1 = pow(xm, q1) /q1;
  Ctor = a1* pow(r_in,-q1) - xm/(r_in-rgToR0);
  Kbar = (Ctor + xm/(xm-rgToR0) - 1./q1) /(nAd+1);


#endif

#ifdef SOLOV
  Kbar_Sol =(1.5 + Br0_Sol)/(1 + nAd);
  Kbar = Kbar_Sol;
#endif

  //    	printf( " %e %e %e %e\n", Dsc, nc0, M_MW, M_U);
  printf("Ctor, Kbar= %e %e \n", Ctor, Kbar );
}

static void printProblemParameters(){
  //	#ifdef SOLOV
  //	  printf("solov is not implemented"); getchar();
  //	#endif

  printf("parameters \n");

  printf("F2Fedd = %e \n", F2Fedd);
  printf("nc0 = %e \n", nc0);
  printf("R0 = %e * R_g \n", r0);
  printf("Ctor, Kbar= %e %e \n", Ctor, Kbar );

  //	getchar();
}



/* ============================================================================== */
/*                           Boundary conditions */
/* ============================================================================== */


void bvals_tau_init(MeshS *pM){

 GridS *pG;
 DomainS *pD;
 int l,m,n,nx1t,nx2t,nx3t,size;
 int x1cnt=0, x2cnt=0, x3cnt=0; /* Number of words passed in x1/x2/x3-dir. */
 
/* Figure out largest size needed for send/receive buffers with MPI ----------*/ 

 if (pM->NLevels > 1)
   ath_error("[bval_rray_init]: xray module does not support SMR\n");
#ifdef MPI_PARALLEL

 pD = &(pM->Domain[0][0]);  /* ptr to Domain */
 pG = pD->Grid; /* ptr to Grid */
 
    for (n=0; n<(pD->NGrid[2]); n++){
      for (m=0; m<(pD->NGrid[1]); m++){
	for (l=0; l<(pD->NGrid[0]); l++){
/* x1cnt is surface area of x1 faces */
	if(pD->NGrid[0] > 1){
	  nx2t = pD->GData[n][m][l].Nx[1];
	  if(nx2t > 1) nx2t += 1;

	  nx3t = pD->GData[n][m][l].Nx[2];
	  if(nx3t > 1) nx3t += 1;
          if(nx2t*nx3t > x1cnt) x1cnt = nx2t*nx3t;
	}
	
/* x2cnt is surface area of x2 faces */
	if(pD->NGrid[1] > 1){
	  nx1t = pD->GData[n][m][l].Nx[0];
	  if(nx1t > 1) nx1t += 2*nghost;
	  
	  nx3t = pD->GData[n][m][l].Nx[2];
	  if(nx3t > 1) nx3t += 1;
          if(nx1t*nx3t > x2cnt) x2cnt = nx1t*nx3t;
	}

/* x3cnt is surface area of x3 faces */
	if(pD->NGrid[2] > 1){
	  nx1t = pD->GData[n][m][l].Nx[0];
	  if(nx1t > 1) nx1t += 2*nghost;

	  nx2t = pD->GData[n][m][l].Nx[1];
	  if(nx2t > 1) nx2t += 2*nghost;
          if(nx1t*nx2t > x3cnt) x3cnt = nx1t*nx2t;
	}
	}
      }
    }
    size = x1cnt > x2cnt ? x1cnt : x2cnt;
    size = x3cnt >  size ? x3cnt : size;

    size *= nghost; /* Multiply by the third dimension */
    if (size > 0) {
    if((send_buf = (double**)calloc_2d_array(2,size,sizeof(double))) == NULL)
      ath_error("[bvals_init]: Failed to allocate send buffer\n");

    if((recv_buf = (double**)calloc_2d_array(2,size,sizeof(double))) == NULL)
      ath_error("[bvals_init]: Failed to allocate recv buffer\n");
  }

  if((recv_rq = (MPI_Request*) calloc_1d_array(2,sizeof(MPI_Request))) == NULL)
    ath_error("[bvals_init]: Failed to allocate recv MPI_Request array\n");
  if((send_rq = (MPI_Request*) calloc_1d_array(2,sizeof(MPI_Request))) == NULL)
    ath_error("[bvals_init]: Failed to allocate send MPI_Request array\n");

 #endif /* MPI_PARALLEL */

  if(pG->Nx[0] > 1) {
    if(apply_ix1 == NULL){
      apply_ix1 = boundCondOptDepthLike_ix1;      
    }
    else {
      ath_perr(-1,"[bvals_xray_init: ix1] error \n");      
      exit(EXIT_FAILURE);
    }
    if(apply_ox1 == NULL){
      apply_ox1 = boundCondOptDepthLike_ox1;      
    }
    else {
      ath_perr(-1,"[bvals_xray_init: ox1] error \n");      
      exit(EXIT_FAILURE);
    }
  }

  if(pG->Nx[1] > 1) {
    if(apply_ix2 == NULL){
      apply_ix2 = boundCondOptDepthLike_ix2;      
    }
    else {
      ath_perr(-1,"[bvals_xray_init: ix2] error \n");      
      exit(EXIT_FAILURE);
    }
    if(apply_ox2 == NULL){
      apply_ox2 = boundCondOptDepthLike_ox2;      
    }
    else {
      ath_perr(-1,"[bvals_xray_init: ox2] error \n");      
      exit(EXIT_FAILURE);
    }
  }
 
  if(pG->Nx[2] > 1) {
    if(apply_ix3 == NULL){
      apply_ix3 = boundCondOptDepthLike_ix3;      
    }
    else {
      ath_perr(-1,"[bvals_xray_init: ix3] error \n");      
      exit(EXIT_FAILURE);
    }
    if(apply_ox3 == NULL){
      apply_ox3 = boundCondOptDepthLike_ox3;      
    }
    else {
      ath_perr(-1,"[bvals_xray_init: ox3] error \n");      
      exit(EXIT_FAILURE);
    }
  }

  return;
}



void bvals_tau(DomainS *pD)
{
  GridS *pGrid = (pD->Grid);

  
#ifdef MPI_PARALLEL
  //  int cnt1, cnt2, cnt3, cnt, ierr, ;
  int cnt,ierr,mIndex;
#endif /* MPI_PARALLEL */

/*--- Step 1. ------------------------------------------------------------------
 * Boundary Conditions in x1-direction */

int  mpi1=1;
 while( mpi1==0 );
  
  if (pGrid->Nx[0] > 1){

#ifdef MPI_PARALLEL

    cnt = nghost*(pGrid->Nx[1])*(pGrid->Nx[2]);

/* MPI blocks to both left and right */
    if (pGrid->rx1_id >= 0 && pGrid->lx1_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx1_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx1_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_tau_ix1(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx1_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_tau_ox1(pGrid);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx1_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_tau_ix1(pGrid);
      if (mIndex == 1) unpack_tau_ox1(pGrid);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_tau_ix1(pGrid);
      if (mIndex == 1) unpack_tau_ox1(pGrid);

    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx1_id >= 0 && pGrid->lx1_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx1_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));
 
      /* pack and send data R */
      pack_tau_ox1(pGrid);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx1_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*apply_ix1)(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_tau_ox1(pGrid);

    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx1_id < 0 && pGrid->lx1_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx1_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_tau_ix1(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx1_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*apply_ox1)(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_tau_ix1(pGrid);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx1_id < 0 && pGrid->lx1_id < 0) {
      (*apply_ix1)(pGrid);
      (*apply_ox1)(pGrid);

    }

  }

/*--- Step 2. ------------------------------------------------------------------
 * Boundary Conditions in x2-direction */

  if (pGrid->Nx[1] > 1){

#ifdef MPI_PARALLEL

    cnt = (pGrid->Nx[0] + 2*nghost)*nghost*(pGrid->Nx[2]);

/* MPI blocks to both left and right */
    if (pGrid->rx2_id >= 0 && pGrid->lx2_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx2_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx2_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_tau_ix2(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx2_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_tau_ox2(pGrid);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx2_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_tau_ix2(pGrid);
      if (mIndex == 1) unpack_tau_ox2(pGrid);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_tau_ix2(pGrid);
      if (mIndex == 1) unpack_tau_ox2(pGrid);

    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx2_id >= 0 && pGrid->lx2_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx2_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_tau_ox2(pGrid);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx2_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*apply_ix2)(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_tau_ox2(pGrid);

    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx2_id < 0 && pGrid->lx2_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx2_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_tau_ix2(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx2_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*apply_ox2)(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_tau_ix2(pGrid);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx2_id < 0 && pGrid->lx2_id < 0) {
      (*apply_ix2)(pGrid);
      (*apply_ox2)(pGrid);
    }

/* shearing sheet BCs; function defined in problem generator */
#ifdef SHEARING_BOX
    get_myGridIndex(pD, myID_Comm_world, &myL, &myM, &myN);
    if (myL == 0) {
      ShearingSheet_grav_ix1(pD);
    }
    if (myL == (pD->NGrid[0]-1)) {
      ShearingSheet_grav_ox1(pD);
    }
#endif

  }

/*--- Step 3. ------------------------------------------------------------------
 * Boundary Conditions in x3-direction */

  if (pGrid->Nx[2] > 1){

#ifdef MPI_PARALLEL

    cnt = (pGrid->Nx[0] + 2*nghost)*(pGrid->Nx[1] + 2*nghost)*nghost;

/* MPI blocks to both left and right */
    if (pGrid->rx3_id >= 0 && pGrid->lx3_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx3_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx3_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_tau_ix3(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx3_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_tau_ox3(pGrid);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx3_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_tau_ix3(pGrid);
      if (mIndex == 1) unpack_tau_ox3(pGrid);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_tau_ix3(pGrid);
      if (mIndex == 1) unpack_tau_ox3(pGrid);

    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx3_id >= 0 && pGrid->lx3_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx3_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_tau_ox3(pGrid);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx3_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*apply_ix3)(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_tau_ox3(pGrid);

    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx3_id < 0 && pGrid->lx3_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx3_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_tau_ix3(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx3_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*apply_ox3)(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_tau_ix3(pGrid);
    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx3_id < 0 && pGrid->lx3_id < 0) {
      (*apply_ix3)(pGrid);
      (*apply_ox3)(pGrid);
    }

  }

  return;
}





#ifdef MPI_PARALLEL  

/*! \fn static void pack_tau_ix1(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Inner x1 boundary */

static void pack_tau_ix1(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

/* Pack only tau into send buffer */
  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=is+(nghost-1); i++){
        *(pSnd++) = pG->tau_e[k][j][i];
      }
    }
  }

  return;
}
/*----------------------------------------------------------------------------*/
/*! \fn static void pack_tau_ox1(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Outer x1 boundary */

static void pack_tau_ox1(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

/* Pack only tau into send buffer */
  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=ie-(nghost-1); i<=ie; i++){
        *(pSnd++) = pG->tau_e[k][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_tau_ix2(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Inner x2 boundary */

static void pack_tau_ix2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

/* Pack only tau into send buffer */
  for (k=ks; k<=ke; k++){
    for (j=js; j<=js+(nghost-1); j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        *(pSnd++) = pG->tau_e[k][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_tau_ox2(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Outer x2 boundary */

static void pack_tau_ox2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

/* Pack only tau into send buffer */

  for (k=ks; k<=ke; k++){
    for (j=je-(nghost-1); j<=je; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        *(pSnd++) = pG->tau_e[k][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_tau_ix3(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Inner x3 boundary */

static void pack_tau_ix3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

/* Pack only tau into send buffer */

  for (k=ks; k<=ks+(nghost-1); k++){
    for (j=js-nghost; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        *(pSnd++) = pG->tau_e[k][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_tau_ox3(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Outer x3 boundary */

static void pack_tau_ox3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

/* Pack only tau into send buffer */

  for (k=ke-(nghost-1); k<=ke; k++){
    for (j=js-nghost; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        *(pSnd++) = pG->tau_e[k][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_tau_ix1(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Inner x1 boundary */

static void unpack_tau_ix1(GridS *pG)
{
  int is = pG->is;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);

/* Manually unpack the data from the receive buffer */

  for (k=ks; k<=ke; k++){
//    for (j=js; j<=js; j++){
//MODIFIED BY HAO GONG
    for (j=js; j<=je; j++){
      for (i=is-nghost; i<=is-1; i++){
        pG->tau_e[k][j][i] = *(pRcv++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_tau_ox1(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Outer x1 boundary */

static void unpack_tau_ox1(GridS *pG)
{
  int ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);

/* Manually unpack the data from the receive buffer */

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=ie+1; i<=ie+nghost; i++){
        pG->tau_e[k][j][i] = *(pRcv++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_tau_ix2(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Inner x2 boundary */

static void unpack_tau_ix2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);

/* Manually unpack the data from the receive buffer */

  for (k=ks; k<=ke; k++){
    for (j=js-nghost; j<=js-1; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        pG->tau_e[k][j][i] = *(pRcv++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_tau_ox2(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Outer x2 boundary */

static void unpack_tau_ox2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);

/* Manually unpack the data from the receive buffer */

  for (k=ks; k<=ke; k++){
    for (j=je+1; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        pG->tau_e[k][j][i] = *(pRcv++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_tau_ix3(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Inner x3 boundary */

static void unpack_tau_ix3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks;
  int i,j,k;
  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);

/* Manually unpack the data from the receive buffer */

  for (k=ks-nghost; k<=ks-1; k++){
    for (j=js-nghost; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        pG->tau_e[k][j][i] = *(pRcv++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_tau_ox3(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Outer x3 boundary */

static void unpack_tau_ox3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ke = pG->ke;
  int i,j,k;
  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);

/* Manually unpack the data from the receive buffer */

  for (k=ke+1; k<=ke+nghost; k++){
    for (j=js-nghost; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        pG->tau_e[k][j][i] = *(pRcv++);
      }
    }
  }

  return;
}
#endif /* MPI_PARALLEL */
static void boundCondOptDepthLike_ix1(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef RADIAL_RAY_TEST
  Real x1,x2,x3, rsph;
#endif 
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
	
	pGrid->tau_e[k][j][is-i] = 0.;

#ifdef RADIAL_RAY_TEST	 	 
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	rsph =  sqrt(pow(x1 - rRadSrc[0], 2) + pow(x3 - rRadSrc[1] ,2));	
	pGrid->tau_e[k][j][is-i] = rsph;	
#endif
	
	/* pGrid->tau_e[k][j][is-i] = pGrid->tau_e[k][j][is+(i-1)]; */

      }
    }
  }

  return; 
}
static void boundCondOptDepthLike_ox1(GridS *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->tau_e[k][j][ie+i] = pGrid->tau_e[k][j][ie-(i-1)];
      }
    }
  }

  return;

}

static void boundCondOptDepthLike_ix2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->tau_e[k][js-j][i]    =  pGrid->tau_e[k][js+(j-1)][i];
      }
    }
  }

  return;
}
static void boundCondOptDepthLike_ox2(GridS *pGrid)
{
   int is = pGrid->is, ie = pGrid->ie;
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->tau_e[k][je+j][i] = pGrid->tau_e[k][je-(j-1)][i];
      }
    }
  }

  return;
}


static void boundCondOptDepthLike_ix3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->tau_e[ks-k][j][i] = pGrid->tau_e[ks+(k-1)][j][i];
      }
    }
  }

  return;
}
static void boundCondOptDepthLike_ox3(GridS *pGrid)
{
   int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->tau_e[ke+k][j][i] = pGrid->tau_e[ke-(k-1)][j][i];
      }
    }
  }

  return;
}


/* ============================================================================= */
/*                      end of boundary conditions */
/* ============================================================================= */


void testRayTracings(GridS *pG){
  Real r, t, z, l,dl, ri, tj, zk,
    tau=0,dtau=0,tau_ghost,xi=0,x1,x2,x3;      
   
    int i,j,k,is,ie,js,je,ks,ke, il, iu, jl,ju,kl,ku,ip,jp,kp,knew,
    my_id=0;

    Real abs_cart_norm, cart_norm[3], cyl_norm[3], xyz_pos[3], rtz_pos[3], xyz_p[3],
      res[1];
    int ijk_cur[3],iter,iter_max, lr, ir=0, i0,j0,k0;
    short nroot;

    Real xyz_in[3], radSrcCyl[3], dist, sint, cost, tmp;

       
  printf("hello from test raytracing \n");
  
  is = pG->is;
  ie = pG->ie;
  js = pG->js;
  je = pG->je;
  ks = pG->ks;
  ke = pG->ke;
  il = is - nghost*(ie > is);
  jl = js - nghost*(je > js);
  kl = ks - nghost*(ke > ks);
  iu = ie + nghost*(ie > is);
  ju = je + nghost*(je > js);
  ku = ke + nghost*(ke > ks);

  /* start point */
  radSrcCyl[0]=1.;
  radSrcCyl[1]=-1;
  radSrcCyl[2]=4.;
  

  lr=celli(pG,radSrcCyl[0], 1./pG->dx1, &i0, &ir);
  lr=cellj(pG,radSrcCyl[1], 1./pG->dx2, &j0, &ir);
  lr=cellk(pG,radSrcCyl[2], 1./pG->dx3, &k0, &ir);

  ijk_cur[0] = i0;
  ijk_cur[1] = j0;
  ijk_cur[2] = k0;
  cc_pos(pG,i0,j0,k0,&rtz_pos[0],&rtz_pos[1],&rtz_pos[2]);
  for(i=0;i<=2;i++) radSrcCyl[i]=rtz_pos[i];
  coordCylToCart(xyz_in, radSrcCyl, cos(radSrcCyl[1]), sin(radSrcCyl[1]) );
  
  /* end point */
  x1 = 6.;
  x2 = 1.;
  x3 = -4.;
  lr=celli(pG, x1, 1./pG->dx1, &ip, &ir);
  lr=cellj(pG, x2, 1./pG->dx2, &jp, &ir);
  lr=cellk(pG, x3, 1./pG->dx3, &kp, &ir);    

  
  cc_pos(pG,ip,jp,kp,&rtz_pos[0],&rtz_pos[1],&rtz_pos[2]);
  sint = sin(rtz_pos[1]);
  cost = cos(rtz_pos[1]);
  coordCylToCart(xyz_p, rtz_pos, cost, sint);

  for(i=0;i<=2;i++) cart_norm[i]= xyz_p[i]-xyz_in[i];

  
  dist = absValueVector(cart_norm);
  abs_cart_norm = sqrt(pow(cart_norm[0], 2)+pow(cart_norm[1], 2)+pow(cart_norm[2], 2));
  for(i=0;i<=2;i++) cart_norm[i] = cart_norm[i]/ abs_cart_norm;  

  cartVectorToCylVector(cyl_norm, cart_norm, cos(radSrcCyl[1]), sin(radSrcCyl[1]));		    	     
 

  /* starting point for ray-tracing */
  rtz_pos[0]=radSrcCyl[0];
  rtz_pos[1]=radSrcCyl[1];
  rtz_pos[2]=radSrcCyl[2];
  sint = sin(rtz_pos[1]);
  cost = cos(rtz_pos[1]);
  
  for(i=0;i<=2;i++) xyz_pos[i]=xyz_in[i];

  iter_max = ke*je*ke;
  iter_max=50;
  
  for (iter=0; iter<=iter_max; iter++){
    
    traceGridCell(pG, res, ijk_cur, xyz_pos, rtz_pos,
		  cart_norm, cyl_norm, &nroot);

    /* /\* get a more accurate aim, re-target *\/ */
    /* normVectBToA(cart_norm, xyz_p, xyz_pos, &tmp); */
    
    /* find new cyl coordinates: */
    cc_pos(pG,ijk_cur[0],ijk_cur[1],ijk_cur[2],&rtz_pos[0],&rtz_pos[1],&rtz_pos[2]);

    ri = sqrt(pow(xyz_pos[0],2)+pow(xyz_pos[1],2));
    tj = atan2(xyz_pos[1], xyz_pos[0]);
    zk =xyz_pos[2];

    /* we need only sign(s) of cyl_norm */
    cyl_norm[0] = cart_norm[0]*cos(tj) + cart_norm[1]*sin(tj);
    cyl_norm[1] = -cart_norm[0]*sin(tj) + cart_norm[1]*cos(tj);

    /* get a non-normilized direction in cyl coordinates */    
    //cartVectorToCylVector(cyl_norm, cart_norm, cos(rtz_pos[1]),sin(rtz_pos[1]));
    //printf("\n %f %f \n", cyl_norm[1], res[1]);
    
    // cartVectorToCylVector(cyl_norm, cart_norm, cos(rtz_pos[1]),sin(rtz_pos[1]));
    
    
    dl = res[0];
    l += dl;
    
    printf("ip,jp,kp, ijk,l,dl: %d %d %d %d %d %d  %f %f \n",
	   ip,jp,kp, ijk_cur[0], ijk_cur[1], ijk_cur[2],l, dl);

    printf("xyz_dest= %f %f %f xyz_pos= %f %f %f \n", xyz_p[0],xyz_p[1],xyz_p[2],
	   xyz_pos[0],xyz_pos[1],xyz_pos[2]);
    
    /* s = sqrt(pow(A[0]-B[0], 2)+pow(A[1]-B[1], 2)+pow(A[2]-B[2], 2)); */
 
    if( ijk_cur[0]==ip &&  ijk_cur[1]==jp &&  ijk_cur[2]==kp ){  /* || s > s_prev){ */

      cc_pos(pG,ijk_cur[0], ijk_cur[1], ijk_cur[2], &x1,&x2,&x3);
      printf("compare: %f %f \n", l, dist);
      break;
    }
    /* s_prev = s; */    
    printf("iteration: %d \n", iter);       
   
  }
  
  printf(" test raytracing done ");
  getchar();
  return;
}


void traceGridCell(GridS *pG, Real *res, int *ijk_cur, Real *xyz_pos,
		   Real* rtz_pos, const Real *cart_norm, const Real *cyl_norm,
		   short* nroot) {
  /* ijk_cur is the 3d index of the current cell */
  /*  xpos is the (x,y,z) exact cart. position, should be located at one of the cell's boundaries; */
  /*  x1x2x3 is the 3d - r,t,z -center coordinates of the cell to which the 
      above boundaries belong */
  // traces a cingle cell; cnorm is in Cart. coordinates;

  int change_indx, icur, jcur, kcur;

  Real a,b,c, rfc, zfc, d,tfc,sint,cost,
    L=HUGE_NUMBER, L1=HUGE_NUMBER, L2=HUGE_NUMBER,
     xcur, ycur, zcur, one_over_a;
     int is,ie,js,je,ks,ke;
     *nroot=0;
     
     is = pG->is;
     ie = pG->ie;
     js = pG->js;
     je = pG->je;
     ks = pG->ks;
     ke = pG->ke;
  
     icur= ijk_cur[0]; /* --- icur etc. are on cyl. grid ----*/
     jcur= ijk_cur[1];
     kcur= ijk_cur[2];

     xcur = xyz_pos[0];   /*  are on Cart. grid */
     ycur=  xyz_pos[1];
     zcur=  xyz_pos[2];
  
      /*  1)  intersection with cylinders  */     
     
     if(cyl_norm[0]>0){       
       rfc = pG->MinX[0] + ((Real)(icur - pG->is)+1.)*pG->dx1;
     }
     else if (cyl_norm[0]<0.){
       rfc = pG->MinX[0] + ((Real)(icur - pG->is))*pG->dx1;
     }
       else{ /* nr==0: almost never happens*/
	 goto z_planes;
     }
     change_indx = 0;
       
     a = pow(cart_norm[0],2) + pow(cart_norm[1],2);
     b = 2.*(xyz_pos[0]*cart_norm[0]+xyz_pos[1]*cart_norm[1]); 
     c = pow(xyz_pos[0],2)  + pow(xyz_pos[1],2) - pow(rfc,2);
     d = pow(b,2)-4.*a*c;

     if (d<0){ /* possibly at a grazing angle to a cylinder */
       printf("Discr<0. in traceACell\n");
       goto z_planes;
     }

      one_over_a =copysign(1./fmax(fabs(a), tiny), a);
      L1 = fabs(-b+sqrt(d))/2.*one_over_a;
      L2 = fabs(-b-sqrt(d))/2.*one_over_a;
      L = fmin(fabs(L1),fabs(L2));
      nroot++; 
      printf("Lr= %f nr= %f \n", L, cyl_norm[0]);
      
      /*  3)  intersection with planes z_k=const  */
    z_planes:
      if(cyl_norm[2]>0){
	zfc = pG->MinX[2] + ((Real)(kcur - pG->ks) +1.)*pG->dx3;
      }
      else if (cyl_norm[2]<0.){
	zfc = pG->MinX[2] + ((Real)(kcur - pG->ks))*pG->dx3;	
      }
      else{  /* nz==0 */
	goto phi_planes;	 
      }
      L1 = fabs( (zfc - xyz_pos[2])/cart_norm[2]);

     

      if (L1 < L){ /*min of Lr or Lz */
       change_indx = 2;
       L=L1;
       nroot++;

       printf("Lz= %f , nr= %f \n", L1, cyl_norm[1]);
      }
   
      /*  3)  intersection with planes phi_j=const */
    phi_planes:
      if(cyl_norm[1]>0){
       tfc = pG->MinX[1] + ((Real)(jcur - pG->js)+1)*pG->dx2;
      }
      else if(cyl_norm[1]<0.){
       tfc = pG->MinX[1] + ((Real)(jcur - pG->js))*pG->dx2;
      }
      else{
      	return;
      }
      sint = sin(tfc);
      cost = cos(tfc);
      d = cart_norm[1]*cost - cart_norm[0]*sint;
      
      if(d != 0.){
      	L1  = fabs((xyz_pos[0]*sint -xyz_pos[1]*cost)/d);
      	if (L1 < L){ /*min of Lr or Lz */
      	  change_indx = 1;
      	  L=L1;
          nroot++;	  
	  printf("Lt= %f \n", L1);
      	}
      }

      res[0] = L;

      xyz_pos[0] += cart_norm[0]*L; 
      xyz_pos[1] += cart_norm[1]*L;
      xyz_pos[2] += cart_norm[2]*L;

      ijk_cur[change_indx] += (int)copysign(1, cyl_norm[change_indx]);

} //end  traceCell

void coordCylToCart(Real *xyz, const Real *rtz,
			   const Real cost, const Real sint ){ 
    xyz[0] = rtz[0]*cost;
    xyz[1] = rtz[0]*sint;
    xyz[2] = rtz[2];
}

void cartVectorToCylVector(Real *cyl, const Real *cart,
			   const Real cost, const Real sint ){
      cyl[0] = cart[0]*cost + cart[1]*sint;
      cyl[1] = -cart[0]*sint + cart[1]*cost;
      cyl[2] = cart[2];
}

void vectBToA(Real* v, const Real* A, const Real*B){
    v[0]= A[0]-B[0];
    v[1]= A[1]-B[1];
    v[2]= A[2]-B[2];      
}

Real absValueVector(const Real* A){
  return( (Real)(sqrt(pow(A[0], 2)+pow(A[1], 2)+pow(A[2], 2))) );
}

void normVectBToA(Real* v, const Real* A, const Real*B, Real* tmp){
    v[0]= A[0]-B[0];
    v[1]= A[1]-B[1];
    v[2]= A[2]-B[2];
    *tmp = absValueVector(v);
    v[0] /= *tmp;
    v[1] /= *tmp;
    v[2] /= *tmp;
}


