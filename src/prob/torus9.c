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

//#define HAWL

#define SOLOV

//#if !defined(HAWL)
//	#define SOLOV
//#endif


//****************************************

static void inX1BoundCond(GridS *pGrid);
static void diode_outflow_ix3(GridS *pGrid);
static void diode_outflow_ox3(GridS *pGrid);

void plot(MeshS *pM, char name[16]);
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



static Real
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

	a1_Sol =1,
	a2_Sol = 1,
	b1_Sol =3.,

	w0_Sol =2.,
	w1_Sol =0.;
#endif


//----------------------------

//extern x1GravAcc;

/*! \fn static Real grav_pot(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential */


void getIndexOfCellByZ(GridS *pG, float z, int  *k){
	Real xcur;
//	Real dx = pG->dx3;
//	xcur =  (z - pG->MinX[2]  + 0.5*pG->dx3)/pG->dx3+ pG->ks;
	xcur =  (z - pG->MinX[2])/pG->dx3+ pG->ks;
	*k = (int)(xcur);
//	*k = round(xcur);
}


////	int unsigned long n, float x, unsigned long *j
//	unsigned long ju,jm,jl;
////	float xx[],
//	int ascnd;

//
//	jl=0;
//	ju=n+1;
//	ascnd=(xx[n] > xx[1]);
//
//	while (ju-jl > 1) {
//		jm=(ju+jl) >> 1;
//
////		xcur = pG->MinX[0] + ((Real)(i - pG->is) + 0.5)*pG->dx1;
//
//		xcur  = pG->MinX[numOfCoord] + ((Real)(jm - pG->startIndex))*pG->dx3;
//
//		if (x > xcur == ascnd)
//			jl=jm;
//		else
//			ju=jm;
//	}
//	*j=jl;

//	*px1 = pG->MinX[0] + ((Real)(i - pG->is) + 0.5)*pG->dx1;
//	*px2 = pG->MinX[1] + ((Real)(j - pG->js) + 0.5)*pG->dx2;
//	*px3 = pG->MinX[2] + ((Real)(k - pG->ks) + 0.5)*pG->dx3;


void apause(){
	printf("pause.. press any key ss... ");
	getchar();
}

#ifdef XRAYS

//void optDepthFunctions(MeshS *pM){

void optDepthFunctions(GridS *pG){
	// it is assumed that the source is located on the axis of symmetry

//	GridS *pG=pM->Domain[0][0].Grid;

	Real r, t, z, dl, sToX, sToZ,tau,dtau, xi,x1,x2,x3, ro,rad, nnorm,norm[2],rCur[2],rNextX[2],rNextZ[2],colDens, \
			den, xBndMax,zBndMin,zBndMax,x_is,z_is,res, sToX1, sToZ1, rsph;
	int i,j,k,is,ie,js,je,ks,ke, il, iu, jl,ju,kl,ku,m,ip,jp,kp,knew;
	float sgnf=1.;

	is = pG->is; ie = pG->ie;
    js = pG->js; je = pG->je;
    ks = pG->ks; ke = pG->ke;
    il = is - nghost*(ie > is);
    jl = js - nghost*(je > js);
    kl = ks - nghost*(ke > ks);
    iu = ie + nghost*(ie > is);
    ju = je + nghost*(je > js);
    ku = ke + nghost*(ke > ks);

//    getIndexOfCellByZ(pG, z, &k);
//    cc_pos(pG,is,js,k,&x1,&x2,&x3);
//    printf("%i %f %f \n", k, x3,pG->dx3); apause();


    for (kp=ks; kp<=ke; kp++) {

      for (jp=js; jp<=je; jp++) {

    	   for (ip= is+1;     ip<=ie; ip++) {

    		pG->tau_e[kp][jp][ip]=0.;

    		cc_pos(pG,ip,jp,kp,&x1,&x2,&x3);

        	nnorm = sqrt(pow(x1-rRadSrc[0], 2)+pow( x3-rRadSrc[1], 2));
        	norm[0]=(x1-rRadSrc[0])/nnorm;
        	norm[1]=(x3-rRadSrc[1])/nnorm;
        	colDens =0.;

        	x_is = pG->MinX[0] + 0.5*pG->dx1;
        	z_is = norm[1]/fmax(norm[0],tiny)*(x_is - rRadSrc[0]) + rRadSrc[1];
        	rCur[0]=x_is;
        	rCur[1]=z_is;
 	 	xBndMax = pG->MaxX[0];
      	zBndMin = pG->MinX[2];
      	zBndMax = pG->MaxX[2];
      	 rsph =  sqrt(pow(rCur[0] - rRadSrc[0], 2) + pow(rCur[1] - rRadSrc[1] ,2));
      	 tau = 0.;
      	 k = (int) ((z_is - pG->MinX[2])/pG->dx3 + pG->ks);
       	 i=is;

    	 	 while (i < ie)	{

        		// crossing face, xf,i+1
        		rNextX[0] = pG->MinX[0] +    (i+1 - pG->is)*pG->dx1;
        		rNextX[1]= norm[1]/norm[0]*(rNextX[0]-rRadSrc[0])+rRadSrc[1];

        		if( rNextX[0]>pG->MaxX[0] ||  rNextX[1]>pG->MaxX[2] || rNextX[1]<pG->MinX[2] ){
        			break;
        		}

        		sToX = sqrt(pow(rNextX[0] - rCur[0], 2) + pow(rNextX[1] - rCur[1] ,2));
        		sToX1 = sqrt(pow(rNextX[0] - rRadSrc[0], 2) + pow(rNextX[1] - rRadSrc[1] ,2));

        		knew = k +( copysign(1., norm[1]) );
        		rNextZ[1] = pG->MinX[2] + (  (Real) (knew - pG->ks)  )*pG->dx3;

        		if (norm[1] != 0){
        			rNextZ[0] =(norm[0]/norm[1]) * ( rNextZ[1]-rRadSrc[1] )+rRadSrc[0];
        		} else{
        			rNextZ[0] =norm[0]/fmax(norm[1],tiny)*(rNextZ[1]-rRadSrc[1])+rRadSrc[0];
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

        		dtau = dl * den ;

        		tau  += dtau;

        		if ((k >=ke) || (k<ks) || (i==ip)){
        			break;
        		}

        	}	/* =========== end of ray tracing =========== */

//    	 	 printf("%e %e \n", Rsc, Nsc); apause();


    	 	pG->tau_e[kp][jp][ip]= Rsc*KPE*Dsc*tau;

    	 	rsph =  sqrt(pow(x1 - rRadSrc[0], 2) + pow(x3 - rRadSrc[1] ,2));

    	 	if ((pG->U[k][jp][i].d > tiny) || (rsph > tiny)){

    	 		xi = Lx / (Nsc* fmax(pG->U[ kp ][ jp ][ ip ].d , tiny))/ pow( Rsc*rsph, 2 );

    	 	}
    	 	else{
    	 		printf("quite possibly something is negative in opticalProps ");
    	 		apause();
    	 	}
    	 	pG->xi[kp][jp][ip]  = xi;

    	 	pG->xi[kp][jp][ip]  *=  exp(- (pG->tau_e[kp][jp][ip]) );

//    	printf("%e \n", 	exp(- (pG->tau_e[kp][jp][ip]) ));
//    	 	printf("%e %e \n", Rsc, Dsc); apause();

//    	 	printf("tau = %f  %f \n",  (pG->tau_e[kp][jp][ip]), rsph*Rsc*KPE*Dsc*pG->U[kp][jp][ip].d );

        } //i = x
//    	   apause();
//    	   if (pG->xi[kp][jp][ip]<0.){
//    		   printf("xi = %f %d %d %d\n",  pG->xi[kp][jp][ip], kp, jp, ip );
//    	   }
//    	   printf("%f %d %d %d\n",  pG->tau_e[kp][jp][ip], kp, jp, ip);
//    	   	   apause();
      }	// j = phi
    }  //k = z

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
	int i,j,k,is,ie,js,je,ks,ke,nx1,nx2,nx3,nr,nz;
	char cwd[1024];
	getcwd(cwd, sizeof(cwd));
	GridS *pG=pM->Domain[0][0].Grid;

	printf("WORKING DIR: %s \n" ,  getcwd(cwd, sizeof(cwd) ) );
	FILE *f = fopen("athena_plot.tmp", "w");

	if (f == NULL){

		printf("Error opening file!\n");
	    exit(1);
	}
	is = pG->is; ie = pG->ie;
	js = pG->js; je = pG->je;
	ks = pG->ks; ke = pG->ke;
	nr= (ie-is)+1;
	nz=(ke-ks)+1;
	nx1 = (ie-is)+1 + 2*nghost;
	nx2 = (je-js)+1 + 2*nghost;
	nx3 = (ke-ks)+1 + 2*nghost;

//	printf("%d %d %d %d %d %d %d\n", nx1,nx2,nx3,nghost, (ie-is)+1 ,(je-js)+1,ie); apause();

	j=js;
	fprintf(f, "%d  %d\n", nr, nz );

	for (k=ks; k<=ke; k++) {
	   for (i=is; i<=ie; i++) {
			 if (strcmp(name, "d") == 0){
				 fprintf(f, "%f\n", pG->U[k][j][i].d );
			 }
			 if (strcmp(name, "E") == 0){
				 fprintf(f, "%f\n", pG->U[k][j][i].E );
			 }
#ifdef XRAYS
			 if (strcmp(name, "tau") == 0){
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

#ifdef SOLOV
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
  
#ifdef RESISTIVITY

  eta_Ohm = 0.0;
  Q_AD    = 0.0;
  Q_Hall  = par_getd("problem","Q_H");
  //d_ind   = 1.0;


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

          pG->U[k][j][i].B1c = 0.5*((x1-0.5*pG->dx1)*pG->B1i[k][j][i] + (x1+0.5*pG->dx1)*pG->B1i[k][j][i+1])/x1;


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

  MPI_Reduce(&Pgas, &TotPgas, 1, MPI_Real, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Pb, &TotPb, 1, MPI_Real, MPI_SUM, 0, MPI_COMM_WORLD);

  if (my_id == 0) {

    printf("Total gas pressure = %f\n", TotPgas);
    printf("Total magnetic pressure = %f\n", TotPb);
    apause();
  }

  MPI_Bcast(&TotPgas, 1, MPI_Real, 0, MPI_COMM_WORLD);
  MPI_Bcast(&TotPb, 1, MPI_Real, 0, MPI_COMM_WORLD);

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
                + SQR(pG->U[k][j][i].M3))  /  (pG->U[k][j][i].d);

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

  bvals_mhd_fun(pDomain, left_x1,  inX1BoundCond );
  bvals_mhd_fun(pDomain, left_x3,  diode_outflow_ix3 );
  bvals_mhd_fun(pDomain, right_x3,  diode_outflow_ox3);

//  plot(pG, "tau");

#ifdef XRAYS

  optDepthFunctions(pG);

#endif

//plot(pDomain, "tau");

//plot(pDomain, "xi");


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
	Real IntE, KinE, MagE=0.0, x1, x2, x3, DivB,Pgas,rad, dns,as_lim;

	Real di, v1, v2, v3, qsq,p, asq,b1,b2,b3, bsq,rho_max;


	static Real TotMass=0.0;

  //A.D.
  GridS *pG=pM->Domain[0][0].Grid;


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
optDepthFunctions(pM);
#endif

//plot(pM, "d");
//plot(pM, "tau");
//plot(pM, "xi");
//plot(pM, "E");


  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {

    	  	cc_pos(pG,i,j,k,&x1,&x2,&x3);

//        rad = sqrt(SQR(x1) + SQR(x3));
//        if (rad< R_ib){
//        		pG->U[k][j][i].d = rho0;
//        		 pG->U[k][j][i].E = e0;
//        		pG->U[k][j][i].M1= 0.;
//        		pG->U[k][j][i].M2= 0.;
//        		pG->U[k][j][i].M3= 0.;
//        }


        KinE = 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2)
                + SQR(pG->U[k][j][i].M3))/(pG->U[k][j][i].d);

#ifdef MHD

        MagE = 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c)
                + SQR(pG->U[k][j][i].B3c));
#endif

        IntE = pG->U[k][j][i].E - KinE - MagE;
       
        if (isnan(pG->U[k][j][i].d) || isnan(pG->U[k][j][i].E)) {
        				printf("At pos isnan: (%d,%d,%d) (%f,%f,%f), Den1 = %f, E1 = %f\n", i, j, k, x1,x2,x3,pG->U[k][j][i].d, pG->U[k][j][i].E);
        				printf("KinE1, MagE1, IntE (%f,%f,%f), \n", KinE, MagE, IntE);
        				apause();

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

	  if (my_id == 0) {

		  printf(" my_id =%d\n",  my_id);

//		  apause();
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

//void traceCell(MeshS *pM, int indxpos[3], Real xpos[3], Real cnorm[3], Real cylDir[3]) {
//	// traces a cingle cell; cnorm is in Cart. coordinates;
////      use Grids, only: x1b,x2b,x3b,is,ie,js,je, ks,ke
////      use paramStruct
////      use physconst_mod, only: tiny, huge
////      implicit none
//
////      Real cnorm(3), cylDir(3)
////      integer,intent(inout) ::indxpos(3)
////      real(knd), intent(inout) ::  xpos(3)
//	  GridS *pG=pM->Domain[0][0].Grid;
//      int inext, jnext, knext, icur, jcur, kcur, ind[1];
//
//      Real a,b,c, rcyl,d,L1,L2,Lfin,L[3],phnext,sinph,cosph,
//                          xnext,ynext,znext,tanph,cotph,nzbar,den,xcur, ycur, zcur,
//                          swtch[3];
//  	  int is = pG->is,
//  	  ie = pG->ie,
//  	  js = pG->js,
//	  je = pG->je;
//
//      icur= indxpos[0]; /* --- icur etc. are on cyl. grid ----*/
//      jcur= indxpos[1];
//      kcur=indxpos[2];
//
//      xcur = xpos[0];   /*  are on Cart. grid */
//      ycur=xpos[1];
//      zcur=xpos[2];
//      swtch[0]=0.;
//      swtch[1]=0.;
//      swtch[2]=0.;
//
//       /*  1)  intersection with cylinders  */
//
////       inext = icur +int(sign(1._knd, cylDir(1)))
//
//       inext = icur + copysign(1, cylDir[0] );
//
//       if (inext<is) inext = is;
//       if (inext>ie) inext = ie;
//
//
//       rcyl = fmax(inRadOfMesh, x1b(inext));
//
//       a = pow(cnorm[0],2) + pow(cnorm[1],2);
//
//       b = 2.*(xcur*cnorm[0]+ycur*cnorm[1]);
//
//       c = pow(xcur,2)  + pow(ycur,2) - pow(rcyl,2);
//
//       d = pow(b,2)-4.*a*c;
//
//       if (d<0){
////          print*, a,b,c, inext,icur,d,rcyl
////          pause 'Discr<0. in traceACell'
//    	   	   apause();
//       }
//
////        a = sign(max(tiny, abs(a)), a);
//
//       a = copysign(fmax(tiny, fabs(a)), a);
//
//        L1 = abs(-b+sqrt(d))/2./a;
//        L2 = abs(-b-sqrt(d))/2./a;
//
//        L[0] = fmin(L1, L2);
//        if (L[0]==0) L[0]=HUGE_NUMBER;
//
//      /*  2)  intersection with planes phi_j=const */
//
//          jnext= jcur;
//          L[1]=HUGE_NUMBER;
//
//          if (cylDir[1] /= 0.) {
//              jnext= jcur + ( copysign(1, cylDir[1]) );
//          }
//
//          if (jnext<js) jnext = je;
//          if (jnext>je) jnext = js;
//
//
//          phnext = x2b[jnext];
//          cosph = cos(phnext);
//          sinph = sin(phnext);
//
//          if ( fabs(cosph)>0.1) {
//              tanph=sinph/cosph;
//              den = cnorm[1] - tanph*cnorm[0];
//              den = copysign(fmax(tiny, fabs(den)), den);
//              L[1] = fabs((xcur*tanph-ycur)/den);
//
//          else
//              cotph= cosph/sinph
//              den = (cotph*cnorm[1] - cnorm[0])
//              den = sign(max(tiny, abs(den)), den)
//              L[1] = abs((xcur - cotph*ycur)/den)
//          }
//          if (L[1]==0) L[1]=HUGE_NUMBER;
//
//
//
//     /*  3)  intersection with planes z_k=const */
//
//     knext=kcur + copysign(1, cylDir[2]);
//
//      nzbar = fmax(fabs(cnorm[2]),tiny);
//      nzbar = copysign(nzbar, cnorm[2]);
//      L[2] = fabs((x3b[knext] - zcur)/nzbar);
//      if (L[2]==0)  L[2]=HUGE_NUMBER;
//
//        ind = minloc(L);
//        Lfin =  L(ind(1))
//        swtch(ind(1)) =1.
//
//        indxpos =[icur, jcur, kcur] + swtch*([inext, jnext, knext] - [icur, jcur, kcur])
//
//        xpos = [xcur, ycur, zcur] + cnorm*Lfin
//
//        /* special case of corner values */
//        if( (L(1) /= L(2)) .and. (L(2) /= L(3)) .and. (L(1) /= L(3)) ) then
//          return
//        else if (  (L(1)==L(2))  .and. (L(2) < L(3))) then
//          indxpos(1) = inext
//          indxpos(2) = jnext
//
//        else if ( (L(1)==L(3))  .and. (L(3) < L(2))) then
//          indxpos(1) = inext
//          indxpos(3) = knext
//
//        else if (  (L(1)==L(2)) .and. (L(2)==L(3)) .and. L(3)>0) then
//            indxpos = [inext, jnext, knext]
//        end if
//
//} //end  traceCell

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
