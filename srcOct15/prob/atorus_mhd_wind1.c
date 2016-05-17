#include "copyright.h"
/*============================================================================*/
/*! \file hkdisk.c
 *  \brief Problem generator for Hawley Krolik disk (Specifically, GT4)
 */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
//#include "cyl.h"


// Input parameters
static Real q, r0, rhomax, r_in, rho0, e0, dcut, beta, seed;

// Derived quantities
//static Real f, C, Kbar, n;


static Real Ctor, Kbar, nAd, q1,rgToR0,Tx,Rg,MBH_in_gram,Tgmin,
Dsc,Rsc,Usc,Psc, Time0,Tvir0,Lamd0,Evir0, nc0,Ledd_e,
Lx, Luv;


static Real F2Fedd = 0.01, fx =0.5, fuv =0.5,
	nc0 = 1E8,
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
	M_MW = 1.; //mean mol. weight
//	PI = 3.14159265359;



//extern x1GravAcc;



/*! \fn static Real grav_pot(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential */




void pause(){
	printf("pause.. press any key ss... ");
	getchar();
}

Real xRayHeatCool(const Real dens, const Real Press, const Real dt)//, const Real xi)
{

	Real Gx,Tg;

	Real res=0.,
    c1 = 8.9e-36,
	c2 = 1.5e-21,
	al0 = 3.3e-27,
	al1 = 1.7e-18,
	al2 = 1.3e5,
	bta = -1.;



//       egas = Evir0 * eij
//
//
//       Tg = abs( (gam -1.)*m_mw *egas/ (nd* nc0*m_u) /Rgas )

	 Tg = abs( Press * Psc * M_MW) / (dens*Dsc*RGAS );

	//
//       taucool=0.
//       if((Lbr_exp_Lalpha_tau .eqv.  .true.)) then
//               kpx = kappa_e
//               ro = Dsc* sum(gas1%dns(i-1:i+1, j-1:j+1) ) /27.
//
//               dtz = kpx * ro * Rsc* ( grid1%x1a(i+1) - grid1%x1a(i-1))/2.
//               dtr = kpx * ro * Rsc* ( grid1%x2a(j+1) - grid1%x2a(j-1))/2.
//               taucool = min(dtz, dtr)
//       endif
//
//       dlt = deltaXCool*exp( - taucool)
//       Gc = 0.
//       Gx =0.
//       Lbr =0.
//


//			 here
//		Gx= c2*pow(xi,0.25) /pow(Tg,0.5)*(1. - Tg/Tx)


//       Cal=0.101
//       N22 = log10(gas1% col_dens(i,j))/22.
//
//       Fx = Lx/4./pi/rcm**2 !erg cm^-2 s^-1
//
//
//       tau_x = gas1% tau_e(i,j)
//
//
//       if (optical_depth_attenuate1) tau_x=0.
//
//
//       Fx= max( Fx * exp( - tau_x ),  tiny)
//       S1 = 1.
//       Gx_mol = 4.2e-23 * Cal * Fx*S1 / N22**0.9 !erg s^-1
//       Gx_mol = Gx_mol/(nd*nc0)
//
//
//       if ( log10(gas1% col_dens(i,j)) >  lg10_col_dns_cr ) then
//!                     Gx = min(Gx_mol,Gx)
//
//!                     Gx = Gx_mol
//
//!                    print*, "Gc=", Gc, "Gx=", Gx,"Gx_mol=", Gx_mol, "Lbr=", Lbr, "xi= ", xi,"n=",nd*nc0
//       endif
//
//       if (XRHC(3)==1) then
//           Gc = c1* xi*(Tx - 4.*Tg)
//       endif
//
//       lbr1= al0*sqrt(Tg)
//       lbr2= ( al1*exp(-al2/Tg)* (xi**bta) /sqrt(Tg)+1.e-24)*dlt
//       Lbr=lbr1 + lbr2 !(erg cm^3/s)
//
//
//       if ( Lbr_T4_Lbr ) limitfactor=exp(-10.**4 /Tg)
//       Lbr = Lbr* limitfactor
//
//!                Lbr  = Lbr* exp(-taucool)
//
//!ro = Dsc*gas1%dns(i,j)
//!Lbr = ro*Lbr
//
//          dGcdT = - 4.*c1*xi
//          dGxdT = - 0.5*c2/Tx * (xi**0.25) *(Tx/Tg+1.)* (Tg**(-0.5))
//          dLbrdT = 0.5* al0 *Tg**(-0.5) + dlt*al1* (xi**bta) *exp(-al2/Tg)* &
//                      ( al2 -  0.5* Tg  )  * (Tg**-2.5)
//
//
//       if ( (XRHC(2)==1) .and. (XRHC(1)==1 )) then
//           Hx = Gc +Gx - Lbr
//           dHxde = dGcdT + dGxdT - dLbrdT
//      else if ( (XRHC(1) ==1 ) .and. (XRHC(2) ==0 )  ) then
//           Hx = - Lbr
//           dHxde = - dLbrdT
//       else
//           Hx=0.
//           dHxde = 0.
//
//       endif
//
//       Hx = Hx/Lamd0 * (nd*nc0)**2
//
//!print*,  Hx, Lamd0 , nd, nc0; pause
//
//!                print*, this%A2, Lamd0; pause
//
//
//       fac1 = (gam-1.)*m_mw*Evir0/Dsc/Rgas
//       dHxde = dHxde *  fac1  /Lamd0 * nd * ( nc0 )**2
	Tg=0.;
    return(Tg);
}

void plot(const Real *A22){

	FILE *f = fopen("file.txt", "w");

	if (f == NULL){

		printf("Error opening file!\n");
	    exit(1);
	}

//    status= getcwd( wrkdir )
//
//	print*, wrkdir(1:index(wrkdir,' ')-1)//'/plot_from_fort.py'
//
//    fname1= wrkdir(1:index(wrkdir,' ')-1)//'/plot_from_fort.py'
//
//    status =  system(  fname1   )
//	write(11, *) N1, N2
//	do  j = 1, N2
//		do  i = 1, N1
// 			write(11, *)  X(i,j)
//		enddo
//		write(11, *)  ' '
//	enddo


//	 call system( ffile )
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

Real density(Real x1, Real x2, Real x3) {
  Real rad, tmp, d;

//  x1=1.; x2=0.; x3=0.;

  rad = sqrt( SQR(x1) + SQR(x3));

  tmp = (Ctor + (1.0/(rad-rgToR0)) - pow(x1,-q1)/q1)/(nAd-1)/Kbar;

  d = pow(tmp,nAd)*( x1>=r_in );



//  printf("in density()  %f %f %f %f %f %f\n", r_in, nAd, x1, rad, tmp, d); pause();

//  if ( x1>=r_in  && tmp>0. ){
//	  d = pow(tmp,nAd);
//	  printf("%f %f %f %f %f %f\n", r_in, nAd, x1, rad, tmp, d); pause();
//  }

  d = MAX(d, rho0);


//  temp = C + (1.0/(rad-2.0)) - f*pow(x1,-2.0*(q-1.0));
//  temp = MAX(temp,0.0);
//  d = pow(temp/Kbar,n)*(x1>=r_in);
//  d = MAX(d,rho0);
//   d = x1 < r_in ? rho0 : d;
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
//#define VP(R) ((sqrt(r0)/(r0-2))*pow(R/r0,-q+1))

#define VP(R) pow(R,-q+1)


void problem(DomainS *pDomain)
{


  GridS *pG=(pDomain->Grid);


  int i,j,k;
  int is,ie,js,je,ks,ke,nx1,nx2,nx3;
  int il,iu,jl,ju,kl,ku;
  Real rad, IntE, KinE, MagE;
  Real x1,x2,x3, x1i, x2i, x3i,rhoa;
  //  , T, Ta, rhoa,q1,r02;

  Real Pgas, Pb, TotPgas, TotPb;

  Real scl, ***Ap, ***Bri, ***Bzi;
  Real divB=0.0, maxdivB=0.0;




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
  q      = par_getd("problem","q");
  r0     = par_getd("problem","r0");
  rhomax = par_getd("problem","rhomax");
  r_in   = par_getd("problem","r_in");
  rho0   = par_getd("problem","rho0");
  e0     = par_getd("problem","e0");
  seed   = par_getd("problem","seed");

#ifdef MHD
  dcut = par_getd("problem","dcut");
  beta = par_getd("problem","beta");
#endif


// Set up physical parameters
//  n = 1.0/Gamma_1;
//  q1 = 2.0*(q-1.0);
//  r02 = (r0-2.0);
//  f = pow(r0, 2.0*q-1.0) / (q1*pow(r02,2.0));
//  C = f*pow(r_in,-1.0*q1) - 1.0/(r_in-2.0);
//  Kbar = (C + (1.0/r02) - f*pow(r0,-1.0*q1))/pow(rhomax, 1.0/n);

  //  A.D. calculate local model parameters
    Real rTorBnd, tmp, mue;

    MSOLYR = 1.989e33/YR;
    SGMB = ARAD*CL/4;

    MBH_in_gram = M2Msun*MSUN;
    Rg = 2*GRV*MBH_in_gram/pow(CL,2);
    r0*=Rg;
    rgToR0 = Rg/r0;

    Tx=1.0e8; //Compton temperature
    	Tgmin = 2.3;
    	Dsc = nc0*M_MW*	M_U;
//    	printf( " %e %e %e %e\n", Dsc, nc0, M_MW, M_U);

   	Rsc = r0;
    	Usc = sqrt(GRV*MBH_in_gram/Rsc);
    	Time0=sqrt( pow(Rsc,3) /GRV/MBH_in_gram);
    	Evir0=Dsc*GRV*MBH_in_gram/Rsc;  //(erg/cm^3)
    	Psc = Evir0;
    	mue = M_MW;

	Tvir0=GRV*MBH_in_gram/RGAS/Rsc*mue;
    	Lamd0=Evir0* Time0;  // (erg/cm^3/s)

    	Ledd_e = 4.*PI *CL*GRV*MBH_in_gram/KPE;

    	Lx = fx * F2Fedd*Ledd_e;

    	Luv= fuv *F2Fedd*Ledd_e;


    nAd = 1.0/Gamma_1;
    q1 = 2.0*(q-1.0);

//    rgToR0 = 0.;

    Ctor = pow(r_in,-1.0*q1)/q1 - 1.0/(r_in-rgToR0);

    Kbar = (Ctor + 1.0/(1.-rgToR0) - pow(1.,-q1)/q1) /(nAd-1);


    printf(" %e %e %e %e %e %e\n", Rsc, Dsc, Usc, Time0,Tvir0,Lamd0); pause();


//    tmp = (Ctor + (1.0/(1.-rgToR0)) - pow(1.,-q1)/q1)/(nAd-1)/Kbar;
//    tmp = pow(tmp,nAd);
//    printf("1 %f\n",  tmp);
//    printf("Ctor =  %f   Kbar = %f\n", Ctor, Kbar); pause();
//  ************************************************


  // Initialize data structures, set density, energy, velocity and Ap
  // Loop over all cells (including ghosts)
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {

//    	    printf("%d\n", ju);

    	    cc_pos(pG,i,j,k,&x1,&x2,&x3);

        rad = sqrt(SQR(x1) + SQR(x3));

        x1i = x1 - 0.5*pG->dx1;
        x2i = x2 - 0.5*pG->dx2;
        x3i = x3 - 0.5*pG->dx3;

        pG->U[k][j][i].d = rho0;
        pG->U[k][j][i].M1 = 0.0;
        pG->U[k][j][i].M2 = VP(x1)*rho0;
        pG->U[k][j][i].M3 = 0.0;
        pG->U[k][j][i].E   = e0;
        Ap[k][j][i] = 0.0;

#ifdef MHD
        pG->B1i[k][j][i] = 0.0;
        pG->B2i[k][j][i] = 0.0;
        pG->B3i[k][j][i] = 0.0;
        pG->U[k][j][i].B1c = 0.0;
        pG->U[k][j][i].B2c = 0.0;
        pG->U[k][j][i].B3c = 0.0;

#endif

        // Set up torus

        rTorBnd =  1/( pow(x1, -1.0*q1 ) -  Ctor ) + rgToR0;



//        if ( (x1 >= r_in) && (rad <= rTorBnd) ) { //Checks to see if cell is in torus

        if ( x1 >= r_in ) { //Checks to see if cell is in torus

        	rhoa = density(x1, x2, x3);


        	IntE = pow(rhoa,GAM53)/Gamma_1;

		// Add pressure fluctuations to seed instability
		IntE = IntE*(1.0 - seed*sin(x2));

		if ((IntE >= e0) && (rhoa >= rho0)) {
						//If the values are above cutoff, set up the cell
            pG->U[k][j][i].d = rhoa;
            pG->U[k][j][i].M2 = VP(x1)*pG->U[k][j][i].d;
            pG->U[k][j][i].E = IntE;

            //Note, at this point, E holds only the internal energy.  This must
            //be fixed later
          }



        }



      }
    }
  }


  // Calculate density at corner and set up Ap if appropriate
  for (k=kl; k<=ku+1; k++) {
    for (j=jl; j<=ju+1; j++) {
      for (i=il; i<=iu+1; i++) {

    	    cc_pos(pG,i,j,k,&x1,&x2,&x3);

    	    rad = sqrt(SQR(x1) + SQR(x3));

    	    x1i = x1 - 0.5*pG->dx1;
        x2i = x2 - 0.5*pG->dx2;
        x3i = x3 - 0.5*pG->dx3;

        rhoa = density(x1i,x2i,x3i);

        if (rhoa >= dcut) {
          // Ap = (density-cutoff)^2
          Ap[k][j][i] = SQR(rhoa-dcut);
        }
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

      }
    }
  }

  // Calculate cell centered fields
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {

        cc_pos(pG,i,j,k,&x1,&x2,&x3);

        if (i==iu)
          pG->U[k][j][i].B1c = pG->B1i[k][j][i];
        else

          pG->U[k][j][i].B1c = 0.5*((x1-0.5*pG->dx1)*pG->B1i[k][j][i] + (x1+0.5*pG->dx1)*pG->B1i[k][j][i+1])/x1;

        if (k==ku)
          pG->U[k][j][i].B3c = pG->B3i[k][j][i];
        else

          pG->U[k][j][i].B3c = 0.5*(pG->B3i[k+1][j][i] + pG->B3i[k][j][i]);
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
  MPI_Reduce(&Pgas, &TotPgas, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Pb, &TotPb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (pG->my_id == 0) {
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

  CoolingFunc = xRayHeatCool;

//  x1GravAcc = grav_acc;
//  set_bvals_fun(left_x1,  disk_ir_bc);
//  set_bvals_fun(right_x1,  disk_or_bc);
//  set_bvals_fun(left_x3,  diode_outflow_ix3);
//  set_bvals_fun(right_x3,  diode_outflow_ox3);



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
  q = par_getd("problem","q");
  r0 = par_getd("problem","r0");
  rhomax = par_getd("problem","rhomax");
  r_in = par_getd("problem","r_in");
  rho0 = par_getd("problem","rho0");
  e0 = par_getd("problem","e0");
	seed = par_getd("problem","seed");

#ifdef MHD
  dcut = par_getd("problem","dcut");
  beta = par_getd("problem","beta");
#endif

  /* Enroll the gravitational function and radial BC */
  StaticGravPot = grav_pot;

//  x1GravAcc = grav_acc;


//  set_bvals_fun(left_x1,disk_ir_bc);
//  set_bvals_fun(right_x1,disk_or_bc);
//  set_bvals_fun(left_x3,diode_outflow_ix3);
//  set_bvals_fun(right_x3,diode_outflow_ox3);

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
	int prote, protd;
  Real IntE, KinE, MagE=0.0, x1, x2, x3, DivB;
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


  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);

        if (isnan(pG->U[k][j][i].d) || isnan(pG->U[k][j][i].E)) {
          printf("At pos (%f,%f,%f), Den = %f, E = %f\n", x1,x2,x3,pG->U[k][j][i].d, pG->U[k][j][i].E);
        }

        KinE = 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2)
                + SQR(pG->U[k][j][i].M3))/(pG->U[k][j][i].d);

#ifdef MHD

        MagE = 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c)
                + SQR(pG->U[k][j][i].B3c));
#endif
        IntE = pG->U[k][j][i].E - KinE - MagE;

        if (pG->U[k][j][i].d < rho0) {

          pG->U[k][j][i].d = rho0;

          KinE = 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2)
                  + SQR(pG->U[k][j][i].M3))/(pG->U[k][j][i].d);

          IntE = e0;  // Set this as well to keep c_s reasonable

          pG->U[k][j][i].E = IntE + KinE + MagE;

          protd++;
		  prote++;

        }

        else if (IntE < e0) {

          pG->U[k][j][i].E = e0 + KinE + MagE;
					prote++;
        }
        IntE = pG->U[k][j][i].E - KinE - MagE;

      }
    }
  }

#ifdef MPI_PARALLEL
	if (pG->my_id == 0) {
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


}

void Userwork_after_loop(MeshS *pM)
{
}
