/*---------------------------------------------------------------------------------------------*/
/*  Abed Zadehgol               d2z8-CP (plain version)                 DATE : 2016/09/16      */
/*  email:   a.zadehgol@me.iut.ac.ir   &    abed.zadehgol@gmail.com                            */
/*                                                                                             */
/*  This is a sample C code, based on the generalized BGK-Boltzmann Equation proposed in       */
/*  Zadehgol [Phys. Rev. E (94), 023326 (2016)] (with complex-valued probability density)      */      
/*  to simulate incompressible fluid flows. The present simulation uses a single relaxation    */
/*  time method.  Note that the apparent jumps, are given by d = +/- (0.5*(sqrt(2.0)-1.0)),    */
/*  while the radius of the collision circle is given by v=(0.5*(sqrt(2.0)+1.0)).              */
/*                                                                                             */
/*  Additional notes:                                                                          */
/*                                                                                             */
/*    (1) The present code is a complete and self-sufficient version (additional codes are     */
/*        not needed to complile and run this code). This program has successfuly been         */
/*        complied in Microsoft Visual Studio 2008 environment.                                */
/*                                                                                             */
/*    (2) The program generates an output data file, "FIELD.DAT", in "Tecplot" format.         */
/*         The output file contains the information of the following fields:                   */
/*                [a] velocity  (ux, uy)                                                       */
/*				  [b] density   (ro)                                                           */
/*				  [c] pressure  (P)                                                            */
/*                [d] si        (for the steamlines)                                           */ 
/*                                                                                             */
/*    (3)   To see the contour plots, load the file FIELD.DAT into Tecplot, go to the          */
/*           contour mode, switch from "flood" to "lines", and choose the desired field.       */
/*                                                                                             */
/*    (4)   Adjust the adjustable parameters below.                                            */
/*                                                                                             */
/*---------------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <stdlib.h>
#include <complex>

// ---------------------------------  the adjustable parameters --------------------------------

#define MAX_T    500.0       // Duration of simulation (in dimensionless time)

#define Re       400.0       // Reynolds nunber 
#define MACH     0.10        // Mach number (try the 0.025 - 0.10 range)

#define y_dim    128         // hight of the cavity (in lattice units)
#define x_dim    128         // width of the cavity (in lattice units)

#define d_iter_01          50
#define d_iter_02          500
#define d_iter_03          20

// ---------------------------------   other parameters  ---------------------------------------
#define V       (0.5*(sqrt(2.0)+1.0))   // magitude of the microscopic velocities
#define delta   (0.5*(sqrt(2.0)-1.0))   // magnitude of the apparent jumps (in streaming stage)

#define cs                 (V*sqrt(0.5))  // speed of sound 

#define DENSITY            1.0  
#define U                  (MACH*cs)    // velocity of the upper lid

#define pi                 (4.0*atan(1.0))

#define Ndirs              9         

// ------ data structures  ---------
using namespace std;
complex <double> ttau[Ndirs];
typedef struct{complex <double> DIR[Ndirs+1];} EQU;
typedef struct{double ro1; double ro2; double ux; double uy;} STATE;

//-------------------------  declarations ------------------------
double H, visc, COEF;
double si[x_dim][y_dim];
int    ex[Ndirs], ey[Ndirs]; 
int    iter,n_iter_01, n_iter_02,n_iter_03,x_pos,y_pos;
EQU    *F_in, *F_out;
int    *bdr_state;

STATE get_state(int i,int j);

/*---------------------------------------------------------------------------------------------*/
/*------------------>                 allocate_memory                       <------------------*/
/*---------------------------------------------------------------------------------------------*/
void allocate_memory(void){

	int size_a = (x_dim*y_dim)*sizeof(EQU);	
	int size_b = (x_dim*y_dim)*sizeof(int);	
	int size_c = (x_dim*y_dim)*sizeof(double);	

	// allocate memory for all transferable and resident cpu and gpu veriables 
	F_in  = (EQU*) malloc(size_a);
	F_out = (EQU*) malloc(size_a);
	bdr_state = (int*) malloc(size_b);

}

/*---------------------------------------------------------------------------------------------*/
/*------------------>     the equilibrium distribution of the d2z8-CP mode  <------------------*/
/*---------------------------------------------------------------------------------------------*/
EQU get_equ(STATE state){

  using namespace std;

  EQU res;
  int i;
  double v1, v2, ang, del_ang, ux, uy, ro, w;
  complex <double> u, v;

  w = 1.0/(Ndirs - 1.0);

  ro = state.ro1;
  ux = state.ux;
  uy = state.uy;

  u  = complex <double> (ux,  uy);
  del_ang = 2.0*pi/(Ndirs-1.0);

  for(i=1; i<= Ndirs-1; ++i){	 	       
	 ang = (i-1.0)*del_ang;	 
	 v1 = V*cos(ang);
	 v2 = V*sin(ang);
	 v  = complex <double> (v1, v2);
	 res.DIR[i] = w*ro*v/(v-u);
  }

  res.DIR[0] = complex <double> (0.0, 0.0);

  return(res);

}

/*---------------------------------------------------------------------------------------------*/
/*--------------->      implementing boundary conditions for the upper lid        <------------*/
/*---------------------------------------------------------------------------------------------*/
void reset_F_in(void){

  int id, i, j;
  STATE state;

  EQU equ_dist;

  // --------------------  upper lid -----------------
  state.ro1 = DENSITY;
  state.ro2 = 0.0;
  state.ux = U;   
  state.uy = 0.0;
  equ_dist = get_equ(state);
  j = y_dim-1;
  for(i=0;i<=x_dim-1;++i){
	id = j*x_dim+i;

	// implementing a simple boundary condition at the upper lid
	//F_in[id] = equ_dist;
	
	// implementing the Zuo-He boundary condition (bounce back of the non-equilib. dists.) 
	F_in[id].DIR[6] = F_in[id].DIR[2] +(equ_dist.DIR[6] - equ_dist.DIR[2]);
	F_in[id].DIR[7] = F_in[id].DIR[3] +(equ_dist.DIR[7] - equ_dist.DIR[3]);
	F_in[id].DIR[8] = F_in[id].DIR[4] +(equ_dist.DIR[8] - equ_dist.DIR[4]);

  }  

}
/*---------------------------------------------------------------------------------------------*/
/*----------------->              initialize the main parameters           <-------------------*/
/*---------------------------------------------------------------------------------------------*/
void init_params(void){

  COEF = 20.0;  
  H = (y_dim + 0.0);
  visc = H*(0.5*U)/Re;

  // --- code parameters ---
  n_iter_01 = d_iter_01;
  n_iter_02 = d_iter_02;
  n_iter_03 = d_iter_03;

  x_pos = x_dim/2;
  y_pos = y_dim/2;

  // ------ intervals for the neighbour nodes -------
  ex[1] =  1;  ey[1] =  0;
  ex[2] =  1;  ey[2] =  1;
  ex[3] =  0;  ey[3] =  1;
  ex[4] = -1;  ey[4] =  1;
  ex[5] = -1;  ey[5] =  0;
  ex[6] = -1;  ey[6] = -1;
  ex[7] =  0;  ey[7] = -1;
  ex[8] =  1;  ey[8] = -1;
  ex[0] =  0;  ey[0] =  0;
}


/*---------------------------------------------------------------------------------------------*/
/*-------------------->    specify the boundry types of the nodes   <---------------------------*/
/*---------------------------------------------------------------------------------------------*/
void init_bdrs(void){
  int i,j, id;
  
  // -------------- initializing the fluid nodes ----------------
  for(i=0;i<=(x_dim-1);++i){
    for(j=0;j<=(y_dim-1);++j){
	  id = j*x_dim + i;
      bdr_state[id] = 0;
    }
  }


  // -------------- left wall ----------------
  i = 0;
  for(j=0;j<=(y_dim-1);++j){
    id = j*x_dim + i;
    bdr_state[id] = 1;
  }


  // -------------- right wall ----------------
  i = x_dim-1;
  for(j=0;j<=(y_dim-1);++j){
    id = j*x_dim + i;
    bdr_state[id] = 1;
  }


  // -------------- buttom wall ----------------
  j = 0;
  for(i=0;i<=(x_dim-1);++i){
    id = j*x_dim + i;
    bdr_state[id] = 1;
  }

}

/*---------------------------------------------------------------------------------------------*/
/*----->                           initialize all nodes                                 <-----*/
/*---------------------------------------------------------------------------------------------*/
void init_nodes(void){
  int id, i,j,k;
  STATE state;

  EQU equ_dist;

  state.ux = 0.0;
  state.uy  = 0.0;
  state.ro1 = DENSITY;
  state.ro2 = 0.0;

  equ_dist = get_equ(state);

  /* initializing the general nodes */
  for(j=0;j<=y_dim-1;++j){
    for(i=0;i<=x_dim-1;++i){
	  id = j*x_dim+i;
      for(k=1;k<=Ndirs-1;++k){
		  F_in[id].DIR[k] = equ_dist.DIR[k];
		  F_out[id].DIR[k] = equ_dist.DIR[k];
	  }
    }
  }
}

/*---------------------------------------------------------------------------------------------*/
/*--------------->                  initialize the variables                      <------------*/
/*---------------------------------------------------------------------------------------------*/
void initialize(void){
  init_params();
  init_bdrs();
  init_nodes();
}

/*---------------------------------------------------------------------------------------------*/
/*                              calculate the macroscopic state properties                     */
/*---------------------------------------------------------------------------------------------*/
STATE get_state(int i,int j){

  double v1, v2, ang, del_ang;
  STATE res;
  int id, k;

  using namespace std;
  complex <double> f, ro, u, v;

  ro = complex <double> (0.0, 0.0);
  u  = complex <double> (0.0, 0.0);

  id = j*x_dim+i;

  for(k=1;k<=Ndirs-1;++k){
     del_ang = 2.0*pi/(Ndirs-1.0);
	 ang = (k-1.0)*del_ang;
	 
	 v1 = V*cos(ang);
	 v2 = V*sin(ang);
	 v  = complex <double> (v1, v2);

	 f = F_in[id].DIR[k];

	 ro +=  f;
	 u  +=  f*v;
  }

  res.ro1 = ro.real();

  if(abs(ro) != 0.0){
    res.ux = (u/res.ro1).real();
    res.uy = (u/res.ro1).imag();
  }else{
    res.ux = 0.0;
    res.uy = 0.0;
  }

  return(res);
}

/*---------------------------------------------------------------------------------------------*/
/*                              calculate the pressure                                         */
/*---------------------------------------------------------------------------------------------*/
double calc_P(int i,int j){
  STATE state;
  double res, ro, ux, uy, u__2, v__2;

  state = get_state(i,j);
  ro = state.ro1;
  ux = state.ux;
  uy = state.uy;

  u__2 = ux*ux + uy*uy;
  v__2 = V*V; 

  res = 0.5*ro*(v__2 - u__2);

  return(res);

}
/*---------------------------------------------------------------------------------------------*/
/* --------------->               non-boundary collisions                      <---------------*/
/*---------------------------------------------------------------------------------------------*/
void ord_col(int i,int j){

  using namespace std;

  int id, k;
  double ro,u,v, tau, u__2, cs2;
  complex <double> ttau;

  STATE state;
  EQU equ_dist;

  id = j*x_dim+i;

  /*-- calculating the local properties ---*/

  state = get_state(i,j);

  ro = state.ro1;
  u  = state.ux;
  v  = state.uy;

  u__2 = u*u + v*v;

  cs2 = 0.5*(V*V - u__2);
  tau = 4.0*visc/(V*V) + 0.5;
  ttau = complex <double> (tau, 0.0);

  /*-- calculating the local post-colision out-states ---*/
  equ_dist = get_equ(state);

  for(k=1; k<=Ndirs-1; ++k){
    F_out[id].DIR[k]  = F_in[id].DIR[k]-(F_in[id].DIR[k] - equ_dist.DIR[k])/ttau;
  }

}

/*---------------------------------------------------------------------------------------------*/
/*--------------->                    boundary collisions                  <-------------------*/
/*---------------------------------------------------------------------------------------------*/
void bdr_col(int i,int j){
   int id = j*x_dim+i;
   F_out[id].DIR[1] = F_in[id].DIR[5];
   F_out[id].DIR[2] = F_in[id].DIR[6];
   F_out[id].DIR[3] = F_in[id].DIR[7];
   F_out[id].DIR[4] = F_in[id].DIR[8];
   F_out[id].DIR[5] = F_in[id].DIR[1];
   F_out[id].DIR[6] = F_in[id].DIR[2];
   F_out[id].DIR[7] = F_in[id].DIR[3];
   F_out[id].DIR[8] = F_in[id].DIR[4];
   F_out[id].DIR[0] = F_in[id].DIR[0];
}

/*---------------------------------------------------------------------------------------------*/
/*                Post-Collision streaming effect on the neighbouring nodes.                   */
/*---------------------------------------------------------------------------------------------*/
void stream(int i,int j){
  
  int id, nid, k,ni,nj;
  id = j*x_dim + i;
  
  for(k=1;k<=Ndirs-1;++k){
    ni = i + ex[k];
    nj = j + ey[k];
	nid = nj*x_dim + ni;
	if(ni>=0 && ni<=x_dim-1 && nj>=0 && nj<=y_dim-1){
		F_in[nid].DIR[k] = F_out[id].DIR[k];
	}
  }

}

/*---------------------------------------------------------------------------------------------*/
/*--------------                  Calculate the Si field                  ---------------------*/
/*---------------------------------------------------------------------------------------------*/
void rec_si(void){
  int id, i,j;

  double dx = 1.0;

  STATE state;

  /*--------- calculating the si field ----------*/
  for(i=0;i<=x_dim-1;++i){ si[i][0]=0.0; si[i][y_dim-1]=0.0;}
  for(j=0;j<=y_dim-1;++j){ si[0][j]=0.0; si[x_dim-1][j]=0.0;}

  for(i=1;i<=x_dim-2;++i){
    for(j=y_dim-2;j>=1;--j){
		state = get_state(i,j);
	  si[i][j] = si[i-1][j] - state.uy*dx;	
	}
  }
  
  for(i=0;i<=x_dim-1;++i){
    for(j=y_dim-1;j>=0;--j){
		id = j*x_dim +i;
	  if(bdr_state[id] == 1){si[i][j] = 0.0;}
	}
  }
  
}

/*---------------------------------------------------------------------------------------------*/
/*---                     Preparing and recording the field data                           ----*/
/*---------------------------------------------------------------------------------------------*/
void rec_field(void){
  FILE *out;
  int i,j, k;

  STATE state;
  double T, xx, yy, ro, u, v, P;
  double dx, dy,  sai, u__2;

  int GRID_SIZE = x_dim - y_dim/2;

  dx = 0.5;
  dy = 0.5*sqrt(3.0);

  T = (iter + 0.0)*U/H;

  out=fopen("FIELD.dat","w");
  
  fprintf(out,"title=\"Re=%2.0f iter=%d  T=%2.2f  GRID=%dx%d\"\n ",
	  Re,iter, T, GRID_SIZE, GRID_SIZE);

  fprintf(out,"variables=X, Y, si, u, v, ro, P\n");
  fprintf(out,"zone t=\"d2z6\",i=%d,j=%d\n",x_dim-2,y_dim-2);
  

  // ---------------- calculate the field variables ------------------
  k = 0;
  for(j=y_dim-2;j>=0;--j){
    for(i=1;i<=x_dim-2;++i){
        sai = si[i][j];

        state = get_state(i,j);

        ro = state.ro1;
        u  = state.ux;
        v  = state.uy;
		u__2 = u*u + v*v;

		xx = (i + 0.0)/(y_dim - 1.0);
		yy = (j + 0.0)/(y_dim - 1.0);

		P = calc_P(i, j);

		fprintf(out,"%f   %f   %f  %1.12f  %1.12f   %1.12f   %1.12f\n",
        xx,  yy,  sai, u/U,  v/U,  ro, P);

    }
  }
  fclose(out);

}


/*---------------------------------------------------------------------------------------------*/
/*---------------------->       cpu_run forward one time increment        <-----------------------*/
/*---------------------------------------------------------------------------------------------*/
void simulate(void){

  int i,j,id, b_state;


  // ----------   x direction gradient --------
   reset_F_in();

  // --- colision stage ---
  for(i=0;i<=x_dim-1;++i){
    for(j=0;j<=y_dim-1;++j){
		id = j*x_dim + i;
        b_state = bdr_state[id];
        if(b_state == 0){ord_col(i,j);}
        if(b_state == 1){bdr_col(i,j);}
    }
  }

  // --- streaming stage ---
  for(i=0;i<=x_dim-1;++i){
    for(j=0;j<=y_dim-1;++j){
        stream(i,j);
    }
  }

}

// ************************************************************************************
//                                     Display function
// ************************************************************************************
void display(void){

	  double y0, T;

      y0 = 0.5*sqrt(3.0)*(y_dim-1.0);
      T = (iter + 0.0)*U/y0;

      simulate(); 

	  if(iter > n_iter_01){
	    printf(" Iter = %d     T=%2.2f\n", iter, T);
	    n_iter_01 = n_iter_01+d_iter_01;
      }

	  if(iter > n_iter_02){
	    printf("\n\n updating the output data files ... \n\n");
        rec_si();
        rec_field();
	    n_iter_02=n_iter_02+d_iter_02;
      }

	  if(T >= MAX_T){
        rec_field();
		exit(0);
	  }

	  ++iter;
}

/*---------------------------------------------------------------------------------------------*/
/* ----------------->                 The main program                     <-------------------*/
/*---------------------------------------------------------------------------------------------*/
int main(int argc, char **argv){    
	 iter = 0;
	 allocate_memory();
	 initialize();

	 while ( 1==1){display();}


}

