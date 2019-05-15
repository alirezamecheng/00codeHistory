/* 
	* Version 2
		1- Adding file Output.
		2- Adding some print to stdout fore more information.
		3- Adding error criteria. 







*/


/*---------------------------------------------------------------------------------------------*/
/*  Abed Zadehgol              Plain version (without graphical display)         2016/09/24    */
/*  email:  a.zadehgol@me.iut.ac.ir     &    abed.zdehgol@gmail.com              version: 2    */
/*                                                                                             */
/*  This is a sample C code for the d2z6 model proposed in [J. Comp. Phys. 274, 803 (2014)]    */
/*  & [Phys. Rev. E 91, 063311 (2015)]. This sample code simulates the flow of a visous and    */
/*  nearly incompressible fluid inside a 2D square cavity.                                     */
/*                                                                                             */
/*  Additional notes:                                                                          */
/*                                                                                             */
/*    (1) The present code is a complete and self-sufficient version (additional codes are     */
/*        not needed to complile and run this code). This program has successfuly been         */
/*        complied in Microsoft Visual Studio 2008 environment.                                */
/*                                                                                             */
/*    (2) The program generates an output data file, "FIELD.DAT", in "Tecplot" format.         */
/*         The output file contains the information of the following fields:                   */
/*                [a] si        (for the steamlines)                                           */ 
/*                [b] velocity  (ux, uy)                                                       */
/*				  [c] density   (ro)                                                           */
/*				  [d] pressure  (P)                                                            */
/*                [e] entropy   (S)                                                            */ 
/*                                                                                             */
/*    (3)   To see the contour plots, load the file FIELD.DAT into Tecplot, go to the          */
/*           contour mode, switch from "flood" to "lines", and choose the desired field.       */
/*                                                                                             */
/*    (4)   Adjust the adjustable parameters below.                                            */
/*                                                                                             */

// ---------------------------------  the adjustable parameters --------------------------------

#define BS                 52              //  grid size (increase BS for finer grids) 
#define Re                 30000.0          //  Reynolds nunber (try the 100 - 10,000 range)
//#define MACH               0.10           //  Mach number (try the 0.025 - 0.10 range)
#define U                  (0.1)
#define MAX_T              500.0          //  Duration of simulation (in dimensionless time)
#define MAX_ITTERATION     500000         // Max ittertation 
// ---------------------------------   the rest of the program  -------------------------------
#include <stdio.h>
//#include <conio.h>
#include <math.h>
#include <stdlib.h>

char title[] = "A sample CSKM code (graphics included)";  
int windowHeight = 400;                  // Windowed mode's height
int windowWidth  = 400;                  // Windowed mode's width
int windowPosX   = 400;                  // Windowed mode's top-left corner x
int windowPosY   = 100;                  // Windowed mode's top-left corner y

#define V                  1.0

#define ALF                0.0  // rest particle fraction 
#define gama               ALF  //(ALF - u__2)    

#define cs2                0.5*(1.0 - ALF)*V*V
#define cs                 (sqrt(cs2))

#define DENSITY            1.0
//#define U                  (MACH*cs)

#define y_dim              BS*3	// originally 77
#define x_dim              BS*4	// originally 105

#define pi                 (4.0*atan(1.0))

#define d_iter_01          100  // Showing ittertations intervals
#define d_iter_02          1000 // File Outputs intervals
#define d_iter_03          100	// does nothing here.
#define NAMEFRACTION	   1000000

#define MAX_DIR            7

//-------------------------  variable declarations ------------------------
double H, visc, COEF;
int   iter, n_iter_01, n_iter_02, n_iter_03, x_pos, y_pos;
double si[x_dim][y_dim];

double Vold,Vnew;

// ---------------- structure for the constant parameters -----------------
int    ex[MAX_DIR]; int ey[MAX_DIR]; 
double cx[MAX_DIR]; double cy[MAX_DIR]; 
  
// ------ structure for the transferable grid dependent variables ---------
typedef struct{ double DIR[MAX_DIR]; } EQU;

EQU    *F_in, *F_out;
int    *bdr_state;

/*---------------------------------------------------------------------------------------------*/
/*------------------>                 allocate_memory                       <------------------*/
/*---------------------------------------------------------------------------------------------*/
void allocate_memory(void){

	int size_a = (x_dim*y_dim)*sizeof(EQU);	
	int size_b = (x_dim*y_dim)*sizeof(int);	

	// allocate memory for all transferable and resident cpu and gpu veriables 
	F_in  = (EQU*) malloc(size_a);
	F_out = (EQU*) malloc(size_a);

	bdr_state = (int*) malloc(size_b);

}

/*---------------------------------------------------------------------------------------------*/
/*------------------>                  d2z6 formulation                     <------------------*/
/*---------------------------------------------------------------------------------------------*/
EQU get_equ(double u1, double u2, double ro){

   EQU res;
   double e;
   double v1, v2, uv, u__2, v__2;
   int i;

   u__2 = u1*u1 + u2*u2;

   // ----------------- calc the streaming vector  --------------
   for(i=1; i<=6; ++i){

     v1 = cx[i];
     v2 = cy[i];
	 v__2 = v1*v1 + v2*v2;

     uv = v1*u1 + v2*u2;
	 e = (v__2 - u__2)/(v__2 + u__2 - 2.0*uv);

     res.DIR[i] = ro*(e - gama)/6.0;
   }

   res.DIR[0] = ro*gama;

   return(res);

}

/*---------------------------------------------------------------------------------------------*/
/*--------------->             setting the boundar nodes (at the upper lid)       <------------*/
/*---------------------------------------------------------------------------------------------*/
void reset_driving_nodes(void){
  int i,j, id;

  EQU equ;
  equ = get_equ(U, 0.0, DENSITY);

  //  -------------  the left and right triangles ----------------
  j=y_dim-1;
  for(i=0;i<=(x_dim-1);++i){
	   id = j*x_dim+i;
	   F_out[id] = equ;
  }
}

/*---------------------------------------------------------------------------------------------*/
/*----------------->              initialize the main parameters           <-------------------*/
/*---------------------------------------------------------------------------------------------*/
void init_params(void){

  int k;
  double ang;

  COEF = 15.0;
  
  // ---- geometric parameters ----
  H=(y_dim+0.0)*sqrt(3.0)/2.0;
  
  // --- fluid parameters ----
  visc = H*U/Re;

  // --- code parameters ---
  n_iter_01 = d_iter_01;
  n_iter_02 = d_iter_02;
  n_iter_03 = d_iter_03;

  x_pos=0.5*(x_dim+0.0);
  y_pos=0.5*(y_dim+0.0);

  // ---- intervals for the neighbour nodes ---
  for(k=1;k<=6;++k){
	  ang = (k-1)*pi/3.0;
      cx[k] = V*cos(ang); 
      cy[k] = V*sin(ang); 
  }

  cx[0]=  0.0;  
  cy[0]=  0.0;

  // ------ intervals for the neighbour nodes -------
  ex[1] =  1;  ey[1] =  0;
  ex[2] =  0;  ey[2] =  1;
  ex[3] = -1;  ey[3] =  1;
  ex[4] = -1;  ey[4] =  0;
  ex[5] =  0;  ey[5] = -1;
  ex[6] =  1;  ey[6] = -1;
  ex[0] =  0;  ey[0] =  0;
}


/*---------------------------------------------------------------------------------------------*/
/*-------------------->    specify the boundry type of all nodes   <---------------------------*/
/*---------------------------------------------------------------------------------------------*/
void init_bdrs(void){
  double x0, x1, x;
  int i,j, id;
  
  // -------------- initializing the fluid nodes ----------------
  for(i=0;i<=(x_dim-1);++i){
    for(j=0;j<=(y_dim-1);++j){
	  id = j*x_dim + i;
      bdr_state[id] = 0;
    }
  }

  //  -------------  the left and right triangles ----------------
  for(i=0;i<=(x_dim-1);++i){
    for(j=0;j<=(y_dim-1);++j){
	  id = j*x_dim+i;
	  x = (i+0.0) + (j+0.0)*cos(pi/3.0);
	  x0 = (y_dim+0.0)*cos(pi/3.0) + 1.0;
	  x1 = (x_dim+0.0) - 1.0;
	  if(x<x0 || x>x1){bdr_state[id] = 1;}
	}
  }

   //  -------------  defining the bottom wall  ----------------
  j = 0;
  for(i=0;i<=(x_dim-1);++i){
	  id = j*x_dim+i;
	  bdr_state[id] = 1;
  }

   //  -------------  defining the left wall  ----------------
  i = 0;
  for(j=0;j<=(y_dim-1);++j){
	  id = j*x_dim+i;
	  bdr_state[id] = 1;
  }

   //  -------------  defining the right wall  ----------------
  i = x_dim-1;
  for(j=0;j<=(y_dim-1);++j){
	  id = j*x_dim+i;
	  bdr_state[id] = 1;
  }

}

/*---------------------------------------------------------------------------------------------*/
/*----->                           inittialize all nodes                                 <-----*/
/*---------------------------------------------------------------------------------------------*/
void init_nodes(void){
  int id, i,j,k;

  EQU equ;
  equ = get_equ(0.0, 0.0, DENSITY);

  /* initializing the general nodes */
  for(j=0;j<=y_dim-1;++j){
    for(i=0;i<=x_dim-1;++i){
	  id = j*x_dim+i;
      for(k=0;k<=6;++k){
	    F_in[id].DIR[k]  = equ.DIR[k];
	    F_out[id].DIR[k] = equ.DIR[k];
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
/*                              calculate the local density                                    */
/*---------------------------------------------------------------------------------------------*/
double calc_ro(int i,int j){
  double res;
  int id, k;
 
  id = j*x_dim+i;

  res=0.0;
  for(k=0;k<=6;++k){res +=  F_in[id].DIR[k];}
  return(res);
}

/*---------------------------------------------------------------------------------------------*/
/*                       calculate the u-component of the velocity                             */
/*---------------------------------------------------------------------------------------------*/
double calc_u(int i,int j){
  int id, k;
  double res,den, v1;

  id = j*x_dim+i;

  res=0.0;
  den=calc_ro(i,j);
  for(k=1;k<=6;++k){
	  v1 = cx[k];
	  res += F_in[id].DIR[k]*v1;
  }
  if(den != 0.0){res=res/den;}else{res=0.0;}
  return(res);
}

/*---------------------------------------------------------------------------------------------*/
/*                        calculate the v-component of the velocity                            */
/*---------------------------------------------------------------------------------------------*/
double calc_v(int i,int j){
  int id, k;
  double res,den, v2;

  id = j*x_dim+i;

  res=0.0;
  den=calc_ro(i,j);
  for(k=1;k<=6;++k){
	  v2 = cy[k];
	  res += F_in[id].DIR[k]*v2;
  }
  if(den != 0.0){res=res/den;}else{res=0.0;}
  return(res);
}


/*---------------------------------------------------------------------------------------------*/
/*                              calculate the local density                                    */
/*---------------------------------------------------------------------------------------------*/
double calc_S(int i,int j){
  double res, m;
  int id, k;
 
  id = j*x_dim+i;

  res=0.0;
  for(k=1;k<=6;++k){
	  m = F_in[id].DIR[k];
	  res +=  log(m);
  }

  res = res/6.0;

  return(res);
}
/*---------------------------------------------------------------------------------------------*/
/*                              calculate the local density                                    */
/*---------------------------------------------------------------------------------------------*/
double calc_P(int i,int j){
  double res, ro, u, v, u__2, v__2;

  ro = calc_ro(i,j);
  u = calc_u(i,j);
  v = calc_v(i,j);

  v__2 = V*V;
  u__2 = u*u + v*v;

  res = 0.5*ro*(v__2 - u__2 - gama);
 

  return(res);
}


/*---------------------------------------------------------------------------------------------*/
/* --------------->               non-boundary collisions                      <---------------*/
/*---------------------------------------------------------------------------------------------*/
void ord_col(int i,int j){

  int id, k;
  double ro, u1, u2, u__2, tau;

  EQU equ_dist;

  id = j*x_dim+i;

  /*-- calculating the local properties ---*/
  ro=calc_ro(i,j);
  u1=calc_u(i,j);
  u2=calc_v(i,j);
  u__2 = u1*u1 + u2*u2;

  tau = 4.0*visc/(V*V - u__2) + 0.5;

  equ_dist = get_equ(u1, u2, ro);

  for(k=0;k<=6;++k){
    F_out[id].DIR[k]  = F_in[id].DIR[k]-(F_in[id].DIR[k] - equ_dist.DIR[k])/tau;
  }

}

/*---------------------------------------------------------------------------------------------*/
/*--------------->                    boundary collisions                  <-------------------*/
/*---------------------------------------------------------------------------------------------*/
void bdr_col(int i,int j){
   int id = j*x_dim+i;

   F_out[id].DIR[1] = F_in[id].DIR[4];
   F_out[id].DIR[2] = F_in[id].DIR[5];
   F_out[id].DIR[3] = F_in[id].DIR[6];
   F_out[id].DIR[4] = F_in[id].DIR[1];
   F_out[id].DIR[5] = F_in[id].DIR[2];
   F_out[id].DIR[6] = F_in[id].DIR[3];
   F_out[id].DIR[0] = F_in[id].DIR[0];

}

/*---------------------------------------------------------------------------------------------*/
/*                Post-Collision streaming effect on the neighbouring nodes.                   */
/*---------------------------------------------------------------------------------------------*/
void stream(int i,int j){
  
  int id, nid, k,ni,nj;
  id = j*x_dim + i;
  
  for(k=0;k<=6;++k){
    ni = i+ex[k];
    nj = j+ey[k];
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
  double dx, dy;

  dx = 0.5;
  dy = 0.5*sqrt(3.0);

  /*--------- calculating the si field ----------*/
  for(i=0;i<=x_dim-1;++i){ si[i][0]=0.0; si[i][y_dim-1]=0.0;}
  for(j=0;j<=y_dim-1;++j){ si[0][j]=0.0; si[x_dim-1][j]=0.0;}

  for(i=1;i<=x_dim-2;++i){
    for(j=y_dim-2;j>=1;--j){
	  si[i][j] = si[i-1][j] - (calc_v(i-1,j)*dx - calc_v(i-1,j)*dy);	
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
/*---                        Preparing and recording the field data                        ----*/
/*---------------------------------------------------------------------------------------------*/
void rec_params(void){

	char fileNameString[32]; // creating file name buffer.
	// put "file" name in the string  it makes file name in the string each time.
	snprintf(fileNameString,sizeof(char)*32, "v_d2z6_%lf.dat",(double)iter/NAMEFRACTION);
				
  FILE *out;

  int i,j, ss;
  double xx,yy,u,v;
  double L_xx = (x_dim+0.0) - 0.5*(y_dim+0.0);
  double L_yy = sin(pi/3.0)*(y_dim+0.0);

  /*------  v/U values at horizontal centerline --------*/
  out=fopen(fileNameString,"w");
  fprintf(out,"title=\"d2z6      Re=%2.0f     iterations=%d \"\n ",Re,iter);
  fprintf(out,"zone t=\"d2z6\"");  
  fprintf(out,"variables=X,v/U\n");

  j = y_dim/2;
  for(i=0;i<=x_dim-1;++i){
    xx = ((i + 0.0) + 0.5*(j+0.0) - 0.5*(y_dim+0.0))/L_xx;	  
	if(xx >= 0.0 && xx <= 1.0){
	   	v=calc_v(i,y_dim/2)/U;
        fprintf(out,"%f %f\n",xx,v);
	}
  }

  fclose(out);

  char fileNameString2[32]; // creating file name buffer.
	// put "file" name in the string  it makes file name in the string each time.
	snprintf(fileNameString2,sizeof(char)*32, "u_d2z6_%lf.dat",(double)iter/NAMEFRACTION);
	
  /*------  u/U values at vertical centerline --------*/
  out=fopen(fileNameString2,"w");
  fprintf(out,"title=\"d2z6      Re=%2.0f     iterations=%d \"\n ",Re,iter);
  fprintf(out,"zone t=\"d2z6\"");  
  fprintf(out,"variables=u/U,Y\n");
  for(j=0;j<=y_dim-2;++j){
     yy =  (j+0.0)*sin(pi/3.0)/L_yy;
	 ss = 0;
     for(i=0;i<=x_dim-1  && ss==0 ;++i){
	  xx = ((i + 0.0) + 0.5*(j+0.0) - 0.5*y_dim)/L_xx;
	  if(xx>0.5 ){
		 ss = 1;
	     u=calc_u(i,j)/U;
         fprintf(out,"%f %f\n",u,yy);
	  }
    }
  }
  fclose(out);
}


/*---------------------------------------------------------------------------------------------*/
/*---                     Preparing and recording the field data                           ----*/
/*---------------------------------------------------------------------------------------------*/
void rec_field(void){
	char fileNameString3[32]; // creating file name buffer.
	// put "file" name in the string  it makes file name in the string each time.
	snprintf(fileNameString3,sizeof(char)*32, "Results_%lf.dat",(double)iter/NAMEFRACTION);
				
  FILE *out;
  int i,j;

  double xx, yy, ro, u, v, P, S, T, h0;
  double x0, x1, y0,  sai;

  int GRID_SIZE = x_dim - y_dim/2;

  h0 = 0.5*sqrt(3.0)*(y_dim-1.0);
  T = (iter + 0.0)*U/h0;

  x0 = 0.5*(y_dim-1.0);
  x1 = x_dim - x0;
  y0 = 0.5*sqrt(3.0)*(y_dim-1.0);
  T = (iter + 0.0)*U/y0;

  out=fopen(fileNameString3,"w");
  
  fprintf(out,"title=\"Re=%2.0f iter=%d  T=%2.2f  GRID=%dx%d\"\n ",
	  Re,iter, T, GRID_SIZE, GRID_SIZE);

  fprintf(out,"variables=X, Y, si, u, v, ro, P, S\n");
  fprintf(out,"zone t=\"d2z6\",i=%d,j=%d\n",x_dim,y_dim);
  

  // ---------------- calculate the field variables ------------------
  for(j=y_dim-1;j>=0;--j){
    for(i=0;i<=x_dim-1;++i){
        sai = si[i][j];

		S = calc_S(i,j);

		ro = calc_ro(i,j);
		u  = calc_u(i,j);
        v  = calc_v(i,j);

		xx=(i-0+0.0) + cos(pi/3.0)*(j-0+0.0);
        yy =  sin(pi/3.0)*(j-0+0.0);
	    P = calc_P(i, j);

		xx = (xx - x0)/x1;
		yy = yy/y0;
		if(xx < 0.0 || xx > 1.0 || yy < 0.0 || yy> 1.0){ 

			sai = 0.0;
			P = 0.0;
			u = 0.0;
			v = 0.0; 
			ro = 0.0;
			S = 0.0;
		}

		fprintf(out,"%f   %f   %f  %1.12f  %1.12f   %1.12f   %1.12f   %1.12f\n",
        xx,  yy,  sai, u/U,  v/U,  ro, P, S);
    }
  }
  fclose(out);
     
	 printf("Files has been saved at %s Re(%1.1f) U(%1.4f) G(%d)\n \n \n",fileNameString3,Re,U,GRID_SIZE);
	 fflush(stdout);
}



/*---------------------------------------------------------------------------------------------*/
/*---------------------->       cpu_run forward one time increment        <-----------------------*/
/*---------------------------------------------------------------------------------------------*/
void simulate(void){

  int i,j,id, b_state;

//  reset_driving_nodes();

  // --- colision stage ---
  for(i=0;i<=x_dim-1;++i){
    for(j=0;j<=y_dim-1;++j){
		id = j*x_dim + i;
      b_state = bdr_state[id];
      if(b_state == 0){ord_col(i,j);}
      if(b_state == 1){bdr_col(i,j);}
    }
  }

  reset_driving_nodes();
  // --- streaming stage ---
  for(i=0;i<=x_dim-1;++i){
    for(j=0;j<=y_dim-1;++j){
      stream(i,j);
    }
  }

  // --- streaming stage ---
  for(i=0;i<=x_dim-1;++i){
    for(j=0;j<=y_dim-1;++j){
		id = j*x_dim + i;
    }
  }
}

// ************************************************************************************
//                                     Display function
// ************************************************************************************
void display(void){

	  double y0, T, MyError;
		int iii = x_dim/2;
		int jjj = y_dim/2;
      y0 = 0.5*sqrt(3.0)*(y_dim-1.0);
      T = (iter + 0.0)*U/y0;

      simulate(); 
	   
	   Vnew  = calc_v(iii,jjj);
	   MyError = fabs((Vold - Vnew)/Vold);
	   Vold = Vnew;
	   
	  if(iter > n_iter_01){
	    printf(" Iter = %d | T=%2.2f | Center Velocity = %1.12f | Relative Error = %1.12f\n", iter, T, Vnew,MyError);
	    n_iter_01 = n_iter_01+d_iter_01;
      }

	  if(iter > n_iter_02){
	    printf("\n\n updating the output data files ... \n\n");
        rec_si();
        rec_params();
        rec_field();
	    n_iter_02=n_iter_02+d_iter_02;
      }


	  //if(T >= MAX_T){
	  if(iter >= MAX_ITTERATION){
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
	 
	Vold = 0.000;
	
	 allocate_memory();
	 initialize();

		int GRID_SIZE = x_dim - y_dim/2;
		printf(" Grid Size = %d \n", GRID_SIZE);
		
	 while ( 1 == 1){
	 display();
	 fflush(stdout);
	 }
}

	 
