/*---------------------------------------------------------------------------------------------*/
/*  Abed Zadehgol            Graphical version (with graphical display)            2016/09/24  */
/*  email:  a.zadehgol@me.iut.ac.ir     &    abed.zdehgol@gmail.com                version: 2  */
/*                                                                                             */
/*  This is a sample C code for the d2z6 model proposed in [J. Comp. Phys. 274, 803 (2014)]    */
/*  [Phys. Rev. E 91, 063311 (2015)]. This sample code simulates the flow of a visous and      */
/*  nearly incompressible fluid inside a 2D square cavity.                                     */
/*                                                                                             */
/*  Additional notes:                                                                          */
/*    (1) To run this version of the program, you first need to install GLUT. Next, you need   */
/*         to set the compiler, properly. The present code is a complete and self-sufficient   */  
/*         version (additional codes are not needed to complile and run this code). This       */
/*         program has successfuly been complied in Microsoft Visual Studio 2008 environment.  */
/*                                                                                             */
/*    (2) The program generates an output data file, "FIELD.DAT", in "Tecplot" format.         */
/*         The output file contains the information of the following fields:                   */
/*                [a] si        (for the steamlines)                                           */ 
/*                [b] velocity  (ux, uy)                                                       */
/*				  [c] density   (ro)                                                           */
/*				  [d] pressure  (P)                                                            */
/*                [e] entropy   (S)                                                            */ 
/*                                                                                             */
/*    (3) You can use the following hotkeys (in the graphical version only)                    */
/*                [a] Press "r" to update the output data file (in Tecplot format)             */
/*                [b] Press "e" to exit the program                                            */
/*                [c] Press "+" or "-" to change the intensity of the colors                   */
/*                                                                                             */
/*    (4)   To see the contour plots, load the file FIELD.DAT into Tecplot, go to the          */
/*           contour mode, switch from "flood" to "lines", and choose the desired field.       */
/*                                                                                             */
/*    (5)   Adjust the adjustable parameters below.                                            */
/*                                                                                             */

// ---------------------------------  the adjustable parameters --------------------------------

#define BS                 2              //  grid size (increase BS for finer grids) 
#define Re                 400.0        //  Reynolds nunber (try the 100 - 10,000 range)
#define MACH               0.10           //  Mach number (try the 0.025 - 0.10 range)
#define MAX_T              500.0          //  Duration of simulation (in dimensionless time)

// ---------------------------------   the rest of the program  -------------------------------
#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <stdlib.h>
#include "glut.h"

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
#define U                  (MACH*cs)

#define y_dim              BS*77
#define x_dim              BS*105

#define pi                 (4.0*atan(1.0))

#define d_iter_01          100
#define d_iter_02          2000
#define d_iter_03          100

#define MAX_DIR            7

//-------------------------  variable declarations ------------------------
double H, visc, COEF;
int   iter, n_iter_01, n_iter_02, n_iter_03, x_pos, y_pos;
double si[x_dim][y_dim];

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

  FILE *out;

  int i,j, ss;
  double xx,yy,u,v;
  double L_xx = (x_dim+0.0) - 0.5*(y_dim+0.0);
  double L_yy = sin(pi/3.0)*(y_dim+0.0);

  /*------  v/U values at horizontal centerline --------*/
  out=fopen("v_d2z6.dat","w");
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

  /*------  u/U values at vertical centerline --------*/
  out=fopen("u_d2z6.dat","w");
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

  out=fopen("FIELD.dat","w");
  
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

/*---------------------------------------------------------------------------------------------*/
/* --------------                           show_results                      -----------------*/
/*---------------------------------------------------------------------------------------------*/
void show_results(void){

	int i,j,id, b_state;
	double xp1,xp2,xp3,xp4, xx;
	double yp1,yp2,yp3,yp4, yy;
	double del_x, del_y, VELOCITY;
	
	double xx_ctr, yy_ctr;

    double u,v;
	double pc0, pc1, pc2,pc3;

	del_x = 1.25/(x_dim+0.0);
	del_y = 1.25/(y_dim+0.0);

	// calculate and display the field values every n_iter_01 iterations 
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
    glLoadIdentity(); 

    for(j=0; j<y_dim-1;++j){
	    for(i=0; i<x_dim-1;++i){

		  yy = (j + 0.0)*sqrt(3.0)/2.0;
		  xx = (i + 0.0) + 0.5*(j+0.0);

		  xx_ctr = -1.6 + 1.75*xx/(y_dim+0.0);
		  yy_ctr = -0.8 + 1.75*yy/(y_dim+0.0);

		  u = calc_u(i,j);
		  v = calc_v(i,j);
		  VELOCITY = sqrt(u*u + v*v);

          xp1 = xx_ctr - del_x;
		  yp1 = yy_ctr - del_y;

		  xp2 = xx_ctr + del_x;
		  yp2 = yy_ctr - del_y;

		  xp3 = xx_ctr + del_x;
		  yp3 = yy_ctr + del_y;

		  xp4 = xx_ctr - del_x;
		  yp4 = yy_ctr + del_y;

		  id = j*x_dim + i;

          b_state = bdr_state[id];
		  if(b_state == 1){pc0 = 0.0;}else{pc0 = VELOCITY*COEF;}

		  pc1 = 1.0*pc0 + 0.1*sin(10*pc0); 
		  pc2 = 1.5*pc0 + 0.2*sin(20*pc0);  
		  pc3 = 2.5*pc0 + 0.1*sin(30*pc0);

          // draw a solid rectangle 
          glBegin(GL_QUADS);//Start drawing quads
          glColor3f(pc1, pc2 , pc3);  
          glVertex2f(xp1,yp1);//first coordinate
          glVertex2f(xp2,yp2);//second coordinate
          glVertex2f(xp3,yp3);//third coordinate (now blue)
          glVertex2f(xp4,yp4);//last coordinate
          glEnd();//Stop drawing quads
	    }
	}

    glutSwapBuffers();
    glutPostRedisplay();	  
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
        rec_params();
        rec_field();
	    n_iter_02=n_iter_02+d_iter_02;
      }

	  if(iter > n_iter_03){
	    show_results();
	    n_iter_03 = n_iter_03 + d_iter_03;
      }

	  if(T >= MAX_T){
        rec_field();
		exit(0);
	  }

	  ++iter;
}


/*-------------------------------   define some usefull hotkeys  ------------------------------*/
/*---------------------------------------------------------------------------------------------*/
static void Key(unsigned char key, int x, int y){ 
  if(key == 'e'){	  exit(0);  }

  if(key == '+'){	  
	  COEF *= 1.05;  
	  printf("\n COEF = %2.2f <----- \n\n", COEF);
  }

  if(key == '-'){	  
	  COEF /= 1.05;  
	  printf("\n COEF = %2.2f <----- \n\n", COEF);
  }

  if(key == 'r'){
	  printf("\n recording the field data ...\n\n");
	  rec_si();
	  rec_field();
  }

}

/*---------------------------------------------------------------------------------------------*/
/* ----------------->                 The main program                     <-------------------*/
/*---------------------------------------------------------------------------------------------*/
int main(int argc, char **argv){ 

	 iter = 0;

	 allocate_memory();
	 initialize();

	 // ---------------------   prepare for graphics  -------------------
     glutInit(&argc, argv);     
     glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH); 
     glutInitWindowSize(windowWidth, windowHeight);  
     glutInitWindowPosition(windowPosX, windowPosY); 
     glutCreateWindow(title);  
	 glMatrixMode(GL_MODELVIEW);
	 glLoadIdentity();
     glClearColor (0.0, 0.0, 0.0, 0.0);
     glutDisplayFunc(display);    
     glutKeyboardFunc(Key);  
     glutIdleFunc(display);  
     glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
     glutMainLoop();   

}

