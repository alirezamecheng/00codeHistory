# In the name of Allah.
# here is the code for D2Z6 CSKM. written by Abed Zadehgol in C++
# and translated to Julia language by Alireza Ghavami nia.
workspace();


# ---------------------------------  the adjustable parameters --------------------------------

const BS =                2              #  grid size (increase BS for finer grids)
const Re  =               400.0          #  Reynolds nunber (try the 100 - 10,000 range)
const MACH =              0.10           #  Mach number (try the 0.025 - 0.10 range)
const MAX_T =             500.0          #  Duration of simulation (in dimensionless time)

global title = "A sample CSKM code (graphics included)";
global windowHeight = 400;                  # Windowed
global windowWidth  = 400;                  # Windowed mode's width
global windowPosX   = 400;                  # Windowed mode's top-left corner x
global windowPosY   = 100;                  # Windowed mode's top-left corner y

const V  =                1.0

const ALF    =            0.0  # rest particle fraction
global gama    =           ALF  #(ALF - u__2)

global cs2       =         0.5*(1.0 - ALF)*V*V
global cs   =              (sqrt(cs2))

const DENSITY =            1.0
global U        =          (MACH*cs)

global y_dim =             BS*77
global x_dim   =           BS*105

# const pi =                (4.0*atan(1.0))

const d_iter_01   =       100
const d_iter_02    =      2000
const d_iter_03     =     100

const MAX_DIR        =    7

#=
Uses of non-constant globals can be optimized by annotating their types at the point of use:

global x
y = f(x::Int + 1)
=#

#-------------------------  variable declarations ------------------------
# double H, visc, COEF;
# int   iter, n_iter_01, n_iter_02, n_iter_03, x_pos, y_pos;
# double si[x_dim][y_dim];

# ---------------- structure for the constant parameters -----------------
# int    ex[MAX_DIR]; int ey[MAX_DIR];
# double cx[MAX_DIR]; double cy[MAX_DIR];

# ------ structure for the transferable grid dependent variables ---------
# typedef struct{ double DIR[MAX_DIR]; } EQU;
#
# EQU    *F_in, *F_out;
# int    *bdr_state;


#=---------------------------------------------------------------------------------------------*/
/*------------------>                 allocate_memory                       <------------------*/
/*---------------------------------------------------------------------------------------------=#

function allocate_memory()

	size_a::Int64 = (x_dim*y_dim)*sizeof(EQU);
	size_b::Int64 = (x_dim*y_dim)*sizeof(Int64);

	// allocate memory for all transferable and resident cpu and gpu veriables
	F_in  = (EQU*) malloc(size_a);
	F_out = (EQU*) malloc(size_a);

	bdr_state = (int*) malloc(size_b);

end # end function [allocate_memory]

#=---------------------------------------------------------------------------------------------*/
/*------------------>                  d2z6 formulation                     <------------------*/
/*---------------------------------------------------------------------------------------------=#
function get_equ(double u1, double u2, double ro) # خروجی این فانکشن یک دیتا تایپ ای کیو یو است

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
