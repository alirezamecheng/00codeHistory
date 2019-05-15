workspace()






const BS = 30              #  grid size (increase BS for finer grids)
const Re = 100.0           #

const MACH = 0.1           #  Mach number (try the 0.025 - 0.10 range)
# const TestVisc = 0.000837520938023436






const StabilityRun	= "Yes" # if Yes will break at MAX_UnSteady_Itteration
# const StabilityRun	= "No" # if defined will break at MAX_UnSteady_Itteration

const MAX_UnSteady_Itteration  =   100001         # Max ittertation
const DesiredError           =     0.00001         # 1*10^-5


const MakeErrorCSVFile = "Yes"       # if defined will make if comment will not make
#const MakeErrorCSVFile = "No"       # if defined will make if comment will not make
#const MakeMultipleFile        # if defined will make if comment will not make


const d_iter_01 = 1  # error calc intervals
const d_iter_02 = 10  # File Outputs intervals
const d_iter_03 = 50	 # showing error intervals




const NAMEFRACTION	 =  1000000


#const MAX_T = 500.0          #  Duration of simulation (in dimensionless time)

const V = 1.0

const ALF = 0.0  # rest particle fraction
const gama = ALF  #(ALF - u__2)

const cs2 = 0.5*(1.0 - ALF)*V*V
const cs = (sqrt(cs2))

const DENSITY = 1.0
const U = (MACH*cs)

const y_dim = BS*3	# originally 77
const x_dim = BS*4	# originally 105



const MAX_DIR = 7

#const U                  (0.122474487) # equal mack
# for k=2:7
#     ang = (k-2)*pi/3.0;
#     cx[k] = V*cos(ang);
#     cy[k] = V*sin(ang);
# end
# cx[1]=  0.0;
# cy[1]=  0.0;




const COEF = 15.0;

# ---- geometric parameters ----
const H=(y_dim+0.0)*sqrt(3.0)/2.0;

# --- fluid parameters ----
const visc = H*U/Re;
# visc = TestVisc;


const x_pos=floor(Int64,0.5*(x_dim+0.0));
const y_pos=floor(Int64,0.5*(y_dim+0.0));



const cx = [0,V*cos(0),V*cos(pi/3),V*cos(2*pi/3),V*cos(3*pi/3),V*cos(4*pi/3),V*cos(5*pi/3)];
const cy = [0,V*sin(0),V*sin(pi/3),V*sin(2*pi/3),V*sin(3*pi/3),V*sin(4*pi/3),V*sin(5*pi/3)];
const ex = [0,1,0,-1,-1,0,1];
const ey = [0,0,1,1,0,-1,-1];



const MakeMultipleFile = "Yes"
const resultsAddress = "../Results/CSKM_$(Dates.format(now(),"YY.mm.dd.HH.MM"))_Re($(Re))_Gr($(x_dim)_$(y_dim))_Uw($(U))"; # for using in cluster.

# const resultsAddress = "../Results/$(MethodName)_$(Dates.format(now(),"YY.mm.dd.HH.MM"))_Re($(Re))_Gr($(Nx1)_$(Ny1))_Uw($(uw_top))"; # for using in my own pc.

mkpath(resultsAddress); # to creat specified path


function get_equ(u1,u2,ro)
    res = Array{Float64}(7);
    u__2 = u1*u1 + u2*u2;
    for i = 2:7
        v1 = cx[i]
        v2 = cy[i]
        v__2 = v1*v1 + v2*v2;
        uv = v1*u1 + v2*u2;
	    e1 = (v__2 - u__2)/(v__2 + u__2 - 2.0*uv);
        res[i] = ro*(e1 - gama)/6.0;

    end
    res[1] = ro*gama;
    return res;
end

function reset_driving_nodes(F_out)
    equ = get_equ(U,0.0,DENSITY);

    # -------------  the left and right triangles ----------------

    j = y_dim;
    for i = 1:x_dim
        for k = 1:MAX_DIR
        F_out[i,j,k] = equ[k];
    end
    end
end

# function init_params()
#     # int k;
#     #   double ang;
#
#     #   COEF = 15.0;
#     #
#     #   # ---- geometric parameters ----
#     #   H=(y_dim+0.0)*sqrt(3.0)/2.0;
#     #
#     #   # --- fluid parameters ----
#     #  visc = H*U/Re;
#     # # visc = TestVisc;
#     #
#     #   # --- code parameters ---
#     #   n_iter_01 = d_iter_01;
#     #   n_iter_02 = d_iter_02;
#     #   n_iter_03 = d_iter_03;
#     #
#     #   x_pos=floor(Int64,0.5*(x_dim+0.0));
#     #   y_pos=floor(Int64,0.5*(y_dim+0.0));
#
#       # ---- intervals for the neighbour nodes ---
#     #   for k=2:7
#     # 	  ang = (k-2)*pi/3.0;
#     #       cx[k] = V*cos(ang);
#     #       cy[k] = V*sin(ang);
#     #   end
#       #
#     #   cx[1]=  0.0;
#     #   cy[1]=  0.0;
#
#       # ------ intervals for the neighbour nodes -------
#     #   ex[2] =  1;  ey[2] =  0;
#     #   ex[3] =  0;  ey[3] =  1;
#     #   ex[4] = -1;  ey[4] =  1;
#     #   ex[5] = -1;  ey[5] =  0;
#     #   ex[6] =  0;  ey[6] = -1;
#     #   ex[7] =  1;  ey[7] = -1;
#     #   ex[1] =  0;  ey[1] =  0;
#
# end


function init_bdrs(bdr_state)

    # double x0, x1, x;
    #   int i,j, id;

      # -------------- initializing the fluid nodes ----------------
      for i=1:x_dim
        for j=1:y_dim
          bdr_state[i,j] = 0;
        end
    end

      #  -------------  the left and right triangles ----------------
      for i=1:x_dim
        for j=1:y_dim
    	  x = (i+0.0) + (j+0.0)*cos(pi/3.0);
    	  x0 = (y_dim+0.0)*cos(pi/3.0) + 1.0;
    	  x1 = (x_dim+0.0) - 1.0;
    	  if (x<x0 || x>x1)
              bdr_state[i,j] = 1;
              end
          end
      end

       #  -------------  defining the bottom wall  ----------------
      j = 1;
      for i=1:x_dim
    	  bdr_state[i,j] = 1;
      end

       #  -------------  defining the left wall  ----------------
      i = 1;
      for j=1:y_dim
    	  bdr_state[i,j] = 1;
      end

       #  -------------  defining the right wall  ----------------
      i = x_dim;
      for j=1:y_dim
    	  bdr_state[i,j] = 1;
      end

end


function init_nodes(F_in,F_out)
    # int id, i,j,k;

      # EQU equ;

      equ = get_equ(0.0, 0.0, DENSITY);

      # initializing the general nodes
      for i=1:x_dim
        for j=1:y_dim
            for k=1:7
    	        F_in[i,j,k] = equ[k];
    	        F_out[i,j,k] = equ[k];
            end
        end
    end
end

function initialize(bdr_state,F_in,F_out)
  # init_params();
  init_bdrs(bdr_state);
  init_nodes(F_in,F_out);
end

function calc_ro(i,j,F_in,F_out)
    # double res;
    #   int id, k;

      res=0.0;
      for k=1:7
          res = res +  F_in[i,j,k];
      end
      return res;
end

function calc_u(i,j,F_in,F_out)
  # int id, k;
  # double res,den, v1;

  # id = j*x_dim+i;

  res=0.0;
  den=calc_ro(i,j,F_in,F_out);

  for k=2:7
	  v1 = cx[k];
	  res =  res + F_in[i,j,k]*v1;
  end
  if den != 0.0
      res=res/den;
    else
        res=0.0;
    end

  return res;
end


function calc_v(i,j,F_in,F_out)

    # int id, k;
    #   double res,den, v2;

    #   id = j*x_dim+i;

      res=0.0;
      den = calc_ro(i,j,F_in,F_out);
      for k=2:7
    	  v2 = cy[k];
    	  res = res +  F_in[i,j,k]*v2;
      end
      if den != 0.0
          res=res/den;
      else
          res=0.0;
      end

      return res;

end


function calc_S(i,j,F_in,F_out)
  # double res, m;
  # int id, k;

  # id = j*x_dim+i;

  res=0.0;
  for k=:7
	  m = F_in[i,j,k];
	  res = res + log(m);
  end


  res = res/6.0;

  return res;
end



function calc_P(i,j,F_in,F_out)
  # double res, ro, u, v, u__2, v__2;

  ro = calc_ro(i,j,F_in,F_out);
  u = calc_u(i,j,F_in,F_out);
  v = calc_v(i,j,F_in,F_out);

  v__2 = V*V;
  u__2 = u*u + v*v;

  res = 0.5*ro*(v__2 - u__2 - gama);


  return res;
end


function ord_col(i,j,F_in,F_out)

  # int id, k;
  # double ro, u1, u2, u__2, tau;

  # EQU equ_dist;

  # id = j*x_dim+i;

  #-- calculating the local properties ---*/
  ro=calc_ro(i,j,F_in,F_out);
  u1=calc_u(i,j,F_in,F_out);
  u2=calc_v(i,j,F_in,F_out);
  u__2 = u1*u1 + u2*u2;

  tau = 4.0*visc/(V*V - u__2) + 0.5;

  equ_dist = get_equ(u1, u2, ro);

  for k=:7
    F_out[i,j,k]  = F_in[i,j,k]-(F_in[i,j,k] - equ_dist[k])/tau;
end

end


function bdr_col(i,j,F_in,F_out)
   # int id = j*x_dim+i;

   F_out[i,j,2] = F_in[i,j,5];
   F_out[i,j,3] = F_in[i,j,6];
   F_out[i,j,4] = F_in[i,j,7];
   F_out[i,j,5] = F_in[i,j,2];
   F_out[i,j,6] = F_in[i,j,3];
   F_out[i,j,7] = F_in[i,j,4];
   F_out[i,j,1] = F_in[i,j,1];

end



function stream(i,j,F_in,F_out)
  # int id, nid, k,ni,nj;
  # id = j*x_dim + i;

  for k=1:7
    ni = i + ex[k];
    nj = j + ey[k];
	nid = nj*x_dim + ni;
	if ni>0 && ni<=x_dim && nj>0 && nj<=y_dim
		F_in[ni,nj,k] = F_out[i,j,k];
	end
end

end


function rec_si(iter,F_in,F_out,si)
  # int id, i,j;
  # double dx, dy;

  dx = 0.5;
  dy = 0.5*sqrt(3.0);

  # /*--------- calculating the si field ----------*/
  for i=1:x_dim
      si[i,1]=0.0;
      si[i,y_dim]=0.0;
  end
  for j=1:y_dim
      si[1,j]=0.0;
      si[x_dim,j]=0.0;
  end

  for i=2:(x_dim-1)
    for j=(y_dim-1):2
	  si[i,j] = si[i-1,j] - (calc_v(i-1,j,F_in,F_out)*dx - calc_v(i-1,j,F_in,F_out)*dy);
  end
  end


  for i=1:x_dim
    for j=y_dim:1
		# id = j*x_dim +i;
	  if bdr_state[i,j] == 1
          si[i,j] = 0.0;
      end
    end
  end

end



function rec_params(iter,F_in,F_out)

	# char fileNameString[32]; # creating file name buffer.
	@static if  MakeMultipleFile == "Yes"
	# put "file" name in the string  it makes file name in the string each time.
	fileNameString = @sprintf("v_d2z6_%lf.dat",iter/NAMEFRACTION);
	else
	fileNameString = @sprintf("v_d2z6.dat");
	end

  # FILE *out;

  # int i,j, ss;
  # double xx,yy,u,v;
   L_xx = (x_dim+0.0) - 0.5*(y_dim+0.0);
   L_yy = sin(pi/3.0)*(y_dim+0.0);

  # /*------  v/U values at horizontal centerline --------*/
  out=open(fileNameString,"w");
  @printf(out,"title=\"d2z6      Re=%2.0f     iterations=%d \"\n ",Re,iter);
  @printf(out,"zone t=\"d2z6\"");
  @printf(out,"variables=X,v/U\n");

  j = floor(Int64,y_dim/2);
  for i=1:x_dim
    xx = ((i + 0.0) + 0.5*(j+0.0) - 0.5*(y_dim+0.0))/L_xx;
	if xx >= 0.0 && xx <= 1.0
	   	v = calc_v(i,floor(Int64,y_dim/2),F_in,F_out)/U;
        @printf(out,"%f %f\n",xx,v);
	end
end

  close(out);

  # char fileNameString2[32]; # creating file name buffer.
	@static if  MakeMultipleFile == "Yes"
	# put "file" name in the string  it makes file name in the string each time.
	fileNameString2 = @sprintf("u_d2z6_%lf.dat",iter/NAMEFRACTION);
	else
	fileNameString2 = @sprintf("u_d2z6.dat");
	end

  # /*------  u/U values at vertical centerline --------*/
  out=open(fileNameString2,"w");
  @printf(out,"title=\"d2z6      Re=%2.0f     iterations=%d \"\n ",Re,iter);
  @printf(out,"zone t=\"d2z6\"");
  @printf(out,"variables=u/U,Y\n");
  for j=1:y_dim-1
     yy =  (j+0.0)*sin(pi/3.0)/L_yy;
	 ss = 0;
     for i=1:x_dim
        if ss==0
            xx = ((i + 0.0) + 0.5*(j+0.0) - 0.5*y_dim)/L_xx;
            if xx>0.5
	               ss = 1;
                u=calc_u(i,j,F_in,F_out)/U;
                @printf(out,"%f %f\n",u,yy);
            end
        end
     end
  end
  close(out);
end


function rec_field(iter,F_in,F_out,si)
	# char fileNameString3[32]; # creating file name buffer.
	@static if MakeMultipleFile == "Yes"
	# put "file" name in the string  it makes file name in the string each time.
	fileNameString3 = @sprintf("Results_%lf.dat",iter/NAMEFRACTION);
	else
	fileNameString3 = @sprintf("Results.dat");
    end


  # FILE *out;
  # int i,j;

  # double xx, yy, ro, u, v, P, S, T, h0;
  # double x0, x1, y0,  sai;

  GRID_SIZE = floor(Int64,x_dim - y_dim/2);

  h0 = 0.5*sqrt(3.0)*(y_dim-1.0);
  T = (iter + 0.0)*U/h0;

  x0 = 0.5*(y_dim-1.0);
  x1 = x_dim - x0;
  y0 = 0.5*sqrt(3.0)*(y_dim-1.0);
  T = (iter + 0.0)*U/y0;

  out=open(fileNameString3,"w");

  @printf(out,"title=\"Re=%2.0f iter=%d  T=%2.2f  GRID=%dx%d\"\n ",
	  Re,iter, T, GRID_SIZE, GRID_SIZE);

  @printf(out,"variables=X, Y, si, u, v, ro, P, S\n");
  @printf(out,"zone t=\"d2z6\",i=%d,j=%d\n",x_dim,y_dim);


  # ---------------- calculate the field variables ------------------
  for j=y_dim:1
    for i=1:x_dim
        sai = si[i,j];

	S = calc_S(i,j,F_in,F_out);
	ro = calc_ro(i,j,F_in,F_out);
	u  = calc_u(i,j,F_in,F_out);
        v  = calc_v(i,j,F_in,F_out);

        xx=(i-0+0.0) + cos(pi/3.0)*(j-0+0.0);
        yy =  sin(pi/3.0)*(j-0+0.0);
        P = calc_P(i, j,F_in,F_out);

	xx = (xx - x0)/x1;
	yy = yy/y0;
	if xx < 0.0 || xx > 1.0 || yy < 0.0 || yy> 1.0

		sai = 0.0;
		P = 0.0;
		u = 0.0;
		v = 0.0;
		ro = 0.0;
		S = 0.0;
	end

	@printf(out,"%f   %f   %f  %1.12f  %1.12f   %1.12f   %1.12f   %1.12f\n",
        xx,  yy,  sai, u/U,  v/U,  ro, P, S);
    end
end
  close(out);

	 @printf(" Files has been saved at (%s) Re(%2.0f) U(%1.4f) G(%d,%d)\n\n",fileNameString3,Re,U,GRID_SIZE,GRID_SIZE);
	 flush(STDOUT);
end

function simulate(bdr_state,F_in,F_out)

  # int i,j,id, b_state;

# #  reset_driving_nodes();

  # # --- colision stage ---
  for i=1:x_dim
    for j=1:y_dim
		# id = j*x_dim + i;
      b_state = bdr_state[i,j];
      if b_state == 0
          ord_col(i,j,F_in,F_out);
      end
      if b_state == 1
          bdr_col(i,j,F_in,F_out);
      end
  end
  end


  reset_driving_nodes(F_out);


  # # --- streaming stage ---
  for i=1:x_dim
    for j=1:y_dim
      stream(i,j,F_in,F_out);
  end
  end

#
# #   # --- streaming stage --- I did comment
# #   for(i=0;i<=x_dim-1;++i){
# #     for(j=0;j<=y_dim-1;++j){
# # 		id = j*x_dim + i;
# #     }
# #   }

end

function fieldError()

    # #ux0[x_dim,y_dim], uy0[x_dim,y_dim],ux[x_dim,y_dim],uy[x_dim,y_dim]
    numberOfNodes = 0;
    # int i,j;
    # double e1, e2, e3, xx, yy, x0, x1, y0,u,v;
    esum = 0.0;
    RMS_Error = 0.0;
    MAX_ABS_U_Error = 0.0;
    MAX_ABS_U_Error_X_position  = 0;
    MAX_ABS_U_Error_Y_position  = 0;
    MAX_ABS_V_Error  = 0.0;
    MAX_ABS_V_Error_X_position  = 0;
    MAX_ABS_V_Error_Y_position  = 0;
    MAX_RS_Error  = 0.0;
    MAX_RS_Error_X_Position  = 0;
    MAX_RS_Error_Y_Position  = 0;
    Bulk_Velocity_ABS_Error  = 0.0;
    Bulk_Velocity_ABS_Error_X_Position  = 0;
    Bulk_Velocity_ABS_Error_Y_Position  = 0;


        x0 = 0.5*(y_dim-1.0);
        x1 = x_dim - x0;
        y0 = 0.5*sqrt(3.0)*(y_dim-1.0);

      for j=y_dim:0
        for i=1:x_dim

            u  = calc_u(i,j,F_in,F_out);
            v  = calc_v(i,j,F_in,F_out);

            xx=(i-0+0.0) + cos(pi/3.0)*(j-0+0.0);
            yy =  sin(pi/3.0)*(j-0+0.0);

                xx = (xx - x0)/x1;
                yy = yy/y0;

                if (xx < 0.0 || xx > 1.0 || yy < 0.0 || yy> 1.0)

                    u = 0.0;
                    v = 0.0;

                else

                    numberOfNodes = numberOfNodes + 1;
                    ux[i,j] = u;
                    uy[i,j] = v;
                    e1 = (ux[i,j]-ux0[i,j])*(ux[i,j]-ux0[i,j]) + (uy[i,j] - uy0[i,j])*(uy[i,j] - uy0[i,j]);
                    e2 = ux[i,j]*ux[i,j] + uy[i,j]*uy[i,j];
                    e3 = ux0[i,j]*ux0[i,j] + uy0[i,j]*uy0[i,j];
                    esum = esum + e1/e2;

                    if (fabs(ux[i,j] - ux0[i,j]) > MAX_ABS_U_Error)

                        MAX_ABS_U_Error = fabs(ux[i,j] - ux0[i,j]);
                        MAX_ABS_U_Error_X_position  = i;
                        MAX_ABS_U_Error_Y_position  = j;
                    end

                    if (fabs(uy[i,j] - uy0[i,j]) > MAX_ABS_V_Error)

                        MAX_ABS_V_Error = fabs(uy[i,j] - uy0[i,j]);
                        MAX_ABS_V_Error_X_position  = i;
                        MAX_ABS_V_Error_Y_position  = j;
                    end

                    if (sqrt(e1/e2) > MAX_RS_Error)

                        MAX_RS_Error = sqrt(e1/e2);
                        MAX_RS_Error_X_Position  = i;
                        MAX_RS_Error_Y_Position  = j;
                    end

                    if (Bulk_Velocity_ABS_Error < fabs(e2 - e3))

                        Bulk_Velocity_ABS_Error = fabs(e2 - e3);
                        Bulk_Velocity_ABS_Error_X_Position  = i;
                        Bulk_Velocity_ABS_Error_Y_Position  = j;
                    end

                    ux0[i,j] = ux[i,j];
                    uy0[i,j] = uy[i,j];

                end # end of else

                RMS_Error = sqrt(esum/numberOfNodes);
        end
    end
    return RMS_Error,MAX_ABS_U_Error,MAX_ABS_U_Error_X_position,MAX_ABS_U_Error_Y_position,MAX_ABS_V_Error,MAX_ABS_V_Error_X_position,MAX_ABS_V_Error_Y_position,MAX_RS_Error,MAX_RS_Error_X_Position,MAX_RS_Error_Y_Position,Bulk_Velocity_ABS_Error,Bulk_Velocity_ABS_Error_X_Position,Bulk_Velocity_ABS_Error_Y_Position
end

function main()

    iter = 0;
    # int i,j;

    for j=y_dim:1
        for i=1:x_dim
            ux0[i,j] = 0.0;
            uy0[i,j] = 0.0;
        end
    end

    F_in = Array{Float64}(x_dim,y_dim,MAX_DIR);
    F_out = Array{Float64}(x_dim,y_dim,MAX_DIR);
    bdr_state = Array{Int64}(x_dim,y_dim);
    si = Array{Float64}(x_dim,y_dim);


    n_iter_01 = d_iter_01;
    n_iter_02 = d_iter_02;
    n_iter_03 = d_iter_03;


    initialize(bdr_state,F_in,F_out);





    GRID_SIZE = floor(Int64,x_dim - y_dim/2);
    @static if MakeErrorCSVFile == "Yes"
        # char fileNameString4[64]; # creating file name buffer.
        # put "file" name in the string  it makes file name in the string each time.
        fileNameString4 = @sprintf("1_Errors_G(%d)_U(%1.4f)_R(%2.0f).CSV",GRID_SIZE,U,Re);
        # FILE *ErrorCSV;
        ErrorCSV = open(fileNameString4,"w");
	print(ErrorCSV,"Itterate ,Time ,RMS_Error ,MAX_ABS_U_Error ,MAX_ABS_U_Error_position_X ,MAX_ABS_U_Error_position_Y ,MAX_ABS_V_Error ,MAX_ABS_V_Error_position_X, MAX_ABS_V_Error_position_Y ,MAX_RS_Error ,MAX_RS_Error_Position_X ,MAX_RS_Error_Position_Y ,Bulk_Velocity_ABS_Error ,Bulk_Velocity_ABS_Error_Position_X ,Bulk_Velocity_ABS_Error_Position_Y\n");
    end #end static if



    @printf(" CSKM D2Z6 Gr(%d,%d) U(%1.4f) Re(%2.0f) \n", GRID_SIZE,GRID_SIZE,U,Re);

    #while (RMS_Error > DesiredError || iter <= (d_iter_01*3 + 1))
	while(iter < MAX_UnSteady_Itteration)
		# double y0, T;

        y0 = 0.5*sqrt(3.0)*(y_dim-1.0);
        T = (iter + 0.0)*U/y0;

        simulate(bdr_state,F_in,F_out);


        if (iter > n_iter_01)

            RMS_Error,MAX_ABS_U_Error,MAX_ABS_U_Error_X_position,MAX_ABS_U_Error_Y_position,MAX_ABS_V_Error,MAX_ABS_V_Error_X_position,MAX_ABS_V_Error_Y_position,MAX_RS_Error,MAX_RS_Error_X_Position,MAX_RS_Error_Y_Position,Bulk_Velocity_ABS_Error,Bulk_Velocity_ABS_Error_X_Position,Bulk_Velocity_ABS_Error_Y_Position = fieldError();

			@static if  MakeErrorCSVFile == "Yes"
			@printf(ErrorCSV,"%d,%g,%g,%g,%d,%d,%g,%d,%d,%g,%d,%d,%g,%d,%d,\n",iter,T,RMS_Error,MAX_ABS_U_Error,MAX_ABS_U_Error_X_position,MAX_ABS_U_Error_Y_position,MAX_ABS_V_Error,MAX_ABS_V_Error_X_position,MAX_ABS_V_Error_Y_position,MAX_RS_Error,MAX_RS_Error_X_Position,MAX_RS_Error_Y_Position,Bulk_Velocity_ABS_Error,Bulk_Velocity_ABS_Error_X_Position,Bulk_Velocity_ABS_Error_Y_Position);
			flush(ErrorCSV);
        end #end static if
			#@printf(" iter = %d | T=%2.2f \n", iter, T);
			n_iter_01 = n_iter_01+d_iter_01;

        end # end if

		if (iter > n_iter_03)

			# Showing errors in terminal
                            println("======================================================================");
                            @printf(" Itterate = %d \t T=%2.2f \n RMS_Error       = %1.12f \n MAX_ABS_U_Error = %1.12f \t\t @ -> (%d,%d) \n MAX_ABS_V_Error = %1.12f \t\t @ -> (%d,%d) \n MAX_RS_Error    = %1.12f \t\t @ -> (%d,%d) \n Bulk_Velocity_ABS_Error = %1.12f \t @ -> (%d,%d)\n",iter,T,RMS_Error,MAX_ABS_U_Error,MAX_ABS_U_Error_X_position,MAX_ABS_U_Error_Y_position,MAX_ABS_V_Error,MAX_ABS_V_Error_X_position,MAX_ABS_V_Error_Y_position,MAX_RS_Error,MAX_RS_Error_X_Position,MAX_RS_Error_Y_Position,Bulk_Velocity_ABS_Error,Bulk_Velocity_ABS_Error_X_Position,Bulk_Velocity_ABS_Error_Y_Position);
                            println("======================================================================");

                            n_iter_03 = n_iter_03+d_iter_03;
		end # end if


		if (iter > n_iter_02)

            println(" Updating the output data files ...");
            rec_si(iter,F_in,F_out,si);
            rec_params(iter,F_in,F_out);
            rec_field(iter,F_in,F_out,si);
            n_iter_02=n_iter_02+d_iter_02;
        end # end if


	  #if(T >= MAX_T){
    	@static if  StabilityRun =="Yes"
    	  if (((Re> 7500) && (iter > MAX_UnSteady_Itteration)) || (iter > MAX_UnSteady_Itteration))

                println("Reached MAX_UnSteady_Itteration: converged.");
                rec_si(iter,F_in,F_out,si);
                rec_params(iter,F_in,F_out);
                rec_field(iter,F_in,F_out,si);
                flush(STDOUT);
                break;
            end # end if
        end # end static if

	  iter = iter + 1;
      flush(STDOUT);
  end # end while

	@static if  MakeErrorCSVFile == "Yes"
        close(ErrorCSV);
	end # end static if

    println(" CONVERGED by Error criteria.");

end # end function

main()
