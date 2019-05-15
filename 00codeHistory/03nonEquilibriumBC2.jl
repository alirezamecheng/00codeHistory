# in the name of allah
# cleaning workspace
workspace();
#=
Coordinates

Ny
r
|
|
1---------> c Nx

=#



workspace();

const Nx = 99        # number of cells in the x-direction -1
const Ny = 99    # number of cells in the y-direction -1
const Re = 200.0
const errorlimitation = 1.0e-5
const ErrCalculationSteps = 10
const MaxItteration = 1000
#------------------  Printing And I/O Parameters  ------------------------------
const PrintDecimal = 1000000 # number of Decimals for file name.

const PrintOutput = "yes"
# const PrintOutput = "No"
const PrintSteps = 10

# const CSVoutput = "yes"
const CSVoutput = "No"
const CSVoutputFrequency = 1000

const CenterLinePrint = "yes"
# const CenterLinePrint = "No"
const CenterlinePrintFrequency = 20000

# const DistributionFunctions_Out = "yes"
const DistributionFunctions_Out = "No"
const DistributionFunctionsFrequency = 100

const DistributionFunctions_sum_Out = "yes"
# const DistributionFunctions_sum_Out = "No"
const DistributionFunctionsSumFrequency = 10

const ErrorPrintingInInFilePermission = "yes" # error calculation steps is equal with "const ErrCalculationSteps"
# const ErrorPrintingInInFilePermission = "No"

#/==============================================================================
const Nx1 = Nx + 1
const Ny1 = Ny + 1
const L = Ny + 1     # width of the cavity
const Q = 9          # number of discrete velocities
const rho0 = 1.0     # initial density
const ux0 = 0.0      # initial velocity component in x direction
const uy0 = 0.0      # initial velocity component in y direction
#/=========================== BOUNDARY WALLS VELOCITIES ========================
const uw_top    = 0.1 # Lid U. Velocity or top X velocity
const vw_top    = 0.0 # Lid V. Velocity or top Y velocity

const uw_left   = 0.0 # Left wall U
const vw_left   = 0.0 # Left wall V

const uw_right  = 0.0 # Right wall U
const vw_right  = 0.0 # Right wall V

const uw_bottom = 0.0 # Bottom wall U
const vw_bottom = 0.0 # Bottom wall V

#/==============================================================================
const cx = [1, 0, -1, 0, 1, -1, -1, 1, 0]; # lattice velocities changed 0-->9
const cy = [0, 1, 0, -1, 1, 1, -1, -1, 0]; # lattice velocities changed 0-->9
const w  = [1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36, 4.0/9]; # lattice weights changed 0-->9
#/======================= MRT method constants =================================
const D  = [9, 36, 36, 6, 12, 6, 12, 4, 4]; # D = M*M^T
const rc = [0, 3, 4, 1, 2, 7, 8, 5, 6];  # index of reversed velocity
#/==============================================================================

#/==============================================================================

const resultsAddress = "/home/aghnia/Results/SRT_Semi_Non_equilibrum_boundary_$(Dates.format(now(),"YY.mm.dd.HH.MM"))_SRT_Re($(Re))_Gr($(Nx1)_$(Ny1))_Uw($(uw_top))" # for Server
# const resultsAddress = "/home/ghavami/mscThesis/runZone/SRT_Non_equilibrum_boundary_$(Dates.format(now(),"YY.mm.dd.HH.MM"))_SRT_Re($(Re))_Gr($(Nx1)_$(Ny1))_Uw($(uw_top))" # for Local
mkpath(resultsAddress) # to creat specified path



# Calculation the equilibrium distribution
function feq(RHO,U,V,k)  # RHO,U,V,k in this function are "Numbers".
    cu = cx[k]*U + cy[k]*V # c K*u
    U2 = U*U + V*V         # u*u
    return w[k]*RHO*(1.0 + 3.0*cu + 4.5*cu*cu - 1.5*U2)
end # end function [feq]



# Initialization
function Init_Eq!(rho, ux, uy, f)
# چون یک بار اجرا می شود ارزش این که بهینه تعریف شود را ندارد
    for c = 1:Nx1
      for r = 1:Ny1
        rho[r,c] = rho0;
        ux[r,c] = ux0;
        uy[r,c] = uy0;
        for k = 1:Q
            f[r,c,k] = feq(rho[r,c], ux[r,c], uy[r,c], k); # recall feq function
        end #end for K
      end #end for r
    end #end for c
    return rho, ux, uy, f
end #end function Init_Eq.



# BGK collision
function Coll_BGK!(rho, ux, uy, f, f_post, tau)
    #double FEQ
    for k = 1:Q
        for c = 1:Nx1
            for r = 1:Ny1
                FEQ = feq(rho[r,c], ux[r,c], uy[r,c], k) #EDF
                f_post[r,c,k] = f[r,c,k] - (f[r,c,k] - FEQ)/tau; #Post-collision DFs
            end # for r
        end # for c
    end #for k
    return f_post
end #function [Coll_BGK]



# Streaming
function Streaming!(f,f_post)
    for k = 1:Q
        for c = 1:Nx1
            for r = 1:Ny1
                cd = c - cx[k] # upwind node
                rd = r - cy[k]
                if cd>0 && rd>0 && cd<=Nx1 && rd<=Ny1
                    f[r,c,k] = f_post[rd,cd,k] # streaming
                end # if
            end # for r
        end # for c
    end # for k
    return f
end # function [Streaming]



# Streaming For Non_equilibrium_Boundarycondition
function Non_Equilibrium_Streaming!(f,f_post)
    for k = 1:Q
        for c = 2:Nx1-1
            for r = 2:Ny1-1
                cd = c - cx[k] # upwind node
                rd = r - cy[k]
                if cd>0 && rd>0 && cd<=Nx1 && rd<=Ny1
                    f[r,c,k] = f_post[rd,cd,k] # streaming
                end # if
            end # for r
        end # for c
    end # for k
    return f
end # function [Streaming]






# Bounce-back scheme #برای مساله درب متحرک
function Bounce_back!(f,f_post,rho)

     #  C = 1: left wall
     for r = 1:Ny1
       f[r,1,1] = f_post[r,1,3];
       f[r,1,5] = f_post[r,1,7]; # 1 singular node
       f[r,1,8] = f_post[r,1,6]; # 1 singular node
     end

     # c = Nx: right wall
     for r = 1:Ny1
       f[r,Nx1,3] = f_post[r,Nx1,1];
       f[r,Nx1,7] = f_post[r,Nx1,5]; # 1 singular node
       f[r,Nx1,6] = f_post[r,Nx1,8]; # 1 singular node
     end

     # r = Ny top plate
     for c = 1:Nx1
         f[Ny1,c,4] = f_post[Ny1,c,2];
         f[Ny1,c,7] = f_post[Ny1,c,5]+6*rho[Ny1,c]*w[7]*cx[7]*uw_top;
         f[Ny1,c,8] = f_post[Ny1,c,6]+6*rho[Ny1,c]*w[8]*cx[8]*uw_top;
     end

     # r = 1: bottom plate
     for c = 1:Nx1
       f[1,c,2] = f_post[1,c,4];
       f[1,c,5] = f_post[1,c,7];
       f[1,c,6] = f_post[1,c,8];
     end

     return f
end # function [Bounce_back!]



# Non-equilibrium boundary condition
function Non_equilibrium_Boundarycondition!(f,rho,ux,uy)

    # r = Ny top plate
  for c = 1:Nx1
    f[Ny1,c,1] = feq(rho[Ny,c], uw_top, vw_top , 1) + f[Ny,c,1] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 1);
    f[Ny1,c,2] = feq(rho[Ny,c], uw_top, vw_top , 2) + f[Ny,c,2] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 2);
    f[Ny1,c,3] = feq(rho[Ny,c], uw_top, vw_top , 3) + f[Ny,c,3] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 3);
    f[Ny1,c,4] = feq(rho[Ny,c], uw_top, vw_top , 4) + f[Ny,c,4] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 4);
    f[Ny1,c,5] = feq(rho[Ny,c], uw_top, vw_top , 5) + f[Ny,c,5] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 5);
    f[Ny1,c,6] = feq(rho[Ny,c], uw_top, vw_top , 6) + f[Ny,c,6] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 6);
    f[Ny1,c,7] = feq(rho[Ny,c], uw_top, vw_top , 7) + f[Ny,c,7] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 7);
    f[Ny1,c,8] = feq(rho[Ny,c], uw_top, vw_top , 8) + f[Ny,c,8] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 8);
    f[Ny1,c,9] = feq(rho[Ny,c], uw_top, vw_top , 9) + f[Ny,c,9] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 9);
  end



    #   Feq help : feq(rho[r,c], ux[r,c], uy[r,c], k)
    #   C = 1: left wall
  for r = 1:Ny1
   f[r,1,1] = feq(rho[r,2], uw_left, vw_left , 1) + f[r,2,1] - feq(rho[r,2], ux[r,2], uy[r,2], 1);
   f[r,1,2] = feq(rho[r,2], uw_left, vw_left , 2) + f[r,2,2] - feq(rho[r,2], ux[r,2], uy[r,2], 2);
   f[r,1,3] = feq(rho[r,2], uw_left, vw_left , 3) + f[r,2,3] - feq(rho[r,2], ux[r,2], uy[r,2], 3);
   f[r,1,4] = feq(rho[r,2], uw_left, vw_left , 4) + f[r,2,4] - feq(rho[r,2], ux[r,2], uy[r,2], 4);
   f[r,1,5] = feq(rho[r,2], uw_left, vw_left , 5) + f[r,2,5] - feq(rho[r,2], ux[r,2], uy[r,2], 5);
   f[r,1,6] = feq(rho[r,2], uw_left, vw_left , 6) + f[r,2,6] - feq(rho[r,2], ux[r,2], uy[r,2], 6);
   f[r,1,7] = feq(rho[r,2], uw_left, vw_left , 7) + f[r,2,7] - feq(rho[r,2], ux[r,2], uy[r,2], 7);
   f[r,1,8] = feq(rho[r,2], uw_left, vw_left , 8) + f[r,2,8] - feq(rho[r,2], ux[r,2], uy[r,2], 8);
   f[r,1,9] = feq(rho[r,2], uw_left, vw_left , 9) + f[r,2,9] - feq(rho[r,2], ux[r,2], uy[r,2], 9);
  end

    #   c = Nx: right wall
 for r = 1:Ny1
   f[r,Nx1,1] = feq(rho[r,Nx], uw_right, vw_right , 1) + f[r,Nx,1] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 1);
   f[r,Nx1,2] = feq(rho[r,Nx], uw_right, vw_right , 2) + f[r,Nx,2] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 2);
   f[r,Nx1,3] = feq(rho[r,Nx], uw_right, vw_right , 3) + f[r,Nx,3] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 3);
   f[r,Nx1,4] = feq(rho[r,Nx], uw_right, vw_right , 4) + f[r,Nx,4] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 4);
   f[r,Nx1,5] = feq(rho[r,Nx], uw_right, vw_right , 5) + f[r,Nx,5] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 5);
   f[r,Nx1,6] = feq(rho[r,Nx], uw_right, vw_right , 6) + f[r,Nx,6] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 6);
   f[r,Nx1,7] = feq(rho[r,Nx], uw_right, vw_right , 7) + f[r,Nx,7] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 7);
   f[r,Nx1,8] = feq(rho[r,Nx], uw_right, vw_right , 8) + f[r,Nx,8] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 8);
   f[r,Nx1,9] = feq(rho[r,Nx], uw_right, vw_right , 9) + f[r,Nx,9] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 9);
 end



    #  r = 1: bottom plate
 for c = 1:Nx1
   f[1,c,1] = feq(rho[2,c], uw_bottom, vw_bottom , 1) + f[2,c,1] - feq(rho[2,c], ux[2,c], uy[2,c], 1);
   f[1,c,2] = feq(rho[2,c], uw_bottom, vw_bottom , 2) + f[2,c,2] - feq(rho[2,c], ux[2,c], uy[2,c], 2);
   f[1,c,3] = feq(rho[2,c], uw_bottom, vw_bottom , 3) + f[2,c,3] - feq(rho[2,c], ux[2,c], uy[2,c], 3);
   f[1,c,4] = feq(rho[2,c], uw_bottom, vw_bottom , 4) + f[2,c,4] - feq(rho[2,c], ux[2,c], uy[2,c], 4);
   f[1,c,5] = feq(rho[2,c], uw_bottom, vw_bottom , 5) + f[2,c,5] - feq(rho[2,c], ux[2,c], uy[2,c], 5);
   f[1,c,6] = feq(rho[2,c], uw_bottom, vw_bottom , 6) + f[2,c,6] - feq(rho[2,c], ux[2,c], uy[2,c], 6);
   f[1,c,7] = feq(rho[2,c], uw_bottom, vw_bottom , 7) + f[2,c,7] - feq(rho[2,c], ux[2,c], uy[2,c], 7);
   f[1,c,8] = feq(rho[2,c], uw_bottom, vw_bottom , 8) + f[2,c,8] - feq(rho[2,c], ux[2,c], uy[2,c], 8);
   f[1,c,9] = feq(rho[2,c], uw_bottom, vw_bottom , 9) + f[2,c,9] - feq(rho[2,c], ux[2,c], uy[2,c], 9);
 end

 return f
end # function [Bounce_back!]



# Non-equilibrium Modified boundary condition
function Non_equilibrium_Boundarycondition_modified!(f,rho,ux,uy)

  # r = Ny top plate
  for c = 1:Nx1
    f[Ny1,c,1] = feq(rho[Ny,c], uw_top, vw_top , 1) + f[Ny,c,1] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 1);
    f[Ny1,c,2] = feq(rho[Ny,c], uw_top, vw_top , 2) + f[Ny,c,2] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 2);
    f[Ny1,c,3] = feq(rho[Ny,c], uw_top, vw_top , 3) + f[Ny,c,3] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 3);
    f[Ny1,c,4] = feq(rho[Ny,c], uw_top, vw_top , 4) + f[Ny,c,4] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 4);
    f[Ny1,c,5] = feq(rho[Ny,c], uw_top, vw_top , 5) + f[Ny,c,5] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 5);
    f[Ny1,c,6] = feq(rho[Ny,c], uw_top, vw_top , 6) + f[Ny,c,6] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 6);
    f[Ny1,c,7] = feq(rho[Ny,c], uw_top, vw_top , 7) + f[Ny,c,7] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 7);
    f[Ny1,c,8] = feq(rho[Ny,c], uw_top, vw_top , 8) + f[Ny,c,8] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 8);
    # f[Ny1,c,9] = feq(rho[Ny,c], uw_top, vw_top , 9) + f[Ny,c,9] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 9);
    f[Ny1,c,9] = 1 - (f[Ny1,c,1] + f[Ny1,c,2] + f[Ny1,c,3] + f[Ny1,c,4] + f[Ny1,c,5] + f[Ny1,c,6] + f[Ny1,c,7] + f[Ny1,c,8]); 
  end



  # Feq help : feq(rho[r,c], ux[r,c], uy[r,c], k)
  #  C = 1: left wall
  for r = 1:Ny1
   f[r,1,1] = feq(rho[r,2], uw_left, vw_left , 1) + f[r,2,1] - feq(rho[r,2], ux[r,2], uy[r,2], 1);
   f[r,1,2] = feq(rho[r,2], uw_left, vw_left , 2) + f[r,2,2] - feq(rho[r,2], ux[r,2], uy[r,2], 2);
   f[r,1,3] = feq(rho[r,2], uw_left, vw_left , 3) + f[r,2,3] - feq(rho[r,2], ux[r,2], uy[r,2], 3);
   f[r,1,4] = feq(rho[r,2], uw_left, vw_left , 4) + f[r,2,4] - feq(rho[r,2], ux[r,2], uy[r,2], 4);
   f[r,1,5] = feq(rho[r,2], uw_left, vw_left , 5) + f[r,2,5] - feq(rho[r,2], ux[r,2], uy[r,2], 5);
   f[r,1,6] = feq(rho[r,2], uw_left, vw_left , 6) + f[r,2,6] - feq(rho[r,2], ux[r,2], uy[r,2], 6);
   f[r,1,7] = feq(rho[r,2], uw_left, vw_left , 7) + f[r,2,7] - feq(rho[r,2], ux[r,2], uy[r,2], 7);
   f[r,1,8] = feq(rho[r,2], uw_left, vw_left , 8) + f[r,2,8] - feq(rho[r,2], ux[r,2], uy[r,2], 8);
#    f[r,1,9] = feq(rho[r,2], uw_left, vw_left , 9) + f[r,2,9] - feq(rho[r,2], ux[r,2], uy[r,2], 9);
    f[r,1,9] = 1 - (f[r,1,1] + f[r,1,2] + f[r,1,3] + f[r,1,4] + f[r,1,5] + f[r,1,6] + f[r,1,7] + f[r,1,8])
  end

  # c = Nx: right wall
 for r = 1:Ny1
   f[r,Nx1,1] = feq(rho[r,Nx], uw_right, vw_right , 1) + f[r,Nx,1] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 1);
   f[r,Nx1,2] = feq(rho[r,Nx], uw_right, vw_right , 2) + f[r,Nx,2] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 2);
   f[r,Nx1,3] = feq(rho[r,Nx], uw_right, vw_right , 3) + f[r,Nx,3] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 3);
   f[r,Nx1,4] = feq(rho[r,Nx], uw_right, vw_right , 4) + f[r,Nx,4] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 4);
   f[r,Nx1,5] = feq(rho[r,Nx], uw_right, vw_right , 5) + f[r,Nx,5] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 5);
   f[r,Nx1,6] = feq(rho[r,Nx], uw_right, vw_right , 6) + f[r,Nx,6] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 6);
   f[r,Nx1,7] = feq(rho[r,Nx], uw_right, vw_right , 7) + f[r,Nx,7] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 7);
   f[r,Nx1,8] = feq(rho[r,Nx], uw_right, vw_right , 8) + f[r,Nx,8] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 8);
#    f[r,Nx1,9] = feq(rho[r,Nx], uw_right, vw_right , 9) + f[r,Nx,9] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 9);
    f[r,Nx1,9] = 1 - (f[r,Nx1,1] + f[r,Nx1,2] + f[r,Nx1,3] + f[r,Nx1,4] + f[r,Nx1,5] + f[r,Nx1,6] + f[r,Nx1,7] + f[r,Nx1,8])
 end



 # r = 1: bottom plate
 for c = 1:Nx1
   f[1,c,1] = feq(rho[2,c], uw_bottom, vw_bottom , 1) + f[2,c,1] - feq(rho[2,c], ux[2,c], uy[2,c], 1);
   f[1,c,2] = feq(rho[2,c], uw_bottom, vw_bottom , 2) + f[2,c,2] - feq(rho[2,c], ux[2,c], uy[2,c], 2);
   f[1,c,3] = feq(rho[2,c], uw_bottom, vw_bottom , 3) + f[2,c,3] - feq(rho[2,c], ux[2,c], uy[2,c], 3);
   f[1,c,4] = feq(rho[2,c], uw_bottom, vw_bottom , 4) + f[2,c,4] - feq(rho[2,c], ux[2,c], uy[2,c], 4);
   f[1,c,5] = feq(rho[2,c], uw_bottom, vw_bottom , 5) + f[2,c,5] - feq(rho[2,c], ux[2,c], uy[2,c], 5);
   f[1,c,6] = feq(rho[2,c], uw_bottom, vw_bottom , 6) + f[2,c,6] - feq(rho[2,c], ux[2,c], uy[2,c], 6);
   f[1,c,7] = feq(rho[2,c], uw_bottom, vw_bottom , 7) + f[2,c,7] - feq(rho[2,c], ux[2,c], uy[2,c], 7);
   f[1,c,8] = feq(rho[2,c], uw_bottom, vw_bottom , 8) + f[2,c,8] - feq(rho[2,c], ux[2,c], uy[2,c], 8);
#    f[1,c,9] = feq(rho[2,c], uw_bottom, vw_bottom , 9) + f[2,c,9] - feq(rho[2,c], ux[2,c], uy[2,c], 9);
    f[1,c,9] = 1 - (f[1,c,1] + f[1,c,2] + f[1,c,3] + f[1,c,4] + f[1,c,5] + f[1,c,6] + f[1,c,7] + f[1,c,8])
 end

 return f
end # function [Bounce_back!]




# Semi_Non-equilibrium_boundary_condition
function Semi_Non_equilibrium_Boundarycondition!(f,rho,ux,uy)

  # Feq help : feq(rho[r,c], ux[r,c], uy[r,c], k)
  #  C = 1: left wall
  for r = 1:Ny1
   f[r,1,1] = feq(rho[r,2], uw_left, vw_left , 1) + f[r,2,1] - feq(rho[r,2], ux[r,2], uy[r,2], 1);
   f[r,1,5] = feq(rho[r,2], uw_left, vw_left , 5) + f[r,2,5] - feq(rho[r,2], ux[r,2], uy[r,2], 5);
   f[r,1,8] = feq(rho[r,2], uw_left, vw_left , 8) + f[r,2,8] - feq(rho[r,2], ux[r,2], uy[r,2], 8);
  end

  # c = Nx: right wall
 for r = 1:Ny1
   f[r,Nx1,3] = feq(rho[r,Nx], uw_right, vw_right , 3) + f[r,Nx,3] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 3);
   f[r,Nx1,6] = feq(rho[r,Nx], uw_right, vw_right , 6) + f[r,Nx,6] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 6);
   f[r,Nx1,7] = feq(rho[r,Nx], uw_right, vw_right , 7) + f[r,Nx,7] - feq(rho[r,Nx], ux[r,Nx], uy[r,Nx], 7);
 end

  # r = Ny top plate
  for c = 1:Nx1
    f[Ny1,c,4] = feq(rho[Ny,c], uw_top, vw_top , 4) + f[Ny,c,4] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 4);
    f[Ny1,c,7] = feq(rho[Ny,c], uw_top, vw_top , 7) + f[Ny,c,7] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 7);
    f[Ny1,c,8] = feq(rho[Ny,c], uw_top, vw_top , 8) + f[Ny,c,8] - feq(rho[Ny,c], ux[Ny,c], uy[Ny,c], 8);
  end

 # r = 1: bottom plate
 for c = 1:Nx1
   f[1,c,2] = feq(rho[2,c], uw_bottom, vw_bottom , 2) + f[2,c,2] - feq(rho[2,c], ux[2,c], uy[2,c], 2);
   f[1,c,5] = feq(rho[2,c], uw_bottom, vw_bottom , 5) + f[2,c,5] - feq(rho[2,c], ux[2,c], uy[2,c], 5);
   f[1,c,6] = feq(rho[2,c], uw_bottom, vw_bottom , 6) + f[2,c,6] - feq(rho[2,c], ux[2,c], uy[2,c], 6);
 end

 return f
end # function [Bounce_back!]




# Fluid variables (density and velocity)
function Den_vel!(rho,ux,uy,f)
    for c = 1:Nx1
        for r = 1:Ny1
            rho[r,c]=f[r,c,9]+f[r,c,1]+f[r,c,2]+f[r,c,3]+f[r,c,4]+f[r,c,5]+f[r,c,6]+f[r,c,7]+f[r,c,8];
            ux[r,c]=(f[r,c,1]+f[r,c,5]+f[r,c,8]-f[r,c,3]-f[r,c,6]-f[r,c,7])/rho[r,c];
            uy[r,c]=(f[r,c,5]+f[r,c,6]+f[r,c,2]-f[r,c,7]-f[r,c,8]-f[r,c,4])/rho[r,c];
        end # for r
    end # for c
    return rho , ux , uy
end # function  [Den_vel!]



# function Err!(ux, uy, u0, v0)
#      e1 = 0.0;
#      e2 = 0.0;
#     for c = 2:Nx
#         for r = 2:Ny # حدود محاسبه نود هایی هست که روی آنها عملیات مرزی انجام نشده است
#             e1 = e1 + sqrt((ux[r,c]-u0[r,c])^2 + (uy[r,c]-v0[r,c])^2);
#             e2 = e2 + sqrt(ux[r,c]*ux[r,c] + uy[r,c]*uy[r,c]);
#             u0[r,c] = ux[r,c];
#             v0[r,c] = uy[r,c];
#         end # for r
#     end # for c
#     return e1/e2 ,u0 , v0
# end # function [Err]


# Error RMS
function Err!(ux, uy, u0, v0)
     e1 = 0.0;
     e2 = 0.0;
     s  = 0.0;
    for c = 2:Nx
        for r = 2:Ny # حدود محاسبه نود هایی هست که روی آنها عملیات مرزی انجام نشده است
            e1 = (ux[r,c]-u0[r,c])^2 + (uy[r,c]-v0[r,c])^2;
            e2 = (ux[r,c]        )^2 + (uy[r,c]        )^2;
            s = s + e1/e2;
            u0[r,c] = ux[r,c];
            v0[r,c] = uy[r,c];
        end # for r
    end # for c
    return sqrt(s/((Nx-1)*(Ny-1))) ,u0 , v0
end # function [Err]



# function Err2!(ux, uy, u0, v0)
#   ERR = abs(ux - u0)
#   err = maximum(ERR)
#   u0 = deepcopy(ux)
#   v0 = deepcopy(uy)
#   return err , u0 , v0
# end

function Data_output(LoopItterateCounter,rho,ux,uy)
    FileNameAndAddress = @sprintf("%s/Results%lf.plt",resultsAddress,LoopItterateCounter/PrintDecimal)
    File = open(FileNameAndAddress,"w")
      write(File,"TITLE=\"SRT Re = $(Re)\"\n")
      # write(File,"SOLUTIONTIME = $(LoopItterateCounter)\n")
      # write(File,"T = \"Time=$(LoopItterateCounter)\"\n")
      write(File,"VARIABLES = X,Y,Rho,U,V\n")
      write(File,"ZONE\n")
      write(File,"I = $(Nx1)\t J = $(Ny1)\n")
      write(File,"F=POINT\n")
      for c = 1:Nx1
          for r = 1:Ny1
              write(File,@sprintf("%d,%d,%lf,%lf,%lf\n",c,r,rho[r,c],ux[r,c],uy[r,c]))
          end
      end
    close(File)
    println("Results saved at : ",FileNameAndAddress)
flush(STDOUT);
end

function Data_output_CSV(LoopItterateCounter,rho,ux,uy)
    UxFileName = @sprintf("%s/Ux%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal)
    File = open(UxFileName,"w")
    writecsv(File,ux)
    close(File)

    VyFileName = @sprintf("%s/Vy%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal)
    File = open(VyFileName,"w")
    writecsv(File,uy)
    close(File)

    RhoFileName = @sprintf("%s/Rho%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal)
    File = open(RhoFileName,"w")
    writecsv(File,rho)
    close(File)
    println("Ux , Vy and Rho has been saved in \"CSV\" file format.")
flush(STDOUT);
end #function[Data_output_CSV]

function CenterlineOutput(LoopItterateCounter,rho,ux,uy)

    Xcenter = floor(Int32,Nx1/2);
    Ycenter = floor(Int32,Ny1/2);

    UxFileName = @sprintf("%s/CenterUx%lf.plt",resultsAddress,LoopItterateCounter/PrintDecimal);
    File = open(UxFileName,"w");
    writecsv(File,ux[:,Xcenter]);
    close(File);

    VyFileName = @sprintf("%s/CenterVy%lf.plt",resultsAddress,LoopItterateCounter/PrintDecimal);
    File = open(VyFileName,"w");
    writecsv(File,uy[Ycenter,:]);
    close(File);

    println("Center line values saved\n");
flush(STDOUT);

end # end function[centerlineOutput]





function DistributionFunctions_output_CSV(LoopItterateCounter,f)
    Distribution_Functions = @sprintf("%s/D_F[1]_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal) # file name creation
    File = open(Distribution_Functions,"w") # file open
    writecsv(File,f[:,:,1]) # writing to file
    close(File) # closing the file

    Distribution_Functions = @sprintf("%s/D_F[2]_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal) # file name creation
    File = open(Distribution_Functions,"w") # file open
    writecsv(File,f[:,:,2]) # writing to file
    close(File) # closing the file

    Distribution_Functions = @sprintf("%s/D_F[3]_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal) # file name creation
    File = open(Distribution_Functions,"w") # file open
    writecsv(File,f[:,:,3]) # writing to file
    close(File) # closing the file

    Distribution_Functions = @sprintf("%s/D_F[4]_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal) # file name creation
    File = open(Distribution_Functions,"w") # file open
    writecsv(File,f[:,:,4]) # writing to file
    close(File) # closing the file

    Distribution_Functions = @sprintf("%s/D_F[5]_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal) # file name creation
    File = open(Distribution_Functions,"w") # file open
    writecsv(File,f[:,:,5]) # writing to file
    close(File) # closing the file

    Distribution_Functions = @sprintf("%s/D_F[6]_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal) # file name creation
    File = open(Distribution_Functions,"w") # file open
    writecsv(File,f[:,:,6]) # writing to file
    close(File) # closing the file

    Distribution_Functions = @sprintf("%s/D_F[7]_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal) # file name creation
    File = open(Distribution_Functions,"w") # file open
    writecsv(File,f[:,:,7]) # writing to file
    close(File) # closing the file

    Distribution_Functions = @sprintf("%s/D_F[8]_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal) # file name creation
    File = open(Distribution_Functions,"w") # file open
    writecsv(File,f[:,:,8]) # writing to file
    close(File) # closing the file

    Distribution_Functions = @sprintf("%s/D_F[9]_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal) # file name creation
    File = open(Distribution_Functions,"w") # file open
    writecsv(File,f[:,:,9]) # writing to file
    close(File) # closing the file

    println("Distibution functions has been saved in \"CSV\" file format.\n")
    flush(STDOUT);
end #function[Data_output_CSV]




#///////////////////////////////////////////////////////////////////////////////
#
#                              Function main
#
#///////////////////////////////////////////////////////////////////////////////

function main()
    println("Start...\n");
    flush(STDOUT);

    #/================================ Allocation ==========================
    # double tau , err
    f = Array{Float64}(Ny1,Nx1,Q);  # array of the distribution functions (DFs)

    # f_post = similar(f);  # array of the post-collision DFs
    f_post = Array{Float64}(Ny1,Nx1,Q)

    # arrays of fluid density and velocity
    rho = Array{Float64}(Ny1,Nx1);
    ux = Array{Float64}(Ny1,Nx1);
    uy = Array{Float64}(Ny1,Nx1);
    # s = Array{Float64}(Q) # relaxation rates for MRT model

    M2 = floor(Int32,Ny/2);
    N2 = floor(Int32,Nx/2);

    LoopItterateCounter = 0;

    err = 1.0000000001;
    u0 = Array{Float64}(Ny1,Nx1); # for calculation on Err function
    v0 = Array{Float64}(Ny1,Nx1);

    #/================= Constant Variable calculation ==========================

    tau = 3*L*uw_top/Re + 0.5; # relaxation time for BGK


    #/================================ Initialization ==========================
    rho, ux, uy, f = Init_Eq!(rho, ux, uy, f);
    # Debuging print
    #println("rho \n"); show(rho); println("ux \n"); show(ux); println("uy \n"); show(uy); println("f \n"); show(f);

    if DistributionFunctions_sum_Out == "yes"
        Distribution_Functions_sum = @sprintf("%s/A_Distribution_functions_Sum_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal) # file name creation
        Filesum = open(Distribution_Functions_sum,"w") # file open
    end

    if ErrorPrintingInInFilePermission == "yes"
      ErrorFileName = @sprintf("%s/A_Error_File_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal) # file name creation
      ErrorFilePointer = open(ErrorFileName,"w") # file open
    end

    #/==========================================================================
    #/============================= START OF MAIN LOOP =========================
    #/==========================================================================

    # starting of solving main loop
    while (err > errorlimitation) || (LoopItterateCounter <= (ErrCalculationSteps*2 + 1)) # method 1
    # while ((err > errorlimitation) || (LoopItterateCounter <= MaxItteration)) && (LoopItterateCounter <= (ErrCalculationSteps*2 + 1)) # method 2
    # for it = 1:MaxItteration # method 3
        LoopItterateCounter = LoopItterateCounter + 1;


        #/==========================================================================
        #/============================= solving algorithem =========================
        #/==========================================================================

        
        
        # BGK Simple Boundary ======================================================

#        f_post = Coll_BGK!(rho, ux, uy, f, f_post, tau);              # BGK collision
#        f = Streaming!(f,f_post)                                      # Streaming for BGK
#        f = Bounce_back!(f,f_post,rho)                                # No-Slip boundary condition
#        rho, ux, uy = Den_vel!(rho,ux,uy,f)                           # Fluid variables calculation for all models
        
        # BGK Non_equilibrium_Boundarycondition ====================================

#         f_post = Coll_BGK!(rho, ux, uy, f, f_post, tau);              # BGK collision
#         f = Non_Equilibrium_Streaming!(f,f_post)                      # Streaming For Non_equilibrium_Boundarycondition
#         rho, ux, uy = Den_vel!(rho,ux,uy,f)                           # Fluid variables calculation for all models
#         f = Non_equilibrium_Boundarycondition!(f,rho,ux,uy)           # Non_equilibrium_Boundarycondition

        # BGK Non_equilibrium_Modified_Boundarycondition ============================

#         f_post = Coll_BGK!(rho, ux, uy, f, f_post, tau);              # BGK collision
#         f = Non_Equilibrium_Streaming!(f,f_post)                      # Streaming For Non_equilibrium_Boundarycondition
#         rho, ux, uy = Den_vel!(rho,ux,uy,f)                           # Fluid variables calculation for all models
#         f = Non_equilibrium_Boundarycondition_modified!(f,rho,ux,uy)  # Non_equilibrium_Modified_Boundarycondition

        # BGK Semi_Non_equilibrium_Boundarycondition =================================

#           ************ Method is NOT WORKING Code works. ************
#           Non_Equilibrium_Streaming is Better than Simple Streaming. But the boundary at the bottom at which the velocities are Zero Does not have any good respond.

#          f_post = Coll_BGK!(rho, ux, uy, f, f_post, tau);              # BGK collision
# #          f = Streaming!(f,f_post)                                      # Streaming for BGK 
#          f = Non_Equilibrium_Streaming!(f,f_post)                      # Streaming For Non_equilibrium_Boundarycondition
#          rho, ux, uy = Den_vel!(rho,ux,uy,f)                           # Fluid variables calculation for all models
#          f = Semi_Non_equilibrium_Boundarycondition!(f,rho,ux,uy)      # Non_equilibrium_Boundarycondition
#         
        
        
        #/==========================================================================
        #/=========================== Post-processing ==============================
        #/==========================================================================
        if LoopItterateCounter % ErrCalculationSteps == 0
            err ,u0 ,v0 = Err!(ux, uy, u0, v0) # Velocity differences between two successive 1000 steps
            # println("Err1 = ",err) # Error with difference method
            # err ,u0 ,v0 = Err2!(ux, uy, u0, v0)
            println("err = ", err , " ux_center = ", ux[M2,N2] , " uy_center = ", uy[M2,N2] , " itt = ", LoopItterateCounter);  # Display some results
            flush(STDOUT);

            if ErrorPrintingInInFilePermission == "yes"
              write(ErrorFilePointer,"$(LoopItterateCounter)\t$(err)\n") # writing to Filesum
            end # end if[ErrorPrintingInInFilePermission]
            
        end # end if[ErrCalculationSteps]

        if LoopItterateCounter % PrintSteps == 0 && PrintOutput == "yes"
            Data_output(LoopItterateCounter,rho,ux,uy)
        end # end if [PrintSteps]

        if LoopItterateCounter % CSVoutputFrequency == 0 && CSVoutput=="yes"
            Data_output_CSV(LoopItterateCounter,rho,ux,uy)
        end # end if[LoopItterateCounter % CSVoutputFrequency]

        if LoopItterateCounter % CenterlinePrintFrequency == 0 && CenterLinePrint == "yes"
            CenterlineOutput(LoopItterateCounter,rho,ux,uy)
        end # end if[LoopItterateCounter % CenterlinePrintFrequency]

        if LoopItterateCounter % DistributionFunctionsFrequency == 0 && DistributionFunctions_Out == "yes"
            DistributionFunctions_output_CSV(LoopItterateCounter,f)
        end # end of if [LoopItterateCounter % DistributionFunctionsFrequency]

        if LoopItterateCounter % DistributionFunctionsSumFrequency == 0 && DistributionFunctions_sum_Out == "yes"
          sum = 0;
          for i = 1:Nx1
            for j = 1:Ny1
              for k = 1:Q
                sum = sum + f[j,i,k];
              end
            end
          end

          write(Filesum,"$(LoopItterateCounter)\t$(sum)\n") # writing to Filesum
          flush(STDOUT);
        end # end of if [LoopItterateCounter % DistributionFunctionsFrequency]

        # println(LoopItterateCounter) # printing the LoopItterateCounter in the STDOUT






    end # end while [main loop]







    

    if DistributionFunctions_sum_Out == "yes"
      close(Filesum) # closing the file
    end

    if ErrorPrintingInInFilePermission == "yes"
      close(ErrorFilePointer) # closing the file
    end

    Data_output(LoopItterateCounter,rho,ux,uy) # Output last simulation data
    Data_output_CSV(LoopItterateCounter,rho,ux,uy) # Output last simulation data in CSV
    CenterlineOutput(LoopItterateCounter,rho,ux,uy) # output last simulation centerline data in csv
    return rho,ux,uy
end # end function [main]

# rho,ux,uy = main()
@time main()
