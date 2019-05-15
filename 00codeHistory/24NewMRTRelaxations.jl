#=
                    In the name of Allah.

By: Alireza Ghavami nia
=#







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

#/==============================================================================
#                                   Changes!
#/==============================================================================
#=  
    * version 24
        1- New MRT relaxations.

    * version 23 
        1- Back to Old MRT Meq
    
    * version 22 based on version 21
        1- fixing the cascaded central moments space LBM. (God bless Geier for his mistakes) 
        2- Adding Grid const to avoiding mistakes on gird size.
        3- fixing a bug in printing grid size on the start line of the code. 
    
    * version 21 based on version 19
        1- adding new criteria for Reynolds numbers more than 7500;

    * version 19 
        1- equations of m_eq has been fixed. (God bless Guo for his mistakes).
        
    * version 18
        1- adding some other relaxation factors.

    * version 17
        1- adding static if for functions.

    * version 16
        1- fixing viscosity of central moments space.
        2- speed-up and parallelism of central moments space.

    * version 15
        2- semi parallelism with simd macro.
        3- speed-up with eliminationg of bounds cheque with inbounds macro.

    * version 14
        1- Adding () to solving memory allocation bug.

    * version 09 and 10 :
        1- Adding some more Error calculation methods.
        2- Adding compiler directives.

    * version 08, code name: 08 IOCleanUp
        1- I/O flushing bug solved.
        2- I/O string interpertation in write function solved.

    * version 07, code name: 07nonOrthogonalCMs
        1- Implimentation of Centeral moments space method.

    * version 06, code name: 06CleanupAndComment
        1- adding some forgotten ";".
        2- more comments in code.

    * version 04, code name: 04MRTGuo
        1- All Int32 variable types changed to Int64.
        2- MRT code with Guo variables added.
        3- Adding a variable to determine the name of saving folder.

=#




#/==============================================================================
#               Definition of constants and global variables
#/==============================================================================

const Grid = 201
const Nx = Grid - 1       # number of cells in the x-direction -1
const Ny = Grid - 1       # number of cells in the y-direction -1
const Re = 11000.0
const errorlimitation = 1.0e-5;
const ErrCalculationSteps = 1000;
const MAX_UnSteady_Itteration = 500000;
const MaxItteration = 1000;

#/==============================================================================
#                             Choosing Method!
#/==============================================================================
# const MethodName = "BGK";                 # BGK Simple Boundary

# const MethodName = "BGK_NEBC";            # BGK Non_equilibrium_Boundarycondition
# const MethodName = "BGK_MNEBC";           # BGK Non_equilibrium_Modified_Boundarycondition
# const MethodName = "BGK_SNEBC";           # BGK Semi_Non_equilibrium_Boundarycondition

# const MethodName = "MRT_Guo";             # MRT From LBM01.c Guo book
# const MethodName = "MRT_Lalle"            # MRT From lallemand2000 article
# const MethodName = "Pradipto"             # MRT From Pradipto2017

# const MethodName = "MRT_Khatoon"          # proposed by Khatoonabadi
# const MethodName = "MRT_EAslan";          # MRT From E.AslanI.Taymaz2014



# const MethodName = "CMs_NON_Orthogonal";  # CMs. Non-Orthogonal Central Moments space

# const MethodName = "Cas_Geier"            # cascaded Central moments space
# const MethodName = "Cas_Lalle";           # CascadedLBM From lallemand
# const MethodName = "Cas_Chavez";          # CascadedLBM From Chavez-Madena






# test method
# const MethodName = "MRT2"                 # different algorithm different rho

# Not working Methods.

println("Your Selected Method for solver is :",MethodName);
#/==============================================================================
#                       Printing And I/O Parameters
#/==============================================================================

const PrintDecimal = 1000000; # number of Decimals for file names.
#const STDOUT_Error_Report_Size = "Short";
const STDOUT_Error_Report_Size = "Long";
# const STDOUT_Error_Report_Size = "NON";

const PrintOutput = "yes";                      # Main output function that prints results for tecplot. in plt format.
# const PrintOutput = "No";
const PrintSteps = 100000;                        # onley controls Results file output

# const CSVoutput = "yes";
const CSVoutput = "No";
const CSVoutputFrequency = 1000;

const CenterLinePrint = "yes";
# const CenterLinePrint = "No";
const CenterlinePrintFrequency = 100000;

# const DistributionFunctions_Out = "yes";     # Saves populations of each direction in a seperate CSV file.(Makes too much file)
const DistributionFunctions_Out = "No";
const DistributionFunctionsFrequency = 1000;

# const DistributionFunctions_sum_Out = "yes";   # Makes a file and saves the sum of population to cheque the continuity
const DistributionFunctions_sum_Out = "No";
const DistributionFunctionsSumFrequency = 1000;

const ErrorPrintingInInFilePermission = "yes"; # error calculation steps is equal with "const ErrCalculationSteps"
# const ErrorPrintingInInFilePermission = "No";

#/==============================================================================
#                     Geometry VARIABLES and definition
#/==============================================================================
# WARNING: DO NOT Change ANY line below unless otherwise you are 200% sure what
# you are doing.

const Nx1 = Nx + 1;   # actual number of nodes in directions.
const Ny1 = Ny + 1;
const L = Ny + 1;     # width of the cavity
const Q = 9;          # number of discrete velocities
const rho0 = 1.0;     # initial density
const ux0 = 0.0;      # initial velocity component in x direction
const uy0 = 0.0;      # initial velocity component in y direction

#/==============================================================================
#                        BOUNDARY WALLS VELOCITIES
#/==============================================================================

const uw_top    = 0.1; # Lid U. Velocity or top X velocity
const vw_top    = 0.0; # Lid V. Velocity or top Y velocity

const uw_left   = 0.0; # Left wall U
const vw_left   = 0.0; # Left wall V

const uw_right  = 0.0; # Right wall U
const vw_right  = 0.0; # Right wall V

const uw_bottom = 0.0; # Bottom wall U
const vw_bottom = 0.0; # Bottom wall V

#/==============================================================================
#                             Methodes Constants
#/==============================================================================
#/=========================== public variables =================================

const cx = [1, 0, -1, 0, 1, -1, -1, 1, 0]; # lattice velocities changed 0-->9
const cy = [0, 1, 0, -1, 1, 1, -1, -1, 0]; # lattice velocities changed 0-->9

#/======================= BGK method constants =================================
const w  = [1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36, 4.0/9]; # lattice weights changed 0-->9

#/======================= MRT method constants =================================

const D  = [36, 36, 6, 12, 6, 12, 4, 4, 9]; # D = M*M^T changed 0-->9

#/==============================================================================
#                           Results files Address
#/==============================================================================
# making Results adress path: It will be used for all Output savings.

const resultsAddress = "../Results/$(MethodName)_$(Dates.format(now(),"YY.mm.dd.HH.MM"))_Re($(Re))_Gr($(Nx1)_$(Ny1))_Uw($(uw_top))"; # for using in cluster.

# const resultsAddress = "../Results/$(MethodName)_$(Dates.format(now(),"YY.mm.dd.HH.MM"))_Re($(Re))_Gr($(Nx1)_$(Ny1))_Uw($(uw_top))"; # for using in my own pc.

mkpath(resultsAddress); # to creat specified path


println("Re = ",Re," Ux = ",uw_top," (Nx1 , Ny1) ==> (",Nx1,",",Ny1,")");
flush
#/==============================================================================
#                             Functions definition
#/==============================================================================



# Calculation the equilibrium distribution
function feq(RHO::Float64,U::Float64,V::Float64,k::Int64)  # RHO,U,V,k in this function are "Numbers".
    cu = cx[k]*U + cy[k]*V; # c K*u
    U2 = U*U + V*V;         # u*u
    return w[k]*RHO*(1.0 + 3.0*cu + 4.5*cu*cu - 1.5*U2);
end # end function [feq]


# calculation the equilibrium moment for MRT method
@static if MethodName == "MRT_Guo" || MethodName == "MRT_EAslan" || MethodName == "MRT_Lalle" || MethodName == "Pradipto" ||  MethodName == "MRT_Khatoon"
function meq(RHO::Float64, U::Float64, V::Float64, k::Int64)

#with Rho^2

    if k == 9
        x = RHO;
    elseif k == 1
        x = RHO*(-2+3*(U*U+V*V));
    elseif k == 2
        x = RHO - 3*RHO*(U*U+V*V);
    elseif k == 3
        x = RHO*U;
    elseif k == 4
        x = -RHO*U;
    elseif k == 5
        x = RHO*V;
    elseif k == 6
        x = -RHO*V;
    elseif k == 7
        x = RHO*(U*U-V*V);
    elseif k == 8
        x = RHO*U*V;
    else
        x = 0.0;
    end
  return x;
end # end function[meq]
end # end static if [MRT]

# Initialization
function Init_Eq!(rho::Array{Float64}, ux::Array{Float64}, uy::Array{Float64}, f::Array{Float64})
# چون یک بار اجرا می شود ارزش این که بهینه تعریف شود را ندارد
    @inbounds for c = 1:Nx1
      for r = 1:Ny1
        rho[r,c] = rho0;
        ux[r,c] = ux0;
        uy[r,c] = uy0;
       @simd for k = 1:Q
            f[r,c,k] = feq(rho[r,c], ux[r,c], uy[r,c], k); # recall feq function
        end #end for K
      end #end for r
    end #end for c
    return rho, ux, uy, f;
end #end function Init_Eq.



# BGK collision
@static if MethodName == "BGK" || MethodName == "BGK_SNEBC" || MethodName == "BGK_MNEBC" || MethodName == "BGK_NEBC"
function Coll_BGK!(rho::Array{Float64}, ux::Array{Float64}, uy::Array{Float64}, f::Array{Float64}, f_post::Array{Float64}, tau::Float64)
    #double FEQ
 @inbounds for k = 1:Q
        for c = 1:Nx1
           @simd for r = 1:Ny1
                FEQ = feq(rho[r,c], ux[r,c], uy[r,c], k); #EDF
                f_post[r,c,k] = f[r,c,k] - (f[r,c,k] - FEQ)/tau; #Post-collision DFs
            end # for r
        end # for c
    end #for k
    return f_post;
end #function [Coll_BGK]
end # end static if [BGK]


# MRT collision
@static if MethodName == "MRT_Guo" || MethodName == "MRT_EAslan" || MethodName == "MRT_Lalle" || MethodName == "Pradipto" ||  MethodName == "MRT_Khatoon"
function Coll_MRT!(rho::Array{Float64}, ux::Array{Float64}, uy::Array{Float64}, f::Array{Float64}, f_post::Array{Float64}, m::Array{Float64}, s::Array{Float64})

    MEQ::Float64 = 0.0;

  @inbounds  for c = 1:Nx1
        for r = 1:Ny1

            # Transformation from velocity space to moment space:
            m[9] =     (f[r,c,9] +    f[r,c,1]  + f[r,c,2]) + (f[r,c,3] + f[r,c,4]  +    f[r,c,5]) + (f[r,c,6] + f[r,c,7] + f[r,c,8]);
            m[1] =  -4*f[r,c,9] -    f[r,c,1]  - f[r,c,2] - f[r,c,3] - f[r,c,4]  + 2*(f[r,c,5] + f[r,c,6] + f[r,c,7] + f[r,c,8]);
            m[2] =   4*f[r,c,9] - 2*(f[r,c,1]  + f[r,c,2] + f[r,c,3] + f[r,c,4]) +    f[r,c,5] + f[r,c,6] + f[r,c,7] + f[r,c,8];
            m[3] =     f[r,c,1] -    f[r,c,3]  + f[r,c,5] - f[r,c,6] - f[r,c,7]  +    f[r,c,8];
            m[4] = -2*(f[r,c,1] -    f[r,c,3]) + f[r,c,5] - f[r,c,6] - f[r,c,7]  +    f[r,c,8];
            m[5] =     f[r,c,2] -    f[r,c,4]  + f[r,c,5] + f[r,c,6] - f[r,c,7]  -    f[r,c,8];
            m[6] = -2*(f[r,c,2] -    f[r,c,4]) + f[r,c,5] + f[r,c,6] - f[r,c,7]  -    f[r,c,8];
            m[7] =     f[r,c,1] -    f[r,c,2]  + f[r,c,3] - f[r,c,4];
            m[8] =     f[r,c,5] -    f[r,c,6]  + f[r,c,7] - f[r,c,8];


            # Relaxation in moment space:
         @simd   for k = 1:Q

                MEQ  = meq(rho[r,c], ux[r,c], uy[r,c], k);
                m[k] = m[k] - s[k]*(m[k] - MEQ);            # relaxation
                m[k] = m[k] / D[k];                         # rescaling

            end #end[for k]

            # Transforming back to the velocity space:

            f_post[r,c,9] = m[9] -  4*m[1] + 4* m[2];
            f_post[r,c,1] = m[9] -    m[1] - 2*(m[2] + m[4]) + m[3] + m[7];
            f_post[r,c,2] = m[9] -    m[1] - 2*(m[2] + m[6]) + m[5] - m[7];
            f_post[r,c,3] = m[9] -    m[1] - 2*(m[2] - m[4]) - m[3] + m[7];
            f_post[r,c,4] = m[9] -    m[1] - 2*(m[2] - m[6]) - m[5] - m[7];
            f_post[r,c,5] =(m[9] +  2*m[1])+   (m[2]  + m[3] + m[4])+(m[5] + m[6] + m[8]);
            f_post[r,c,6] = m[9] +  2*m[1] +    m[2]  - m[3] - m[4] + m[5] + m[6] - m[8];
            f_post[r,c,7] = m[9] +  2*m[1] +    m[2]  - m[3] - m[4] - m[5] - m[6] + m[8];
            f_post[r,c,8] = m[9] +  2*m[1] +    m[2]  + m[3] + m[4] - m[5] - m[6] - m[8];

        end #end [for r]
    end #end [for c]

    return f_post;
end # function [Coll_MRT]
end # end static if [MRT]




@static if MethodName == "MRT2"
function Coll_MRT2!(rho::Array{Float64}, ux::Array{Float64}, uy::Array{Float64}, f::Array{Float64}, f_post::Array{Float64})

        Lenght::Int64    = L # = N + 1;              # For halfway Bondary Condition

        tau::Float64  = 3*Lenght*uw_top/Re + 0.5; # relaxation time for BGK


@inbounds for c = 1:Nx1
        @simd for r = 1:Ny1



            den_d = rho[r,c];
            u_d = ux[r,c];

            v_d = uy[r,c];

            fin0_d = f[r,c,9];
            fin1_d = f[r,c,1];
            fin2_d = f[r,c,2];
            fin3_d = f[r,c,3];
            fin4_d = f[r,c,4];
            fin5_d = f[r,c,5];
            fin6_d = f[r,c,6];
            fin7_d = f[r,c,7];
            fin8_d = f[r,c,8];

            M0=(     fin0_d)+(     fin1_d)+(     fin2_d)+(     fin3_d)+(     fin4_d)+(     fin5_d)+(     fin6_d)+(     fin7_d)+(     fin8_d);

            M1=(-4.0*fin0_d)+(-1.0*fin1_d)+(-1.0*fin2_d)+(-1.0*fin3_d)+(-1.0*fin4_d)+( 2.0*fin5_d)+( 2.0*fin6_d)+( 2.0*fin7_d)+( 2.0*fin8_d);

            M2=( 4.0*fin0_d)+(-2.0*fin1_d)+(-2.0*fin2_d)+(-2.0*fin3_d)+(-2.0*fin4_d)+(     fin5_d)+(     fin6_d)+(     fin7_d)+(     fin8_d);

            M3=( 0.0*fin0_d)+(     fin1_d)+( 0.0*fin2_d)+(-1.0*fin3_d)+( 0.0*fin4_d)+(     fin5_d)+(-1.0*fin6_d)+(-1.0*fin7_d)+(     fin8_d);

            M4=( 0.0*fin0_d)+(-2.0*fin1_d)+( 0.0*fin2_d)+( 2.0*fin3_d)+( 0.0*fin4_d)+(     fin5_d)+(-1.0*fin6_d)+(-1.0*fin7_d)+(     fin8_d);

            M5=( 0.0*fin0_d)+( 0.0*fin1_d)+(     fin2_d)+( 0.0*fin3_d)+(-1.0*fin4_d)+(     fin5_d)+(     fin6_d)+(-1.0*fin7_d)+(-1.0*fin8_d);

            M6=( 0.0*fin0_d)+( 0.0*fin1_d)+(-2.0*fin2_d)+( 0.0*fin3_d)+( 2.0*fin4_d)+(     fin5_d)+(     fin6_d)+(-1.0*fin7_d)+(-1.0*fin8_d);

            M7=( 0.0*fin0_d)+(     fin1_d)+(-1.0*fin2_d)+(     fin3_d)+(-1.0*fin4_d)+( 0.0*fin5_d)+( 0.0*fin6_d)+( 0.0*fin7_d)+( 0.0*fin8_d);

            M8=( 0.0*fin0_d)+( 0.0*fin1_d)+( 0.0*fin2_d)+( 0.0*fin3_d)+( 0.0*fin4_d)+(     fin5_d)+(-1.0*fin6_d)+(     fin7_d)+(-1.0*fin8_d);




            Meq0=den_d;

            Meq1=-2.0*den_d+den_d*3.0*(((u_d)*(u_d))+(v_d*v_d));

            Meq2=den_d-den_d*3.0*(((u_d)*(u_d))+(v_d*v_d));

            Meq3=den_d*(u_d);

            Meq4=-den_d*(u_d);

            Meq5=den_d*v_d;

            Meq6=-den_d*v_d;

            Meq7=den_d*(((u_d)*(u_d))-(v_d*v_d));

            Meq8=den_d*((u_d)*(v_d));


            # 

            # 	//	S0=1.0;

            # 	//	S1=1.64;

            # 	//	S2=1.54;

            # 	//	S3=1.0;

            # 	//	S4=1.9;

            # 	//	S5=1.0;

            # 	//	S6=1.9;

            # 	//	S7=1.0/TAU;

            # 	//	S8=1.0/TAU;


            S0=0.0;

            S1=1.6;

            S2=1.8;

            S3=0.0;

            S4=8.0*(2.0-(1.0/tau))/(8.0-(1.0/tau));

            S5=0.0;

            S6=8.0*(2.0-(1.0/tau))/(8.0-(1.0/tau));

            S7=1.0/tau;

            S8=1.0/tau;


            # 		//S0=1.0;

            # 		//S1=1.4;

            # 		//S2=1.4;

            # 		//S3=1.0;

            # 		//S4=1.2;

            # 		//S5=1.0;

            # 		//S6=1.2;

            # 		//S7=2.0/(3.0*tau-0.5);

            # 		//S8=2.0/(3.0*tau-0.5);

            # 

            # 		//S0=1.0;

            # 		//S1=2.0/(3.0*tau-0.5);

            # 		//S2=2.0/(3.0*tau-0.5);

            # 		//S3=1.0;

            # 		//S4=8.0*(2.0-(1.0/tau))/(8.0-(1.0/tau));

            # 		//S5=1.0;

            # 		//S6=8.0*(2.0-(1.0/tau))/(8.0-(1.0/tau));

            # 		//S7=2.0/(3.0*tau-0.5);

            # 		//S8=2.0/(3.0*tau-0.5);


            M0_d=M0+(S0*(Meq0-M0));

            M1_d=M1+(S1*(Meq1-M1));

            M2_d=M2+(S2*(Meq2-M2));

            M3_d=M3+(S3*(Meq3-M3));

            M4_d=M4+(S4*(Meq4-M4));

            M5_d=M5+(S5*(Meq5-M5));

            M6_d=M6+(S6*(Meq6-M6));

            M7_d=M7+(S7*(Meq7-M7));

            M8_d=M8+(S8*(Meq8-M8));


                                        
            f_post[r,c,9] =((1.0/9.0)*M0_d)+((-4.0/36.0)*M1_d)+(( 4.0/36.0)*M2_d)+(( 0.0/6.0)*M3_d)+(( 0.0/12.0)*M4_d)+(( 0.0/6.0)*M5_d)+(( 0.0/12.0)*M6_d)+(( 0.0/4.0)*M7_d)+(( 0.0/4.0)*M8_d);

            f_post[r,c,1] =((1.0/9.0)*M0_d)+((-1.0/36.0)*M1_d)+((-2.0/36.0)*M2_d)+(( 1.0/6.0)*M3_d)+((-2.0/12.0)*M4_d)+(( 0.0/6.0)*M5_d)+(( 0.0/12.0)*M6_d)+(( 1.0/4.0)*M7_d)+(( 0.0/4.0)*M8_d);

            f_post[r,c,2] =((1.0/9.0)*M0_d)+((-1.0/36.0)*M1_d)+((-2.0/36.0)*M2_d)+(( 0.0/6.0)*M3_d)+(( 0.0/12.0)*M4_d)+(( 1.0/6.0)*M5_d)+((-2.0/12.0)*M6_d)+((-1.0/4.0)*M7_d)+(( 0.0/4.0)*M8_d);

            f_post[r,c,3] =((1.0/9.0)*M0_d)+((-1.0/36.0)*M1_d)+((-2.0/36.0)*M2_d)+((-1.0/6.0)*M3_d)+(( 2.0/12.0)*M4_d)+(( 0.0/6.0)*M5_d)+(( 0.0/12.0)*M6_d)+(( 1.0/4.0)*M7_d)+(( 0.0/4.0)*M8_d);				

            f_post[r,c,4] =((1.0/9.0)*M0_d)+((-1.0/36.0)*M1_d)+((-2.0/36.0)*M2_d)+(( 0.0/6.0)*M3_d)+(( 0.0/12.0)*M4_d)+((-1.0/6.0)*M5_d)+(( 2.0/12.0)*M6_d)+((-1.0/4.0)*M7_d)+(( 0.0/4.0)*M8_d);				

            f_post[r,c,5] =((1.0/9.0)*M0_d)+(( 2.0/36.0)*M1_d)+(( 1.0/36.0)*M2_d)+(( 1.0/6.0)*M3_d)+(( 1.0/12.0)*M4_d)+(( 1.0/6.0)*M5_d)+(( 1.0/12.0)*M6_d)+(( 0.0/4.0)*M7_d)+(( 1.0/4.0)*M8_d);				

            f_post[r,c,6] =((1.0/9.0)*M0_d)+(( 2.0/36.0)*M1_d)+(( 1.0/36.0)*M2_d)+((-1.0/6.0)*M3_d)+((-1.0/12.0)*M4_d)+(( 1.0/6.0)*M5_d)+(( 1.0/12.0)*M6_d)+(( 0.0/4.0)*M7_d)+((-1.0/4.0)*M8_d);				

            f_post[r,c,7] =((1.0/9.0)*M0_d)+(( 2.0/36.0)*M1_d)+(( 1.0/36.0)*M2_d)+((-1.0/6.0)*M3_d)+((-1.0/12.0)*M4_d)+((-1.0/6.0)*M5_d)+((-1.0/12.0)*M6_d)+(( 0.0/4.0)*M7_d)+(( 1.0/4.0)*M8_d);				

            f_post[r,c,8] =((1.0/9.0)*M0_d)+(( 2.0/36.0)*M1_d)+(( 1.0/36.0)*M2_d)+(( 1.0/6.0)*M3_d)+(( 1.0/12.0)*M4_d)+((-1.0/6.0)*M5_d)+((-1.0/12.0)*M6_d)+(( 0.0/4.0)*M7_d)+((-1.0/4.0)*M8_d);	

        end #end [for r]
    end #end [for c]

    return f_post;
end # function [Coll_MRT2]
end # end static if [MRT2]




# function Non-Orthogonal Central Moments Space Collision
@static if MethodName == "CMs_NON_Orthogonal"
function Coll_CMs_Non_Orthogonal!(rho::Array{Float64}, ux::Array{Float64}, uy::Array{Float64}, f::Array{Float64}, f_post::Array{Float64}, tau::Float64)

    W0::Float64 = 1.0;
    W1::Float64 = 1.0;
    W2::Float64 = 1.0;
    W3::Float64 = 1.0;
    W4::Float64 = 1.0/tau;
    W5::Float64 = 1.0/tau;
    W6::Float64 = 1.0;
    W7::Float64 = 1.0;
    W8::Float64 = 1.0;


   @inbounds for c = 1:Nx1
       @simd for r = 1:Ny1
            cx0=cx[9]-ux[r,c];
            cx1=cx[1]-ux[r,c];
            cx2=cx[2]-ux[r,c];
            cx3=cx[3]-ux[r,c];
            cx4=cx[4]-ux[r,c];
            cx5=cx[5]-ux[r,c];
            cx6=cx[6]-ux[r,c];
            cx7=cx[7]-ux[r,c];
            cx8=cx[8]-ux[r,c];

            cy0=cy[9]-uy[r,c];
            cy1=cy[1]-uy[r,c];
            cy2=cy[2]-uy[r,c];
            cy3=cy[3]-uy[r,c];
            cy4=cy[4]-uy[r,c];
            cy5=cy[5]-uy[r,c];
            cy6=cy[6]-uy[r,c];
            cy7=cy[7]-uy[r,c];
            cy8=cy[8]-uy[r,c];

            K0 = (f[r,c,9] + f[r,c,1] + f[r,c,2]) + (f[r,c,3] + f[r,c,4] + f[r,c,5]) + (f[r,c,6] + f[r,c,7] + f[r,c,8]);

            K1 = ((cx0)*f[r,c,9] + (cx1)*f[r,c,1] + (cx2)*f[r,c,2]) + ((cx3)*f[r,c,3] + (cx4)*f[r,c,4] + (cx5)*f[r,c,5]) + ((cx6)*f[r,c,6] + (cx7)*f[r,c,7] + (cx8)*f[r,c,8]);

            K2 = ((cy0)*f[r,c,9] + (cy1)*f[r,c,1] + (cy2)*f[r,c,2]) + ((cy3)*f[r,c,3] + (cy4)*f[r,c,4] + (cy5)*f[r,c,5]) + ((cy6)*f[r,c,6] + (cy7)*f[r,c,7] + (cy8)*f[r,c,8]);

            K3 = ((cx0*cx0+cy0*cy0)*f[r,c,9] + (cx1*cx1+cy1*cy1)*f[r,c,1] + (cx2*cx2+cy2*cy2)*f[r,c,2]) + ((cx3*cx3+cy3*cy3)*f[r,c,3] + (cx4*cx4+cy4*cy4)*f[r,c,4] + (cx5*cx5+cy5*cy5)*f[r,c,5]) + ((cx6*cx6+cy6*cy6)*f[r,c,6] + (cx7*cx7+cy7*cy7)*f[r,c,7] + (cx8*cx8+cy8*cy8)*f[r,c,8]);

            K4 = ((cx0*cx0-cy0*cy0)*f[r,c,9] + (cx1*cx1-cy1*cy1)*f[r,c,1] + (cx2*cx2-cy2*cy2)*f[r,c,2]) + ((cx3*cx3-cy3*cy3)*f[r,c,3] + (cx4*cx4-cy4*cy4)*f[r,c,4] + (cx5*cx5-cy5*cy5)*f[r,c,5]) + ((cx6*cx6-cy6*cy6)*f[r,c,6] + (cx7*cx7-cy7*cy7)*f[r,c,7] + (cx8*cx8-cy8*cy8)*f[r,c,8]);

            K5 = ((cx0*cy0)*f[r,c,9] + (cx1*cy1)*f[r,c,1] + (cx2*cy2)*f[r,c,2]) + ((cx3*cy3)*f[r,c,3] + (cx4*cy4)*f[r,c,4] + (cx5*cy5)*f[r,c,5]) + ((cx6*cy6)*f[r,c,6] + (cx7*cy7)*f[r,c,7] + (cx8*cy8)*f[r,c,8]);

            K6 = ((cx0*cx0*cy0)*f[r,c,9] + (cx1*cx1*cy1)*f[r,c,1] + (cx2*cx2*cy2)*f[r,c,2]) + ((cx3*cx3*cy3)*f[r,c,3] + (cx4*cx4*cy4)*f[r,c,4] + (cx5*cx5*cy5)*f[r,c,5]) + ((cx6*cx6*cy6)*f[r,c,6] + (cx7*cx7*cy7)*f[r,c,7] + (cx8*cx8*cy8)*f[r,c,8]);

            K7 = ((cx0*cy0*cy0)*f[r,c,9] + (cx1*cy1*cy1)*f[r,c,1] + (cx2*cy2*cy2)*f[r,c,2]) + ((cx3*cy3*cy3)*f[r,c,3] + (cx4*cy4*cy4)*f[r,c,4] + (cx5*cy5*cy5)*f[r,c,5]) + ((cx6*cy6*cy6)*f[r,c,6] + (cx7*cy7*cy7)*f[r,c,7] + (cx8*cy8*cy8)*f[r,c,8]);

            K8 = ((cx0*cx0*cy0*cy0)*f[r,c,9] + (cx1*cx1*cy1*cy1)*f[r,c,1] + (cx2*cx2*cy2*cy2)*f[r,c,2] + (cx3*cx3*cy3*cy3)*f[r,c,3]) + ((cx4*cx4*cy4*cy4)*f[r,c,4] + (cx5*cx5*cy5*cy5)*f[r,c,5] + (cx6*cx6*cy6*cy6)*f[r,c,6]) + ((cx7*cx7*cy7*cy7)*f[r,c,7] + (cx8*cx8*cy8*cy8)*f[r,c,8]);



            Keq0=rho[r,c];
            Keq1=0.0;
            Keq2=0.0;
            Keq3=(2.0/3.0)*rho[r,c];
            Keq4=0.0;
            Keq5=0.0;
            Keq6=-rho[r,c]*(ux[r,c]^2)*uy[r,c];
            Keq7=-rho[r,c]*ux[r,c]*(uy[r,c]^2);
            Keq8=(1.0/9.0)*rho[r,c]*(27.0*(ux[r,c]^2)*(uy[r,c]^2)+1.0);


            Kout0=rho[r,c] ;
            Kout1=0.0 ;
            Kout2=0.0 ;
            Kout3=K3+W3*(Keq3-K3) ;
            Kout4=K4+W4*(Keq4-K4) ;
            Kout5=K5+W5*(Keq5-K5) ;
            Kout6=K6+W6*(Keq6-K6) ;
            Kout7=K7+W7*(Keq7-K7) ;
            Kout8=K8+W8*(Keq8-K8) ;


            f_post[r,c,9] =(Kout3/2.0)*(-2.0+(ux[r,c]^2)+(uy[r,c]^2))+(Kout4/2.0)*(-(ux[r,c]^2)+(uy[r,c]^2))+4.0*Kout5*ux[r,c]*uy[r,c]+2.0*Kout6*uy[r,c]+2.0*Kout7*ux[r,c]+Kout8+rho[r,c]*(1.0-(ux[r,c]^2)-(uy[r,c]^2)+(ux[r,c]^2)*(uy[r,c]^2));

            f_post[r,c,1] =(Kout3/4.0)*(1.0-ux[r,c]-(ux[r,c]^2)-(uy[r,c]^2))+(Kout4/4.0)*(1.0+ux[r,c]+(ux[r,c]^2)-(uy[r,c]^2))+Kout5*( -uy[r,c]-2.0*ux[r,c]*uy[r,c])-Kout6*uy[r,c]+(Kout7/2.0)*(-1.0-2.0*ux[r,c])-(Kout8/2.0)+(rho[r,c]/2.0)*(ux[r,c]+(ux[r,c]^2)-(ux[r,c]*uy[r,c]^2)-(ux[r,c]^2)*(uy[r,c]^2));

            f_post[r,c,2] =(Kout3/4.0)*(1.0-uy[r,c]-(ux[r,c]^2)-(uy[r,c]^2))+(Kout4/4.0)*(-1.0-uy[r,c]+(ux[r,c]^2)-(uy[r,c]^2))+Kout5*( -ux[r,c]-2.0*ux[r,c]*uy[r,c])+(Kout6/2.0)*(-1.0-2.0*uy[r,c])-Kout7*ux[r,c]-(Kout8/2.0)+(rho[r,c]/2.0)*(uy[r,c]+(uy[r,c]^2)-(ux[r,c]^2*uy[r,c])-(ux[r,c]^2)*(uy[r,c]^2));

            f_post[r,c,3] =(Kout3/4.0)*(1.0+ux[r,c]-(ux[r,c]^2)-(uy[r,c]^2))+(Kout4/4.0)*(1.0-ux[r,c]+(ux[r,c]^2)-(uy[r,c]^2))+Kout5*( uy[r,c]-2.0*ux[r,c]*uy[r,c])-Kout6*uy[r,c]+(Kout7/2.0)*(1.0-2.0*ux[r,c])-(Kout8/2.0)+(rho[r,c]/2.0)*(-ux[r,c]+(ux[r,c]^2)+(ux[r,c]*uy[r,c]^2)-(ux[r,c]^2)*(uy[r,c]^2));

            f_post[r,c,4] =(Kout3/4.0)*(1.0+uy[r,c]-(ux[r,c]^2)-(uy[r,c]^2))+(Kout4/4.0)*(-1.0+uy[r,c]+(ux[r,c]^2)-(uy[r,c]^2))+Kout5*( ux[r,c]-2.0*ux[r,c]*uy[r,c])+(Kout6/2.0)*(1.0-2.0*uy[r,c])-Kout7*ux[r,c]-(Kout8/2.0)+(rho[r,c]/2.0)*(-uy[r,c]+(uy[r,c]^2)+(ux[r,c]^2*uy[r,c])-(ux[r,c]^2)*(uy[r,c]^2));

            f_post[r,c,5] =(Kout3/8.0)*(ux[r,c]+uy[r,c]+(ux[r,c]^2)+(uy[r,c]^2))+(Kout4/8.0)*(-ux[r,c]+uy[r,c]-(ux[r,c]^2)+(uy[r,c]^2))+(Kout5/4.0)*( 1.0+2.0*ux[r,c]+2.0*uy[r,c]+4.0*ux[r,c]*uy[r,c])+(Kout6/4.0)*(1.0+2.0*uy[r,c])+(Kout7/4.0)*(1.0+2.0*ux[r,c])+(Kout8/4.0)+(rho[r,c]/4.0)*((ux[r,c]*uy[r,c])+(ux[r,c]^2*uy[r,c])+(ux[r,c]*uy[r,c]^2)+(ux[r,c]^2)*(uy[r,c]^2));

            f_post[r,c,6] =(Kout3/8.0)*(-ux[r,c]+uy[r,c]+(ux[r,c]^2)+(uy[r,c]^2))+(Kout4/8.0)*(ux[r,c]+uy[r,c]-(ux[r,c]^2)+(uy[r,c]^2))+(Kout5/4.0)*( -1.0+2.0*ux[r,c]-2.0*uy[r,c]+4.0*ux[r,c]*uy[r,c])+(Kout6/4.0)*(1.0+2.0*uy[r,c])+(Kout7/4.0)*(-1.0+2.0*ux[r,c])+(Kout8/4.0)+(rho[r,c]/4.0)*(-(ux[r,c]*uy[r,c])+(ux[r,c]^2*uy[r,c])-(ux[r,c]*uy[r,c]^2)+(ux[r,c]^2)*(uy[r,c]^2));

            f_post[r,c,7] =(Kout3/8.0)*(-ux[r,c]-uy[r,c]+(ux[r,c]^2)+(uy[r,c]^2))+(Kout4/8.0)*(ux[r,c]-uy[r,c]-(ux[r,c]^2)+(uy[r,c]^2))+(Kout5/4.0)*( 1.0-2.0*ux[r,c]-2.0*uy[r,c]+4.0*ux[r,c]*uy[r,c])+(Kout6/4.0)*(-1.0+2.0*uy[r,c])+(Kout7/4.0)*(-1.0+2.0*ux[r,c])+(Kout8/4.0)+(rho[r,c]/4.0)*((ux[r,c]*uy[r,c])-(ux[r,c]^2*uy[r,c])-(ux[r,c]*uy[r,c]^2)+(ux[r,c]^2)*(uy[r,c]^2));

            f_post[r,c,8] =(Kout3/8.0)*(ux[r,c]-uy[r,c]+(ux[r,c]^2)+(uy[r,c]^2))+(Kout4/8.0)*(-ux[r,c]-uy[r,c]-(ux[r,c]^2)+(uy[r,c]^2))+(Kout5/4.0)*( -1.0-2.0*ux[r,c]+2.0*uy[r,c]+4.0*ux[r,c]*uy[r,c])+(Kout6/4.0)*(-1.0+2.0*uy[r,c])+(Kout7/4.0)*(1.0+2.0*ux[r,c])+(Kout8/4.0)+(rho[r,c]/4.0)*(-(ux[r,c]*uy[r,c])-(ux[r,c]^2*uy[r,c])+(ux[r,c]*uy[r,c]^2)+(ux[r,c]^2)*(uy[r,c]^2));

        end # end for r
    end # end for c

    return f_post;

end # end function [Coll_CMs_Non_Orthogonal]
end # static if [CMs_NON_Orthogonal]


@static if MethodName == "Cas_Geier" || MethodName == "Cas_Lalle" || MethodName == "Cas_Chavez"
    # Cascaded LBM Transfer Function
function Ktransform!(RHO::Float64,

                        Ux::Float64,

                        Uy::Float64,

                        omega_bulk::Float64,

                        omega_shear::Float64,

                        omega6::Float64,

                        omega7::Float64,

                        omega8::Float64,

                        r::Float64,

                        et::Float64,

                        n::Float64,

                        w::Float64,

                        s::Float64,

                        ne::Float64,

                        nw::Float64,

                        sw::Float64,

                        se::Float64)



                        omega4 = omega_shear;

                        omega5 = omega_shear;

                        omega3 = omega_bulk;

                        

                        f9 = r;

                        f1 = et;

                        f2 = n;

                        f3 = w;

                        f4 = s;

                        f5 = ne;

                        f6 = nw;

                        f7 = sw;

                        f8 = se;



                        # k3 = kxxpyy omega_bulk

                        # k4 = kxxmyy omega_shear

                        # k5 = kxy    omega_shear

                        # k6 = kxxy

                        # k7 = kxyy

                        # k8 = kxxyy

# based on Ning2016 Numerical study of the properties of the central moment lattice Boltzmann method



kxx = f1 + f3 + f5 + f6 + f7 + f8;

kyy = f2 + f4 + f5 + f6 + f7 + f8;

kxy = f5 + f7 - f6 - f8;

kxxy = f5 + f6 - f7 - f8;

kxyy = f5 + f8 - f6 - f7;

kxxyy = f5 + f6 + f7 + f8;





g3 = (omega3/12)*((2*RHO/3) + RHO*((Ux^2) + (Uy^2)) - (kxx + kyy));

g4 = (omega4/4)*(RHO*((Ux^2) - (Uy^2)) - (kxx - kyy));

g5 = (omega5/4)*(RHO*Ux*Uy - kxy);

g6 = (omega6/4)*(2*RHO*(Ux^2)*Uy + kxxy - 2*Ux*kxy - Uy*kxx) - 0.5*Uy*(3*g3 + g4) - 2*Ux*g5;

g7 = (omega7/4)*(2*RHO*Ux*(Uy^2) + kxyy - 2*Uy*kxy - Ux*kyy) - 0.5*Ux*(3*g3 - g4) - 2*Uy*g5;

g8 = (omega8/4)*(RHO/9 + 3*RHO*(Ux^2)*(Uy^2) - (kxxyy - 2*Ux*kxyy - 2*Uy*kxxy + (Ux^2)*kyy + (Uy^2)*kxx + 4*Ux*Uy*kxy)) -2*g3 -0.5*(Uy^2)*(3*g3 + g4) - 0.5*(Ux^2)*(3g3 - g4) - 4*Ux*Uy*g5 - 2*Uy*g6 - 2*Ux*g7;





        r  = r  - 4(g3 - g8);

        et = et - g3 + g4 + 2*(g7 - g8);

        n  = n  - g3 - g4 + 2*(g6 - g8); 

        w  = w  - g3 + g4 - 2*(g7 + g8);

        s  = s  - g3 - g4 - 2*(g6 + g8);

        ne = ne + 2*g3 + g5 - g6 - g7 + g8;

        nw = nw + 2*g3 - g5 - g6 + g7 + g8;

        sw = sw + 2*g3 + g5 + g6 + g7 + g8;

        se = se + 2*g3 - g5 + g6 - g7 + g8;

        return r, et, n, w, s, ne, nw, sw, se;



    end # end of function [Ktransform!]
end # end of static if



@static if MethodName == "Cas_Geier" || MethodName == "Cas_Lalle" || MethodName == "Cas_Chavez"
    # cascaded LBM Collision
    function CascadedLBMCollition!(rho::Array{Float64},
                                    ux::Array{Float64},
                                    uy::Array{Float64},
                                    f::Array{Float64},
                                    f_post::Array{Float64},
                                    omega_bulk::Float64,
                                    omega_shear::Float64,
                                    omega6::Float64,
                                    omega7::Float64,
                                    omega8::Float64)


       @inbounds for c = 1:Nx1
           @simd for r = 1:Ny1
                f_post[r,c,9],
                f_post[r,c,1],
                f_post[r,c,2],
                f_post[r,c,3],
                f_post[r,c,4],
                f_post[r,c,5],
                f_post[r,c,6],
                f_post[r,c,7],
                f_post[r,c,8] = Ktransform!(rho[r,c],
                             ux[r,c],
                             uy[r,c],
                             omega_bulk,
                             omega_shear,
                             omega6,
                             omega7,
                             omega8,
                             f[r,c,9],
                             f[r,c,1],
                             f[r,c,2],
                             f[r,c,3],
                             f[r,c,4],
                             f[r,c,5],
                             f[r,c,6],
                             f[r,c,7],
                             f[r,c,8]);
            end
        end
        return f_post;
    end # end of function [cascadedLBMCollition!]

end # end static if


# Streaming
function Streaming!(f::Array{Float64},f_post::Array{Float64})
  @inbounds for k = 1:Q
        for c = 1:Nx1
            for r = 1:Ny1
                cd = c - cx[k] # upwind node
                rd = r - cy[k]
                if cd>0 && rd>0 && cd<=Nx1 && rd<=Ny1
                    f[r,c,k] = f_post[rd,cd,k]; # streaming
                end # if
            end # for r
        end # for c
    end # for k
    return f;
end # function [Streaming]



# Streaming For Non_equilibrium_Boundarycondition
@static if MethodName == "BGK_NEBC"
function Non_Equilibrium_Streaming!(f::Array{Float64},f_post::Array{Float64})
  @inbounds  for k = 1:Q
        for c = 2:Nx1-1
            for r = 2:Ny1-1
                cd = c - cx[k] # upwind node
                rd = r - cy[k]
                if cd>0 && rd>0 && cd<=Nx1 && rd<=Ny1
                    f[r,c,k] = f_post[rd,cd,k]; # streaming
                end # if
            end # for r
        end # for c
    end # for k
    return f;
end # function [Non_Equilibrium_Streaming!]
end # end static if [BGK_NEBC]


# Bounce-back scheme
function Bounce_back!(f::Array{Float64},f_post::Array{Float64},rho::Array{Float64})

     #  C = 1: left wall
    @inbounds for r = 1:Ny1
       f[r,1,1] = f_post[r,1,3];
       f[r,1,5] = f_post[r,1,7]; # 1 singular node
       f[r,1,8] = f_post[r,1,6]; # 1 singular node
     end

     # c = Nx: right wall
    @inbounds for r = 1:Ny1
       f[r,Nx1,3] = f_post[r,Nx1,1];
       f[r,Nx1,7] = f_post[r,Nx1,5]; # 1 singular node
       f[r,Nx1,6] = f_post[r,Nx1,8]; # 1 singular node
     end

     # r = Ny top plate
    @inbounds for c = 1:Nx1
         f[Ny1,c,4] = f_post[Ny1,c,2];
         f[Ny1,c,7] = f_post[Ny1,c,5]+6*rho[Ny1,c]*w[7]*cx[7]*uw_top;
         f[Ny1,c,8] = f_post[Ny1,c,6]+6*rho[Ny1,c]*w[8]*cx[8]*uw_top;
     end

     # r = 1: bottom plate
    @inbounds for c = 1:Nx1
       f[1,c,2] = f_post[1,c,4];
       f[1,c,5] = f_post[1,c,7];
       f[1,c,6] = f_post[1,c,8];
     end

     return f;
end # function [Bounce_back!]



# Non-equilibrium boundary condition
@static if MethodName == "BGK_NEBC"
function Non_equilibrium_Boundarycondition!(f::Array{Float64},rho::Array{Float64},ux::Array{Float64},uy::Array{Float64})

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

 return f;
end # function [Non_equilibrium_Boundarycondition!]
end # static if [BGK_NEBC]


# Non-equilibrium Modified boundary condition
@static if MethodName == "BGK_MNEBC"
function Non_equilibrium_Boundarycondition_modified!(f::Array{Float64},rho::Array{Float64},ux::Array{Float64},uy::Array{Float64})

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
   f[r,1,9] = 1 - (f[r,1,1] + f[r,1,2] + f[r,1,3] + f[r,1,4] + f[r,1,5] + f[r,1,6] + f[r,1,7] + f[r,1,8]);
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
   f[r,Nx1,9] = 1 - (f[r,Nx1,1] + f[r,Nx1,2] + f[r,Nx1,3] + f[r,Nx1,4] + f[r,Nx1,5] + f[r,Nx1,6] + f[r,Nx1,7] + f[r,Nx1,8]);
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
   f[1,c,9] = 1 - (f[1,c,1] + f[1,c,2] + f[1,c,3] + f[1,c,4] + f[1,c,5] + f[1,c,6] + f[1,c,7] + f[1,c,8]);
 end

 return f;
end # function [Non_equilibrium_Boundarycondition_modified!]
end # static if [BGK_MNEBC]



# Semi_Non-equilibrium_boundary_condition
@static if MethodName == "BGK_SNEBC"
function Semi_Non_equilibrium_Boundarycondition!(f::Array{Float64},rho::Array{Float64},ux::Array{Float64},uy::Array{Float64})

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

 return f;
end # function [Bounce_back!]
end # static if [BGK_SNEBC]



# Fluid variables (density and velocity)
function Den_vel!(rho::Array{Float64},ux::Array{Float64},uy::Array{Float64},f::Array{Float64})
  @inbounds for c = 1:Nx1
     @simd for r = 1:Ny1
            rho[r,c]=(f[r,c,9]+f[r,c,1]+f[r,c,2])+(f[r,c,3]+f[r,c,4]+f[r,c,5])+(f[r,c,6]+f[r,c,7]+f[r,c,8]);
            ux[r,c]=(f[r,c,1]+f[r,c,5]+f[r,c,8]-f[r,c,3]-f[r,c,6]-f[r,c,7])/rho[r,c];
            uy[r,c]=(f[r,c,5]+f[r,c,6]+f[r,c,2]-f[r,c,7]-f[r,c,8]-f[r,c,4])/rho[r,c];
        end # for r
    end # for c
    return rho , ux , uy;
end # function  [Den_vel!]



# function Err!(ux::Array{Float64}, uy::Array{Float64}, u0::Array{Float64}, v0::Array{Float64})
#      e1 = 0.0;
#      e2 = 0.0;
#     for c = 2:Nx
#         for r = 2:Ny # حدود محاسبه نود هایی هست که روی آنها عملیات مرزی انجام نشده است
#             e1 = e1 + sqrt((ux[r,c]-u0[r,c])^2 + (uy[r,c]-v0[r,c])^2);
#             e2 = e2 + sqrt(ux[r,c]^2 + uy[r,c]^2);
#             u0[r,c] = ux[r,c];
#             v0[r,c] = uy[r,c];
#         end # for r
#     end # for c
#     return e1/e2 ,u0 , v0;
# end # function [Err]




# Error Methods Function
function ErrorManager!(ux::Array{Float64},
                       uy::Array{Float64},
                       u0::Array{Float64},
                       v0::Array{Float64},
                       RMS_Error::Float64,
                       MAX_ABS_U_Error::Float64,
                       MAX_ABS_U_Error_position::Tuple{Int64,Int64},
                       MAX_ABS_V_Error::Float64,
                       MAX_ABS_V_Error_position::Tuple{Int64,Int64},
                       MAX_RS_Error::Float64,
                       MAX_RS_Error_Position::Tuple{Int64,Int64},
                       Bulk_Velocity_ABS_Error::Float64,
                       Bulk_Velocity_ABS_Error_Position::Tuple{Int64,Int64})

    # ================= Allocations for Error function =========================
    RMS_Error = 0.0;
    MAX_ABS_U_Error = 0.0;
    MAX_ABS_U_Error_position  = (0,0);
    MAX_ABS_V_Error  = 0.0;
    MAX_ABS_V_Error_position  = (0,0);
    MAX_RS_Error  = 0.0;
    MAX_RS_Error_Position  = (0,0);
    Bulk_Velocity_ABS_Error  = 0.0;
    Bulk_Velocity_ABS_Error_Position  = (0,0);
    e1::Float64 = 0.0;
    e2::Float64 = 0.0;
    e3::Float64 = 0.0;
    s::Float64  = 0.0;


   @inbounds for c = 2:Nx
        for r = 2:Ny # حدود محاسبه نود هایی هست که روی آنها عملیات مرزی انجام نشده است

            e1 = (ux[r,c]-u0[r,c])^2 + (uy[r,c]-v0[r,c])^2;
            e2 = (ux[r,c]        )^2 + (uy[r,c]        )^2;
            e3 = (u0[r,c]        )^2 + (v0[r,c]        )^2;
            s = s + e1/e2;

            if abs(ux[r,c] - u0[r,c]) > MAX_ABS_U_Error
                MAX_ABS_U_Error = abs(ux[r,c] - u0[r,c]);
                MAX_ABS_U_Error_position = (r,c);
            end

            if abs(uy[r,c] - v0[r,c]) > MAX_ABS_V_Error
                MAX_ABS_V_Error = abs(uy[r,c] - v0[r,c]);
                MAX_ABS_V_Error_position = (r,c);
            end

            # Calculation of MAX_RS_Error and position.
            if sqrt(e1/e2) > MAX_RS_Error
                MAX_RS_Error = sqrt(e1/e2);
                MAX_RS_Error_Position = (r,c);
            end

            # Calculation of Bulk velocity Abs error and position.
            if Bulk_Velocity_ABS_Error < abs(e2 - e3)
                Bulk_Velocity_ABS_Error = abs(e2 - e3);
                Bulk_Velocity_ABS_Error_Position = (r,c);
            end

            # Replacing the preceading time.
            u0[r,c] = ux[r,c];
            v0[r,c] = uy[r,c];

        end # for r
    end # for c

    RMS_Error = sqrt(s/((Nx-1)*(Ny-1)));

    return RMS_Error, MAX_ABS_U_Error, MAX_ABS_U_Error_position, MAX_ABS_V_Error, MAX_ABS_V_Error_position, MAX_RS_Error, MAX_RS_Error_Position, Bulk_Velocity_ABS_Error, Bulk_Velocity_ABS_Error_Position, u0, v0;
end # function [Err]



# Error RMS
# function Err_RMS!(ux::Array{Float64}, uy::Array{Float64}, u0::Array{Float64}, v0::Array{Float64})
#      e1 = 0.0;
#      e2 = 0.0;
#      s  = 0.0;
#     for c = 2:Nx
#         for r = 2:Ny # حدود محاسبه نود هایی هست که روی آنها عملیات مرزی انجام نشده است
#             e1 = (ux[r,c]-u0[r,c])^2 + (uy[r,c]-v0[r,c])^2;
#             e2 = (ux[r,c]        )^2 + (uy[r,c]        )^2;
#             s = s + e1/e2;
#             u0[r,c] = ux[r,c];
#             v0[r,c] = uy[r,c];
#         end # for r
#     end # for c
#     return sqrt(s/((Nx-1)*(Ny-1))) ,u0 , v0;
# end # function [Err]



# function Err2!(ux::Array{Float64}, uy::Array{Float64}, u0::Array{Float64}, v0::Array{Float64})
#   ERR = abs(ux - u0);
#   err = maximum(ERR);
#   u0 = deepcopy(ux);
#   v0 = deepcopy(uy);
#   return err , u0 , v0;
# end




function Data_output(LoopItterateCounter::Int64, rho::Array{Float64}, ux::Array{Float64}, uy::Array{Float64})
    FileNameAndAddress = @sprintf("%s/Results%lf.plt",resultsAddress,LoopItterateCounter/PrintDecimal)
    File = open(FileNameAndAddress,"w");
      println(File,"TITLE=\"",MethodName," Re = ",Re,"\"");
      # write(File,"SOLUTIONTIME = $(LoopItterateCounter)\n");
      # write(File,"T = \"Time=$(LoopItterateCounter)\"\n");
      write(File,"VARIABLES = X,Y,Rho,U,V\n");
      write(File,"ZONE\n");
      println(File,"I = ",Nx1,"\t","J = ",Ny1);
      write(File,"F=POINT\n");
    @inbounds  for c = 1:Nx1
          for r = 1:Ny1
              # write(File,@sprintf("%d,%d,%lf,%lf,%lf\n",c,r,rho[r,c],ux[r,c],uy[r,c]));
              println(File,c,",",r,",",rho[r,c],",",ux[r,c],",",uy[r,c]);
          end
      end
    close(File);
    println("Results saved at : ",FileNameAndAddress);
    flush(STDOUT);
end



function Data_output_CSV(LoopItterateCounter::Int64, rho::Array{Float64}, ux::Array{Float64}, uy::Array{Float64})
    UxFileName = @sprintf("%s/Ux%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal);
    File = open(UxFileName,"w");
    writecsv(File,ux);
    close(File);

    VyFileName = @sprintf("%s/Vy%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal);
    File = open(VyFileName,"w");
    writecsv(File,uy);
    close(File);

    RhoFileName = @sprintf("%s/Rho%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal);
    File = open(RhoFileName,"w");
    writecsv(File,rho);
    close(File);
    println("Ux , Vy and Rho have been saved in \"CSV\" file format.");
    flush(STDOUT);
end #function[Data_output_CSV]



function CenterlineOutput(LoopItterateCounter::Int64, rho::Array{Float64}, ux::Array{Float64}, uy::Array{Float64})
    if isodd(Nx1)
        Xcenter = floor(Int64,Nx1/2) + 1;
    else
        Xcenter = floor(Int64,Nx1/2);
    end
    if isodd(Ny1)
        Ycenter = floor(Int64,Ny1/2) + 1;
    else
        Ycenter = floor(Int64,Ny1/2);
    end
    UxFileName = @sprintf("%s/CenterUx%lf.plt",resultsAddress,LoopItterateCounter/PrintDecimal);
    File = open(UxFileName,"w");
    writecsv(File,ux[:,Xcenter]);
    close(File);

    VyFileName = @sprintf("%s/CenterVy%lf.plt",resultsAddress,LoopItterateCounter/PrintDecimal);
    File = open(VyFileName,"w");
    writecsv(File,uy[Ycenter,:]);
    close(File);

    println("Center line values have been saved.\n");
    flush(STDOUT);

end # end function[centerlineOutput]



function DistributionFunctions_output_CSV(LoopItterateCounter::Int64, f::Array{Float64})
    Distribution_Functions = @sprintf("%s/D_F[1]_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal); # file name creation
    File = open(Distribution_Functions,"w"); # file open
    writecsv(File,f[:,:,1]); # writing to file
    close(File) # closing the file

    Distribution_Functions = @sprintf("%s/D_F[2]_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal); # file name creation
    File = open(Distribution_Functions,"w"); # file open
    writecsv(File,f[:,:,2]); # writing to file
    close(File) # closing the file

    Distribution_Functions = @sprintf("%s/D_F[3]_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal); # file name creation
    File = open(Distribution_Functions,"w"); # file open
    writecsv(File,f[:,:,3]); # writing to file
    close(File); # closing the file

    Distribution_Functions = @sprintf("%s/D_F[4]_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal); # file name creation
    File = open(Distribution_Functions,"w"); # file open
    writecsv(File,f[:,:,4]); # writing to file
    close(File); # closing the file

    Distribution_Functions = @sprintf("%s/D_F[5]_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal); # file name creation
    File = open(Distribution_Functions,"w"); # file open
    writecsv(File,f[:,:,5]); # writing to file
    close(File); # closing the file

    Distribution_Functions = @sprintf("%s/D_F[6]_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal); # file name creation
    File = open(Distribution_Functions,"w"); # file open
    writecsv(File,f[:,:,6]); # writing to file
    close(File); # closing the file

    Distribution_Functions = @sprintf("%s/D_F[7]_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal); # file name creation
    File = open(Distribution_Functions,"w"); # file open
    writecsv(File,f[:,:,7]); # writing to file
    close(File); # closing the file

    Distribution_Functions = @sprintf("%s/D_F[8]_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal); # file name creation
    File = open(Distribution_Functions,"w"); # file open
    writecsv(File,f[:,:,8]); # writing to file
    close(File); # closing the file

    Distribution_Functions = @sprintf("%s/D_F[9]_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal); # file name creation
    File = open(Distribution_Functions,"w"); # file open
    writecsv(File,f[:,:,9]); # writing to file
    close(File); # closing the file

    println("Distibution functions has been saved in \"CSV\" file format.\n");
    flush(STDOUT);
end #function[Data_output_CSV]



#///////////////////////////////////////////////////////////////////////////////
#
##                              Function main
#
#///////////////////////////////////////////////////////////////////////////////

function main()
    # println("Start...\n");
    # flush(STDOUT);
    #/==============================================================================
    #                             Loop variables
    #/==============================================================================

    LoopItterateCounter::Int64 = 0;
    err::Float64 = 1.0000000001;


    #/==============================================================================
    #                files opennings for saving data in the loop
    #/==============================================================================
    @static if DistributionFunctions_sum_Out == "yes"
        Distribution_Functions_sum = @sprintf("%s/A_Distribution_functions_Sum_for_%s_Re%1.1lf_N%d_%d.csv",resultsAddress,MethodName,Re,Nx1,Ny1); # file name creation
        Filesum = open(Distribution_Functions_sum,"w"); # file open
    end # end static if

    @static if ErrorPrintingInInFilePermission == "yes"
      ErrorFileName = @sprintf("%s/A_Error_Values_for_%s_Re%1.1lf_N%d_%d.csv",resultsAddress,MethodName,Re,Nx1,Ny1); # file name creation

      ErrorFilePointer = open(ErrorFileName,"w"); # file open

      println(ErrorFilePointer,"Itterate ,RMS_Error ,MAX_ABS_U_Error ,MAX_ABS_U_Error_position_R ,MAX_ABS_U_Error_position_C ,MAX_ABS_V_Error ,MAX_ABS_V_Error_position_R ,MAX_ABS_V_Error_position_C ,MAX_RS_Error ,MAX_RS_Error_Position_R ,MAX_RS_Error_Position_C ,Bulk_Velocity_ABS_Error ,Bulk_Velocity_ABS_Error_Position_R ,Bulk_Velocity_ABS_Error_Position_C");
      flush(ErrorFilePointer);
    end # end static of



    #/==============================================================================
    #                              Allocations
    #/==============================================================================

    f       = Array{Float64}(Ny1,Nx1,Q); # array of the distribution functions (DFs)
    f_post  = Array{Float64}(Ny1,Nx1,Q);

    rho     = Array{Float64}(Ny1,Nx1);   # arrays of fluid density and velocity
    ux      = Array{Float64}(Ny1,Nx1);
    uy      = Array{Float64}(Ny1,Nx1);

    u0      = Array{Float64}(Ny1,Nx1);   # for calculations in Err function
    v0      = Array{Float64}(Ny1,Nx1);

#     M2      = floor(Int64,Ny/2);         # Center points for centerline
#     N2      = floor(Int64,Nx/2);         # If nodes are even.

# ================= Allocations for Error function =============================

    RMS_Error::Float64 = 0.0;
    MAX_ABS_U_Error::Float64 = 0.0;
    MAX_ABS_U_Error_position::Tuple{Int64,Int64} = (0,0);
    MAX_ABS_V_Error::Float64 = 0.0;
    MAX_ABS_V_Error_position::Tuple{Int64,Int64} = (0,0);
    MAX_RS_Error::Float64 = 0.0;
    MAX_RS_Error_Position::Tuple{Int64,Int64} = (0,0);
    Bulk_Velocity_ABS_Error::Float64 = 0.0;
    Bulk_Velocity_ABS_Error_Position::Tuple{Int64,Int64} = (0,0);

    #/================= SRT Special Variables and allocations ======================

    tau::Float64     = 3*L*uw_top/Re + 0.5;      # relaxation time for BGK



    # ================= MRT_Guo Special Variables and allocations ======================

    @static if MethodName == "MRT_Guo"  # static if MRT_Guo Method

    m       = Array{Float64}(Q);                 # moments array
    s       = Array{Float64}(Q);                 # Relaxation rates array

    Lenght::Int64    = L # = N + 1;              # For halfway Bondary Condition

    tauMRT::Float64  = 3*Lenght*uw_top/Re + 0.5; # relaxation time for BGK

    # S = (0,Se,S_epsilon,0,Sq,0,Sq,Sv,Sv)
    # Sv = viscosity shear 
    # Se = viscosity shear = viscosity bulk

    s[7]    = 1.0/tauMRT;                        # relaxation rates for MRT
    s[8]    = 1.0/tauMRT;
    s[9]    = 0.0;
    s[1]    = 1.6;
    s[2]    = 1.8;
    s[3]    = 0.0;
    s[4]    = 8*(2-s[7])/(8-s[7]);
    s[5]    = 0.0;
    s[6]    = 8*(2-s[7])/(8-s[7]);
    





    end # end [ static if MRT_Guo]
    
    
    
    @static if MethodName == "MRT_Khatoon"  # static if MRT_Guo Method

    m       = Array{Float64}(Q);                 # moments array
    s       = Array{Float64}(Q);                 # Relaxation rates array

    Lenght::Int64    = L # = N + 1;              # For halfway Bondary Condition

    tauMRT::Float64  = 3*Lenght*uw_top/Re + 0.5; # relaxation time for BGK

    # S = (0,Se,S_epsilon,0,Sq,0,Sq,Sv,Sv)
    # Sv = viscosity shear 
    # Se = viscosity shear = viscosity bulk

    s[7]    = 1.0/tauMRT;                        # relaxation rates for MRT
    s[8]    = 1.0/tauMRT;
    s[9]    = 1.0;
    s[1]    = 1.2;
    s[2]    = 1.1;
    s[3]    = 1.0;
    s[4]    = 3*(2-s[7])/(3-s[7]);
    s[5]    = 1.0;
    s[6]    = 3*(2-s[7])/(3-s[7]);
    





    end # end [ static if MRT_Guo]

    
    
    # ================= MRT_Kruger Special Variables and allocations ======================

    @static if MethodName == "MRT_EAslan"  # static if MRT_Guo Method

    m       = Array{Float64}(Q);                 # moments array
    s       = Array{Float64}(Q);                 # Relaxation rates array

    Lenght::Int64    = L # = N + 1;              # For halfway Bondary Condition

    tauMRT::Float64  = 3*Lenght*uw_top/Re + 0.5; # relaxation time for BGK

    # S = (0,Se,S_epsilon,0,Sq,0,Sq,Sv,Sv)
    # Sv = viscosity shear 
    # Se = viscosity shear = viscosity bulk

    s[9]    = 1.0;
    s[1]    = 1.4;
    s[2]    = 1.4;
    s[3]    = 1.0;
    s[4]    = 1.2;
    s[5]    = 1.0;
    s[6]    = 1.2;
    s[7]    = 1.0/tauMRT;                        # relaxation rates for MRT
    s[8]    = 1.0/tauMRT;
    

    end # end [ static if MRT-EAslan]
    
    
    

    # ================= MRT_Lalle Special Variables and allocations ======================

    @static if MethodName == "MRT_Lalle"  # static if MRT_Guo Method

    m       = Array{Float64}(Q);                 # moments array
    s       = Array{Float64}(Q);                 # Relaxation rates array

    Lenght::Int64    = L # = N + 1;              # For halfway Bondary Condition

    tauMRT::Float64  = 3*Lenght*uw_top/Re + 0.5; # relaxation time for BGK

    # s = (s9,s1,s2,s3,s4,s5,s6,s7,s8)

    # Based on Pradipto 2017 And EAslani And lalle 
    
    s[1]    = 1.64;                        # relaxation rates for MRT
    s[2]    = 1.14;
    s[3]    = 0.0;
    s[4]    = 1.92;
    s[5]    = 0.0;
    s[6]    = 1.92;
    s[7]    = 1.0/tauMRT;
    s[8]    = 1.0/tauMRT;
    s[9]    = 0.0;

    end # end [ static if MRT_Lalle]


    # ================= MRT_Lalle Special Variables and allocations ======================

    @static if MethodName == "Pradipto"  # static if MRT_Guo Method

    m       = Array{Float64}(Q);                 # moments array
    s       = Array{Float64}(Q);                 # Relaxation rates array

    Lenght::Int64    = L # = N + 1;              # For halfway Bondary Condition

    tauMRT::Float64  = 3*Lenght*uw_top/Re + 0.5; # relaxation time for BGK

    # s = (s9,s1,s2,s3,s4,s5,s6,s7,s8)


    s[9]    = 0.0;
    s[1]    = 1.64;
    s[2]    = 1.54;
    s[3]    = 0.0;
    s[4]    = 1.92;
    s[5]    = 0;
    s[6]    = 1.92;
    s[7]    = 1.0/tauMRT;                        # relaxation rates for MRT
    s[8]    = 1.0/tauMRT;
    

    end # end [ static if MRT_Lalle2]



    # ================= Cas_Geier Special Variables and allocations ===============

    @static if MethodName == "Cas_Geier"

    # s = (1,1,1,omega_bulk,omega_shear,omega_shear,1,1,1)

    # k3 = kxxpyy omega_bulk
    # k4 = kxxmyy omega_shear
    # k5 = kxy    omega_shear
    # k6 = kxxy
    # k7 = kxyy
    # k8 = kxxyy

    Lenght::Int64    = L # = N + 1;              # For halfway Bondary Condition
    tauCascade::Float64  = 3*Lenght*uw_top/Re + 0.5; # relaxation time for BGK

    #omega_bulk::Float64 = 1/tauCascade;
    omega_bulk::Float64 = 1.0;
    omega_shear::Float64 = 1/tauCascade;
    omega6::Float64 = 1.0;
    omega7::Float64 = 1.0;
    omega8::Float64 = 1.0;

    end # end if [MethodName Cas_Geier]



    # ================= Cas_Lalle Special Variables and allocations ===============

    @static if MethodName == "Cas_Lalle"

    # k3 = kxxpyy omega_bulk
    # k4 = kxxmyy omega_shear
    # k5 = kxy    omega_shear
    # k6 = kxxy
    # k7 = kxyy
    # k8 = kxxyy

    # S = (0,0,0,1.64,omega_shear,omega_shear,1.90,1.90,1.54)

    Lenght::Int64    = L # = N + 1;              # For halfway Bondary Condition
    tauCascade::Float64  = 3*Lenght*uw_top/Re + 0.5; # relaxation time for BGK

    #omega_bulk::Float64 = 1/tauCascade;
    omega_bulk::Float64 = 1.64;
    omega_shear::Float64 = 1/tauCascade;
    omega6::Float64 = 1.90;
    omega7::Float64 = 1.90;
    omega8::Float64 = 1.54;

    end # end if [MethodName Cas_Lalle]


    # ================= Cas_Lalle Special Variables and allocations ===============

    @static if MethodName == "Cas_Chavez"

    # k3 = kxxpyy omega_bulk
    # k4 = kxxmyy omega_shear
    # k5 = kxy    omega_shear
    # k6 = kxxy
    # k7 = kxyy
    # k8 = kxxyy

    # S = (0,0,0,1.64,omega_shear,omega_shear,0.90,0.90,1.54)

    Lenght::Int64    = L # = N + 1;              # For halfway Bondary Condition
    tauCascade::Float64  = 3*Lenght*uw_top/Re + 0.5; # relaxation time for BGK

    #omega_bulk::Float64 = 1/tauCascade;
    omega_bulk::Float64 = 1.64;
    omega_shear::Float64 = 1/tauCascade;
    omega6::Float64 = 0.90;
    omega7::Float64 = 0.90;
    omega8::Float64 = 1.54;

    end # end if [MethodName Cas_Chavez]


    #/==============================================================================
    #                   Initialization of rho, ux, uy and f
    #/==============================================================================

    rho, ux, uy, f = Init_Eq!(rho, ux, uy, f);

    # Debuging print
    #println("rho \n"); show(rho); println("ux \n"); show(ux); println("uy \n"); show(uy); println("f \n"); show(f);



    #/==============================================================================

    #/********************* Itterative loop starts here ****************************

    #/==============================================================================

     while (RMS_Error > errorlimitation) || (LoopItterateCounter <= (ErrCalculationSteps*2 + 1)) # method 1
    # for it = 1:MaxItteration                                                            # method 3
    # for it = 1:5000
        LoopItterateCounter = LoopItterateCounter + 1;
     
     #=
    .__________________________________________________________.
    |                        ATTENTION:                        |          
    |             Another criteria for breaking the while loop |
    | for steady and unsteady simulations is defined under     |
    | error calculation cheque criteria [which is placed after |
    | solving algorithem selection]                            |
    .__________________________________________________________.
     =# 

        #/==========================================================================
        #/======================== solving algorithem selection ====================
        #/==========================================================================



        # BGK Simple Boundary ======================================================
        @static if MethodName == "BGK"

        f_post = Coll_BGK!(rho, ux, uy, f, f_post, tau);               # BGK collision
        f = Streaming!(f,f_post);                                      # Streaming for BGK
        f = Bounce_back!(f,f_post,rho);                                # No-Slip boundary condition
        rho, ux, uy = Den_vel!(rho,ux,uy,f);                           # Fluid variables calculation for all models


        elseif MethodName == "BGK_NEBC"   # static elseif
        # BGK Non_equilibrium_Boundarycondition ====================================

         f_post = Coll_BGK!(rho, ux, uy, f, f_post, tau);              # BGK collision
         f = Non_Equilibrium_Streaming!(f,f_post);                     # Streaming For Non_equilibrium_Boundarycondition
         rho, ux, uy = Den_vel!(rho,ux,uy,f);                          # Fluid variables calculation for all models
         f = Non_equilibrium_Boundarycondition!(f,rho,ux,uy);          # Non_equilibrium_Boundarycondition


        elseif MethodName == "BGK_MNEBC"   # static elseif
        # BGK Non_equilibrium_Modified_Boundarycondition ===========================

         f_post = Coll_BGK!(rho, ux, uy, f, f_post, tau);              # BGK collision
         f = Non_Equilibrium_Streaming!(f,f_post);                     # Streaming For Non_equilibrium_Boundarycondition
         rho, ux, uy = Den_vel!(rho,ux,uy,f);                          # Fluid variables calculation for all models
         f = Non_equilibrium_Boundarycondition_modified!(f,rho,ux,uy); # Non_equilibrium_Modified_Boundarycondition


        elseif MethodName == "BGK_SNEBC"   # static elseif
        # BGK Semi_Non_equilibrium_Boundarycondition ===============================

           #_____________________________________________________________
           #|                                                            |
           #|            Method is NOT WORKING Code works.               |
           #| Non_Equilibrium_Streaming is Better than Simple Streaming. |
           #| But the boundary at the bottom at which the velocities are |
           #| zero does not have any good results.                       |
           #|____________________________________________________________|


          f_post = Coll_BGK!(rho, ux, uy, f, f_post, tau);               # BGK collision
          # f = Streaming!(f,f_post);                                      # Streaming for BGK
          f = Non_Equilibrium_Streaming!(f,f_post);                      # Streaming For Non_equilibrium_Boundarycondition
          rho, ux, uy = Den_vel!(rho,ux,uy,f);                           # Fluid variables calculation for all models
          f = Semi_Non_equilibrium_Boundarycondition!(f,rho,ux,uy);      # Non_equilibrium_Boundarycondition



        elseif MethodName == "MRT_Guo" || MethodName == "MRT_EAslan" || MethodName == "MRT_Lalle" || MethodName == "Pradipto" ||  MethodName == "MRT_Khatoon"  # static elseif
        # MRT From LBM01.c Guo book ================================================
         f_post = Coll_MRT!(rho, ux, uy, f, f_post, m, s);             # BGK collision
         f = Streaming!(f,f_post);                                     # Streaming for BGK
         f = Bounce_back!(f,f_post,rho);                               # No-Slip boundary condition
         rho, ux, uy = Den_vel!(rho,ux,uy,f);                          # Fluid variables calculation for all models

        elseif MethodName == "MRT2"
        # MRT From LBM01.c Guo book ================================================
         f_post = Coll_MRT2!(rho, ux, uy, f, f_post);                 # BGK collision
         f = Streaming!(f,f_post);                                     # Streaming for BGK
         f = Bounce_back!(f,f_post,rho);                               # No-Slip boundary condition
         rho, ux, uy = Den_vel!(rho,ux,uy,f);                          # Fluid variables calculation for all models


         
         
        elseif MethodName == "CMs_NON_Orthogonal" # static elseif
        # CMs. Non-Orthogonal Central Moments space ================================
        f_post = Coll_CMs_Non_Orthogonal!(rho, ux, uy, f, f_post, tau); # CMs Non-Orthogonal collision
        f      = Streaming!(f,f_post);                                  # Streaming for BGK
        f      = Bounce_back!(f,f_post,rho);                            # No-Slip boundary condition
        rho, ux, uy = Den_vel!(rho,ux,uy,f);                            # Fluid variables calculation for all models


        elseif MethodName == "Cas_Geier" || MethodName == "Cas_Lalle" || MethodName == "Cas_Chavez" # static if
        # Cascaded LBM ============================ ================================
        f_post = CascadedLBMCollition!(rho,ux,uy,f,f_post,omega_bulk,omega_shear,omega6,omega7,omega8); # Cascaded LBM collision.
        f = Streaming!(f,f_post);                                     # Streaming for BGK
        f = Bounce_back!(f,f_post,rho);                               # No-Slip boundary condition
        rho, ux, uy = Den_vel!(rho,ux,uy,f);                          # Fluid variables calculation for all models


        else # static else
        println("You did not chosen the right method.\n","Your input for the MethodName is (\"",MethodName,"\") which is not a SELECTABLE method.\n");
        end
        #/==========================================================================
        #/=========================== Post-processing ==============================
        #/==========================================================================


        if LoopItterateCounter % ErrCalculationSteps == 0

            RMS_Error, MAX_ABS_U_Error, MAX_ABS_U_Error_position, MAX_ABS_V_Error, MAX_ABS_V_Error_position, MAX_RS_Error, MAX_RS_Error_Position, Bulk_Velocity_ABS_Error, Bulk_Velocity_ABS_Error_Position, u0, v0 = ErrorManager!(ux, uy, u0, v0, RMS_Error, MAX_ABS_U_Error, MAX_ABS_U_Error_position, MAX_ABS_V_Error, MAX_ABS_V_Error_position, MAX_RS_Error, MAX_RS_Error_Position, Bulk_Velocity_ABS_Error, Bulk_Velocity_ABS_Error_Position); # Velocity differences between two successive [ErrCalculationSteps] steps

            # println("Err1 = ",err); # Error with difference method
            # err ,u0 ,v0 = Err2!(ux, uy, u0, v0);
#             println("err = ", err , " ux_center = ", ux[M2,N2] , " uy_center = ", uy[M2,N2] , " itt = ", LoopItterateCounter);  # Display some results


            # NEW Criteria For Re > 7500 
            if ((Re > 7500) && (LoopItterateCounter > MAX_UnSteady_Itteration)) || (LoopItterateCounter > MAX_UnSteady_Itteration)
            println("Reached MAX_UnSteady_Itteration: converged.")
            break;
            end 

            @static if STDOUT_Error_Report_Size == "Short"
            println("RMS_Error = ", err, " itt = ", LoopItterateCounter);

            elseif STDOUT_Error_Report_Size == "Long" # static elseif

            # Showing errors in terminal
            println("======================================================================")
            println("Itterate = ",
                    LoopItterateCounter,
                    "\t RMS_Error = ",
                    RMS_Error,
                    "\n",
                    "MAX_ABS_U_Error = ",
                    MAX_ABS_U_Error,
                    "\t\t @ ->",
                    MAX_ABS_U_Error_position,
                    "\n",
                    "MAX_ABS_V_Error = ",
                    MAX_ABS_V_Error ,
                    "\t\t @ ->",
                    MAX_ABS_V_Error_position,
                    "\n",
                    "MAX_RS_Error    = ",
                    MAX_RS_Error,
                    "\t\t @ ->",
                    MAX_RS_Error_Position,
                    "\n",
                    "Bulk_Velocity_ABS_Error = " ,
                    Bulk_Velocity_ABS_Error,
                    "\t @ ->", Bulk_Velocity_ABS_Error_Position);
            println("======================================================================")
            flush(STDOUT);

            #=
            +--------------------------------------------------+
            |                      Errors                      |
            +--------------------------------------------------+
            |                Iterate = 0.0000000               |
            +--------------------------------------------------+
            | No |    Method    |       Value       | Position |
            +----+--------------+-------------------+----------+
            |  1 | RMS          | 0.000000000000000 |  (00,00) |
            +----+--------------+-------------------+----------+
            |  2 | MAX Abs      | 0.000000000000000 |  (00,00) |
            +----+--------------+-------------------+----------+
            |  3 | MAX Abs U    | 0.000000000000000 |  (00,00) |
            +----+--------------+-------------------+----------+
            |  4 | MAX Abs V    | 0.000000000000000 |  (00,00) |
            +----+--------------+-------------------+----------+
            |  5 | MAX RS       | 0.000000000000000 |  (00,00) |
            +----+--------------+-------------------+----------+
            |  6 | MAX Abs Bulk | 0.000000000000000 |  (00,00) |
            +----+--------------+-------------------+----------+
            =#

            else #static else

            end # end static if [STDOUT_Error_Report_Size]


            # Saving Errors in a file.
            @static if ErrorPrintingInInFilePermission == "yes"

              println(ErrorFilePointer,
                      LoopItterateCounter,
                      ",",
                      RMS_Error,",",
                      MAX_ABS_U_Error,",",
                      MAX_ABS_U_Error_position[1],",",
                      MAX_ABS_U_Error_position[2],",",
                      MAX_ABS_V_Error,",",
                      MAX_ABS_V_Error_position[1],",",
                      MAX_ABS_V_Error_position[2],",",
                      MAX_RS_Error,",",
                      MAX_RS_Error_Position[1],",",
                      MAX_RS_Error_Position[2],",",
                      Bulk_Velocity_ABS_Error,",",
                      Bulk_Velocity_ABS_Error_Position[1],",",
                      Bulk_Velocity_ABS_Error_Position[2]
                    );

              flush(ErrorFilePointer);                                      # flushing the Buffer.
            end # end static if [ErrorPrintingInInFilePermission]

        end # end if[ErrCalculationSteps]


        @static if PrintOutput == "yes"
        if LoopItterateCounter % PrintSteps == 0
            Data_output(LoopItterateCounter,rho,ux,uy);
        end # end if [PrintSteps]
        end # end static if

        @static if CSVoutput=="yes"
        if LoopItterateCounter % CSVoutputFrequency == 0
            Data_output_CSV(LoopItterateCounter,rho,ux,uy);
        end # end if[LoopItterateCounter % CSVoutputFrequency]
        end #end static if

        @static if CenterLinePrint == "yes"
        if LoopItterateCounter % CenterlinePrintFrequency == 0
            CenterlineOutput(LoopItterateCounter,rho,ux,uy);
        end # end if[LoopItterateCounter % CenterlinePrintFrequency]
        end # end static if

        @static if DistributionFunctions_Out == "yes"
        if LoopItterateCounter % DistributionFunctionsFrequency == 0
            DistributionFunctions_output_CSV(LoopItterateCounter,f);
        end # end of if [LoopItterateCounter % DistributionFunctionsFrequency]
        end # end static if

        @static if DistributionFunctions_sum_Out == "yes"
        if LoopItterateCounter % DistributionFunctionsSumFrequency == 0
#           sum = 0;
#           for i = 1:Nx1
#             for j = 1:Ny1
#               for k = 1:Q
#                 sum = sum + f[j,i,k];
#               end
#             end
#           end
#        println(Filesum,LoopItterateCounter,"\t",sum);# writing to Filesum


        # writing sum of populations by a reduction function. to a file pointed at with Filesum.
          println(Filesum,LoopItterateCounter,"\t",reduce(+,f));# writing sum of populations by a reduction function. to a file pointed at with Filesum.
          flush(Filesum);
        end # end of if [LoopItterateCounter % DistributionFunctionsFrequency]
        end #end static if


        # println(LoopItterateCounter) # printing the LoopItterateCounter in the STDOUT



    end # end [while] or [for] (Itterative loop) function main still continues.



    #/==============================================================================
    #                       Closing of Opened file
    #/==============================================================================

    @static if DistributionFunctions_sum_Out == "yes"
      close(Filesum); # closing the file
    end # end static if



    @static if ErrorPrintingInInFilePermission == "yes"
      close(ErrorFilePointer); # closing the file
    end # end static if



    #/==============================================================================
    #                   Recording of last converged data
    #/==============================================================================

    Data_output(LoopItterateCounter,rho,ux,uy);       # Output last simulation data
    Data_output_CSV(LoopItterateCounter,rho,ux,uy);   # Output last simulation data in CSV
    CenterlineOutput(LoopItterateCounter,rho,ux,uy);  # output last simulation centerline data in csv



    return rho,ux,uy;
end # end function [main]



#/==============================================================================
#                   Execution Of main() function
#/==============================================================================

# rho,ux,uy = main();            # to use final results

@time main();                    # to time and measure performance
