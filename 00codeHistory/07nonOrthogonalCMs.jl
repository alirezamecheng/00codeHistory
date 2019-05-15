#=
                    In the name of Allah.
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
    * version 07, code name: 07nonOrthogonalCMs
        1- Implimentation of Centeral moments space method.
        
    * version 06, code name: 06CleanupAndComment
        1- adding some forgotten ";".
        2- more comments in code.
        
    * version 04, code name: 04MRTGue
        1- All Int32 variable types changed to Int64.
        2- MRT code with Gue variables added.
        3- Adding a variable to determine the name of saving folder.




=#

#/==============================================================================
#               Definition of constants and global variables
#/==============================================================================

const Nx = 127        # number of cells in the x-direction -1
const Ny = 127        # number of cells in the y-direction -1
const Re = 400.0
const errorlimitation = 1.0e-5;
const ErrCalculationSteps = 1000;
const MaxItteration = 1000;


#/==============================================================================
#                             Choosing Method!
#/==============================================================================
# const MethodName = "MRT";
const MethodName = "NOCMS";



#/==============================================================================
#                       Printing And I/O Parameters
#/==============================================================================

const PrintDecimal = 1000000; # number of Decimals for file name.

const PrintOutput = "yes";
# const PrintOutput = "No";
const PrintSteps = 1000;

# const CSVoutput = "yes";
const CSVoutput = "No";
const CSVoutputFrequency = 1000;

# const CenterLinePrint = "yes";
const CenterLinePrint = "No";
const CenterlinePrintFrequency = 20000;

# const DistributionFunctions_Out = "yes";
const DistributionFunctions_Out = "No";
const DistributionFunctionsFrequency = 100;

# const DistributionFunctions_sum_Out = "yes";
const DistributionFunctions_sum_Out = "No";
const DistributionFunctionsSumFrequency = 10;

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

const resultsAddress = "/home/aghnia/Results/$(MethodName)_$(Dates.format(now(),"YY.mm.dd.HH.MM"))_Re($(Re))_Gr($(Nx1)_$(Ny1))_Uw($(uw_top))"; # for using in cluster.

# const resultsAddress = "/home/ghavami/mscThesis/runZone/SRT_Non_equilibrum_boundary_$(Dates.format(now(),"YY.mm.dd.HH.MM"))_SRT_Re($(Re))_Gr($(Nx1)_$(Ny1))_Uw($(uw_top))"; # for using in my own pc.

mkpath(resultsAddress); # to creat specified path



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
function meq(RHO::Float64, U::Float64, V::Float64, k::Int64)



    if k == 9
        x = RHO;
    elseif k == 1
        x = RHO*(-2+3*(U*U+V*V));
    elseif k == 2
        x = RHO*(1-3*(U*U+V*V));
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


# Initialization
function Init_Eq!(rho::Array{Float64}, ux::Array{Float64}, uy::Array{Float64}, f::Array{Float64})
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
    return rho, ux, uy, f;
end #end function Init_Eq.



# BGK collision
function Coll_BGK!(rho::Array{Float64}, ux::Array{Float64}, uy::Array{Float64}, f::Array{Float64}, f_post::Array{Float64}, tau::Float64)
    #double FEQ
    for k = 1:Q
        for c = 1:Nx1
            for r = 1:Ny1
                FEQ = feq(rho[r,c], ux[r,c], uy[r,c], k); #EDF
                f_post[r,c,k] = f[r,c,k] - (f[r,c,k] - FEQ)/tau; #Post-collision DFs
            end # for r
        end # for c
    end #for k
    return f_post;
end #function [Coll_BGK]



# MRT collision
function Coll_MRT!(rho::Array{Float64}, ux::Array{Float64}, uy::Array{Float64}, f::Array{Float64}, f_post::Array{Float64}, m::Array{Float64}, s::Array{Float64})

    MEQ::Float64 = 0.0;

    for c = 1:Nx1
        for r = 1:Ny1

            # Transformation from velocity space to moment space:
            m[9] =     f[r,c,9] +    f[r,c,1]  + f[r,c,2] + f[r,c,3] + f[r,c,4]  +    f[r,c,5] + f[r,c,6] + f[r,c,7] + f[r,c,8];
            m[1] =  -4*f[r,c,9] -    f[r,c,1]  - f[r,c,2] - f[r,c,3] - f[r,c,4]  + 2*(f[r,c,5] + f[r,c,6] + f[r,c,7] + f[r,c,8]);
            m[2] =   4*f[r,c,9] - 2*(f[r,c,1]  + f[r,c,2] + f[r,c,3] + f[r,c,4]) +    f[r,c,5] + f[r,c,6] + f[r,c,7] + f[r,c,8];
            m[3] =     f[r,c,1] -    f[r,c,3]  + f[r,c,5] - f[r,c,6] - f[r,c,7]  +    f[r,c,8];
            m[4] = -2*(f[r,c,1] -    f[r,c,3]) + f[r,c,5] - f[r,c,6] - f[r,c,7]  +    f[r,c,8];
            m[5] =     f[r,c,2] -    f[r,c,4]  + f[r,c,5] + f[r,c,6] - f[r,c,7]  -    f[r,c,8];
            m[6] = -2*(f[r,c,2] -    f[r,c,4]) + f[r,c,5] + f[r,c,6] - f[r,c,7]  -    f[r,c,8];
            m[7] =     f[r,c,1] -    f[r,c,2]  + f[r,c,3] - f[r,c,4];
            m[8] =     f[r,c,5] -    f[r,c,6]  + f[r,c,7] - f[r,c,8];


            # Relaxation in moment space:
            for k = 1:Q

                MEQ  = meq(rho[r,c], ux[r,c], uy[r,c], k);
                m[k] = m[k] - s[k]*(m[k] - MEQ);            # relaxation
                m[k] = m[k] / D[k];                         # rescaling

            end #end[for k]

            # Transforming back to the velocity space:

            f_post[r,c,9] = m[9] - 4*(m[1] -    m[2]);
            f_post[r,c,1] = m[9] -    m[1] - 2*(m[2] + m[4]) + m[3] + m[7];
            f_post[r,c,2] = m[9] -    m[1] - 2*(m[2] + m[6]) + m[5] - m[7];
            f_post[r,c,3] = m[9] -    m[1] - 2*(m[2] - m[4]) - m[3] + m[7];
            f_post[r,c,4] = m[9] -    m[1] - 2*(m[2] - m[6]) - m[5] - m[7];
            f_post[r,c,5] = m[9] +    m[1] +    m[1] + m[2]  + m[3] + m[4] + m[5] + m[6] + m[8];
            f_post[r,c,6] = m[9] +    m[1] +    m[1] + m[2]  - m[3] - m[4] + m[5] + m[6] - m[8];
            f_post[r,c,7] = m[9] +    m[1] +    m[1] + m[2]  - m[3] - m[4] - m[5] - m[6] + m[8];
            f_post[r,c,8] = m[9] +    m[1] +    m[1] + m[2]  + m[3] + m[4] - m[5] - m[6] - m[8];

        end #end [for r]
    end #end [for c]

    return f_post;
end # function [Coll_MRT]



# function Central Moments Space Collision
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

    
    
    for c = 1:Nx1
        for r = 1:Ny1    
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
                
            K0 =                 f[r,c,9]+                  f[r,c,1]+                  f[r,c,2]+                  f[r,c,3]+                  f[r,c,4]+                  f[r,c,5]+                  f[r,c,6]+                  f[r,c,7]+                  f[r,c,8];

            K1=(cx0            )*f[r,c,9]+(cx1            )*f[r,c,1]+(cx2            )*f[r,c,2]+(cx3            )*f[r,c,3]+(cx4            )*f[r,c,4]+(cx5            )*f[r,c,5]+(cx6            )*f[r,c,6]+(cx7            )*f[r,c,7]+(cx8            )*f[r,c,8];

            K2=(cy0            )*f[r,c,9]+(cy1            )*f[r,c,1]+(cy2            )*f[r,c,2]+(cy3            )*f[r,c,3]+(cy4            )*f[r,c,4]+(cy5            )*f[r,c,5]+(cy6            )*f[r,c,6]+(cy7            )*f[r,c,7]+(cy8            )*f[r,c,8];

            K3=(cx0*cx0+cy0*cy0)*f[r,c,9]+(cx1*cx1+cy1*cy1)*f[r,c,1]+(cx2*cx2+cy2*cy2)*f[r,c,2]+(cx3*cx3+cy3*cy3)*f[r,c,3]+(cx4*cx4+cy4*cy4)*f[r,c,4]+(cx5*cx5+cy5*cy5)*f[r,c,5]+(cx6*cx6+cy6*cy6)*f[r,c,6]+(cx7*cx7+cy7*cy7)*f[r,c,7]+(cx8*cx8+cy8*cy8)*f[r,c,8];

            K4=(cx0*cx0-cy0*cy0)*f[r,c,9]+(cx1*cx1-cy1*cy1)*f[r,c,1]+(cx2*cx2-cy2*cy2)*f[r,c,2]+(cx3*cx3-cy3*cy3)*f[r,c,3]+(cx4*cx4-cy4*cy4)*f[r,c,4]+(cx5*cx5-cy5*cy5)*f[r,c,5]+(cx6*cx6-cy6*cy6)*f[r,c,6]+(cx7*cx7-cy7*cy7)*f[r,c,7]+(cx8*cx8-cy8*cy8)*f[r,c,8];

            K5=(cx0*cy0        )*f[r,c,9]+(cx1*cy1        )*f[r,c,1]+(cx2*cy2        )*f[r,c,2]+(cx3*cy3        )*f[r,c,3]+(cx4*cy4        )*f[r,c,4]+(cx5*cy5        )*f[r,c,5]+(cx6*cy6        )*f[r,c,6]+(cx7*cy7        )*f[r,c,7]+(cx8*cy8        )*f[r,c,8];

            K6=(cx0*cx0*cy0    )*f[r,c,9]+(cx1*cx1*cy1    )*f[r,c,1]+(cx2*cx2*cy2    )*f[r,c,2]+(cx3*cx3*cy3    )*f[r,c,3]+(cx4*cx4*cy4    )*f[r,c,4]+(cx5*cx5*cy5    )*f[r,c,5]+(cx6*cx6*cy6    )*f[r,c,6]+(cx7*cx7*cy7    )*f[r,c,7]+(cx8*cx8*cy8    )*f[r,c,8];

            K7=(cx0*cy0*cy0    )*f[r,c,9]+(cx1*cy1*cy1    )*f[r,c,1]+(cx2*cy2*cy2    )*f[r,c,2]+(cx3*cy3*cy3    )*f[r,c,3]+(cx4*cy4*cy4    )*f[r,c,4]+(cx5*cy5*cy5    )*f[r,c,5]+(cx6*cy6*cy6    )*f[r,c,6]+(cx7*cy7*cy7    )*f[r,c,7]+(cx8*cy8*cy8    )*f[r,c,8];

            K8=(cx0*cx0*cy0*cy0)*f[r,c,9]+(cx1*cx1*cy1*cy1)*f[r,c,1]+(cx2*cx2*cy2*cy2)*f[r,c,2]+(cx3*cx3*cy3*cy3)*f[r,c,3]+(cx4*cx4*cy4*cy4)*f[r,c,4]+(cx5*cx5*cy5*cy5)*f[r,c,5]+(cx6*cx6*cy6*cy6)*f[r,c,6]+(cx7*cx7*cy7*cy7)*f[r,c,7]+(cx8*cx8*cy8*cy8)*f[r,c,8];
            
            
            
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



# Streaming
function Streaming!(f::Array{Float64},f_post::Array{Float64})
    for k = 1:Q
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
function Non_Equilibrium_Streaming!(f::Array{Float64},f_post::Array{Float64})
    for k = 1:Q
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



# Bounce-back scheme
function Bounce_back!(f::Array{Float64},f_post::Array{Float64},rho::Array{Float64})

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

     return f;
end # function [Bounce_back!]



# Non-equilibrium boundary condition
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



# Non-equilibrium Modified boundary condition
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




# Semi_Non-equilibrium_boundary_condition
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




# Fluid variables (density and velocity)
function Den_vel!(rho::Array{Float64},ux::Array{Float64},uy::Array{Float64},f::Array{Float64})
    for c = 1:Nx1
        for r = 1:Ny1
            rho[r,c]=f[r,c,9]+f[r,c,1]+f[r,c,2]+f[r,c,3]+f[r,c,4]+f[r,c,5]+f[r,c,6]+f[r,c,7]+f[r,c,8];
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



# Error RMS
function Err!(ux::Array{Float64}, uy::Array{Float64}, u0::Array{Float64}, v0::Array{Float64})
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
    return sqrt(s/((Nx-1)*(Ny-1))) ,u0 , v0;
end # function [Err]



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
      write(File,"TITLE=\"SRT Re = $(Re)\"\n");
      # write(File,"SOLUTIONTIME = $(LoopItterateCounter)\n");
      # write(File,"T = \"Time=$(LoopItterateCounter)\"\n");
      write(File,"VARIABLES = X,Y,Rho,U,V\n");
      write(File,"ZONE\n");
      write(File,"I = $(Nx1)\t J = $(Ny1)\n");
      write(File,"F=POINT\n");
      for c = 1:Nx1
          for r = 1:Ny1
              write(File,@sprintf("%d,%d,%lf,%lf,%lf\n",c,r,rho[r,c],ux[r,c],uy[r,c]));
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
    println("Ux , Vy and Rho has been saved in \"CSV\" file format.");
flush(STDOUT);
end #function[Data_output_CSV]



function CenterlineOutput(LoopItterateCounter::Int64, rho::Array{Float64}, ux::Array{Float64}, uy::Array{Float64})

    Xcenter = floor(Int64,Nx1/2);
    Ycenter = floor(Int64,Ny1/2);

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
    println("Start...\n");
    flush(STDOUT);
    #/==============================================================================
    #                             Loop variables
    #/==============================================================================

    LoopItterateCounter::Int64 = 0;
    err::Float64 = 1.0000000001;


    #/==============================================================================
    #                opening files for recording some data
    #/==============================================================================
    if DistributionFunctions_sum_Out == "yes"
        Distribution_Functions_sum = @sprintf("%s/A_Distribution_functions_Sum_%lf.csv",resultsAddress,LoopItterateCounter/PrintDecimal); # file name creation
        Filesum = open(Distribution_Functions_sum,"w"); # file open
    end

    if ErrorPrintingInInFilePermission == "yes"
      ErrorFileName = @sprintf("%s/A_Error_File_for_%s_Re%1.1lf_N%d_%d.csv",resultsAddress,MethodName,Re,Nx1,Ny1); # file name creation
      ErrorFilePointer = open(ErrorFileName,"w"); # file open
    end



    #/==============================================================================
    #                              Allocations
    #/==============================================================================
    # double tau , err
    f       = Array{Float64}(Ny1,Nx1,Q); # array of the distribution functions (DFs)
    f_post  = Array{Float64}(Ny1,Nx1,Q);

    rho     = Array{Float64}(Ny1,Nx1);   # arrays of fluid density and velocity
    ux      = Array{Float64}(Ny1,Nx1);
    uy      = Array{Float64}(Ny1,Nx1);

    u0      = Array{Float64}(Ny1,Nx1);   # for calculations in Err function
    v0      = Array{Float64}(Ny1,Nx1);

    M2      = floor(Int64,Ny/2);         # Center points for centerline
    N2      = floor(Int64,Nx/2);



    #/================= SRT Special Variables and allocations ======================

    tau::Float64     = 3*L*uw_top/Re + 0.5;      # relaxation time for BGK



    # ================= MRT Special Variables and allocations ======================

    m       = Array{Float64}(Q);                 # moments array
    s       = Array{Float64}(Q);                 # Relaxation rates array

    Lenght::Int64    = L # = N + 1;              # For halfway Bondary Condition

    tauMRT::Float64  = 3*Lenght*uw_top/Re + 0.5; # relaxation time for BGK

    s[7]    = 1.0/tauMRT;                        # relaxation rates for MRT
    s[8]    = 1.0/tauMRT;
    s[9]    = 0.0;
    s[3]    = 0.0;
    s[5]    = 0.0;
    s[4]    = 8*(2-s[7])/(8-s[7]);
    s[6]    = 8*(2-s[7])/(8-s[7]);
    s[1]    = 1.6;
    s[2]    = 1.8;


    #/==============================================================================
    #                   Initialization of rho, ux, uy and f
    #/==============================================================================

    rho, ux, uy, f = Init_Eq!(rho, ux, uy, f);

    # Debuging print
    #println("rho \n"); show(rho); println("ux \n"); show(ux); println("uy \n"); show(uy); println("f \n"); show(f);



    #/==============================================================================

    #/********************* Itterative loop starts here ****************************

    #/==============================================================================

    while (err > errorlimitation) || (LoopItterateCounter <= (ErrCalculationSteps*2 + 1)) # method 1
    # for it = 1:MaxItteration                                                            # method 3

        LoopItterateCounter = LoopItterateCounter + 1;



        #/==========================================================================
        #/======================== solving algorithem selection ====================
        #/==========================================================================



        # BGK Simple Boundary ======================================================

#        f_post = Coll_BGK!(rho, ux, uy, f, f_post, tau);               # BGK collision
#        f = Streaming!(f,f_post);                                      # Streaming for BGK
#        f = Bounce_back!(f,f_post,rho);                                # No-Slip boundary condition
#        rho, ux, uy = Den_vel!(rho,ux,uy,f);                           # Fluid variables calculation for all models



        # BGK Non_equilibrium_Boundarycondition ====================================

#         f_post = Coll_BGK!(rho, ux, uy, f, f_post, tau);              # BGK collision
#         f = Non_Equilibrium_Streaming!(f,f_post);                     # Streaming For Non_equilibrium_Boundarycondition
#         rho, ux, uy = Den_vel!(rho,ux,uy,f);                          # Fluid variables calculation for all models
#         f = Non_equilibrium_Boundarycondition!(f,rho,ux,uy);          # Non_equilibrium_Boundarycondition



        # BGK Non_equilibrium_Modified_Boundarycondition ============================

#         f_post = Coll_BGK!(rho, ux, uy, f, f_post, tau);              # BGK collision
#         f = Non_Equilibrium_Streaming!(f,f_post);                     # Streaming For Non_equilibrium_Boundarycondition
#         rho, ux, uy = Den_vel!(rho,ux,uy,f);                          # Fluid variables calculation for all models
#         f = Non_equilibrium_Boundarycondition_modified!(f,rho,ux,uy); # Non_equilibrium_Modified_Boundarycondition



        # BGK Semi_Non_equilibrium_Boundarycondition =================================

           #_____________________________________________________________
           #|                                                            |
           #|            Method is NOT WORKING Code works.               |
           #| Non_Equilibrium_Streaming is Better than Simple Streaming. |
           #| But the boundary at the bottom at which the velocities are |
           #| zero does not have any good results.                       |
           #|____________________________________________________________|


#          f_post = Coll_BGK!(rho, ux, uy, f, f_post, tau);               # BGK collision
# #          f = Streaming!(f,f_post);                                      # Streaming for BGK
#          f = Non_Equilibrium_Streaming!(f,f_post);                      # Streaming For Non_equilibrium_Boundarycondition
#          rho, ux, uy = Den_vel!(rho,ux,uy,f);                           # Fluid variables calculation for all models
#          f = Semi_Non_equilibrium_Boundarycondition!(f,rho,ux,uy);      # Non_equilibrium_Boundarycondition
#



        # MRT From LBM01.c Gue book ==================================================

#         f_post = Coll_MRT!(rho, ux, uy, f, f_post, m, s);             # BGK collision
#         f = Streaming!(f,f_post);                                     # Streaming for BGK
#         f = Bounce_back!(f,f_post,rho);                               # No-Slip boundary condition
#         rho, ux, uy = Den_vel!(rho,ux,uy,f);                          # Fluid variables calculation for all models



        # CMs. Non-Orthogonal Central Moments space ================================
        
        f_post = Coll_CMs_Non_Orthogonal!(rho, ux, uy, f, f_post, tau); # CMs Non-Orthogonal collision
        f      = Streaming!(f,f_post);                                  # Streaming for BGK
        f      = Bounce_back!(f,f_post,rho);                            # No-Slip boundary condition
        rho, ux, uy = Den_vel!(rho,ux,uy,f);                            # Fluid variables calculation for all models
        
        
        
        #/==========================================================================
        #/=========================== Post-processing ==============================
        #/==========================================================================


        if LoopItterateCounter % ErrCalculationSteps == 0
            err ,u0 ,v0 = Err!(ux, uy, u0, v0); # Velocity differences between two successive 1000 steps
            # println("Err1 = ",err); # Error with difference method
            # err ,u0 ,v0 = Err2!(ux, uy, u0, v0);
            println("err = ", err , " ux_center = ", ux[M2,N2] , " uy_center = ", uy[M2,N2] , " itt = ", LoopItterateCounter);  # Display some results
            flush(STDOUT);

            if ErrorPrintingInInFilePermission == "yes"
              write(ErrorFilePointer,"$(LoopItterateCounter)\t$(err)\n"); # writing to Filesum
            end # end if[ErrorPrintingInInFilePermission]

        end # end if[ErrCalculationSteps]



        if LoopItterateCounter % PrintSteps == 0 && PrintOutput == "yes"
            Data_output(LoopItterateCounter,rho,ux,uy);
        end # end if [PrintSteps]



        if LoopItterateCounter % CSVoutputFrequency == 0 && CSVoutput=="yes"
            Data_output_CSV(LoopItterateCounter,rho,ux,uy);
        end # end if[LoopItterateCounter % CSVoutputFrequency]



        if LoopItterateCounter % CenterlinePrintFrequency == 0 && CenterLinePrint == "yes"
            CenterlineOutput(LoopItterateCounter,rho,ux,uy);
        end # end if[LoopItterateCounter % CenterlinePrintFrequency]



        if LoopItterateCounter % DistributionFunctionsFrequency == 0 && DistributionFunctions_Out == "yes"
            DistributionFunctions_output_CSV(LoopItterateCounter,f);
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

          write(Filesum,"$(LoopItterateCounter)\t$(sum)\n"); # writing to Filesum
          flush(STDOUT);
        end # end of if [LoopItterateCounter % DistributionFunctionsFrequency]



        # println(LoopItterateCounter) # printing the LoopItterateCounter in the STDOUT



    end # end [while] or [for] (Itterative loop) function main still continues.



    #/==============================================================================
    #                       Closing of Opened file
    #/==============================================================================

    if DistributionFunctions_sum_Out == "yes"
      close(Filesum); # closing the file
    end



    if ErrorPrintingInInFilePermission == "yes"
      close(ErrorFilePointer); # closing the file
    end



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
