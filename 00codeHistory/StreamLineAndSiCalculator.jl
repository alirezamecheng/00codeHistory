using Pkg
if haskey(Pkg.installed(), "CSV") != true
    @warn("You dont have the \"CSV\" package installed. The program will install the package automatically.")
    Pkg.add("CSV");
end
@info("requirements: CSV package: pass");

if haskey(Pkg.installed(), "Printf") != true
    @warn("You dont have the \"Printf\" package installed. The program will install the package automatically.")
    Pkg.add("Printf");
end
@info("requirements: Printf package: pass");

if haskey(Pkg.installed(), "REPL") != true
    @warn("You dont have the \"REPL\" package installed. The program will install the package automatically.")
    Pkg.add("REPL");
end
@info("requirements: REPL package: pass");

@info("Loading CSV package");
using CSV
@info("Loading Printf package");
using Printf
@info("Loading REPL package");
using REPL.TerminalMenus



function fileSelectionAndInput()

    #uy = CSV.read("Vy0.002002.csv",header=0);
    #ux = CSV.read("Ux0.002002.csv",header=0);
    SizeMismach = 1;
    global ux, uy;
    while (SizeMismach == 1)
        allFilesInTheDirectoryName = readdir();
        selectableFiles = allFilesInTheDirectoryName[occursin.(".csv",allFilesInTheDirectoryName)];
        println("\nAll SELLECTABLE files found in the working directory : ", "\"", pwd(), "\"");
        for i=1:size(selectableFiles,1)
            println(rpad(selectableFiles[i],20,'.'),lpad(i,5,'.'));
        end

        
    #     menu = RadioMenu(selectableFiles, pagesize=5);
        menu = MultiSelectMenu(selectableFiles)
        choiceUx = collect(request("Choose the Ux file:", menu));
        while length(choiceUx) != 1
            choiceUx = 0;
            println("\nPlease select 1 file only.");
            choiceUx = collect(request("Choose the Ux file:", menu));
        end
        
        println("Your selected file is : ", selectableFiles[choiceUx[1]]);
        println("\n");
        
    #     menu = RadioMenu(selectableFiles, pagesize=5)
        menu = MultiSelectMenu(selectableFiles)
        choiceVy = collect(request("Choose the Vy file:", menu));
        while length(choiceVy) != 1
        choiceVy = 0;
            println("\nPlease select 1 file only.");
            choiceUx = collect(request("Choose the Vy file:", menu));
        end
        println("Your selected file is : ", selectableFiles[choiceVy[1]]);
        
        @info("Loading Data");
        global uy = CSV.read(selectableFiles[choiceVy[1]],header=0);
        global ux = CSV.read(selectableFiles[choiceUx[1]],header=0);

        Nuxi = size(ux,1);
        Nuxj = size(ux,2);
        Nuyi = size(uy,1);
        Nuyj = size(uy,2);
        if (Nuxi != Nuyi || Nuxj != Nuyj)
        @error("Size mismach in the input files. Please select again.");
        else
        SizeMismach = 0;
        end # end if 
    end # end while (SizeMismach)
    return ux , uy;
end # [function fileSelectionAndInput]


function si_cal!(si,v,dx = 1)
    @info("Calculating \"si\"\n")
    for i = 2:size(v,1)-1
        for j = 1:size(v,2)
            si[i+1,j] = si[i-1,j] - dx*(v[i-1,j] + 4*v[i,j] + v[i+1,j])/6
        end
    end
    return si
end #[function si_cal!]

function StreamLines!(st,u,dx = 1)
    @info("Calculating \"Stream lines\"\n")
    
#     for j = 1:size(u,2)
# #         st[1,j] = 0 - dx*(0 + 4*u[1,j] + u[1+1,j])/6
#         st[1,j] = st[i-1,j] - dx*(0 + 4*u[1,j] + u[1+1,j])/6
#     end
    for i = 2:size(u,1)-1
        for j = 1:size(u,2)
            st[i+1,j] = st[i-1,j] - dx*(u[i-1,j] + 4*u[i,j] + u[i+1,j])/6
        end
    end
    return st
end #[function StreamLines!]


function Worticity()

end


# # Zadehgol formulation # dosent work!!!!
# function StreamLines!(st,u,dx = 1,dy =1)
#     for i=1:size(u,1)
#         st[i,1] = 0.0;
#         st[i,size(u,2)] = 0.0;
#     end
#     for j = 1:size(u,2)
#         st[1,j]=0.0;
#         st[size(u,1),j] = 0.0;
#     end
#
#     for i = 2:size(u,1)-1
#         for j = size(u,2)-1:-1:2
#             st[i,j] = st[i-1,j] - (u[i-1,j]*dx - u[i-1,j]*dy); # Zadehgol formulation
#         end
#     end
#     return st
# end #[function StreamLines!]




function Data_output(si::Array{Float64},st::Array{Float64}, ux, uy)

    print("Please specify the output file name : ");
    OutputFileName = strip(readline());
    while length(OutputFileName) == 0
        println("The file name can't be empty")
        print("Please specify the output file name : ");
        OutputFileName = strip(readline());
    end
    
    FileNameAndAddress = @sprintf("./%S.tp",OutputFileName)
    File = open(FileNameAndAddress,"w");
      #println(File,"TITLE=\"",MethodName," Re = ",Re,"\"");
      # write(File,"SOLUTIONTIME = $(LoopItterateCounter)\n");
      # write(File,"T = \"Time=$(LoopItterateCounter)\"\n");
      write(File,"VARIABLES = X,Y,StreamLines,si,U,V\n");
      write(File,"ZONE\n");
      println(File,"I = ",size(uy,1),"\t","J = ",size(uy,2));
      write(File,"F=POINT\n");
    for r = 1:size(uy,2)
        for c = 1:size(uy,1)
            # write(File,@sprintf("%d,%d,%lf,%lf,%lf\n",c,r,rho[r,c],ux[r,c],uy[r,c]));
            println(File,c,",",r,",",st[r,c],",",si[r,c],",",ux[r,c],",",uy[r,c]);
        end
    end
    close(File);
    println("Results has been saved at : ",FileNameAndAddress);
#     flush(stdout);
end

function main()



    ux , uy = fileSelectionAndInput();

    # allocation
    si = zeros(size(uy,1),size(uy,2));
    st = zeros(size(uy,1),size(uy,2));

    # calculations
    si = si_cal!(si,uy);
    st = StreamLines!(st,ux);

    # Data output
    Data_output(si,st,ux,uy);

end #[function main]



main();
println("finished");
