
using PyPlot
using StatsBase 
using Optim 
using DataFrames
PyPlot.close("all")
include("$(pwd())/RCLModule.jl")
include("$(pwd())/ForceManager.jl")
  #  Simulate the first encounter time between two ends of a break in a RCL polymer

  # include("RCLSimulationParams.jl")
function RunRCLBreakSimulation(params)
      
    #  # WORKING IN UNITS OF NM  

    params["numSimulationBatches"] = length(params["numRandomConnectors"])
    results = InitializeResultStruct(params)
   
    #---- Main loop-----
    for sbIdx in 1:length(params["numRandomConnectors"])              
           
      for smIdx in 1:params["numSimulations"]
            println("__________________________________")            
            println("Batch $(sbIdx) Simulation $(smIdx) :")
            
            # Initialize the polymer chain   
            curPos   = RCL.InitializeChain(params["numMonomers"],params["connectorsSTD"],params["diffusionConst"],params["dimension"],params["dt"])
             
            # Laplacian of the graph representation of the polymer
            # only consider configuration for which the chain is not split into two            
            laplacian = GetUnsplitLaplacian(params["numMonomers"],params["numRandomConnectors"][sbIdx],params["fetForMonomers"])

            # Determine the number of relaxation steps according to the laplacian            
            DetermineRelaxationTime!(laplacian,params)
            
            #---------- Run until polymer relaxation time-------------------
            RCL.Step!(params["numRelaxationSteps"],curPos,laplacian,params["springConst"],
                      params["dimension"],params["dtRelaxation"],params["diffusionConst"])
            
            # Disconnect the monomers specified in params["monomersToDetach"] after relaxation steps
               laplacian = InduceDSB(laplacian, params)   

                # Record polymer connectivity after DSB
                connectedMonomersDSB                         = RCL.GetConnectivity(laplacian)               
                numConnectorsDSB                             = length(connectedMonomersDSB[1])
                results["numConnectorsRemoved"][smIdx,sbIdx] = params["numRandomConnectors"][sbIdx]-numConnectorsDSB
                println("Number of cross-links removed: $(results["numConnectorsRemoved"][smIdx,sbIdx])")
             

            # Record connected monomers after the break (in case of a break)
            cm     = RCL.GetConnectivity(laplacian)
            results["connectedMonomers"][sbIdx,smIdx]= [cm[1] cm[2]]

            #------ Start simulation after the break to make the ends drift apart --------

            # Preallocation 
            chainPosSimulation = zeros(params["numMonomers"],params["dimension"],params["numSteps"]);            
            maxDist       = 0;
            initialFETPos = zeros(params["numMonomers"],params["dimension"])# initial chain position for FET simulations 

            for stIdx =1:params["numSteps"] 
                        RCL.Step!(1,curPos,laplacian,params["springConst"],
                        params["dimension"],params["dt"],params["diffusionConst"])
                if params["mechanicalSpringPointForce"]
                        ForceManager.MechanicalSpringPointForce!(curPos,curPos,"out",params["dimension"]*params["diffusionConst"]/params["connectorsSTD"]^2,0.0,params["encounterDist"],params["dt"])#-1 in, 1 out                      
                end
                chainPosSimulation[:,:,stIdx]= copy(curPos);
            
              # Record max distance between fetForMonomers 
              mDist  = sqrt(sum((curPos[params["fetForMonomers"][1],:]-curPos[params["fetForMonomers"][2],:]).^2))
              if mDist>maxDist
                 maxDist    = mDist
                 for n in eachindex(curPos)
                        initialFETPos[n] = curPos[n] # record position at maximal dist (after the break)
                 end
              end# if 
            end
            
            # Record max distance
            results["maxMonomerDist"][sbIdx,smIdx] = maxDist
            # Calculate the apparent spring constant
            results["Kc"][sbIdx,smIdx] = CalculateK_c(chainPosSimulation,params["diffusionConst"],params["dt"])            
            # Calculate L_c         
            results["Lc"][sbIdx,smIdx] = CalculateL_C(chainPosSimulation)          
  
            # Calculate radius of gyration at the end of simulation             
            results["radiusOfGyration"][sbIdx,smIdx] = CalculateRadiusOfGyration(curPos)
            println("ROG: $(results["radiusOfGyration"][sbIdx,smIdx]) nm")
                                                                     
             #----- Start FET simulation---------------
            if params["fetSimulation"]
               # Show snapshot of the chain at the begining of FET simulation 
               ShowChainSnapshot(curPos,laplacian,params)

                if params["fetSimulationStartFromMaxDistance"]             
                    # strat from maximal distance betwee monomers
                    # Copy the initial position into curPos
                    for n in eachindex(initialPos)
                        curPos[n]  = initialFETPos[n] # start from the inital pos                
                    end
                end            
                
                runCondition  = true # condition to run the simulation
                stepIdx       = 1 # keep this index defined outside the while loop
                while runCondition                   
                    RCL.Step!(1,curPos,laplacian,params["springConst"], params["dimension"], params["dt"], params["diffusionConst"])

                    if params["mechanicalSpringPointForce"]
                    ForceManager.MechanicalSpringPointForce!(curPos,curPos,"out",params["springConst"],0.0,params["encounterDist"],params["dt"])#-1 in, 1 out                      
                    end
                                

                    runCondition    = norm(curPos[params["fetForMonomers"][1],:]-curPos[params["fetForMonomers"][2],:],2)>params["encounterDist"] 
                    runCondition    = (runCondition & (stepIdx<params["maxFETSteps"]))
                    # Show simulations                              
                    ShowSimulation(curPos,params)

                    stepIdx += 1                     
                end #  FET simulation    

                # Record First encounter time
                totalFETSteps = stepIdx              
                results["firstEncounterTime"][sbIdx,smIdx] = totalFETSteps*params["dt"]
                println("FET: $(results["firstEncounterTime"][sbIdx,smIdx]) seconds.")
                 # Record last chain position
                results["lastChainPosition"][sbIdx,smIdx]= copy(curPos)

                # Show chain configuration at encounter time 
                ShowChainSnapshot(curPos,laplacian, params)
            end

      end# for each simulation 
        
        # Compute the mean first encounter time        
        results["MFET"][sbIdx] = mean(results["firstEncounterTime"][sbIdx,:])  # MFET[sbIdx]

        # computer the mean square radius of gyration        
        results["MSRG"][sbIdx] = mean(results["radiusOfGyration"][sbIdx,:].^2) # MSRG[sbIdx]

     if params["fitHistToData"] # the FET histogram 
        results["histData"][sbIdx,:], results["bins"][sbIdx,:] = GetHistogram(results["firstEncounterTime"][sbIdx,:],params["numHistBins"])
        results["fitParams"][sbIdx,:]= FitExponentToHistData(results["histData"][sbIdx,:],results["bins"][sbIdx,:],params["numFitTrials"])        
     end# if 

      # Plot histogram and fit
        if params["histPlotFlag"]
            lineWidth = 3
            lineColor = "red"
            fontSize  = 20
            figTitle  = " Experiment $(sbIdx), num. connectors= $(params["numRandomConnectors"][sbIdx])"
            PlotHistAndFit(results["histData"][sbIdx,:],results["bins"][sbIdx,:],
                           results["fitParams"][sbIdx,:],fontSize,lineColor,lineWidth,figTitle)
        end# if

      results["averageMaxMonomerDist"][sbIdx] = mean(results["maxMonomerDist"][sbIdx,:])
       

    end# for, each batch 
      CalculateMeanLc!(results)
      

       PlotEncounterTimeVsConnectivity(results,params)
       PlotMaxMonomerDistanceVsConnectivity(results,params)
       # Export all results 
       ExportResults(results)
      
    return results
end# function

function InitializeResultStruct(params)
    # initialize result dicionary 
        results = Dict([("firstEncounterTime",zeros(params["numSimulationBatches"],params["numSimulations"])),
                    ("radiusOfGyration",zeros(params["numSimulationBatches"],params["numSimulations"])),
                    ("MSRG",zeros(params["numSimulationBatches"])),
                    ("MFET",zeros(params["numSimulationBatches"])),
                    ("fitParams",zeros(params["numSimulationBatches"],2)),
                    ("histData",zeros(params["numSimulationBatches"],params["numHistBins"])),
                    ("bins",zeros(params["numSimulationBatches"],params["numHistBins"])),
                    ("maxMonomerDist",zeros(params["numSimulationBatches"], params["numSimulations"])),
                    ("averageMaxMonomerDist", zeros(params["numSimulationBatches"])),
                    ("position", []),  
                    ("numConnectorsRemoved", zeros(params["numSimulations"], params["numSimulationBatches"])),
                    ("connectorsRemoved",Array{Array{Float64,2},2}(params["numSimulationBatches"],params["numSimulations"])),
                    ("lastChainPosition",Array{Array{Float64,2},2}(params["numSimulationBatches"],params["numSimulations"])),#zeros(params["numMonomers"],params["dimension"],params["numSimulationBatches"],params["numSimulations"])), 
                    ("connectedMonomers",Array{Array{Int64,2},2}(params["numSimulationBatches"],params["numSimulations"])),                                      
                    ("Lc", Array{Array{Float64,1},2}(params["numSimulationBatches"],params["numSimulations"])),
                    ("meanLc",zeros(params["numMonomers"],params["numSimulationBatches"])),
                    ("Kc", Array{Array{Float64,2},2}(params["numSimulationBatches"],params["numSimulations"])),                                        
                    ("params",params)])
return results
end

function InduceDSB(laplacian, params)
    # induce DSB between specified monomers 
        if params["simulateBreak"]
                if params["induceUniformDSB"]
                    # override monomersToDetach 
                    params = UniformDistributedPairs(params)
                end

                for m2dIdx =1:size(params["monomersToDetach"],1) # for all pairs
                    println("DSB between monomers $(params["monomersToDetach"][m2dIdx,1]) and $(params["monomersToDetach"][m2dIdx,2])")
                  RCL.DisconnectMonomers!(params["monomersToDetach"][m2dIdx,1],params["monomersToDetach"][m2dIdx,2], laplacian) # detach pairs return laplacian
                    if params["breakAllConnectorsFromDetachedMonomers"] # break all random connectors 
                        println("Breaking all cross-links to monomers $(params["monomersToDetach"][m2dIdx,1]) and $(params["monomersToDetach"][m2dIdx,2])")
                        mtdetach = params["monomersToDetach"][m2dIdx,:]
                        for mtdIdx = 1:length(mtdetach)                            
                            for mIdx in [1:mtdetach[mtdIdx]-2 ; mtdetach[mtdIdx]+2: params["numMonomers"]]
                               RCL.DisconnectMonomers!(mtdetach[mtdIdx],mIdx,laplacian)# detach random connectors return laplacian 
                            end
                        end
                    end
                end
        end

return laplacian 
end

function UniformDistributedPairs(params)
     # Writes monomersToDetach with numberOfDSB uniformly chosen pairs 
     r = randperm(params["numMonomers"]-1);
     params["monomersToDetach"]= zeros(Int64,params["numberOfDSB"],2)
     for rIdx in 1:params["numberOfDSB"]
      params["monomersToDetach"][rIdx,:]= [r[rIdx] r[rIdx]+1]
     end
     println(params["monomersToDetach"])
    return params
end


function GetUnsplitLaplacian(numMonomers::Int64,numRandomConnectors::Int64,monomersToDetach)
    # Generate a laplacian such that the after the induction of a DSB, the laplacian is not splot into two seperate chains 
     laplacian = zeros(numMonomers,numMonomers)
      # Get connectivity fraction
     ξ  = RCL.Connectors2Connectivity(numRandomConnectors,numMonomers)
     lCondition = true
            while lCondition
              laplacian = RCL.GraphSemiRandomLaplacian(numMonomers,ξ)
              cMonomers = RCL.GetConnectivity(laplacian)
              for c in 1:length(cMonomers[1])
                  if cMonomers[1][c]<monomersToDetach[1] & cMonomers[2][c]>monomersToDetach[2]
                      lCondition = false
                  elseif cMonomers[1][c]>monomersToDetach[2] & cMonomers[2][c]<monomersToDetach[1]
                       lCondition = false 
                  end
              end
            end# while            
return laplacian
end    

function DetermineRelaxationTime!(laplacian,params)
        # Set Relaxation steps     
            if params["calculateRelaxationTime"]
                println("Determining num. relaxation steps automatically.")
                # calculate the relaxation time by tau = 1/(springConst λ) with the spring const ant  λ the smallest non vanishing eigen value  of the laplacian of the chain
                relaxationSteps = RCL.RelaxationTime(laplacian, params["springConst"])
                println(relaxationSteps)
                if relaxationSteps !=0
                    params["numRelaxationSteps"] = relaxationSteps
                    else
                        # keep default values 
                end
            end
return params
end

function ShowChainSnapshot(curPos,laplacian,params)
        if params["showChainSnapshot"]
            PyPlot.figure()
            connectedPairs = RCL.GetConnectivity(laplacian)         
            figTitle       = "Simulation" 
            PlotChain(curPos,params["fetForMonomers"],connectedPairs,figTitle,"green","red")
        end# if
return nothing 
end

function ShowSimulation(curPos,params)
     if params["showSimulation"]
        PyPlot.cla()
        if params["dimension"]==3
            PyPlot.plot3D(curPos[:,1],curPos[:,2],curPos[:,3],marker="o")
        elseif params["dimension"]==2
            PyPlot.plot(curPos[:,1],curPos[:,2],marker="o")
            end
        PyPlot.pause(0.000001)
    end
return nothing
end

function CalculateRadiusOfGyration(pos)
     # Radius of gyration (not squared)
      centerMass = mean(pos,1)      
      rog        = sqrt(sum(broadcast(-,pos,centerMass).^2,2)/size(pos,1))
return rog[1]
end

function CalculateL_C(chainPos)# for each simulation 
   # the std of monomer position          
     Lc = sqrt(sum(var(chainPos,3),2))# variance of monomer position
     Lc = Lc[:,1,1]

    return Lc
    
end

function CalculateMeanLc!(results)
    for bIdx = 1:results["params"]["numSimulationBatches"]
        # println(bIdx)
        lc = zeros(results["params"]["numMonomers"],results["params"]["numSimulations"])
          for sIdx = 1:results["params"]["numSimulations"]
               lc[:,sIdx]= results["Lc"][bIdx,sIdx]
          end
          results["meanLc"][:,bIdx] = mean(lc,2)
          #println(results["meanLc"][:,bIdx])
        end
    return results
end

function CalculateK_c(chainPos,D,dt)
    # compute the apparent spring constant   
   # unpack chainPos

    numMonomers = size(chainPos,1)
    dimension   = size(chainPos,2)
    numSteps    = size(chainPos,3)
    cm          = mean(chainPos,3)# average position  
     

    kc  = zeros(numMonomers)
     for sIdx in 1:numSteps-1
       kc+=(1/(D*dt)).*sum((chainPos[:,:,sIdx].-chainPos[:,:,sIdx+1])./(chainPos[:,:,sIdx].- cm),2)

     end

    kc = kc./(dimension*(numSteps-1))
    kc = kc[:,:,1]
    return kc    
end

function ExportResults(results)
 # Save data as DataFrame tables and export as csv

  if results["params"]["saveResults"]
     resultFolder = "$(results["params"]["resultBaseFolder"])$(results["params"]["resultFolderName"])"

    if !isdir(results["params"]["resultBaseFolder"])        
        mkdir(results["params"]["resultBaseFolder"])       
    end

    if !isdir(resultFolder)
         mkdir(resultFolder)
    end

    fetTable               = DataFrame()
    histTable              = DataFrame()
    mfetTable              = DataFrame()
    averageMaxDists        = DataFrame()
    maxDistTable           = DataFrame()
    maxMonomerDistTable    = DataFrame()
    msrgTable              = DataFrame()   
    lcTable                = DataFrame()
    lcSTDTable             = DataFrame()
    kcTable                = DataFrame() 
    connectorsRemovedTable = DataFrame()
    radiusOfGyrationTable  = DataFrame()
    numConnectorsRemovedTable = DataFrame()

    for bIdx in 1:length(results["params"]["numRandomConnectors"])# for each simulation batch
            # for each simulation 
            str2 = "bins $(results["params"]["numRandomConnectors"][bIdx]) connectors "
            histTable[Symbol(str2)] = results["bins"][bIdx,:]        
                
            str3 = "hist $(results["params"]["numRandomConnectors"][bIdx]) connectors "            
            histTable[Symbol(str3)] = results["histData"][bIdx,:]

            str6 = "$(results["params"]["numRandomConnectors"][bIdx]) connectors "
            numConnectorsRemovedTable[Symbol(str6)]= results["numConnectorsRemoved"][:,bIdx]           
               lcSimulation = zeros(results["params"]["numMonomers"],results["params"]["numSimulations"])
               kcSimulation = zeros(results["params"]["numMonomers"],results["params"]["numSimulations"])

                for sIdx in 1:results["params"]["numSimulations"]
                    lcSimulation[:,sIdx] = results["Lc"][bIdx,sIdx]
                    kcSimulation[:,sIdx] = results["Kc"][bIdx,sIdx]
                end
            
                str4 = "mean Lc[nm] $(results["params"]["numRandomConnectors"][bIdx]) connectors "
                lcTable[Symbol(str4)] =  results["meanLc"][:,bIdx] #mean(lcSimulation,2)[:,1] 

                str7 = "STD Lc[nm] $(results["params"]["numRandomConnectors"][bIdx]) connectors "
                lcSTDTable[Symbol(str7)] = std(lcSimulation,2)[:,1]       
 
                str5 = "Kc $(results["params"]["numRandomConnectors"][bIdx]) connectors "            
                kcTable[Symbol(str5)] = mean(kcSimulation,2)[:,1]
               

   end


        str1 = "numConnectors"
        radiusOfGyrationTable[Symbol(str1)] = results["params"]["numRandomConnectors"]
        maxMonomerDistTable[Symbol(str1)]   = results["params"]["numRandomConnectors"]
        fetTable[Symbol(str1)]              = results["params"]["numRandomConnectors"]
        
        for sIdx = 1:results["params"]["numSimulations"]
            str1 = "simulation $sIdx"
            radiusOfGyrationTable[Symbol(str1)] = results["radiusOfGyration"][:,sIdx].^2
            maxMonomerDistTable[Symbol(str1)]   = results["maxMonomerDist"][:,sIdx]
            fetTable[Symbol(str1)]              = results["firstEncounterTime"][:,sIdx]
       end

            mfetTable[Symbol("numConnectors")] = results["params"]["numRandomConnectors"]
            mfetTable[Symbol("MFET(ms)")]      = results["MFET"]

            msrgTable[Symbol("numConnectors")] = results["params"]["numRandomConnectors"]
            msrgTable[Symbol("MSRG")]          = results["MSRG"]

            averageMaxDists[Symbol("numConnectors")] = results["params"]["numRandomConnectors"]
            averageMaxDists[Symbol("MaxDist(nm)")]   = results["averageMaxMonomerDist"]

            # Save csv  files
            writetable("$resultFolder/FirstEncounterTime.csv", fetTable)
            writetable("$resultFolder/Histograms.csv", histTable)
            writetable("$resultFolder/MFET.csv",mfetTable)
            writetable("$resultFolder/averageMaxDistance.csv",averageMaxDists)
            writetable("$resultFolder/MSRG.csv", msrgTable)
            writetable("$resultFolder/maxMonomerDistance.csv", maxMonomerDistTable)
            writetable("$resultFolder/radiusOfGyration.csv", radiusOfGyrationTable)
            writetable("$resultFolder/meanLC.csv", lcTable)
            writetable("$resultFolder/stdLc.csv",lcSTDTable)
            writetable("$resultFolder/meanKC.csv", kcTable)
            fpTable = DataFrame(λ=results["fitParams"][:,2],a=results["fitParams"][:,1])
            writetable("$resultFolder/FitParams.csv", fpTable)
            writetable("$resultFolder/numConnectorsRemoved.csv",numConnectorsRemovedTable)
            
            # export parameters as a txt file 
            open("$resultFolder/parameters.txt","w") do fid
            for kName in keys(results["params"])
                p = results["params"]["$kName"]
                write(fid, "$kName : $p \n" )
           end
  end

    # export characteristic chain configuration as the last chain position of the last simulation in each batch
    # create chain configuration folder 
    chainPositionFolder = "$resultFolder/lastChainConfiguration"
    
        if !isdir(chainPositionFolder)
            mkdir("$chainPositionFolder")
        end

        polymerPDBFolder    = "$(chainPositionFolder)/PDB"

      if !isdir(polymerPDBFolder)
            mkdir("$polymerPDBFolder")
      end

    for bIdx in 1:length(results["params"]["numRandomConnectors"])# for each simulation batch
        chainConfigurationTable              = DataFrame()
        chainConfigurationTable[Symbol("x")] = results["lastChainPosition"][bIdx,end][:,1]
        chainConfigurationTable[Symbol("y")] = results["lastChainPosition"][bIdx,end][:,2]
        chainConfigurationTable[Symbol("z")] = results["lastChainPosition"][bIdx,end][:,3]
        
        # export to csv
        chainPosFileName = "$chainPositionFolder/chainPos_$(results["params"]["numRandomConnectors"][bIdx])_connectors.csv" 
        writetable(chainPosFileName,chainConfigurationTable)
        
        # export connected monomer list for each batch
        cmFileName = "$chainPositionFolder/connectedMonomers_$(results["params"]["numRandomConnectors"][bIdx])_connectors.csv" 
        connectedMonomersTable = DataFrame()  
        connectedMonomersTable[Symbol("monomer1")] = results["connectedMonomers"][bIdx,end][:,1]
        connectedMonomersTable[Symbol("monomer2")] = results["connectedMonomers"][bIdx,end][:,2]  
        writetable(cmFileName,connectedMonomersTable)

        pdbFileName  = "$polymerPDBFolder/chainConfiguration_$(results["params"]["numRandomConnectors"][bIdx])_connectors.pdb" 
        # export a pdb file from last polyer configuration
        numMonomers              = results["params"]["numMonomers"]
        chainConfigurationStruct = Dict([("X",results["lastChainPosition"][bIdx,end][:,1]),
                                    ("Y",results["lastChainPosition"][bIdx,end][:,2]),
                                    ("Z",results["lastChainPosition"][bIdx,end][:,3]),
                                    ("outFile", pdbFileName),
                                    ("recordName",fill("ATOM",numMonomers)),
                                    ("atomNum",1:numMonomers),
                                    ("atomName",fill("N",numMonomers)),
                                    ("altLoc",fill(' ',numMonomers)),
                                    ("resName",fill("HIS",numMonomers)),
                                    ("chainID",fill("A",numMonomers)),
                                    ("resNum",fill(1,numMonomers)),
                                    ("occupancy",fill(1,numMonomers)),
                                    ("betaFactor",fill(0,numMonomers)),
                                    ("element",fill("N",numMonomers)),
                                    ("charge",fill(' ', numMonomers)),
                                    ("connectedMonomers",results["connectedMonomers"][bIdx,end])])
        Polymer2PDB(chainConfigurationStruct)
      end
 end

end# function

function Polymer2PDB(input)
    #  -- Adapted from  mat2PDB --
    # 
    #  This function creates a PDB from coordinate data. 
    #  -- required inputs (3) --
    # 
    #  input value        meaning
    # 
    #  input.X            orthagonal X coordinate data (angstroms)
    #  input.Y            orthagonal Y coordinate data (angstroms)
    #  input.Z            orthagonal Z coordinate data (angstroms)
    # 
    #  -- optional inputs (12): generates defaults when not user-specified --
    # 
    #  input value        meaning                           default value
    # 
    #  input.outfile      output file name                 "mat2PDB.pdb"
    #  input.recordName   output record name of atoms      "ATOM"
    #  input.atomNum      atom serial number                sequential number
    #  input.atomName     name of atoms                    "OW" (water oxygen)
    #  input.altLoc       alt. location indicator          " "
    #  input.resName      name of residue                  "SOL" (water)
    #  
    #  input.chainID      protein chain identifier         "A"
    #  input.resNum       residue sequence number           sequential number
    #  input.occupancy    occupancy factor                 "1.00"
    #  input.betaFactor   beta factor, temperature         "0.00"
    #  input.element      element symbol                   "O" (oxygen)
    #  input.charge       atomic charge                    " "
    # 
    # 

    #  create PDB

    #  open file    
    FILE = open(input["outFile"], "w");

    #  output data
    for n = 1:length(input["atomNum"])
    # 1 -  6        Record name     "ATOM  "                                            
    # 7 - 11        Integer         Atom serial number.                   
    # 13 - 16        Atom            Atom name.                            
    # 17             Character       Alternate location indicator.         
    # 18 - 20        Residue name    Residue name.                         
    # 22             Character       Chain identifier.                     
    # 23 - 26        Integer         Residue sequence number.              
    # 27             AChar           Code for insertion of residues.       
    # 31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.                       
    # 39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.                            
    # 47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.                            
    # 55 - 60        Real(6.2)       Occupancy.                            
    # 61 - 66        Real(6.2)       Temperature factor (Default = 0.0).                   
    # 73 - 76        LString(4)      Segment identifier, left-justified.   
    # 77 - 78        LString(2)      Element symbol, right-justified.      
    # 79 - 80        LString(2)      Charge on the atom.  
    
    # Find if the monomer is connected (random spring connecter)
    monomerElement = input["element"][n]
        for mIdx=1:size(input["connectedMonomers"][:,1])[1]
            if input["connectedMonomers"][mIdx,1]==n 
                input["element"][n] = "O" # Nitrogen
                input["atomName"][n] = "O"
            end

            if input["connectedMonomers"][mIdx,2]==n
                input["element"][n] = "O" # Nitrogen
                input["atomName"][n] ="O"
            end
        end        

        @printf(FILE,"%-6s%5u%3s%3.1s%3s %1.1s%4u%12.3f%8.3f%8.3f%6.2f%6.2f%12s%2s\n",
            input["recordName"][n], input["atomNum"][n], input["atomName"][n], 
            input["altLoc"][n],    input["resName"][n], input["chainID"][n], 
            input["resNum"][n],    input["X"][n]/100, input["Y"][n]/100, input["Z"][n]/100,
            input["occupancy"][n], input["betaFactor"][n], 
            input["element"][n],   input["charge"][n])
    
    end

    #  prepare a string with connectivity information for each monomer
    # Create the linear backbone 
    connectivity = Array{Array{Int64,1},1}(length(input["X"]))
    connectivity[1] = [1, 2]
    for cIdx in 2:length(input["X"])-1
        connectivity[cIdx] = [cIdx, cIdx-1,cIdx+1]
    end

        # COLUMNS         DATA TYPE        FIELD           DEFINITION
        # ---------------------------------------------------------------------------------
        #  1 -  6         Record name      "CONECT"
        #  7 - 11         Integer          serial          Atom serial number
        # 12 - 16         Integer          serial          Serial number of bonded atom
        # 17 - 21         Integer          serial          Serial number of bonded atom
        # 22 - 26         Integer          serial          Serial number of bonded atom
        # 27 - 31         Integer          serial          Serial number of bonded atom
        # 32 - 36         Integer          serial          Serial number of hydrogen bonded atom                                                
        # 37 - 41         Integer          serial          Serial number of hydrogen bonded atom
        # 42 - 46         Integer          serial          Serial number of salt bridged atom
        # 47 - 51         Integer          serial          Serial number of hydrogen bonded atom 
        # 52 - 56         Integer          serial          Serial number of hydrogen bonded atom 
        # 57 - 61         Integer          serial          Serial number of salt bridged atom

        # define the last entry 

    connectivity[length(input["X"])] = [length(input["X"]), length(input["X"])-1,length(input["X"])+1]
    # Add extra connectors
    for cIdx in 1:size(input["connectedMonomers"],1)
        push!(connectivity[input["connectedMonomers"][cIdx,1]],input["connectedMonomers"][cIdx,2])
    end


    for cIdx = 1:length(connectivity)
        @printf(FILE,"CONECT ")
        for sIdx = 1:length(connectivity[cIdx])
            cm   = connectivity[cIdx][sIdx]
            @printf(FILE,"%4.0d",cm)
        end
        @printf(FILE,"\n")
    end

    print( FILE, "END\n")

    #  Close file
    @printf("   %6.2f%%\n    done! closing file...\n", 100)

    close(FILE)

end

function GetHistogram(data,numBins::Int64)
    # calculate the histogram with numBins
    h        = StatsBase.fit(Histogram,data,linspace(0,maximum(data),numBins+1) )
    println(data)
    bins     = linspace(h.edges[1][1],h.edges[1][end],numBins)
    histData = h.weights/sum(h.weights)
    
    return histData, bins

end

function FitExponentToHistData(histData::Array{Float64,1},bins::Array{Float64},numFitTrials::Int64)

 # Fit an exponent function to histogram data histData with bins
 # Take the best fit out of numFitTrials
    optimFunVal = zeros(numFitTrials)
    f           = []
    # firstInd    = 1
    validInds   = histData.>0
    modelFit(x) = 0.5*sum((x[1].*exp(-x[2].*bins[validInds]) -histData[validInds]).^2)

    for nIdx in 1:numFitTrials   
        push!(f,optimize(modelFit,rand(2), BFGS()))
        optimFunVal[nIdx]=  f[end].minimum
    end

  # get the parameters for the best fit
  optimParams =Optim.minimizer(f[findfirst(optimFunVal.==minimum(optimFunVal))])
  return optimParams
end

function PlotHistAndFit(histData::Array{Float64,1},bins::Array{Float64},fitParams::Array{Float64},
                        fontSize::Int64,lineColor::String,lineWidth::Int64,
                        figTitle::String)
        PyPlot.figure()
        t         = linspace(bins[1],bins[end],80)
        fitString = @sprintf "%1.2f exp(-%1.2f t)" fitParams[1] fitParams[2]
        PyPlot.bar(bins,histData./sum(histData),width=(bins[2]-bins[1]),edgecolor=[1 1 1])
        PyPlot.plot(t,fitParams[1].*exp(-t.*fitParams[2]),linewidth=lineWidth,color=lineColor)
        PyPlot.legend([fitString, "data"],fontsize=fontSize,framealpha=0)
        PyPlot.xlabel("First encounter time [sec]",fontsize=fontSize)
        PyPlot.ylabel("Probability",fontsize=fontSize)
        PyPlot.title(figTitle,fontsize=fontSize)
end

function PlotEncounterTimeVsConnectivity(results,params)
    if params["plotEncounterTimeVsConnectivity"]
       figTitle = "MFET Vs. Connectivity "
       PyPlot.figure()
       PyPlot.plot(params["numRandomConnectors"],results["MFET"])
       PyPlot.xlabel("Number of connectors")
       PyPlot.ylabel("MFET [sec]")
    end

end

function PlotMaxMonomerDistanceVsConnectivity(results,params)
    if params["plotMaxMonomerDistanceVsConnectivity"]
        PyPlot.figure()
        PyPlot.plot(params["numRandomConnectors"], results["averageMaxMonomerDist"])
        PyPlot.xlabel("Number of random connectors")
        PyPlot.ylabel("Mean maximal intermonomer distance")
    end
end

function PlotChain(curPos::Array{Float64,2},highlightMonomers::Array{Int64},connectedPairs::Tuple{Array{Int64,1},Array{Int64,1}},
                   figTitle::String,lineColor::String,highlightColor::String)
       

         PyPlot.title(figTitle)
         if size(curPos,2)==2 # for 2D plot
            PyPlot.plot(curPos[:,1],curPos[:,2],color=lineColor,linewidth=3)
            PyPlot.plot(curPos[highlightMonomers[1],1],curPos[highlightMonomers[1],2],marker="o",color=highlightColor,markersize=10)
            PyPlot.plot(curPos[highlightMonomers[2],1],curPos[highlightMonomers[2],2],marker="o",color=highlightColor,markersize=10)

            # Plot extra connectors
            for cIdx  in 1:length(connectedPairs[1])
                x1Ind = connectedPairs[1][cIdx]
                x2Ind = connectedPairs[2][cIdx]
                PyPlot.plot([curPos[x1Ind,1],curPos[x2Ind,1]],
                            [curPos[x1Ind,2],curPos[x2Ind,2]],color="k")
            end

          elseif size(curPos,2)==3 # for 3D plot
           
            PyPlot.plot3D(curPos[:,1],curPos[:,2],curPos[:,3],color=lineColor,linewidth=2)
            PyPlot.plot3D(curPos[highlightMonomers[1],1]*ones(2),
                          curPos[highlightMonomers[1],2]*ones(2),
                          curPos[highlightMonomers[1],3]*ones(2),
                          marker="o",color=highlightColor,markersize=10)
            PyPlot.plot3D(curPos[highlightMonomers[2],1]*ones(2),
                          curPos[highlightMonomers[2],2]*ones(2),
                          curPos[highlightMonomers[2],3]*ones(2),
                          marker="o",color=highlightColor,markersize=10)
            # Plot added connectors
            for cIdx  in 1:length(connectedPairs[1])
                x1Ind = connectedPairs[1][cIdx]
                x2Ind = connectedPairs[2][cIdx]

                PyPlot.plot3D([curPos[x1Ind,1],curPos[x2Ind,1]],
                              [curPos[x1Ind,2],curPos[x2Ind,2]],
                              [curPos[x1Ind,3],curPos[x2Ind,3]],
                              color="g")
            end
        end
        # PyPlot.pause(0.000001)
        # PyPlot.cla()
        # PyPlot.draw()
        
        
end



