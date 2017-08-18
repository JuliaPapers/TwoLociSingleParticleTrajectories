# This script runs simulation of the RCL polymer with DSB induction. the MFET, max monomer dist, L_C, k_C are estimated and exported as csv files in the result folder
# to run this script type include("scrRunRCLBreakSimulation.jl") in the julia REPL, when the pwd is the code folder
# parameters are defined in RCLBreakParams.jl

include("$(pwd())/RCLBreakSimulation.jl");
include("$(pwd())/RCLBreakParams.jl");

rclDSBParams = RCLBreakParams;# load parameters 

# Set parameters 
rclDSBParams.params["numSimulations"]          = 1       # number of simulatiosn in each batch (Int64)
rclDSBParams.params["numRandomConnectors"]     = ([10, 50, 130]) # keep commas between connectors, wrap vector within brackets ()
rclDSBParams.params["numSimulationBatches"]    = length(rclDSBParams.params["numRandomConnectors"]);
rclDSBParams.params["numRelaxationSteps"]      = 10000; # (Int64)
rclDSBParams.params["numSteps"]                = 7500;  # (Int4) in case it is not a transient FET simulation, a total of 60s of simulation numStps*deltaT=60s
rclDSBParams.params["calculateRelaxationTime"] = false  # (boolean) determine relaxatio by the smallest non-vanishing eig of the connectivity matrix
rclDSBParams.params["numMonomers"]             = 100    # number of monomers of the RCL polymer (Int64)
rclDSBParams.params["dimension"]               = 3      # spatial dimension (Int64)
rclDSBParams.params["connectorsSTD"]           = 200.0  # [nm] (Float64)
rclDSBParams.params["diffusionConst"]          = 8000.0 # [nm]^2 /[sec] (float64)
rclDSBParams.params["dtRelaxation"]            = 0.01   # time step for relaxation stage [sec]
rclDSBParams.params["dt"]                      = 0.008  # time steps for post relaxation stages [sec] 
rclDSBParams.params["springConst"]             = rclDSBParams.params["dimension"] *rclDSBParams.params["diffusionConst"]/rclDSBParams.params["connectorsSTD"]^2# normalized
rclDSBParams.params["encounterDist"]           = 200.0  # [nm]
rclDSBParams.params["simulateBreak"]           = true   # simulate a DSB between monomers monomersToDetach
rclDSBParams.params["fetForMonomers"]          = [50 51] # monomer number to check for DSB
rclDSBParams.params["induceUniformDSB"]        = false
rclDSBParams.params["numberOfDSB"]             = 0 # overrides monomersToDetach in case induceUniformDSB=true, cannot exceed numMonomers/2
rclDSBParams.params["monomersToDetach"]        = [50 51] # matrix of indices pairs of successive monomers to detach. Gets overwritten in induceUniformDSB=true
rclDSBParams.params["fetSimulation"]           = true # simulate the first encounter time between selected monomers. If fetSimulation=false, the values of Lc and Kc are computed
rclDSBParams.params["fetSimulationStartFromMaxDistance"]     = false # start measuring fet when the monomers are at their maximal distance of a period of relaxation time 
rclDSBParams.params["fetSimulationStartFromMeanPosition"]    = false;# start FET simulation from the mean position of the chain over the course of numRelaxationSteps
rclDSBParams.params["maxFETSteps"]             = 1e5 # maximal steps for FET simulation
rclDSBParams.params["breakAllConnectorsFromDetachedMonomers"]= true 
rclDSBParams.params["mechanicalSpringPointForce"] = false
rclDSBParams.params["fitHistToData"]              = false # histograms 
rclDSBParams.params["saveResults"]                = true
rclDSBParams.params["resultBaseFolder"]           = "$(pwd())"#../../../OwnCloud/RCL_DSBSimulationResults/DSB/"# changes for each simulation type 
rclDSBParams.params["resultFolderName"]           = "ResultExp01"
rclDSBParams.params["showSimulation"]             = false  # show each step after relaxation time (boolean)
rclDSBParams.params["showChainSnapshot"]          = false  # show only several snapshot of the chain after break and after first encounter (boolean)
rclDSBParams.params["histPlotFlag"]               = false  # plot FET histogram (boolean)
rclDSBParams.params["plotEncounterTimeVsConnectivity"]      = true  # (boolean)
rclDSBParams.params["plotMaxMonomerDistanceVsConnectivity"] = false # (boolean)
rclDSBParams.params["numFitTrials"]               = 10
rclDSBParams.params["numHistBins"]                = 20

# Run simulation non broken loci
rclDSBParams.params["simulateBreak"]    = false
rclDSBParams.params["resultBaseFolder"] = "$(pwd())/nonBroken/"#../../../OwnCloud/RCL_DSBSimulationResults/nonBroken/" 
 @time resultsUnbroken                  =  RunRCLBreakSimulation(rclDSBParams.params)

# Run DSB simulation
rclDSBParams.params["simulateBreak"]    = true
rclDSBParams.params["resultBaseFolder"] = "$(pwd())/DSB/"#../../../OwnCloud/RCL_DSBSimulationResults/DSB/"
@time resultsDSB                        = RunRCLBreakSimulation(rclDSBParams.params)


