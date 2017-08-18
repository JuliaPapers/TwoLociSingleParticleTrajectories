module RCLBreakParams
# This module file contains a dictionary of parameters used for the RCLBreakSimulation.jl

params =       Dict([("numSimulations", 1),              # Number of simulations in each batch (Int64)
                     ("numRandomConnectors", [0]),       # Number of random connectors (column vecotr Int64). keep commas between connectors e.g [10,20,30]
                     ("numSimulationBatches",1),         # Determined by length of numRandomConnectors (Int64)
                     ("numRelaxationSteps",1000),        # Number of steps to polymer relaxation with dtRelaxation (Int64)
                     ("numSteps",100),                   # Number of simulation steps with dt (Int64)
                     ("calculateRelaxationTime",false),  # Relaxation time for each simulation determined by the eigenvalues of the connectivity matrix (boolean)
                     ("numMonomers",100),                # Number of monomers in the polymer chain (Int64)
                     ("dimension",3),                    # Dimension of the space (int64)
                     ("connectorsSTD",200.0),            # [nm]  
                     ("diffusionConst",8000.0),          # [nm]^2 /[sec]
                     ("dtRelaxation",0.1),               # Time step for the relaxation phase [sec]
                     ("dt",0.008),                       # Time steps [sec]. 
                     ("springConst", 3.0*8000.0/200.0^2),# Determined automatically (divided by the friction coefficient)
                     ("encounterDist",200.0),            # [nm]
                     ("simulateBreak",true),             # Induce a break true/false (boolean)
                     ("induceUniformDSB", false),        # Uniformly distributed DSBs (boolean)
                     ("numberOfDSBs",0),                 # Number of DSB to create in case induceUniformDSB=true (Int64)
                     ("monomersToDetach",[50 51]),       # Induce a break (remove spring) between these monomer pairs (Int64 row vector)
                     ("fetForMonomers",[50 51]),         # The FET between these two monomers (row vector Int64)
                     ("fetSimulation",false),            # Simulate the first encounter time between selected monomers (boolean)
                     ("fetSimulationStartFromMaxDistance",true),       # Start measuring fet when the monomers are at their amximal distance of a period of relaxation time (boolean)
                     ("fetSimulationStartFromMeanPosition",false),     # deprecated (boolean)
                     ("maxFETSteps",1e5),                              # Maximal steps for FET simulation (Int64)
                     ("breakAllConnectorsFromDetachedMonomers", true), # Remove all random connectors following DSB (boolean)
                     ("mechanicalSpringPointForce", false),            # deprecated- Harmonic force around each monomer  (boolean)
                     ("fitHistToData",true),                           # fit FET histograms with exponential (boolean)
                     ("saveResults",true),                             # Save results (boolean)
                     ("resultBaseFolder",  "../../../OwnCloud/RCL_DSBSimulationResults/DSB/"), # Path to base result folder (string)
                     ("resultFolderName","ResultExp01"),                                       # Name of result folder (string)
                     ("showSimulation",false),                                                 # Show real time simulation (boolean)
                     ("showChainSnapshot", false),                                             # Show chain after DSB and after first encounter (boolean)
                     ("histPlotFlag",false),                                                   # Plot FET histograms   (boolean)
                     ("plotEncounterTimeVsConnectivity", false),                               # Plot encounter time vs connectivity (boolean)
                     ("plotMaxMonomerDistanceVsConnectivity", false),                          # Plot maximal fetForMonomers distance vs connectivity (boolean)
                     ("numFitTrials",10),                                                      # Number of trials to fit the FET histograms (Int64)
                     ("numHistBins",20)])                                                      # Number of bins in the FET histogram (Int64)
end                     