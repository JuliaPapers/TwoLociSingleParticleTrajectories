# experimental type for RCL parameters
module  RCLParams
# Simulation parameters
module Simulator
        const numSimulationBatches = 1
        const numSimulations       = 10                     
        const numSteps             = 100
        const numRelexationSteps   = 100
        const dimension            = 3
        const Δt                   = 0.1 # simulation time step
        const encounterDist        = 0.1
        const monomersToDetach     = [50,51]    
end

module Chain
       const numMonomers         = 100
       const numRandomConnectors = 0                                                               
       const connectorsSTD       = 1.0
       const κ                   = 1.0 # kappa spring constant
       const diffusionConst      = 1.0                                     
end

module Graphical
        const plotFlag             = false # for debugging
        const histPlotFlag         = true # plot histograms and fit  
        const numHistBins          = 20
        const lineWidth            = 3
        const lineColor            ="red"
        const fontSize             = 10
        const numFitTrials         = 10
end

end
