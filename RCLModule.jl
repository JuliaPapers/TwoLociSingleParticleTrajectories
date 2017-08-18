module RCL
# this module summarizes function related to the steady state properties of the Randomly Cross Linked polymer
 forceManager = include("ForceManager.jl")

 using Distances

function RCLParams()
  #TODO: RCLparams should be turned into a Type, to be able to add two chains together 
  params = Dict([("numMonomers",100), 
                ("b",1.0),                     
                ("ξ" ,0.0),
                ("monomersOnBoundary",[]),
                ("fixedMonomers",[]),
                ("connectedMonomes",[]),
                ("minMonomerDist",0.0),
                ("laplacian", Array{Float64,2}), 
                ("curPos",Array{Float64,2}),
                ("prevPos",Array{Float64,2}),                    
                ("forceParams",ForceManager.ForceParams())])
  return params
end

# Tranlate connectivity fraction to the number of connectors
function Connectivity2Connectors(N::Int64,connectivityFrac::Float64)
 Int64(round(connectivityFrac*(N-1)*(N-2)./2))
end

# Translate number of random connectors to connectivity fraction
function Connectors2Connectivity(Nc::Int64,N::Int64)
   2*Nc./((N-1)*(N-2))
end

# Eigenvalues of the RCL polymer as a function of the connectivityFraction
function Eigenvalue(numMonomers::Int64,connectivityFrac::Float64,modeNum::Int64)
 # the RCL eigen values are defined as N*xi_p +lambda*(1-xi_p), where \lambda_p is the Rouse eigenvalue
  if modeNum==0
    eigenvalue =0
  else
    eigenvalue = numMonomers*connectivityFrac +4*sin(modeNum*pi/(2*numMonomers))^2 *(1-connectivityFrac)
  end

return eigenvalue
end

# The mean squared radius of gyration 
function MSRG(numMonomers::Int64,connectivityFrac::Float64,connectorsSTD::Float64)
  b    = connectorsSTD 
  N    = numMonomers
  c    = connectivityFrac
  y    = 1+N*c/(2.0*(1.0-c))
  z0   = y +sqrt(y^2 -1.0) 
  z1   = y - sqrt(y^2 -1.0)
  x    = connectivityFrac
  
  #msrg = sqrt(2)*(b^2 /(N^2 *(1-c)*(z0-z1))) *((1+2*z0)*N*(1+N)/(2*z0) +N*(2*(1+z0)^2 -z0^3)/(1-z0^2) -z0^3 *(1-1/(z0^(2*N)))/(1-z0^2)^2 +2*(1+z0)*(1-1/z0^N)/(1-z0)^2)

  
  msrg = ((b^2)./((1-x).*(z0-z1)*N^2)).*
                (z0.^(-2*N) .*(-z0.*(1 + N*(-1 + z0) + z0) + 
                z0.^(2*N) .*(z0 + (-1 + N.*(-1 + z0)).*(-N + (-1 + N).*z0.^2))))./((-1 + z0).^2 .*(1 + z0))

  return msrg
end

function TimeDependentVariance(numMonomers::Int64,t::Float64,modeNum::Int64,connectivityFrac::Float64,connectorsSTD::Float64, diffusionConst::Float64)
  ep = Eigenvalue(numMonomers,connectivityFrac,modeNum)
  # variance of the normal coordinates as the sum of variance in all dimensions
  (connectorsSTD^2 /ep) *(1-exp(-2*diffusionConst*ep*t/connectorsSTD^2))
end

# Calculate the variance of the RCL polymer given average connectivity fraction ξ
function Variance(m::Int64, n::Int64,N::Int64,xi::Float64,b::Float64)
 # Auxilary functions to calculate the RCL variance
  xi=xi+ 0.0000001 # for numerical stability
  y  = (1.0+(N*xi)./(2.0*(1.0-xi)))
  z0 = y+sqrt(y^2.0 -1.0)
  z1 = y-sqrt(y^2.0 -1.0)

    # RCLVar = (b^2/(1-ξ)).*(1-ν0(N,ξ).^(-k)).*(2*ν0(N,ξ) +1-ν0(N,ξ).^(-k))./(ν0(N,ξ).^2 -1)

    
    for nIdx in 1:length(n)
        if n[nIdx]>m       
            # non symmetric 
                #  sigma(nIdx) = (1./z0.^(2*n(nIdx)-1))*((z0^(n(nIdx)-m)-1).^2 -2*z0.^(m+n(nIdx)-1))+2;
                  #  symmetric 
                sigma = ((1./z0.^(2*n[nIdx]-1))*((z0^(n[nIdx]-m)-1).^2 -2*z0.^(m+n[nIdx]-1))+2+(1./z0.^(2*(N-m+1)-1))*((z0.^(N-m+1-(N-n[nIdx]+1))-1).^2 -2*z0.^(N-m+1+N-n[nIdx]+1-1))+2)./2
                #MATLAB version
                #sigma(nIdx) = ((1./z0.^(2*n(nIdx)-1))*((z0^(n(nIdx)-m)-1).^2 -2*z0.^(m+n(nIdx)-1))+2+(1./z0.^(2*(N-m+1)-1))*((z0.^(N-m+1-(N-n(nIdx)+1))-1).^2 -2*z0.^(N-m+1+N-n(nIdx)+1-1))+2)./2;

        else              
            # non symmetric 
            # sigma(nIdx) = (1./z0.^(2*m-1))*((z0^(m-n(nIdx))-1).^2 -2*z0.^(m+n(nIdx)-1))+2;
                
                # symmetric 
                sigma = ((1./z0.^(2*m-1))*((z0.^(m-n[nIdx])-1).^2 -2*z0.^(m+n[nIdx]-1))+2+(1./z0.^(2*(N-n[nIdx]+1)-1))*((z0^(N-n[nIdx]+1-(N-m+1))-1).^2 -2*z0.^(N-m+1+N-n[nIdx]+1-1))+2)./2
                #MATLAB version 
                #sigma(nIdx) = ((1./z0.^(2*m-1))*((z0.^(m-n(nIdx))-1).^2 -2*z0.^(m+n(nIdx)-1))+2+(1./z0.^(2*(N-n(nIdx)+1)-1))*((z0^(N-n(nIdx)+1-(N-m+1))-1).^2 -2*z0.^(N-m+1+N-n(nIdx)+1-1))+2)./2
        end
    end

    sigma = sigma.*((b^2) /((1-xi)*(z0-z1)))

return sigma  
end

# Encounter probability between two monomers at a distance k (integer)
# in dimension d with connector std b
function EncounterProbability(N,ξ::Float64,k::Float64,b::Float64,dimension::Int64)
  a = (dimension./(2*π*Variance(N,ξ,k,b))).^(dimension/2)  
  a./sum(a)
end

# construct Rouse connectivity matrix
function GraphSemiRandomLaplacian(N::Int64,connectivityFraction::Float64)
   M                   = -diagm(ones(N-1),1)-diagm(ones(N-1),-1)   
   sm                  = sum(M,2)
   M                   = M-diagm(sm[:])
  # add random connectivity
   Nl                  = (N-1)*(N-2)>>1 # shift one bit (half the number in integer)
   Nc::Int64           = Connectivity2Connectors(N,connectivityFraction)# number of connectors
   randp               = randperm(Nl)
   z::Array{Float64,2} = fill(0.0,(Nl,1))
   rp                  = randperm(Nl)
   z[rp[1:Nc]]         = -1.0
   B                   = fill(0.0,(N,N))
   nIdx::Int64         = 1
   cs::Int64           = 0
  @inbounds for nIdx in 1:(N-3)
              numElements::Int64 = N-(nIdx+2)+1
              B[nIdx,(nIdx+2):N]=z[cs+1:cs+numElements]
              # z= z[numElements:end]
              cs = cs+ numElements
            end

   B += B' # make symmetric
   B += -diagm(sum(B,1)[:]) # add monomer connectivity
return M+B
end

function GraphAverageLaplacian(N::Int64,connectivityFraction::Float64)
 # The system simulated with the average connectivitymatrix <B(\xi)> 
   M = connectivityFraction.*ones(N,N)
   [M[mIdx,mIdx+1]=1 for mIdx in 1:N-1]
   [M[mIdx,mIdx-1]=1 for mIdx in 2:N]   
   [M[mIdx,mIdx] = 0 for mIdx in 1:N]
   M -=diagm(sum(M,2)[:],0)
   return -M
end


function GraphLaplacian(N::Int64,connectedMonomers::Array{Int64})
 # generate a laplacian matrix with fixed connectors between monomers defined in connectedMonomer array
   M  = -diagm(ones(N-1),1)-diagm(ones(N-1),-1)
   M[connectedMonomers[:,1], connectedMonomers[:,2]]= -1.0
   M[connectedMonomers[:,2], connectedMonomers[:,1]]= -1.0
  # for bIdx in 1:size(connectedMonomers,1)
  #   M[connectedMonomers[bIdx,1],connectedMonomers[bIdx,2]]= -1;
  #   M[connectedMonomers[bIdx,2],connectedMonomers[bIdx,1]]= -1;
  # end
return M-diagm(sum(M,1)[:])
end

# Disconnect two monomers of the chain 
function DisconnectMonomers!(monomer1, monomer2, laplacian)
     numMonomers = size(laplacian,1)
     # if monomer1<=numMonomers & monomer2<=numMonomers
     if monomer1 != monomer2     
      if laplacian[monomer1, monomer2] !=0.0
          laplacian[monomer1,monomer2] = 0.0
          laplacian[monomer2,monomer1] = 0.0
          laplacian[monomer1,monomer1] -= 1.0       
          laplacian[monomer2,monomer2] -= 1.0      
      end
    end
    return laplacian
end

# Connect two monomers of the chain 
function ConnectMonomers!(monomer1::Int64, monomer2::Int64,graphLaplacian)
    numMonomers = size(graphLaplacian,1)
    if monomer1<=numMonomers & monomer2<=numMonomers
      graphLaplacian[monomer1,monomer2] = -1.0
      graphLaplacian[monomer2,monomer1] = -1.0
      graphLaplacian[monomer1,monomer1] = graphLaplacian[monomer1,monomer1]+1.0
      graphLaplacian[monomer2,monomer2] = graphLaplacian[monomer2,monomer2]+1.0
    end
    return graphLaplacian
end

function Step!(nSteps::Int64,pos::Array{Float64,2},laplacian::Array{Float64,2},κ::Float64,dimension::Int64, Δt::Float64,diffusionConst::Float64)
  # Advance nSteps of the simulation  
     s = sqrt(2.0*diffusionConst.*Δt)#./dimension)
     @inbounds for sIdx = 1:nSteps   
                  lt = laplacian*(pos*(-κ)*Δt)        
      @inbounds  for n in eachindex(pos)
                   pos[n]  = pos[n]+ lt[n]  + s*randn()
                 end         
              end
return pos
end

function GetConnectivity(laplacian::Array{Float64,2})
 # off-diagonal connected monomers 
 cm = findn(triu(laplacian-diagm(diag(laplacian))-diagm(diag(laplacian,1),1)))
return cm
end

function RelaxationTime(laplacian, springConst)
  # calculate relaxation time
  # TODO: the diffusion constant and the dimension should be the input 
  # the fomula should read tau = D/dimension *lambda_0

             numRelaxationSteps = 0 
                eVals, eVecs = eig(laplacian)
                if all(imag(eVals).==0)
                    # find the minimal non vanishin real eigen value                  
                    lambdaMinInd       = findfirst(real(eVals).>1e-8)
                    lambda             = eVals[lambdaMinInd]                      
                    numRelaxationSteps = convert(Int64,floor(1/(springConst*lambda)))
                else
                    # keep the default value 
                    numRelaxationSteps = 0
                end                 
  return numRelaxationSteps
end

function InitializeChain(numMonomers::Int64,b::Float64,diffusionConst::Float64,dimension::Int64,dt::Float64)
   # Get an  initial chain position 
   initPos      = zeros(numMonomers,dimension)
   initPos[1,:] = randn(1,dimension)
   for mIdx in 2:numMonomers
     r = randn(dimension)
     #r = r./sqrt(sum(r.^2)) # normalize
     initPos[mIdx,:] = initPos[mIdx-1,:]+ r
   end

   return initPos 
end

function DiffuseOnSphere(initialPoint::Array{Float64,2},numSteps::Int64,
                        Δt::Float64,diffConst::Float64,domainCenter::Array{Float64,2},
                        domainRadius::Float64)
  #check initialPoint is on the sphere or not
  # numSteps should be obsolete
  # currently works only with one point
   numParticles::Int64        = size(initialPoint,1)
   dimension::Int64           = size(initialPoint,2)    
   pathsOnBoundary            = zeros(numParticles,dimension,numSteps+1)
   pathsOnBoundary[:,:,1]     = initialPoint
  
      
   pIdx::Int64=1
   @inbounds for pIdx = 1:numSteps
    # generate a path 2D with numStep;
    # Find the phi and theta on the sphere
    ϕ::Array{Float64,1}    = atan(pathsOnBoundary[:,2,pIdx]./pathsOnBoundary[:,1,pIdx])
    phiInds::Array{Bool,1} = ϕ.<0
    ϕ[phiInds]             = ϕ[phiInds].+ 2*π
   
    θ::Array{Float64,1}     = acos(pathsOnBoundary[:,3,pIdx]./domainRadius)
    ipInds::Array{Bool,1}   = pathsOnBoundary[:,1,pIdx].<0 
    θ[ipInds]               = -θ[ipInds]
    paths                   = [forceManager.BrownianNoise(numParticles,diffConst,2,Δt) zeros(numParticles,1)]
    paths[:,3]              = paths[:,3].+ domainRadius

    mIdx::Int64 = 1
  @inbounds for mIdx in 1:numParticles
     # Rotate the paths to the initialPoint 
      
      Ry::Array{Float64,2} = [cos(θ[mIdx]) 0 sin(θ[mIdx]) ; 0 1 0 ; -sin(θ[mIdx]) 0 cos(θ[mIdx])]
      Rz::Array{Float64,2} = [cos(ϕ[mIdx]) -sin(ϕ[mIdx]) 0;sin(ϕ[mIdx]) cos(ϕ[mIdx]) 0; 0 0 1]
      # println((Rz*Ry*(paths[pIdx,:].-domainCenter[1,:])))
      paths[mIdx,:] = (Rz*Ry*(paths[mIdx,:].-domainCenter[1,:]))'
    end
    #  project the paths back to the sphere     
    distToCenter::Array{Float64,2} = Distances.pairwise(Euclidean(),paths',domainCenter')   
   
    t::Array{Float64,2}         = domainRadius'./distToCenter
    pathsOnBoundary[:,:,pIdx+1] = domainCenter.+t.*(paths.-domainCenter)
          
    end
  return pathsOnBoundary[:,:,2:numSteps+1]
end


end
