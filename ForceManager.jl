
module ForceManager
using Distances

function ForceParams()
 params = Dict([("diffusionForce",true),
                    ("springForce",true),
                    ("mechanicalForce",false),
                    ("diffusionConst",1.0),
                    ("diffusionConstBoundary",0.0),
                    ("springConst",1.0),
                    ("mechanicalForceSpringConst",0.0),
                    ("mechanicalForceSpringUpperCutoff",1.0),
                    ("mechanicalForceSpringLowerCutoff",0.0)])
return params
end

function BrownianNoise(numMonomers::Int64,diffusionConst::Float64,dimension::Int64,Δt::Float64)
     sqrt(2.0*diffusionConst*Δt/Float64(dimension)).*randn(numMonomers,dimension)
end

function SpringForce(particlePos::Array{Float64,2},particleDist::Array{Float64,2},springConst::Float64,minParticleDist::Float64)
        
        numParticles    = size(particlePos,1)
        dimension       = size(particlePos,2)
        connectivityMap = diagm(ones(numParticles-1),1)+diagm(ones(numParticles-1),-1)
        L               = (1.0.-minParticleDist./particleDist).*connectivityMap                
        [L[i,i]=0 for i in 1:numParticles]
        L               = springConst*dimension.*L                
        # println(L[3,:])

        # L               = L.*springConst.*dimension      
        sumForces       = sum(L,2)
        # println(size(L))
        # println(size(diagm(sumForces)))
        # force                        = -1.*(bsxfun(@times,sumForces,eye(numMonomers))-L); % set the maindiagonal                
         
        force           = L-diagm(sumForces[:])
                # println(force[2,:])
        force           = force*particlePos
        ##### from Matlab code 
        #             force       = zeros(size(particlePosition));
        #             numMonomers = size(connectivityMap,1);
        #             dimension   = size(particlePosition,2);
        #             if springForce                
        #                 L                            = (1-(minParticleDist./particleDist));%.*connectivityMap;
        #                 L                            = L.*springConst*dimension;
        #                 L(~connectivityMap)          = 0;
        #                 sumForces                    = sum(L,2);
        #                 force                        = -1.*(bsxfun(@times,sumForces,eye(numMonomers))-L); % set the maindiagonal                
        #                 force                        = force*particlePosition;
        #                 force(particlesOnBoundary,:) = 0;                
        #                 force(fixedParticleNum,:)    = 0;% zero out forces for fixed particles
                        
        # %                 force = SpringForce_mex(particleDist,springConst,connectivityMap,minParticleDist);
        # %                 force(fixedParticleNum,:)    = 0;% zero out forces for fixed particles
        # %                 force                        = force*particlePosition;
        # %                 force(particlesOnBoundary,:) = 0; 

                        #     end


return force
end

function ParticleDist(pos)

  particleDist = Distances.pairwise(Euclidean(),pos',pos')
  return particleDist

end
        
function BendingElasticityForce()
 return nothing 
end

function MechanicalSpringPointForce!(particlePosition::Array{Float64,2},pointSourcePosition::Array{Float64,2}, forceDirection::String,
                                     springConst::Float64,lowerCutoff::Float64,upperCutoff::Float64,dt::Float64)

 # Apply an harmonic pushing force for particles around force sources. 

 # Calculate the distance betwen each particle and the point source position
    #  force = zeros(size(particlePosition,1), size(particlePosition,2))
    # collect force contribution from all sources 

     distToSource = Distances.pairwise(Euclidean(),particlePosition',pointSourcePosition')
     d2sInd       = (distToSource.<=upperCutoff) & (distToSource.>=lowerCutoff)
    #  distToSource = distToSource.*d2sInd
     if forceDirection=="in"
       springCont = -springConst
      end

    @inbounds for r in 1:size(particlePosition,1) # for each particle                 
                 @inbounds  for c in 1:size(pointSourcePosition,2) # for each source
                              if d2sInd[r,c]
                                # force[r,:] = force[r,:]+(springConst*dt)*(particlePosition[r,:] - pointSourcePosition[c,:]) # broadcast(-,particlePosition,pointSourcePosition[sIdx,:]')
                                particlePosition[r,:] = particlePosition[r,:]+(springConst*dt)*(particlePosition[r,:] - pointSourcePosition[c,:]) # broadcast(-,particlePosition,pointSourcePosition[sIdx,:]')                                  
                              end
                            end
              end

              
    #  @inbounds for fIdx in f
    #                 force[r[n],:]+= forceDir[fIdx,:]
    #               end
    #           end

            #    if forceDirection=="in"               
            #  @inbounds  for n in eachindex(force)
            #          force[n] += -springConst.*force[n]#broadcast(*,distToSource,forceDir)
            #         end

            #     elseif forceDirection =="out"
            #  @inbounds for n in eachindex(force)
            #          force[n] += springConst.*force[n]#broadcast(*,distToSource,forceDir)
            #         end
            #         else 
            #           error("unsupported force direction ")
            #    end

# update particle position 
        # @inbounds  for n in eachindex(particlePosition)
        #               particlePosition[n] += force[n]
        #           end
                 
return particlePosition
end

end