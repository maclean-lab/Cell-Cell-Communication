
using DifferentialEquations
using MultiScaleArrays
using Distributions
using Plots
using DelimitedFiles

#-----------------------------------
# Number of cells
#-----------------------------------

nCell = 3

#-----------------------------------
# Signaling topology
#-----------------------------------

topology = "chain" # "loop" or "chain"

#-------------------------------------------------
# Wait time distribution and pulse time
#-------------------------------------------------

mean = 50.0

pulse = 5.0

#-----------------------------------
# Number of iterations
#-----------------------------------

nits = 10

#-----------------------------------
# lambda
#-----------------------------------


lambda = 1.0


#-----------------------------------
# Range of A0 values
#-----------------------------------

A0_values = 0.7:0.005:0.72

#-----------------------------------
# Noise and the variance
#-----------------------------------

noiseE = "off" # "on" or "off", extracellular noise
noiseP = "off" # "on" or "off", intrinsic noise

etaE = Normal(0.0, 0.01)
etaP = Normal(0.0, 0.1)

#-----------------------------------
# Internal ODE system
#-----------------------------------

function GataPu1(du,u,p,t)
    a = p[1:5] 
    b = p[6:8]
    g = p[9:17]
    L = p[18:20] 
    du[1] = -b[1]*u[1] + (a[1]*L[1] + a[2]*u[1])/(1 + g[1]*L[1] + g[2]*u[1] + g[3]*u[1]*u[2])
    du[2] = -b[2]*u[2] + (a[3]*L[2] + a[4]*u[2])/(1 + g[4]*L[2] + g[5]*u[2] + g[6]*u[1]*u[2] + g[7]*u[1]*u[3])
    du[3] = -b[3]*u[3] + (a[5]*u[1])/(1 + g[8]*u[1] + g[9]*L[3])
end

#-----------------------------------
# Parameter vector p
#-----------------------------------

p = zeros(0)
a = [1.0,0.25,1.0,0.25,0.01] 
b = [0.01,0.01,0.01] 
g = [1.0,0.25,1.0,1.0,0.25,1.0,0.13,0.01,10]
L = [1.0,0.5,0] #L = [A,B,C] = [A,0.5,0]
append!(p,a)
append!(p,b)
append!(p,g)
append!(p,L)

#-----------------------------------
# Wait time distribution
#-----------------------------------

function rand_exponential(mean)
    if mean <= 0.0
        error("mean must be positive")
    end
    -mean*log(rand())
end


function wait_time(mean)
    round(rand_exponential(mean))
end

#-----------------------------------
#Simulation Code
#-----------------------------------

struct Cell{B} <: AbstractMultiScaleArrayLeaf{B}
    values::Vector{B}
end
struct Tissue{T<:AbstractMultiScaleArray,B<:Number} <: AbstractMultiScaleArray{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end

GHighCounts = zeros(length(A0_values), nCell+1)

for l = 1:length(A0_values)
    A0 = A0_values[l]
    for j = 1:nits

    
        CellValues = [[20.0, 1.0, 20.0, A0] for i=1:nCell, n=1:1]
    
        if noiseE == "on"
            CellValues = [[20.0, 1.0, 20.0, A0 + rand(etaE)] for i=1:nCell, n=1:1]
        end
    
        CellIts = zeros(nCell)
        CellWaitTimes = [wait_time(mean) for i=1:nCell, n=1:1]
        CellState = ["wait" for i=1:nCell, n=1:1]
    
        while min(CellIts...) < 1000 + 40*nCell
            t = min(CellWaitTimes...)
        
            for i = 1:nCell
                u0 = [CellValues[i][1], CellValues[i][2], CellValues[i][3]]
                p[18] = CellValues[i][4]
                tspan = (0.0, t)
                prob = ODEProblem(GataPu1,u0,tspan,p)
                sol = solve(prob)
                CellValues[i][1] = sol[length(sol)][1]
                CellValues[i][2] = sol[length(sol)][2]
                CellValues[i][3] = sol[length(sol)][3]
            end
        
            CellWaitTimes = CellWaitTimes .- t
        
        
            if topology == "loop"
                if CellWaitTimes[1] == 0
                    if CellState[1] == "wait"
                        CellWaitTimes[1] = pulse
                        CellState[1] = "pulse"
                        CellIts[1] = CellIts[1] + 1
                        CellValues[1][4] = A0 * (CellValues[nCell][1]/(lambda + CellValues[nCell][2]))
                        if noiseP == "on"
                            CellValues[1][4] = A0 * (CellValues[nCell][1]/(lambda + CellValues[nCell][2])) + rand(etaP)
                        end
                    else
                        CellWaitTimes[1] = wait_time(mean)
                        CellState[1] = "wait"
                        CellValues[1][4] = A0
                        if noiseE == "on"
                            CellValues[1][4] = A0 + rand(etaE)
                        end
                    end
                end
            else
                CellWaitTimes[1] = 100.0
                CellIts[1] = 2000 + 40*nCell
            end
        
        
            for i = 2:nCell
                if CellWaitTimes[i] == 0
                    if CellState[i] == "wait"
                        CellWaitTimes[i] = pulse
                        CellState[i] = "pulse"
                        CellIts[i] = CellIts[i] + 1
                        CellValues[i][4] = A0 * (CellValues[i-1][1]/(lambda + CellValues[i-1][2]))
                        if noiseP == "on"
                            CellValues[i][4] = A0 * (CellValues[i-1][1]/(lambda + CellValues[i-1][2])) + rand(etaP)
                        end
                    else
                        CellWaitTimes[i] = wait_time(mean) 
                        CellState[i] = "wait"
                        CellValues[i][4] = A0
                        if noiseE == "on"
                            CellValues[i][4] = A0 + rand(etaE)
                        end
                    end
                end
            end
    
        end
    
        for i = 1:nCell
            if CellValues[i][1] > 30
                GHighCounts[l, i+1] = GHighCounts[l, i+1] + 1
                GHighCounts[l, 1] = A0
            end
        end
    end
    for i = 1:nCell
        GHighCounts[l, i+1] = GHighCounts[l, i+1]/nits
    end
end


writedlm(topology*"_$nCell"*"_cells.txt", GHighCounts, ", ") 
