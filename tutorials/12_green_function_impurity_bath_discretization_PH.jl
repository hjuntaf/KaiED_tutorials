ENV["PROJECT_PATH_ED"]="../envs/KED"
include("../src/mybase.jl")

norb        = 1
nspin       = 2
nspinorb    = norb*nspin

nbath       = 5
@inline IndBath( ibath )            = ibath + nspinorb
@inline IndBathSpinup( ibathorb )   = 2*ibathorb + nspinorb - 1
@inline IndBathSpindn( ibathorb )   = 2*ibathorb + nspinorb

ntot    = nspinorb+nbath
dim     = 2^ntot
@show dim


ebathl     = zeros(nbath)
ebathl     = collect( LinRange( -1, 1, nbath ) ) 
@show ebathl

Vil     = zeros(nspinorb,nbath)
for i in 1:nspinorb
    for j in 1:nbath
        Vil[i,j]    = 2 - ebathl[j]^2
    end
end
@show Vil

tij     = zeros(nspinorb,nspinorb)
tij[1,1]= 0


epsilon = 0.01

beta    = 32
NImFreq = 4*beta
ImFreqGridVal   = GetImFreqValGrid( beta, NImFreq )
ImFreqGrid      = ImFreqGridVal * im

G0iw= GetGzBethe.( ImFreqGrid )

Hyb11iw   = 1. / 4 * G0iw   # = t^2 G = 1/4 G


## Testing functions of flattened-parameters ##
BParamPH    = BathParamFlatten( ebathl[1:div(end,2)], Vil )
Diw  = GetDeltaHybDiscGrid( ebathl, Vil, ImFreqGrid )
Diw2 = GetDeltaHybDiscGridFromFlatPH( BParamPH, ImFreqGrid, nspinorb, nbath )
D11iw    = GetijarrayFromVecMat( Diw, 1, 1 ) 
D11iw2   = GetijarrayFromVecMat( Diw2, 1, 1 ) 

Hybiw        = 1. / 4 * GetGzBetheUniformScaling.( ImFreqGrid )

SParam  = ( 
                Hybiw, 
                ImFreqGrid, 
                nspinorb, 
                nbath 
                )

@inline function GetCostFromFlat( bathparam, systemparam ) 
    HybwOrig    = systemparam[1]
    wgrid       = systemparam[2]
    dimorb      = systemparam[3]
    dimbath     = systemparam[4]
    HybwParam   = GetDeltaHybDiscGridFromFlat( bathparam, wgrid, dimorb, dimbath )
    return GetCost( HybwOrig, HybwParam )
end

using BenchmarkTools
@time cost = GetCostFromFlatPH( BParamPH, SParam )
println( "initial cost : $(cost) " )


# using Optimization
# prob = OptimizationProblem(GetCostFromFlatPH, BParamPH, SParam)
# using OptimizationOptimJL
# sol = solve(prob, NelderMead())
# @show sol.original
# BParamPHNew   = [ sol... ]

using Optim
@time res = optimize( x -> GetCostFromFlatPH(x,SParam), BParamPH, LBFGS(), Optim.Options(iterations=4000))
@show res
BParamPHNew   = [ res.minimizer... ]


@time cost = GetCostFromFlatPH( BParamPHNew, SParam )
println( "final cost : $(cost) " )

enew, vnew  = BathParamReshapePH( BParamPHNew, nbath )
@show enew
println("vnew :")
using DelimitedFiles
writedlm(stdout, vnew)

Eorb        = collect(I(nspinorb)*0.0)
DiwNew      = GetDeltaHybDiscGridFromFlatPH( BParamPHNew, ImFreqGrid, nspinorb, nbath )
GiwNew      = GetGreenLocalFromHybGrid( DiwNew, ImFreqGrid, Eorb )

x   = ImFreqGridVal
y1  = Hyb11iw
y2  = GetijarrayFromVecMat( DiwNew, 1, 1 ) 



## ReFreq spectra ##

NReFreq = 100
epsilon = 0.04
ReFreqGridVal   = collect( LinRange( -3, 3, NReFreq ) )
ReFreqGrid      = ReFreqGridVal .+ im * epsilon

Gw          = GetGzBethe.( ReFreqGrid )
DwNew       = GetDeltaHybDiscGridFromFlatPH( BParamPHNew, ReFreqGrid, nspinorb, nbath )
Gwn         = GetGreenLocalFromHybGrid( DwNew, ReFreqGrid, Eorb )

bPlot = true
if bPlot
    using Plots
    xw  = ReFreqGridVal
    gw  = Gw
    gwn = GetijarrayFromVecMat( Gwn, 1, 1 ) 

    plot(  xw, [ real(gw)  imag(gw)  ] )
    plot!( xw, [ real(gwn) imag(gwn) ] )
end


