ENV["PROJECT_PATH_ED"]="../envs/KED"
include("../src/mybase.jl")

norb        = 1
nspin       = 2
nspinorb    = norb*nspin

nbath       = 10
@inline IndBath( ibath )            = ibath + nspinorb
@inline IndBathSpinup( ibathorb )   = 2*ibathorb + nspinorb - 1
@inline IndBathSpindn( ibathorb )   = 2*ibathorb + nspinorb

ntot    = nspinorb+nbath
dim     = 2^ntot
@show dim


nbathHalf   = div(nbath,2)
ebathl     = zeros(nbath)
ebathl     = collect(LinRange( -1, 1, nbath ))
ebathl[1:2:end] = collect(LinRange(-1,1,nbathHalf))
ebathl[2:2:end] = collect(LinRange(-1,1,nbathHalf))
ebathlUp    = ebathl[1:2:end]
ebathlDn    = ebathl[2:2:end]
@show ebathlUp
@show ebathlDn
@show ebathl

Vil     = zeros(nspinorb,nbath)
for i in 1:norb
    for j in 1:nbathHalf
        Vil[2*i-1,2*j-1]    = 2 - ebathl[2*j-1]^2
        Vil[2*i  ,2*j  ]    = 2 - ebathl[2*j  ]^2
    end
end
# @show Vil
println( "Vil : " ) ; writedlm(stdout, Vil)

IndOrbUp    = [ 2*i-1 for i in 1:norb ]
IndOrbDn    = [ 2*i   for i in 1:norb ]
IndBathUp   = [ i for i in 1:2:2*nbathHalf ]
IndBathDn   = [ i for i in 2:2:2*nbathHalf ]
VilUp   = Vil[ IndOrbUp, IndBathUp ]
VilDn   = Vil[ IndOrbDn, IndBathDn ]
println( "Vil Up : " ) ; writedlm(stdout, VilUp)
println( "Vil Dn : " ) ; writedlm(stdout, VilDn)

tij     = zeros(nspinorb,nspinorb)
tij[1,1]= 0


beta    = 32
NImFreq  = 4*beta
ImFreqGridVal   = GetImFreqValGrid( beta, NImFreq )
ImFreqGrid      = ImFreqGridVal * im

G0iw= GetGzBethe.( ImFreqGrid )

Hyb11iw   = 1. / 4 * G0iw


## Testing functions of flattened-parameters (Spin-up) ##
DimBathMid  = div( nbath+1 , 2 )
@show DimBathMid  
BParamUp    = BathParamFlatten( ebathlUp, VilUp )
BParamUpPH  = BathParamFlatten( ebathlUp[1:div(end,2)], VilUp )
@show BParamUpPH
enew, vnew  = BathParamReshapePH( BParamUpPH, nbathHalf )
println( "" )
println( "ebathl after PH : $(enew)" ) 
println( "Vil Up after PH : " ) ; writedlm(stdout, vnew)
println( "" )
BParam      = BParamUpPH
# Duiw      = GetDeltaHybDiscGrid( ebathl, VilUp, ImFreqGrid )
Duiw2     = GetDeltaHybDiscGridFromFlatPH( BParam, ImFreqGrid, norb, nbathHalf )
Du11iw2   = GetijarrayFromVecMat( Duiw2, 1, 1 ) 

# ## Testing functions of flattened-parameters ##
# BParam  = BathParamFlatten( ebathl, Vil )
# Diw  = GetDeltaHybDiscGrid( ebathl, Vil, ImFreqGrid )
# Diw2 = GetDeltaHybDiscGridFromFlat( BParam, ImFreqGrid, nspinorb, nbath )
# D11iw    = GetijarrayFromVecMat( Diw, 1, 1 ) 
# D11iw2   = GetijarrayFromVecMat( Diw2, 1, 1 ) 

Hybiw        = 1. / 4 * GetGzBetheUniformScaling.( ImFreqGrid )

# SParam  = Any[]
# push!( SParam, Hybiw, ImFreqGrid, nspinorb, nbath )
SParam  = ( 
                Hybiw, 
                ImFreqGrid, 
                norb, # nspinorb, 
                nbathHalf
                )

using BenchmarkTools
@time cost = GetCostFromFlatPH( BParam, SParam )
println( "initial cost : $(cost) " )


# using Optimization
# prob = OptimizationProblem(GetCostFromFlatPH, BParam, SParam)
# using OptimizationOptimJL
# sol = solve(prob, NelderMead())
# @show sol.original
# BParamNew   = [ sol... ]

using Optim
@time res = optimize( x -> GetCostFromFlatPH(x,SParam), BParam, LBFGS(), Optim.Options(iterations=4000))
@show res
BParamNew   = [ res.minimizer... ]



enew, vnew  = BathParamReshapePH( BParamNew, nbathHalf )
@show enew
writedlm(stdout, vnew)

# DiwNew      = GetDeltaHybDiscGridFromFlat( BParamNew, ImFreqGrid, nspinorb, nbath )
DuiwNew      = GetDeltaHybDiscGridFromFlatPH( BParamNew, ImFreqGrid, norb, nbathHalf )

x   = ImFreqGridVal
y1  = Hyb11iw
# y2  = GetijarrayFromVecMat( DiwNew, 1, 1 ) 
y2  = GetijarrayFromVecMat( DuiwNew, 1, 1 ) 

NReFreq     = 200
epsilon     = 0.02
ReFreqGrid      = LinRange( -3, 3, NReFreq )
ReFreqGridBroad = ReFreqGrid .+ im * epsilon
xr  = ReFreqGrid
yr1 = GetGzBethe.( ReFreqGridBroad )
yrD2  = GetDeltaHybDiscGrid( enew, vnew, ReFreqGridBroad )
yr2  = GetGreenImpGrid( zeros(1,1), yrD2, ReFreqGridBroad )
yr2  = GetijarrayFromVecMat( yr2, 1, 1 )



println("## Spin-dn optimization ##")

## Setting BParam ##
BParamDn    = BathParamFlatten( ebathlDn, VilDn )
BParamDnPH  = BathParamFlatten( ebathlDn[1:div(end+1,2)], VilDn )
BParam      = BParamUpPH

## Setting SParam ##
Hybiw        = 1. / 4 * GetGzBetheUniformScaling.( ImFreqGrid )
SParam  = ( 
                Hybiw, 
                ImFreqGrid, 
                norb, # nspinorb, 
                nbathHalf
                )
@time cost = GetCostFromFlatPH( BParam, SParam )
println( "initial cost : $(cost) " )

# prob = OptimizationProblem(GetCostFromFlatPH, BParam, SParam)
# sol = solve(prob, NelderMead())
# @show sol.original
# BParamDnNew   = [ sol... ]

@time res = optimize( x -> GetCostFromFlatPH(x,SParam), BParam, LBFGS(), Optim.Options(iterations=4000))
@show res
BParamDnNew   = [ res.minimizer... ]

ednnew, vdnnew  = BathParamReshapePH( BParamDnNew, nbathHalf )

DuiwNew      = GetDeltaHybDiscGridFromFlatPH( BParamNew, ImFreqGrid, norb, nbathHalf )

x   = ImFreqGridVal
y1  = Hyb11iw
# y2  = GetijarrayFromVecMat( DiwNew, 1, 1 ) 
y2  = GetijarrayFromVecMat( DuiwNew, 1, 1 ) 

NFreq   = 200
ReFreqGridVal   = collect( LinRange( -3, 3, NReFreq ) )
ReFreqGrid      = ReFreqGridVal .+ im * epsilon
yr1 = GetGzBethe.( ReFreqGrid )
yrD2  = GetDeltaHybDiscGrid( ednnew, vdnnew, ReFreqGrid )
yr2  = GetGreenImpGrid( zeros(1,1), yrD2, ReFreqGrid )
yr2  = GetijarrayFromVecMat( yr2, 1, 1 )

doplot = true
if doplot
    using Plots
    xr  = ReFreqGridVal
    plot( xr, [ real(yr1) imag(yr1) ] )
    plot!( xr, [ real(yr2) imag(yr2) ] )
end
