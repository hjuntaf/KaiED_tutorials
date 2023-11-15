include("../src/mybase.jl")

### Impurity Setup ###
norb        = 1
nspin       = 2
nspinorb    = norb*nspin

nbath       = 10
@inline IndHamilBath( ibath )            = nspinorb + ibath 
@inline IndHamilBathSpinup( ibathorb )   = nspinorb + 2*ibathorb - 1
@inline IndHamilBathSpindn( ibathorb )   = nspinorb + 2*ibathorb 

ntot    = nspinorb+nbath
dim     = 2^ntot
@show ntot
@show dim


### BathParameters Setup ###
nbathHalf   = div(nbath,2)
ebathl      = zeros(nbath)
if nbath > 2 
    ebathl[1:2:end] = collect(LinRange( -1, 1, nbathHalf ))
    ebathl[2:2:end] = collect(LinRange( -1, 1, nbathHalf ))
else 
    ebathl[1]   = 0
    ebathl[2]   = 0
end
@show ebathl
ebathlUp    = ebathl[1:2:end]
ebathlDn    = ebathl[2:2:end]

Vil     = zeros(nspinorb,nbath)
for i in 1:nspinorb
    for j in 1:nbath
        Vil[i,j]    = 3 - ebathl[j]^2 
    end
end
# @show Vil
println( "Vil : " ) ; writedlm(stdout, Vil)

IndOrbUp    = [ i for i in 1:2:2*norb ]
IndOrbDn    = [ i for i in 2:2:2*norb ]
IndBathUp   = [ i for i in 1:2:2*nbathHalf ]
IndBathDn   = [ i for i in 2:2:2*nbathHalf ]
VilUp   = Vil[ IndOrbUp, IndBathUp ]
VilDn   = Vil[ IndOrbDn, IndBathDn ]
println( "Vil Up : " ) ; writedlm(stdout, VilUp)
println( "Vil Dn : " ) ; writedlm(stdout, VilDn)



#### Green Function Setup ####

beta    = 32
NImFreq  = 4*beta
ImFreqGridVal   = GetImFreqValGrid( beta, NImFreq )
ImFreqGrid      = ImFreqGridVal * im

G011iw= GetGzBethe.( ImFreqGrid )
Hyb11iw   = 1. / 4 * G011iw

Hybiw        = 1. / 4 * GetGzBetheUniformScaling.( ImFreqGrid ) ## Matrix


#### Bath Optimization Setup ###
println("## Spin-up optimization ##")
DimBathHalfMid  = div( nbathHalf+1 , 2 )
@show DimBathHalfMid
BParamUp    = BathParamFlatten( ebathlUp, VilUp )
BParamUpPH  = BathParamFlatten( ebathlUp[1:div(end,2)], VilUp )
enew, vnew  = BathParamReshapePH( BParamUpPH, nbathHalf )
println( "" )
println( "Initial-parameters :")
println( "ebathl after PH : $(enew)" ) 
println( "Vil Up after PH : " ) ; writedlm(stdout, vnew)
println( "" )
BParam      = BParamUpPH
SParam      = ( 
                Hybiw, 
                ImFreqGrid, 
                norb, # nspinorb, 
                nbathHalf
                )
@show BParam


### Optimizer ##
using Optimization
prob = OptimizationProblem(GetCostFromFlatPH, BParam, SParam)
using OptimizationOptimJL
maxiters = 2300
sol = solve(prob, NelderMead() ; maxiters=maxiters)
@show sol.original

BParamNew   = [ sol... ]
enew, vnew  = BathParamReshapePH( BParamNew, nbathHalf )
println("New parameters for Up")
@show enew
writedlm(stdout, vnew)



println("## Spin-dn optimization ##")
## Setting BParam ##
BParamDn    = BathParamFlatten( ebathlDn, VilDn )
BParamDnPH  = BathParamFlatten( ebathlDn[1:div(end,2)], VilDn )
BParam      = BParamDnPH
println("Initial-parameters :")
@show BathParamReshapePH( BParam, nbathHalf )

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

prob = OptimizationProblem(GetCostFromFlatPH, BParam, SParam)
sol = solve(prob, NelderMead())
@show sol.original
BParamDnNew   = [ sol... ]
ednnew, vdnnew  = BathParamReshapePH( BParamDnNew, nbathHalf )
println("New parameters for Dn")
@show ednnew
writedlm(stdout, vdnnew)

println("New parameters in total")
ebathl  = zero(ebathl)
Vil     = zero(Vil)
ebathl[1:2:end]             = deepcopy(enew)
ebathl[2:2:end]             = deepcopy(ednnew)
ebathl[2:2:end]             = deepcopy(enew)        ## Impose time-reversal sym naively
Vil[ IndOrbUp, IndBathUp ]  = deepcopy(vnew)
Vil[ IndOrbDn, IndBathDn ]  = deepcopy(vdnnew)
Vil[ IndOrbDn, IndBathDn ]  = deepcopy(vnew)        ## Impose time-reversal sym naively
@show ebathl
writedlm(stdout, Vil)



### Binary Basis ###
println("## Hashing index to binary-basis ##" )
outputlevel = 0
hashfsec    = HashingIntToSecarrBin(;norb=ntot, outputlevel=outputlevel)
hashfinvsec = HashingInvSecarrBinToAllNparInt( hashfsec ; outputlevel=outputlevel)
hashfinv    = HashingInvSecarrBinToInt( hashfsec ; outputlevel=outputlevel)
hashfall    = HashingIntToBinall(;norb=ntot, outputlevel=outputlevel)
hashfallinv = HashingInvBinallToInt( hashfall ; outputlevel=outputlevel)
println("## end of Hashing index ##" )



## Ordering ##
println( "\n## Ordering hashfall ##")
hashfallOrdered     = vcat( hashfsec... )
# @show hashfallOrdered     
hashfallOrderedInv  = HashingInvBinallToInt( hashfallOrdered ; outputlevel=outputlevel)
dimSub      = length(hashfallOrdered)
HashF       = hashfallOrdered
HashFInv    = hashfallOrderedInv
wfargs      = [ HashFInv ]
# isector     = 7
# HashF       = hashfsec[isector]
# HashFInv    = hashfinv
# dimSub      = length(HashF)
# wfargs      = [ HashFInv, isector, isector ]

H   = Hamil(dimSub)
H.MatSparse = H.MatSparse * 0.


### Many-body Operators Setup ###
IndHamilBath( indbath )   = nspinorb + indbath

U   = 0.
Op_Uijkl   = OpCCAA[]
if U != 0.0
    push!( Op_Uijkl, OpCCAA( U, [1,2,2,1] ) )
end

Op_Vil  = OpCA[]
for i in 1:nspinorb
    for j in 1:nbath
        jbath   = IndHamilBath(j)
        if Vil[i,j] != 0.0 
            push!( Op_Vil, OpCA(      Vil[i,j] , [i,     jbath] ) )
            push!( Op_Vil, OpCA( conj(Vil[i,j]), [jbath, i    ] ) )
        end
    end
end

Op_ebathl  = OpCA[]
for i in 1:nbath
    ibath   = IndHamilBath(i)
    if ebathl[i] != 0.0
        push!( Op_ebathl, OpCA( ebathl[i], [ibath,ibath] ) )
    end
end

Op_tij = OpCA[]
for i in 1:nspinorb
    chem    = -U / 2
    if chem != 0.0
        push!( Op_tij, OpCA( - U / 2 , [ i, i ] ) )
    end
end
## Zeeman Field for TEST ##
# for i in 1:norb
#     push!( Op_tij, OpCA( 100 , [ 2*i, 2*i ] ) )     # for Spin-Dn
# end

@show Op_Uijkl
@show Op_Vil
@show Op_ebathl
@show Op_tij

outputlevel = 0
for opccaa in Op_Uijkl
    ConstructHamil_ij!( H, opccaa, HashF, HashFInv ; outputlevel=outputlevel )
end
for opca in Op_Vil
    ConstructHamil_ij!( H, opca, HashF, HashFInv ; outputlevel=outputlevel)
end
for opca in Op_ebathl
    ConstructHamil_ij!( H, opca, HashF, HashFInv ; outputlevel=outputlevel)
end
for opca in Op_tij
    ConstructHamil_ij!( H, opca, HashF, HashFInv ; outputlevel=outputlevel)
end

@show issparse(H.MatSparse)
@show H.MatSparse
@show ishermitian(H.MatSparse)
@show issymmetric(H.MatSparse)
@show count(!iszero, H.MatSparse)
@show size(H.MatSparse)

using Arpack

nev     = 10
esAr    = eigs( H.MatSparse ; which=:SR , nev = nev )

println("Eigensystem from eigs()") 
@show typeof(esAr)
@show size(esAr[2])
println("Eigenvalues from eigs()") 
@show typeof(esAr[1])
@show real(esAr[1])
println("diff : ", real(esAr[1] .- esAr[1][1]) )

println("")
println("Ground-states : ")
println( format("nconv = {}, niter = {} , nmult = {}" ,esAr[3], esAr[4], esAr[5] ) )
nconvAr = esAr[3]
gsarr   = Vector{Any}(undef,nconvAr)
OpSz        = GetOpSz(ntot)
OpSzOrb     = GetOpSz(nspinorb)
Sztotval    = Vector{ComplexF64}(undef,nconvAr)
Szorbval    = Vector{ComplexF64}(undef,nconvAr)
for i in 1:nconvAr
    gsarr[i]    = esAr[2][:,i]
    prob_i, indBin  = findmax( abs.(gsarr[i]) )
    println( format("GS $i : Maximal basis states = {} [{}] N={} p={}",
                HashF[indBin], getbittaili(HashF[indBin]), OpNtot(HashF[indBin]), prob_i ))
end
println("")
for i in 1:nconvAr
    gsarr[i]    = esAr[2][:,i]
    prob_i, indBin  = findmax( abs.(gsarr[i]) )
    println( format("GS $i : Maximal basis states = {} [{}] N={} p={}",
                HashF[indBin], getbittaili(HashF[indBin]), OpNtot(HashF[indBin]), prob_i ))
    gswf    = WF(dimSub, HashF, gsarr[i], wfargs... )
    gswfres = WF(dimSub, HashF, wfargs... )
    gswf.Probamp    = collect( dropzeros!( sparse(gswf.Probamp) ) )


    for opca in OpSz
        OpCAWFadd( opca, gswf, gswfres ; outputlevel = 0)  
    end
    Sztot_i = WFInner( gswf, gswfres )
    println( "    : <Sz> = ", Sztot_i )

    for opca in OpSzOrb
        OpCAWFadd( opca, gswf, gswfres ; outputlevel = 0)  
    end
    gswfres.Probamp = zero(gswfres.Probamp)
    Szorb_i = WFInner( gswf, gswfres )
    println( "    : <Sz_orb> = ", Szorb_i )

    Sztotval[i]    += Sztot_i
    Szorbval[i]    += Szorb_i
end

@show real(esAr[1])
println("diff : ", real(esAr[1] .- esAr[1][1]) )
println( "Szval = ", real(Sztotval) , " , ", real(Szorbval) )


# esH = eigen(collect(H.MatSparse))


Eimp0   = esAr[1][1]
wf      = WF(dimSub, HashF, gsarr[1], wfargs... )
ciwf    = WF(dimSub, HashF, wfargs... )
cdagiwf    = WF(dimSub, HashF, wfargs... )

# v0  = esH.vectors[1,:]
# wf  = WF( dimSub, HashF, v0, HashFInv )
# ciwf= WF( dimSub, HashF, HashFInv )

# oparr  = OpA[]
# for i in nspinorb
#     push!( oparr, OpA( 1, [i] ) ) 
# end
ci      = OpA( 1, [1] ) 
OpAWFadd( ci, wf, ciwf ) 
cdagi   = OpC( 1, [1] ) 
OpCWFadd( cdagi, wf, cdagiwf ) 

# @show sparse(v0)
# @show sparse(wf.Probamp)
# @show sparse(ciwf.Probamp)
# ShowWFBasis(ciwf)

# normalize!(ciwf.Probamp)
cv1  = ciwf.Probamp
# @show cv1

# normalize!(cdagiwf.Probamp)
cdv1  = cdagiwf.Probamp
# @show cv1


m   = 700
Tc,  V = lanczos_algorithm(H.MatSparse, cv1,  m)
Tcd, V = lanczos_algorithm(H.MatSparse, cdv1, m)
# T_orig  = deepcopy(T)
#
# @show sparse(T)
# @show typeof(T)
# @show size(T)
# @show typeof(V)
# @show size(V)


# Choose an energy value for the Green's function
energy = 0.0  # Replace this with your desired energy value

# Calculate the continued fraction representation of the Green's function
# G = continued_fraction(T, m, energy + 0.001 *im)

# println("Green's function at energy ", energy, ": ", G)


# println("")
# println("Test CF-number convergence")
# xdat = collect(2:m)
# ydat = zeros(ComplexF64,length(xdat))
# for i_m in 1:m-1
#     m_i = xdat[i_m]
#     T   = deepcopy(T_orig)
#     ydat[m_i]   = continued_fraction(T, i_m, energy + 0.001 *im)
#     println( "$(m_i) : $(ydat[m_i]) " )
# end


println("")
println("## Spectrum in energy-window ##")
println("")
println("GS energy : $(Eimp0) " )

ydat1   = GetGreenFromHCFcoeffA(    H.MatSparse, cv1,  ImFreqGrid, Eimp0 )
ydat2   = GetGreenFromHCFcoeffAdag( H.MatSparse, cdv1, ImFreqGrid, Eimp0 )
ydat    = ydat1+ydat2

yrD2  = GetDeltaHybDiscGrid( enew, vnew, ImFreqGrid )
yr2  = GetGreenImpGrid( zeros(1,1), yrD2, ImFreqGrid )
yr2  = GetijarrayFromVecMat( yr2, 1, 1 )

bPlot=true
if bPlot
    using Plots
    # writedlm( "log", [  ImFreqGridVal imag(ydat1) imag(ydat2) imag(ydat1+ydat2) ] ) 
    plot( ImFreqGridVal, [ imag(ydat) imag(yr2) ] )
end


