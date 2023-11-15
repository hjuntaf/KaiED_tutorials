ENV["PROJECT_PATH_ED"]="../envs/KED"
include("../src/mybase.jl")

norb        = 3
nspin       = 2
nspinorb    = norb*nspin

ntot    = nspinorb
dim     = 2^ntot
@show dim

println("## Hashing index to binary-basis ##" )
outputlevel = 2
hashf    = HashingIntToSecarrBin(;norb=ntot, outputlevel=outputlevel)
hashfinv = HashingInvSecarrBinToAllNparInt( hashf ; outputlevel=outputlevel)
hashfall    = HashingIntToBinall(;norb=ntot, outputlevel=outputlevel)
hashfallinv = HashingInvBinallToInt( hashfall ; outputlevel=outputlevel)
println("## end of Hashing index ##" )


Uijkl   = OpCCAA[]
U       = 2.0
JHund   = 0.0
for i in 1:norb 
    push!( Uijkl, OpCCAA( U,             [2*i-1, 2*i  , 2*i,     2*i-1] ) ) ## t2g : intra-orbital density-density  ##
    for j in 1:norb 
        if i!=j
            push!( Uijkl, OpCCAA( U-2.0*JHund,   [2*i-1, 2*j  , 2*j,     2*i-1] ) ) ## t2g : inter-orbital density-density  ## # i != j , opposite-spin
            push!( Uijkl, OpCCAA( -JHund,        [2*i-1, 2*j,   2*j-1,   2*i  ] ) ) ## t2g : pair-exchange                  ## # i != j , up|dn => dn|up
            push!( Uijkl, OpCCAA(  JHund,        [2*i-1, 2*i,   2*j  ,   2*j-1] ) ) ## t2g : pair-hopping                   ## # i != j , i,up|i,dn => j,up|j,dn
        end
        if j<i
            push!( Uijkl, OpCCAA( U-3.0*JHund,   [2*i-1, 2*j-1, 2*j-1,   2*i-1] ) ) ## t2g : inter-orbital density-density  ## # i <  j , up-spin
            push!( Uijkl, OpCCAA( U-3.0*JHund,   [2*i  , 2*j  , 2*j  ,   2*i  ] ) ) ## t2g : inter-orbital density-density  ## # i <  j , dn-spin
        end
    end
end

# for op in Uijkl
#     println( op.i, op.v )
# end

tij = OpCA[]
for i in 1:nspinorb
    chem    = 2.5*U - 5.0*JHund     ##  (0.5*U-1.5*JHund , 1.5*U-4.5*JHund , 2.5*U-5.0*JHund , 3.5*U-5.5*JHund , 4.5*U-8.5*JHund ) : chemical-potentials for integer-fillings
    push!( tij, OpCA( -chem , [ i, i ] ) )  
end


## Ordering ##
println( "\n## Ordering hashfall ##")
hashfallOrdered     = vcat( hashf... )
@show hashfallOrdered     
hashfallOrderedInv  = HashingInvBinallToInt( hashfallOrdered ; outputlevel=outputlevel)
dimSub      = length(hashfallOrdered)
HashF       = hashfallOrdered
HashFInv    = hashfallOrderedInv

H   = Hamil(dimSub)
H.MatSparse = H.MatSparse * 0.
for opccaa in Uijkl
    ConstructHamil_ij!( H, opccaa, HashF, HashFInv ; outputlevel=2 )
end
for opca in tij
    ConstructHamil_ij!( H, opca, HashF, HashFInv ; outputlevel=2 )
end

@show issparse(H.MatSparse)
@show H.MatSparse
@show ishermitian(H.MatSparse)
@show issymmetric(H.MatSparse)
@show count(!iszero, H.MatSparse)

eigensys    = eigen( collect(H.MatSparse) )

for evitem in Iterators.enumerate( eigensys.values ) 
    @show evitem
end
