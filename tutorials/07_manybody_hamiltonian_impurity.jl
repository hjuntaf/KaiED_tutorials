include("../src/mybase.jl")

norb        = 1
nspin       = 2
nspinorb    = norb*nspin

nbath       = 3
@inline IndHamilBath( ibath )            = ibath + nspinorb
@inline IndHamilBathSpinup( ibathorb )   = 2*ibathorb + nspinorb - 1
@inline IndHamilBathSpindn( ibathorb )   = 2*ibathorb + nspinorb

ntot    = nspinorb+nbath
dim     = 2^ntot
@show dim

println("## Hashing index to binary-basis ##" )
outputlevel = 0
hashf    = HashingIntToSecarrBin(;norb=ntot, outputlevel=outputlevel)
hashfinv = HashingInvSecarrBinToAllNparInt( hashf ; outputlevel=outputlevel)
hashfall    = HashingIntToBinall(;norb=ntot, outputlevel=outputlevel)
hashfallinv = HashingInvBinallToInt( hashfall ; outputlevel=outputlevel)
println("## end of Hashing index ##" )


U   = 10.
Op_Uijkl   = OpCCAA[]
push!( Op_Uijkl, OpCCAA( U, [1,2,2,1] ) )

tij     = zeros(norb,norb)
tij[1,1]= 0
Op_tij  = OpCA[]
for i in 1:norb
    for j in 1:norb
        push!( Op_tij, OpCA( tij[i,j], [i,j] ) )
    end
end

Vil     = zeros(norb,nbath)
Vil[1,1]= 3
Op_Vil  = OpCA[]
for i in 1:norb
    for j in 1:nbath
        jbath   = IndHamilBath( j ) 
        push!( Op_Vil, OpCA( Vil[i,j], [i,jbath] ) )
    end
end


ebathl     = zeros(nbath)
ebathl     = LinRange( -1, 1, nbath )
@show ebathl
Op_ebathl  = OpCA[]
for i in 1:nbath
    push!( Op_ebathl, OpCA( ebathl[i], [i,i] ) )
end

iNsector    = 4
dimNsector  = length( hashf[iNsector] )

outputlevel = 1
# H   = Hamil(dim)
H   = Hamil(dimNsector)
H.MatSparse = H.MatSparse * 0.
for opccaa in Op_Uijkl
    ConstructHamil_ij!( H, opccaa, hashf[iNsector], hashfinv ; outputlevel=outputlevel)
end
for opca in Op_Vil
    ConstructHamil_ij!( H, opca, hashf[iNsector], hashfinv ; outputlevel=outputlevel)
end
for opca in Op_ebathl
    ConstructHamil_ij!( H, opca, hashf[iNsector], hashfinv ; outputlevel=outputlevel)
end

@show H.MatSparse
@show issparse(H.MatSparse)
@show ishermitian(H.MatSparse)
@show issymmetric(H.MatSparse)
@show count(!iszero, H.MatSparse)
