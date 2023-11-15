include("../src/mybase.jl")


init_println("PARAMETER SETUP")
norb    = 2
nspin   = 2
nspinorb    = norb*nspin
ntot    = nspinorb
dim     = 2^ntot
@show dim
final_println("PARAMETER SETUP")

init_println("Initialize Hashing Function")
outputlevel = 2
hashfsec    = HashingIntToSecarrBin(;norb=ntot, outputlevel=outputlevel)
hashfinvsec = HashingInvSecarrBinToAllNparInt( hashfsec ; outputlevel=outputlevel)
hashfinv    = HashingInvSecarrBinToInt( hashfsec ; outputlevel=outputlevel)
hashfall    = HashingIntToBinall(;norb=ntot, outputlevel=outputlevel)
hashfallinv = HashingInvBinallToInt( hashfall ; outputlevel=outputlevel)
final_println("Initialize Hashing Function")

## Ordering ##
init_println("Initialize Hashing Function Ordered")
hashfallOrdered     = HashingIntToBinallOrdered(;norb=ntot)
@show hashfallOrdered     
hashfallOrderedInv  = HashingInvBinallToInt( hashfallOrdered ; outputlevel=outputlevel)
dimSub      = length(hashfallOrdered)
HashF       = hashfallOrdered
HashFInv    = hashfallOrderedInv


init_println("Start Operations")
tij = OpCA( -1, [1, 1] )
cdi = OpC( 10, [3] )
ci  = OpA(-10, [2] )
@show tij
@show cdi
@show ci

wfall   = WF(dim, hashfallOrdered) 
v0      = collect(1:dim)
wfall.Probamp = v0
wfall2  = WF( wfall )
wfall.HashfInv   = hashfallOrderedInv
wfall2.HashfInv  = hashfallOrderedInv
# OpCAWFadd( tij, wfall, wfall2 )
# OpCWFadd( cdi, wfall, wfall2 )
# OpAWFadd( ci, wfall, wfall2 )

# @show sparse(wfall.Probamp)
# @show sparse(wfall2.Probamp)

isectorN3   = 4
Nel         = 3
dimN3       = length(hashfsec[isectorN3])
wfN3_1      = WF(dimN3, hashfsec[isectorN3], hashfinv, isectorN3, Nel )
v0          = collect(1:dimN3)
wfN3_1.Probamp = v0
@show sparse(wfN3_1.Probamp)


println( "" )
println( "wfN3_2 :: " )
wfN3_2  = WF(wfN3_1)
@show sparse(wfN3_2.Probamp)

println( "OpCAWFadd( tij, wfN3_1, wfN3_2 ) " )
OpCAWFadd( tij, wfN3_1, wfN3_2 ; outputlevel= 2)
@show sparse(wfN3_2.Probamp)

println( "OpCWFadd( cdi, wfN3_1, wfN3_2 ) " )
OpCWFadd( cdi, wfN3_1, wfN3_2 )
@show sparse(wfN3_2.Probamp)

println( "OpAWFadd(  ci, wfN3_1, wfN3_2 ) " )
OpAWFadd( ci, wfN3_1, wfN3_2 )
@show sparse(wfN3_2.Probamp)

println( "" )
println( "ShowWFBasis(wfN3_1) ") 
ShowWFBasis(wfN3_1)


println( "" )
println( "wfN4 :: "  )
isectorN4   = 5
Nel         = 4
dimN4       = length(hashfsec[isectorN4])
wfN4_1      = WF(dimN4, hashfsec[isectorN4], hashfinv, isectorN4, Nel )

println( "OpCWFadd( cdi, wfN3_1, wfN4_1 ) " )
OpCWFadd( cdi, wfN3_1, wfN4_1 )
@show sparse(wfN4_1.Probamp)

println( "" )
println( "ShowWFBasis(wfN4_1) ") 
ShowWFBasis(wfN4_1)


println( "" )
println( "wfN2 :: "  )
isectorN2   = 3
Nel         = 2
dimN2       = length(hashfsec[isectorN2])
wfN2_1      = WF(dimN2, hashfsec[isectorN2], hashfinv, isectorN2, Nel )

println( "OpAWFadd( ci, wfN3_1, wfN2_1 ) " )
OpAWFadd( ci, wfN3_1, wfN2_1 )
@show sparse(wfN3_1.Probamp)
@show sparse(wfN2_1.Probamp)

println( "" )
println( "ShowWFBasis(wfN2_1) ") 
ShowWFBasis(wfN2_1)
