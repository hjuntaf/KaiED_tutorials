include("../src/mybase.jl")


println("\n######### PARAMETER SETUP ##########")
norb    = 3
dim = 2^norb
@show norb
@show dim
@show typeof(dim)
println("####################################\n")

init_println("(1) Initialize Hashing Function")
hashfsec   = HashingIntToSecarrBin(; norb=norb, outputlevel=2)
println("\nShowHashfsec() ::")
ShowHashfsec(hashfsec)
final_println("(1) Initialize Hashing Function")

init_println("(2) HashingInvSecarrBinToInt")
hashfsecinv    = HashingInvSecarrBinToInt( hashfsec ; outputlevel=2 ) 
println("\nShowHashfInv() ::")
ShowHashfInv(hashfsecinv)
# @show hashfsecinv
final_println("(2) HashingInvSecarrBinToInt")


init_println("(3) HashingInvSecarrBinToAllNparInt")
hashfinvsec  = HashingInvSecarrBinToAllNparInt( hashfsec ; outputlevel=2 )
final_println("(3) HashingInvSecarrBinToAllNparInt")



init_println( "(4) HashingIntToBinall " )
hashfall    = HashingIntToBinall( ; norb=norb, outputlevel=2 ) ## It's just an identity, i.e. x[i] = i
final_println( "(4) HashingIntToBinall " )


init_println("(5) HashingInvBinallToInt")
hashfallinv = HashingInvBinallToInt( hashfall ; outputlevel=2 ) ## It's just an inverse of identity, that is, an identity, x[i] = i
final_println("(5) HashingInvBinallToInt")

