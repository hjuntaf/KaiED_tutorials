include("../src/mybase.jl")


init_println("PARAMETER SETUP")
norb    = 3
dim = 2^norb
@show norb
@show dim
@show typeof(dim)
final_println("PARAMETER SETUP")


init_println("Initialize Hashing Function")
hashfsec   = HashingIntToSecarrBin(;outputlevel=1)
println("\nShowHashfsec() ::")
ShowHashfsec(hashfsec)
final_println("Initialize Hashing Function")


init_println( "HashingInvSecarrBinToAllNparInt" )
hashfinv    = HashingInvSecarrBinToAllNparInt( hashfsec ; outputlevel=1 ) 
println("\nshow hashfinv ::")
@show hashfinv
final_println( "HashingInvSecarrBinToAllNparInt" )


init_println("Initialize Wavefunction")
wf  = WF(dim)
@show wf
final_println("Initialize Wavefunction")



init_println("(1) WF.BasisBinstate : " )
println( wf.BasisBinstate )

init_println( "(2) WF.Probamp : " )
println( wf.Probamp )

init_println( "(3) Another copy of wf : ")
ipar    = 2
npar    = ipar-1
wf  = WF(length(hashfsec[ipar]), hashfsec[ipar])
@show wf
@show hashfsec[ipar]

init_println( "(4) After changing a ibinst-value of hashfsec[ipar] : ")
hashfsec[ipar][1] = 999
@show wf
@show hashfsec[ipar]

init_println( "(5) wf2 :" )
wf2 = WF(wf) 
@show wf2
