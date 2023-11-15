ENV["PROJECT_PATH_ED"]="../envs/KED"
include("../src/mybase.jl")


######## Correlated Orbital Info #########
norb    = 1
nspin   = 2
nspinorb    = norb * nspin

# IndOrbUp, IndOrbDn = GetOrbUpDn( nspinorb )
IndOrbUp    = [ i for i in 1:2:nspinorb ]
IndOrbDn    = [ i for i in 2:2:nspinorb ]


######## Imaginary Frequency Green Function Construction #########
beta    = 32
NImFreq  = 8*beta
ImFreqGridVal   = GetImFreqValGrid( beta, NImFreq )
ImFreqGrid      = ImFreqGridVal * im
# Hybiwup_init  = 1. / 4 * GetGzBetheUniformScaling.( ImFreqGrid )
# Hybiwdn_init  = 1. / 4 * GetGzBetheUniformScaling.( ImFreqGrid )
D   = 1
gbetheiw    = GetGzBetheDim.( ImFreqGrid, nspinorb )
gbetheiwup  = GetGzBetheUniformScaling.( ImFreqGrid )
# Hybiw_init       = 1. / 4 * GetGzBetheDim.( ImFreqGrid, nspinorb )
Hybiw_init       = 1. / 4 * GetGzBetheUniformScaling.( ImFreqGrid )

# for iorb in 1:nspinorb
#     WriteVecMat( real(ImFreqGridVal), GetijarrayFromVecMat( Hybiw_init, iorb, iorb ), "Hybiw_$(iorb)_$(iorb)_init.dat" )
# end


######## Bath Orbital Info #########
nbath   = 14
ntot    = nspinorb + nbath
dim     = 2^ntot
@show dim
nbathHalf   = div(nbath,2)
# IndBathUp, IndBathDn = GetOrbUpDn( nbath )
IndBathUp   = [ i for i in 1:2:nbath ]
IndBathDn   = [ i for i in 2:2:nbath ]

ebathl     = zeros(nbath)
ebathl[1:2:end]     = collect(LinRange( -1, 1, nbathHalf ))
ebathl[2:2:end]     = collect(LinRange( -1, 1, nbathHalf ))
@show ebathl

Vil     = zeros(nspinorb,nbath)
for i in 1:norb
    for j in 1:nbathHalf
        Vil[2*i-1,2*j-1]    = 3 - ebathl[2*i-1]^2
        Vil[2*i  ,2*j  ]    = 3 - ebathl[2*i  ]^2
    end
end
println( "Vil : " ) ; writedlm(stdout, Vil)


######## Bath Discretization Setup (Spin-up/dn) #########

######## Setting Up Hamiltonian Operators #########
U           = 8. # 2
JHund       = 0. # 0.167
chem        = 0.5*U  # 0.5*U for single-band half-filling , 2.5*U-5.0*JHund for t2g-multiband half-filling
opccaavec   = [ GetOpUSlaterKanamori( ; U=U, JHund=JHund, norb=norb ) ]

outputlevel = 0
dirnameData   = "./data/"
run( `mkdir -p $(dirnameData)` )

Hybiw   = deepcopy(Hybiw_init)
ndmftloop   = 1
for idmft in 1:ndmftloop
    init_println( "DMFT loop ($(idmft))" )
    Gimpiw, G0impiw, Gnewiw, Hybiwnew, Selfiwnew, ebathlnew, Vilnew, GSSectorInfo = SolveEDBetheDegen( Hybiw, ebathl, Vil, ImFreqGrid, ntot, nspinorb, chem, opccaavec ; TolEGS=6e-5 )
    global Hybiw    = Hybiwnew
    global ebathl   = deepcopy( ebathlnew )
    global Vil      = deepcopy( Vilnew )

    WriteBathAll( ebathl, Vil ; iiter=idmft )

    # gup = GetijarrayFromVecMat( Gnewiw, 1, 1 )
    # gdn = GetijarrayFromVecMat( Gnewiw, 2, 2 )
    # selfup  = GetijarrayFromVecMat( Selfiwnew, 1, 1 )
    # selfdn  = GetijarrayFromVecMat( Selfiwnew, 2, 2 )

    NReFreq = 500
    epsilon = 0.04
    ReFreqGridVal   = collect( LinRange( -10, 10, NReFreq ) )
    ReFreqGrid      = ReFreqGridVal .+ im * epsilon
    G0imp_w = GetGreenDiscGrid( ebathl, Vil, ReFreqGrid )
    Gimpw, G0impw, Gneww, Hybwnew, Selfwnew, ebathldum, Vildum, GSSectorInfo = SolveEDBetheDegen( Hybiw, ebathl, Vil, ReFreqGrid, ntot, nspinorb, chem, opccaavec ; BathOpt=false, GSSectorInfo=GSSectorInfo )
    SelfwUp  = GetijarrayFromVecMat( Selfwnew, 1, 1 )
    SelfwDn  = GetijarrayFromVecMat( Selfwnew, 2, 2 )
    GbethewUp    = GetGzBetheFromSelf( ReFreqGrid, SelfwUp )
    GbethewDn    = GetGzBetheFromSelf( ReFreqGrid, SelfwDn )

    for iorb in 1:nspinorb
        fnameTailOrb    = "_$(iorb)_$(iorb)_i$(idmft).dat"
        WriteVecMat( real(ImFreqGridVal), GetijarrayFromVecMat( Gimpiw, iorb, iorb ),  dirnameData*"Gimpiw"*fnameTailOrb )
        WriteVecMat( real(ImFreqGridVal), GetijarrayFromVecMat( G0impiw, iorb, iorb ), dirnameData*"G0impiw"*fnameTailOrb )
        WriteVecMat( real(ImFreqGridVal), GetijarrayFromVecMat( Selfiwnew, iorb, iorb ), dirnameData*"Selfiw"*fnameTailOrb )
        WriteVecMat( real(ImFreqGridVal), GetijarrayFromVecMat( Hybiwnew, iorb, iorb ), dirnameData*"Hybnewiw"*fnameTailOrb )
        WriteVecMat( real(ReFreqGridVal),    GetijarrayFromVecMat( Gimpw,  iorb, iorb ), dirnameData*"Gimpw"*fnameTailOrb )
        WriteVecMat( real(ReFreqGridVal),    GetijarrayFromVecMat( G0imp_w, iorb, iorb ), dirnameData*"G0impw"*fnameTailOrb )
        WriteVecMat( real(ReFreqGridVal),    GetijarrayFromVecMat( Selfwnew,  iorb, iorb ), dirnameData*"Selfw"*fnameTailOrb )
        WriteVecMat( real(ReFreqGridVal),    GetijarrayFromVecMat( Hybwnew,  iorb, iorb ), dirnameData*"Hybneww"*fnameTailOrb )
    end
    WriteVecMat( real(ReFreqGrid),    GbethewUp , dirnameData*"Gw_1_1_i$(idmft).dat" )
    WriteVecMat( real(ReFreqGrid),    GbethewDn , dirnameData*"Gw_2_2_i$(idmft).dat" )

    println( "Hybiw[div(end,2)] before ConstrainTrev!()" )
    @show Hybiw[div(end,2)]
    ConstrainTrev!( Hybiw, [ 1, 2 ] )
    println( "Hybiw[div(end,2)] after ConstrainTrev!()" )
    @show Hybiw[div(end,2)]
end

######## Lattice Green Function Construction #########
# For Bethe lattice model,
# self_w   = 
# Gbethe_w = Get( self_w )

# For TB model,
# Glocal_iw    = Get( tb, self_iw )
# G0new_iw     = inv( inv(Glocal_iw) + self_iw )
# Same things with ReFreq for spectral analysis

