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
beta    = 100
NImFreq  = 4*beta
ImFreqGridVal   = GetImFreqValGrid( beta, NImFreq )
ImFreqGrid      = ImFreqGridVal * im
# Hybiwup_init  = 1. / 4 * GetGzBetheUniformScaling.( ImFreqGrid )
# Hybiwdn_init  = 1. / 4 * GetGzBetheUniformScaling.( ImFreqGrid )
D   = 1
gbetheiw    = GetGzBetheDim.( ImFreqGrid, nspinorb )
gbetheiwup  = GetGzBetheUniformScaling.( ImFreqGrid )
Hybiw       = 1. / 4 * GetGzBetheDim.( ImFreqGrid, nspinorb )



######## Bath Orbital Info #########
nbath   = 10
ntot    = nspinorb + nbath
dim     = 2^ntot
@show dim
nbathHalf   = div(nbath,2)
# IndBathUp, IndBathDn = GetOrbUpDn( nbath )
IndBathUp   = [ i for i in 1:2:nbath ]
IndBathDn   = [ i for i in 2:2:nbath ]

ebathl     = zeros(nbath)
ebathl     = collect(LinRange( -1, 1, nbath ))
@show ebathl

Vil     = zeros(nspinorb,nbath)
for i in 1:nspinorb
    for j in 1:nbath
        Vil[i,j]    = ebathl[j]
    end
end
println( "Vil : " ) ; writedlm(stdout, Vil)


######## Bath Discretization Setup (Spin-up/dn) #########
ebathlnew, Vilnew   = BathDiscHybSpinPH( ebathl, Vil, Hybiw, ImFreqGrid )
ShowBathParam( ebathlnew, Vilnew )


######## Setting Up Hamiltonian Operators #########
U           = 8
JHund       = 0
chem        = 0.5*U  # 0.5*U for single-band , 2.5*U-5.0*JHund for t2g-multiband
opcavec     = [ GetOpBathParam(ebathlnew, Vilnew, ibath->ibath+nspinorb),
                GetOpChemPot(chem, nspinorb) ]
opccaavec   = [ GetOpUSlaterKanamori( ; U=U, JHund=JHund, norb=norb ) ]
# @show typeof(opcavec)
# @show typeof(opccaavec)
# @show opcavec
# @show opccaavec

######## Fock/Hilbert Space Construction #########
outputlevel = 0

######## Searching Ground-sector #########
Emin_arr    = SearchGSSector( ntot, opcavec, opccaavec ; outputlevel=1 )

######## Choosing Ground-sectors #########
# IndGSSector = [ [ isector, e0, boltzweight ] ]
iGSSector   = argmin(Emin_arr)

######## Accurate Hamiltonian Diagonalization (for the choosen sectors) #########
esyssec_AR = GetGSFromSector( iGSSector, ntot, opcavec, opccaavec ; outputlevel=1 )

######## Boltzmann Weight #########
# IndGSSector = [ [ isector, e0, boltzweight, evec ] ]
ieval        = 1
IndGSSector = [ [ iGSSector, esyssec_AR[1][ieval], 1.0, esyssec_AR[2][:,ieval] ] ]



######## Impurity Green Function Construction #########
Gimp    = GetGreenImpurityFromGS(  nspinorb, IndGSSector[1], ntot, opcavec, opccaavec, ImFreqGrid )
gimpup  = GetijarrayFromVecMat( Gimp, 1, 1 )
gimpdn  = GetijarrayFromVecMat( Gimp, 2, 2 )

Eorb    = collect( I(nspinorb)*0.0 )
G0imp   = GetGreenDiscGrid( ebathlnew, Vilnew, ImFreqGrid, Eorb )
g0impup = GetijarrayFromVecMat( G0imp, 1, 1 )
g0impdn = GetijarrayFromVecMat( G0imp, 2, 2 )

selfup  = GetSelf( g0impup, gimpup, chem )
selfdn  = GetSelf( g0impdn, gimpdn, chem )

G0newimp    = GetGreenImpGrid( zeros(nspinorb,nspinorb), (D*D / 2. / 2. ) * Gimp, ImFreqGrid )

using Plots
xdat    = ImFreqGridVal
plot(  xdat, imag( GetijarrayFromVecMat(gbetheiw,1,1) ) )
plot!( xdat, imag( GetijarrayFromVecMat(G0newimp,1,1) ) )

######## Lattice Green Function Construction #########
# For Bethe lattice model,
# self_w   = 
# Gbethe_w = Get( self_w )

# For TB model,
# Glocal_iw    = Get( tb, self_iw )
# G0new_iw     = inv( inv(Glocal_iw) + self_iw )
# Same things with ReFreq for spectral analysis

