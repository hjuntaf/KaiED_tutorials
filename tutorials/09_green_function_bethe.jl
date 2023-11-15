include("../src/mybase.jl")

epsilon = 0.01

NFreq   = 30
ReFreqGridVal   = collect( LinRange( -3, 3, NFreq ) ) 
ReFreqGrid      = ReFreqGridVal .+ im * epsilon

G0w = GetGzBethe.( ReFreqGrid )

Hybw    = 1. / 4 * G0w

norb    = 1
tij     = zeros(norb,norb)
tij[1,1]= 0

@show typeof( G0w )
@show typeof( Hybw )
G0wNew  = GetGreenImpGrid( tij[1,1], Hybw, ReFreqGrid )


bPlot = true
if bPlot
    xw  = ReFreqGrid
    gw  = G0w
    gw2 = G0wNew
end


