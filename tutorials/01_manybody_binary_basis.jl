include("../src/mybase.jl")


function TestBinaryBasisCreAnni()
    # println( "(x, i), format('{:s} {:s} {:s}',getbittaili(x), getbittaili(cdagi_zbase(x,i)), getbittaili(ci_zbase(x,i)))" )
    for binary_st in 0:3
        x   = binary_st
        println("")
        for i in 0:3
            # println( (x, i), format("{:s} {:s} {:s}",getbittaili(x), getbittaili(cdagi_zbase(x,i)), getbittaili(ci_zbase(x,i))) )
            cdx = cdagi_zbase(x,i)
            cx  = ci_zbase(x,i)
            println( "INPUT ::" )
            println( ("     |state $x> (for $i-th orbital operations) with isocc_zbase(), isemp_zbase(), NFerSwap(), Ntot() ") )
            println( ("cdagi|state $x> (for $i-th orbital operations) with isocc_zbase(), isemp_zbase(), NFerSwap(), Ntot() ") )
            println( ("   ci|state $x> (for $i-th orbital operations) with isocc_zbase(), isemp_zbase(), NFerSwap(), Ntot() ") )
            println( "OUTPUT ::" )
            println( format("  |{:s}> {:s} {:s} ; nswap={:d} ntot={:d}",getbittaili(x),   isocc_zbase(x,i)  , isemp_zbase(x,i)  , NFerSwap(x,i)   , OpNtot(x)   ))
            println( format("  |{:s}> {:s} {:s} ; nswap={:d} ntot={:d}",getbittaili(cdx), isocc_zbase(cdx,i), isemp_zbase(cdx,i), NFerSwap(cdx,i) , OpNtot(cdx) ))
            println( format("  |{:s}> {:s} {:s} ; nswap={:d} ntot={:d}",getbittaili(cx),  isocc_zbase(cx,i) , isemp_zbase(cx,i) , NFerSwap(cx,i)  , OpNtot(cx)  ))
            println("")
        end
    end
    return 0
end

function TestBinaryBasisNtot()
    # println( "(x, i), format('{:s} {:s} {:s}',getbittaili(x), getbittaili(cdagi_zbase(x,i)), getbittaili(ci_zbase(x,i)))" )
    for binary_st in 0:3
        x   = binary_st
        println("")
        for i in 0:3
            cdx = cdagi_zbase(x,i)
            cx  = ci_zbase(x,i)
            println( "INPUT ::" )
            println( ("     |state $x> (for $i-th orbital operations)", "isocc_zbase() isemp_zbase(), Ntot()") )
            println( ("cdagi|state $x> (for $i-th orbital operations)", "isocc_zbase() isemp_zbase(), Ntot()") )
            println( ("   ci|state $x> (for $i-th orbital operations)", "isocc_zbase() isemp_zbase(), Ntot()") )
            println( "OUTPUT ::" )
            println( format("  |{:s}> {:s} {:s} ; ntot={:d} ",getbittaili(x),   isocc_zbase(x,i)  , isemp_zbase(x,i)  , OpNtot(x)   ))
            println( format("  |{:s}> {:s} {:s} ; ntot={:d} ",getbittaili(cdx), isocc_zbase(cdx,i), isemp_zbase(cdx,i), OpNtot(cdx) ))
            println( format("  |{:s}> {:s} {:s} ; ntot={:d} ",getbittaili(cx),  isocc_zbase(cx,i) , isemp_zbase(cx,i) , OpNtot(cx)  ))
            println("")
        end
    end
    return 0
end

using Test 
@testset "Binary test" begin 
    @test TestBinaryBasisCreAnni() == 0 
    end ;
