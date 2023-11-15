include("../src/mybase.jl")

println("")
println(" Vij operation ( Vij |psi> = |psi'> ) :" )
for x in 1:3 
    println("")
    for i in 1:3 
        for j in 1:3 
            x_new   = act_vij(x,i,j)
            println(format("C({:d}) A({:d}) |{:s}>  = {:s}", i, j, getbittaili(x), getbittaili(x_new) ))
        end
    end
end

println("")
println(" Vijkl operation ( Vijkl |psi> = |psi'> ) :" )
for x in 1:6 
    println("")
    for (i,j,k,l) in [ [1,2,3,4], [1,2,2,1], [1,2,1,2], [1,3,2,1] ]
        x_new   = act_vijkl(x,i,j,k,l)
        println(format("CCAA({}) |{:s}> = {:s}", [i,j,k,l], getbittaili(x), getbittaili(x_new) ))
    end
end
