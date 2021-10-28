
using Random
using RCall

R"source('src/argparse.R')"
@rget J;
@rget r_0;
@rget c;
@rget met;
@rget Kon_;
@rget kR;
@rget kB;
@rget k0;
@rget k1;
@rget Mtot;
@rget X;
@rget lattice;
@rget TIME;
@rget dt;
@rget epsilon;
@rget cores;
@rget rep;
@rget outfol

R"source(paste0('src/', lattice, '.R'))"
@rget n;
@rget L;
@rget adjL


Random.seed!(1)

function initialiseA(n)
    A = rand([0 1], n)
    return A
end


function initialiseM(n, Mtot, met)

    if met == "RB+"
        M = rand(0:Mtot,n)
    else
        M = repeat(Mtot, n)
    end
    return M

end

A = initialiseA(n)
M = initialiseM(n, Mtot, met)

@rput A
@rput M



Random.seed!(42)
randN = zeros(Float64, Int(1e04))
for i in 1:10000
    randN[i] = rand(Float64, 1)[1]
end

@rput randN
