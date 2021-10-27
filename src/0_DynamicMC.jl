### ARRAY DYNAMICS SIMULATION ###
# input:        -
# output:       functions for Gillespie simulation of array dynamics
# description:  fast Julia functions for Gillespie simulation
# author:       HR

using Random
using Distributed
using ProgressMeter
using RCall
using HDF5
import Distributions: Uniform

println("------------------------------------------------")
println("DYNAMIC MONTE CARLO SIMULATION OF ARRAY DYNAMICS")
println("------------------------------------------------")


# ----- import variables from R -----
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
@rget outfol

R"source(paste0('src/', lattice, '.R'))"
@rget n;
@rget L;
@rget adjL


# ----- initialise A and M ------

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



# ----- calculate activity changes -----

function calcΔA(A)
    ΔA = @. 1-2A
    return ΔA
end


# ----- determining coupling energies in a lattice -----

function calcJ(i, A, adjL, n, J)

    if lattice == "Kagome"
        k = collect(skipmissing(adjL[i,:]))
        j = @. 2J*(2A[k]-1)
        return sum(j)

    else

        k = adjL[i,:]
        j = sum(@. 2J[1]*(2A[collect(skipmissing(k[1:2]))]-1)) +
            sum(@. 2J[2]*(2A[collect(skipmissing(k[3:4]))]-1)) +
            sum(@. 2J[3]*(2A[collect(skipmissing(k[5:6]))]-1))
        return j

    end
end



# ----- energy of activity state -----

function calcΔH(ΔA, A, M, c, k0, k1, Kon_, n, J, ln, idx)

    idx = idx == nothing ? collect(1:n) : idx
    ln = ln == nothing ? log((1+c)/(1+c/Kon_)) : ln

    Δf = @. k0-k1*M[idx]+ln
    H_int = [calcJ(i, A, adjL, n, J) for i in idx]
    ΔH = @. ΔA[idx]*(Δf-H_int)

    return ΔH
end


# ----- simulation -----

@everywhere begin
    function SIMULATION(adjL, k0, k1, c, Kon_, J, r_0, n, Mtot, met, rx)

        Random.seed!(rx)

        # --- initialise ---
        A = initialiseA(n)
        M = initialiseM(n, Mtot, met)
        ΔA = calcΔA(A)

        # --- record states ---
        ASimResults = zeros(Int64, (Int(TIME), n))
        ASimResults[1, :] = A

        MSimResults = zeros(Int64, (Int(TIME), n))
        MSimResults[1, :] = M

        Ts = zeros(Float64, Int(TIME))

        # --- calculate initial rates ---
        # activity rates
        ln = log((1 + c) / (1 + c/Kon_))
        ΔH = calcΔH(ΔA, A, M, c, k0, k1, Kon_, n, J, ln, nothing)
        rateA = @. r_0*exp(-0.5ΔH)

        # methylation rates only where methylation / demethylation can happen
        rateM = @. kR*(1-A)*(M < Mtot) - kB*A*(M > 0)

        # --- iteration ---
        pb = Progress(Int(TIME))
        t = 0
        cnt = 2

        while t < TIME

            update!(pb, cnt)

            # --- sum of all rates ---
            r = vcat(rateA, abs.(rateM))
            R = sum(r)

            # --- pick time step ---
            u_1 = rand(Uniform(epsilon, 1-epsilon), 1)[1]
            τ = (1/R)*log(1/u_1)

            # --- pick reaction ---
            u_2 = rand(Uniform(epsilon, 1-epsilon), 1)[1]
            j = findmax(cumsum(r) .> u_2*R)[2]  # speed up!

            # --- perform reaction and increment time ---
            # decide which reaction (activity or methylation change) happens

            if j <= length(rateA)

                # change activity for selected site
                A[j] += ΔA[j]
                ΔA[j] = calcΔA(A[j])

                # --- calculate new rates ---
                # update the rates of the neighbors, too
                k = collect(skipmissing(adjL[j,:]))
                idx = vcat(j,k)

                newΔH = calcΔH(ΔA, A, M, c, k0, k1, Kon_, n, J, ln, idx)
                rateA[idx] = @. r_0*exp(-0.5newΔH)

            else

                j -= n
                # decide whether to methylate or to demethylate
                if rateM[j] < 0
                    M[j] -= 1
                elseif rateM[j] > 0
                    M[j] += 1
                end

                # --- calculate new rates ---
                newΔH = calcΔH(ΔA, A, M, c, k0, k1, Kon_, n, J, ln, j)
                rateA[j] = @. r_0*exp(-0.5newΔH)

            end

            # update metylation rates
            rateM[j] = kR*(1-A[j])*(M[j] < Mtot) - kB*A[j]*(M[j] > 0)

            # --- record the state of the system every N steps
            if t >= cnt*dt

                ASimResults[cnt, :] = A
                MSimResults[cnt, :] = M
                Ts[cnt] = t

                cnt += 1
            end

            # --- increment time ---
            t += τ

        end

        SIMresults = Dict("A" => ASimResults,
                            "M" => MSimResults,
                            "Ts" => Ts)

        return SIMresults

    end
end

# ----- execute in parallelised loop -----

TIME = 100
rep = 2

if cores != rep
    println("NUMBER OF THREADS MUST BE EQUAL TO NUMBER OF REPLICATES")
end


Threads.@threads for rx in 1:cores

    SIMresults = SIMULATION(adjL, k0, k1, c, Kon_, J, r_0, n, Mtot, met, rx)

    outfile = string(outfol,"_c",c,"_met-",met,"_rep",rx,".h5")
    h5open(outfile, "w") do file
        write(file, "SIMresults", SIMresults)
    end

end


# ----- put results back to R -----
@rput ALLRESULTS
