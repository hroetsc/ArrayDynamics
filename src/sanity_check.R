
library(tseries)

load("results/SIMresults/benchmarking_RvsJulia/Julia_singleCore/_c0_met-RB+_rep1.RData")
julia_sim = SIMresults

load("results/SIMresults/benchmarking_RvsJulia/R/_c0_met-RB+_rep1.RData")
r_sim = SIMresults

julia_source('src/sanity_check.jl')
julia_randN = randN

set.seed(1)
r_randN = rep(NA, 1e04)
for (i in 1:10000) {
  r_randN[i] = runif(1)
}

plot(density(julia_randN))
lines(density(r_randN), col='red')

plot(density(diff(julia_randN)))
lines(density(diff(r_randN)), col='red')

sd(julia_randN)
sd(r_randN)


r = 1
leap = 10

diffR = rep(NA, length(r_randN)-leap)
diffJulia = rep(NA, length(r_randN)-leap)
while (r <= length(r_randN)-leap) {
  diffR[r] = r_randN[r] - r_randN[r+leap]
  diffJulia[r] = julia_randN[r] - julia_randN[r+leap]
  r = r + 1
}

plot(density(diffR))
plot(density(diffJulia), col='red')


# autocorrelation
x11()
par(mfrow = c(1,2))
acf(r_randN, lag.max = 100, plot = T)
acf(julia_randN, lag.max = 100, plot = T)


