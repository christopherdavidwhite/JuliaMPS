push!(LOAD_PATH, "./")
using JuliaMPS
using Base.Test

#test rfheis
L = 5
hj = rand(L)
H = rfheis_W(1, hj) |> mpo
σz = convert(Array{Complex{Float64},2}, [1.0 0.0; 0 -1])
σzcpt = [real(trace(H*onsite_mpo(σz, j, L))[1])/2^L for j in 1:L]
σzσzcpt = [real(trace(H*onsite_mpo(σz, j, L)*onsite_mpo(σz, j+1, L))[1])/2^L for j in 1:(L-1)]
@test σzcpt ≈ hj
@test σzσzcpt ≈ ones(L-1)

#take an rfheis and test multiplication
L = 5
hj = rand(L)
H = rfheis_W(1, hj) |> mpo
H = H/10
σzcpt = [real(trace(H*onsite_mpo(σz, j, L))[1])/2^L for j in 1:L]
σzσzcpt = [real(trace(H*onsite_mpo(σz, j, L)*onsite_mpo(σz, j+1, L))[1])/2^L for j in 1:(L-1)]
@test σzcpt ≈ hj/10
@test σzσzcpt ≈ ones(L-1)/10

#take an rfheis and test canonical_form (no trunc, preserve_mag)
L = 5
hj = rand(L)
H = rfheis_W(1, hj) |> mpo
H = canonical_form(H, preserve_mag = true, χmax = 0)
σzcpt = [real(trace(H*onsite_mpo(σz, j, L))[1])/2^L for j in 1:L]
σzσzcpt = [real(trace(H*onsite_mpo(σz, j, L)*onsite_mpo(σz, j+1, L))[1])/2^L for j in 1:(L-1)]
@test σzcpt ≈ hj
@test σzσzcpt ≈ ones(L-1)


#take an rfheis and test canonical_form (no trunc, no preserve_mag)
L = 5
hj = rand(L)
H = rfheis_W(1, hj) |> mpo
f = trace(H*H)[1]
@test f ≈ rfheis_frobenius2(1.0, hj)
H = canonical_form(H, preserve_mag = false, χmax = 0)
σzcpt = [real(trace(H*onsite_mpo(σz, j, L))[1])/2^L for j in 1:L]
σzσzcpt = [real(trace(H*onsite_mpo(σz, j, L)*onsite_mpo(σz, j+1, L))[1])/2^L for j in 1:(L-1)]
@test σzcpt ≈ hj/sqrt(f)
@test σzσzcpt ≈ ones(L-1)/sqrt(f)

#take two rfheis and test (-)
L = 5
hj1 = rand(L)
hj2 = rand(L)
δH = (rfheis_W(1, hj1) |> mpo) - (rfheis_W(1, hj2))
f = trace(H*H)[1]
σzcpt = [real(trace(H*onsite_mpo(σz, j, L))[1])/2^L for j in 1:L]
σzσzcpt = [real(trace(H*onsite_mpo(σz, j, L)*onsite_mpo(σz, j+1, L))[1])/2^L for j in 1:(L-1)]
@test σzcpt ≈ hj1 - hj2
@test σzσzcpt ≈ zeros(L-1)
