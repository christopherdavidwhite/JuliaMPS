push!(LOAD_PATH, "./")
using JuliaMPS
using Base.Test

#test rfheis
L = 5
hj = rand(L)
H = rfheis_W(1, hj) |> mpo
σz = convert(Array{Complex{Float64},2}, [1.0 0.0; 0 -1])
σx = convert(Array{Complex{Float64},2}, [0.0 1.0; 1 0])
σy = convert(Array{Complex{Float64},2}, [0.0 -1.0im; 1im 0])
σzcpt = [real(trace(H*onsite_mpo(σz, j, L))[1])/2^L for j in 1:L]
σzσzcpt = [real(trace(H*onsite_mpo(σz, j, L)*onsite_mpo(σz, j+1, L))[1])/2^L for j in 1:(L-1)]
σxσxcpt = [real(trace(H*onsite_mpo(σx, j, L)*onsite_mpo(σx, j+1, L))[1])/2^L for j in 1:(L-1)]
σyσycpt = [real(trace(H*onsite_mpo(σy, j, L)*onsite_mpo(σy, j+1, L))[1])/2^L for j in 1:(L-1)]
@test σzcpt ≈ hj
@test σzσzcpt ≈ ones(L-1)
@test σxσxcpt ≈ ones(L-1)
@test σyσycpt ≈ ones(L-1)

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
δH = (rfheis_W(1, hj1) |> mpo) - (rfheis_W(1, hj2) |> mpo)
σzcpt = [real(trace(δH*onsite_mpo(σz, j, L))[1])/2^L for j in 1:L]
σzσzcpt = [real(trace(δH*onsite_mpo(σz, j, L)*onsite_mpo(σz, j+1, L))[1])/2^L for j in 1:(L-1)]
@test σzcpt ≈ hj1 - hj2
@test σzσzcpt ≈ zeros(L-1)


#test chebyshev space
L = 5
h = 1
H = rfheis_W(1.0, h*(2*rand(L)-1)) |> mpo
H = H/(3*(L-1) + L*h)
E = H |> full |> eigvals |> real |> sort
n = 5
T = chebyshev_space(H, n)
for m in 1:n
    v = zeros(Complex{Float64},n)
    v[m] = 1
    el = element(v,T)
    fel = full(el)
    d = fel |> eigvals |> real |> sort
    @test isapprox( d,  sort(cos.((m-1)*acos.(E))) , rtol=1e-12)
end

#test chebyshev space
L = 5
h = 1
H = rfheis_W(1.0, h*(2*rand(L)-1)) |> mpo
H = H/(3*(L-1) + L*h)
E = H |> full |> eigvals |> real |> sort
n = 5
traces = chebyshev_traces(H, n)
for m in 0:n-1
    @test abs( sum(traces[m+1]) - sum(cos.(m*acos.(E))) ) < 1e-12
end


# test matrix_element
L = 5
σy = [0 -1.0im; 1.0im 0]
φst = rand([0,1], L)
φ = φst |> prod_σz_eigstate
ψst = copy(φst)
ψst[3] = 1 - ψst[3]
ψ = ψst |> prod_σz_eigstate
@test( (ψst[3] - φst[3])*im == matrix_element(φ, onsite_mpo(σy, 3, L), ψ)[] )


L = 5
H1 = rfheis_W(1, rand(L)) |> mpo
H2 = rfheis_W(1, rand(L)) |> mpo
H3 = rfheis_W(1, rand(L)) |> mpo
H4 = rfheis_W(1, rand(L)) |> mpo

H1f = full(H1)
H2f = full(H2)
H3f = full(H3)
H4f = full(H4)

I = mpoeye(L,2)
H1mps = convert(MPS, H1)
H4mps = convert(MPS, H4)
@test( matrix_element(H1mps, kron(H2,I), H4mps) ≈ trace(H1 * H2 * H4) )

@test( matrix_element(H1mps, kron(H2, H3), H4mps) ≈ trace(H3*H1*H2*H4) )
