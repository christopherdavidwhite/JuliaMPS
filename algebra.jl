Algebra = MPO #this is really a vector space (ish)

function element(v :: Array{Complex{Float64}}, A :: Algebra)
    el = deepcopy(A)

    assert(A.χ[0] == length(v))

    v = reshape(v, (1, length(v)))
    W1 = A.W[1]
    W1p = zeros(Complex{Float64}, (A.d,A.d, 1, A.χ[1]))
    @tensor W1p[s,sp, al, ar] = v[al, g]*W1[s,sp,g,ar]
    
    A.W[1] = W1p
    A.χ[0] = 1
    return A
end

function ⊕(A :: Algebra, B:: MPO)
    # almost total overlap with cApdB
    # should be able to re-write
    
    assert(A.L == B.L)
    assert(A.d == B.d)
    L = A.L
    d = A.d 
    
    C = mpo(L,d)
    for j in 1:L
        C.χ[j-1] = A.χ[j-1] + B.χ[j-1]
        C.χ[j]   = A.χ[j]   + B.χ[j]
        
        W = zeros(Complex{Float64},(d,d,C.χ[j-1],C.χ[j],))

        W[:,:,1:A.χ[j-1],1:A.χ[j]] = A.W[j]
        W[:,:,A.χ[j-1]+1:end,A.χ[j] + 1:end] = B.W[j]
        AWj = A.W[j]
        BWj = B.W[j]
        C.W[j] = W
    end
    
    rbc = reshape([1 1],2,1)
    
    WL = C.W[L]
    WLp = zeros(Complex{Float64},(d,d,C.χ[L-1],1))
    @tensor WLp[s,sp,al,ar] = WL[s,sp,al,g]*rbc[g,ar]
    C.W[L] = WLp

    C.χ[0] = 1 + A.χ[0]
    C.χ[L] = 1
    return C
end


# To check:
#
# 1. Take rf heis length 10
# 2. compute spectrum
# 3. compute chebyshev_space (call this function)
# 4. compute spectra
# 5. check that spectra(chebyshev(H)) = chebyshev(spectra(H))
#

function chebyshev_space(H :: MPO, n :: Int, χmax = 0, verbose :: Bool = false)
    assert(n > 2)

    L = H.L
    d = H.d
    
    Trec = [mpoeye(L, d),H]
    T = mpoeye(L,d)⊕H
    
    for j = 1:(n-2)
        Tnext = 2*H*Trec[end] - Trec[end-1]
        sanity_check(Tnext)
        T = T⊕Tnext
        T     = canonical_form(T,     true, χmax)
        Tnext = canonical_form(Tnext, true, χmax)
        Trec[1] = Trec[2]
        Trec[2] = Tnext

        if verbose
            @show T.χ
            @show Tnext.χ
        end
    end

    return T
end

function trace(A :: Algebra)
    I = eye(A.d)
    C = ones(Complex{Float64},1)
    for j in A.L:-1:1
        Cp = zeros(Complex{Float64}, A.χ[j-1])
        Wj = A.W[j]
        @tensor Cp[al] = I[s,sp] * Wj[s,sp,al,ar] * C[ar]
        C = Cp
    end
    return C
end
