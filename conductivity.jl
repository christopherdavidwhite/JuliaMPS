#functions for measuring conductivity

function single_μ(Tn :: MPO, Tm :: MPO, jj :: MPO)
    Tnmps = convert(MPS, Tn)#does this reallocate?
    Tmmps = convert(MPS, Tm)#does this reallocate?
    return matrix_element(Tnmps, jj, Tmmps)[]
end

# note that these are different! One takes j⊗j, the other takes just j!
# This is not so good. I do it so that I can run through the same code
# using full arrays or MPOs.

function single_μ(Tn :: Array{Complex{Float64},2}, Tm :: Array{Complex{Float64},2}, j :: Array{Complex{Float64},2})
    return trace(Tn*j*Tm*j)
end

function next_T(Trec :: Array{MPO,1}, H :: MPO, χmax)
    T = canonical_form(2.0*H*Trec[end],   preserve_mag = true, χmax = χmax)
    T = canonical_form(T - Trec[end-1], preserve_mag = true, χmax = χmax)
    Trec[1] = Trec[2]
    Trec[2] = T
    return T, Trec
end

function next_T(Trec :: Array{Array{Complex{Float64},2},1}, H :: Array{Complex{Float64},2}, χmax)
    T = 2*H*Trec[end] - Trec[end-1]
    Trec = [Trec[2], T]
    return T, Trec
end

# note: if H is an MPO, jj should be an mpo current operator
# outer product itself j⊗j^T. 
# If H is an Array{Complex{Float64},2}, it should be a like array just current operator.
# This is confusing and weird. It gets passed into single_μ.
function all_μ(H, jj, N :: Int, Trec, χmax :: Int, prog_per :: Int = 10)
    assert(H==Trec[2])
    μ = OffsetArray(zeros(Complex{Float64}, (N,N)), (0:N-1, 0:N-1))
    I = Trec[1]
    μ[0,0] = single_μ(I,I, jj)
    μ[1,0] = μ[0,1] = single_μ(I, H, jj)
    μ[1,1] = single_μ(H, H, jj)

    if prog_per != 0
        tic()
    end

    for n in 2:N-1
        Tn, Trec = next_T(Trec, H, χmax)
        μ[n, n] = single_μ(Tn, Tn, jj)

        inner_Trec = copy(Trec) #unclear how much, exactly, this copies
        for m in n+1:N-1
            Tm, inner_Trec = next_T(inner_Trec, H, χmax)
            μ[n,m] = μ[m,n] = single_μ(Tn, Tm, jj)
            μ[n,m] = μ[m,n] = single_μ(Tn, Tm, jj)
        end

        if prog_per != 0 && n % prog_per == 0
            @show n, toq()
            tic()
        end
    end
    
    return μ
end
