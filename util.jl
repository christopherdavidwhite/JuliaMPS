import Base.squeeze


squeeze{T,N}(A :: Array{T,N})  = squeeze(A, (find(x -> x == 1, size(A))...))

#squeeze to scalar. Separate function from squeeze because if not,
#return type is indefinite. This seems bad.
function ssqueeze{T,N}(A :: Array{T,N})
    B = squeeze(A, (find(x -> x == 1, size(A))...))
    assert(size(B) == ())
    return B[]
end

function tittums(n :: Int)
    assert(0 == n % 2)
    m = Int(n/2)
    t = zeros(Int, n)
    for j in 1:m
        t[2j-1] = j
        t[2j]   = j+m
    end
    return t
end

function rfheis_frobenius2(J :: Float64, h :: Array{Float64,1})
    #assumes open BC
    L = size(h,1)
    return 2^L*(3*J*(L-1) + sum(h.^2))
end
                    
#this is kind of ridiculous
#can do better either iteratively or (for certain x) by FFT method in Weiss et al.
function chebyshev_polys(x :: Array{Float64,1}, n :: Int)
    assert(all(abs(x .< 1)))
    
    N_pts = length(x)
    
    polys = zeros(N_pts, n)
    
    for m = 0:(n-1)
        polys[:,m+1] = cos(m * acos(x))
    end
    return polys
end

function jackson_kernel(N :: Int)
    m = 0:(N-1)
    return ( (N - m + 1).*cos(π*m/(N+1)) + sin(π*m/(N+1))*cot(π/(N+1)) )/ (N+1)
end


#sparse Pauli matrices
function pauli_matrices_sparse(L :: Int64)
    sigx = sparse([0 1; 1 0])
    sigy = sparse([0 -im; im 0])
    sigz = sparse([1 0; 0 -1])
    sigp = sparse([0 1; 0 0])
    sigm = sparse([0 0; 1 0])
    
    X = [reduce(kron, (speye(2^(j-1)), sigx, speye(2^(L - j)))) for j in 1:L]
    Y = [reduce(kron, (speye(2^(j-1)), sigy, speye(2^(L - j)))) for j in 1:L]
    Z = [reduce(kron, (speye(2^(j-1)), sigz, speye(2^(L - j)))) for j in 1:L]
    
    P = [reduce(kron, (speye(2^(j-1)), sigp, speye(2^(L - j)))) for j in 1:L]
    M = [reduce(kron, (speye(2^(j-1)), sigm, speye(2^(L - j)))) for j in 1:L]
    
    X = convert(Array{SparseMatrixCSC{Float64,Int64}, 1}, X)
    Y = convert(Array{SparseMatrixCSC{Complex{Float64},Int64}, 1}, Y)
    Z = convert(Array{SparseMatrixCSC{Float64,Int64}, 1}, Z)
    P = convert(Array{SparseMatrixCSC{Float64,Int64}, 1}, P)
    M = convert(Array{SparseMatrixCSC{Float64,Int64}, 1}, M)
    return (X,Y,Z,P,M)
end

#sparse random field Heisenberg with specified fields
function rfheis_sparse(J :: Number, hj :: Array{Float64, 1})
    L = length(hj)
    
    (X,Y,Z,P,M) = pauli_matrices_sparse(L)
    H = sum([X[l]*X[l+1] + Y[l]*Y[l+1] + Z[l]*Z[l+1] for l in 1:L-1])
    H += sum([hj[l]*Z[l] for l in 1:L])
    return H
end
