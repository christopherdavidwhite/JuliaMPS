import Base.*
import Base.'
import Base./
import Base.+
import Base.-

#right now always in B-form
#s is informational, no need to keep it up to date
type MPO
    W :: Array{Array{Complex{Float64},4}, 1}
    χ :: OffsetArrays.OffsetArray{Int64,1,Array{Int64,1}}
    d :: Int
    L :: Int
    s :: Array{Array{Float64,1},1}
end

function convert(::Type{MPS}, A :: MPO)
    d = A.d
    Wp = [reshape(A.W[j], (d^2,A.χ[j-1],A.χ[j])) for j in 1:A.L]
    MPS(Wp, A.χ, d^2, A.L, A.s)
end

function convert(::Type{MPO}, ψ :: MPS)
    d  = Int(sqrt(ψ.d))
    Wp = [reshape(ψ.W[j], (d,d,ψ.χ[j-1],ψ.χ[j])) for j in 1:ψ.L]
    return MPO(Wp, ψ.χ, d, ψ.L, ψ.s)
end

function canonical_form(A :: MPO; preserve_mag :: Bool = false, χmax :: Int = 0, runtime_check = false)
    ψ = convert(MPS, A)
    ψ = canonical_form!(ψ, preserve_mag = preserve_mag, χmax = χmax, runtime_check = runtime_check)
    A = convert(MPO, ψ)
    return A
end

function mpo(Ws :: Array{Array{Complex{Float64},4}, 1})
    L = length(Ws)
    
    χ = OffsetArray(zeros(Int64, L+1), 0:L)
    χ[0]   = size(Ws[1],3)
    χ[1:L] = [size(W,4) for W in Ws]
    
    d = size(Ws[1], 1)
    for W = Ws
        assert(size(W,1) == d)
        assert(size(W,2) == d)
    end
    s = Array{Array{Float64,1}}(L)
    return MPO(Ws,χ,d,L,s)
end

function mpo(L :: Int64, d :: Int64)
    return MPO(Array{Array{Complex{Float64}, 4}}(L),
               OffsetArray(Array{Int64}(L+1), 0:L),
               d,L,
               Array{Array{Float64,1}}(L))
end

#physical, physical, bond, bond
function sanity_check(op :: MPO, verbose :: Bool = false)
    for j in 1:op.L
        if verbose
            @show j, size(op.W[j])
        end
        assert(size(op.W[j]) == (op.d, op.d, op.χ[j-1], op.χ[j]))
    end
end

function (*)(A :: MPO, c :: Number)
    B = deepcopy(A)
    L = B.L
    for j in 1:L
        B.W[j] = B.W[j]*c^(1/L)
    end
    return B
end

(*)(c :: Number, A :: MPO) = A*c
(/)(A :: MPO, c :: Number) = A*(1/c)


function rfheis_W(J :: Number, h :: Array{Float64,1})
    σ0 = [1 0; 0 1]
    σx = [0 1; 1 0]
    σy = [0 -im; im 0]
    σz = [1 0;0 -1]
    σ = [σ0, σx, σy, σz]
    d = 2
    rbdry = [1,0,0,0]
    lbdry = [0,0,0,1]
    L = length(h)
    Ws = Array{Array{Complex{Float64},4}}(L)
    for l in 1:length(h)
        W = zeros(Complex{Float64}, 2,2,5,5)
        for j in 1:4
            W[:,:,j,1]   = sqrt(J)*σ[j]
        end
        for j in 2:4
            W[:,:,5,j] = sqrt(J)*σ[j]
        end
        
        W[:,:,5,1] = h[l]*σz
        W[:,:,5,5] = σ0
        W[:,:,1,1] = σ0
        Ws[l] = W
    end
    rbc = reshape([1 0 0 0 0;],5,1)
    lbc = reshape([0,0,0,0,1;], (1,5))
    W1 = Ws[1]
    W1p = zeros(Complex{Float64}, (d,d,1,size(W1,4)))
    @tensor W1p[s,sp,al,ar] = lbc[al,g]*W1[s,sp,g,ar]
    Ws[1] = W1p
        
    WL = Ws[L]
    WLp = zeros(Complex{Float64}, (d,d,size(WL,4),1))
    @tensor WLp[s,sp,al,ar] = WL[s,sp,al,g]*rbc[g,ar]
    Ws[L] = WLp
    return Ws
end

function sum_charge_current_W(L :: Int64)
    σ0 = [1 0; 0 1]
    σx = [0 1; 1 0]
    σy = [0 -im; im 0]
    d = 2
    Ws = Array{Array{Complex{Float64},4}}(L)
    for l in 1:L
        W = zeros(Complex{Float64}, 2,2,4,4)

        W[:,:,2,1] = σx
        W[:,:,4,2] = 2*σy
        W[:,:,3,1] = -σy
        W[:,:,4,3] = 2*σx
        
        W[:,:,4,4] = σ0
        W[:,:,1,1] = σ0
        Ws[l] = W
    end
    rbc = reshape([1 0 0 0;],4,1)
    lbc = reshape([0,0,0,1;], (1,4))
    W1 = Ws[1]
    @tensor W1p[s,sp,al,ar] := lbc[al,g]*W1[s,sp,g,ar]
    Ws[1] = W1p
        
    WL = Ws[L]
    @tensor WLp[s,sp,al,ar] := WL[s,sp,al,g]*rbc[g,ar]
    Ws[L] = WLp
    return Ws
end

function cApdB(c :: Number, A :: MPO, d :: Number, B :: MPO)
    return element([c+0.0im,d], A⊕B)
                   
end

(+)(A :: MPO, B :: MPO) = cApdB(1,A,1,B)
(-)(A :: MPO, B :: MPO) = cApdB(1,A,-1,B)

function mpoeye(L :: Int, d :: Int)
    I = eye(Complex{Float64}, d)
    Ws = [reshape(I, (d,d,1,1)) for j in 1:L]
    return mpo(Ws)
end


function (*)(A :: MPO, B :: MPO)
    assert(A.L == B.L)
    assert(A.d == B.d)
    
    C = mpo(A.L,A.d)
    d = A.d
    for j in 1:A.L
        C.χ[j-1] = A.χ[j-1] * B.χ[j-1]
        C.χ[j]   = A.χ[j]   * B.χ[j]
        
        AWj = A.W[j]
        BWj = B.W[j]
        @tensor W[s,sp,al1,al2,ar1,ar2] := AWj[s,z,al1,ar1]*BWj[z,sp,al2,ar2]
        W = reshape(W, (d,d,C.χ[j-1],C.χ[j]))
        C.W[j] = W
    end
    return C
end

function full(A :: MPO)
    L = A.L
    d = A.d
    χ = A.χ
    C = ones(Complex{Float64},(1,χ[0]))
    for j in 1:L
        M = size(C,1)
        Wj = A.W[j]
        Cp = zeros(Complex{Float64}, (M, d,d, χ[j]))
        @tensor Cp[m, s, sp, ar] = C[m,g]*Wj[s,sp,g,ar]
        C = reshape(Cp, (M*d*d, χ[j]))
    end
    D = reshape(C, tuple([d for j in 1:2L]...))
    t = tittums(2L)
    tinv = sortperm(t)
    E = permutedims(D, tinv)
    F = reshape(E, (d^L, d^L))
    return F
end

function onsite_mpo(A :: Array{Complex{Float64},2}, j :: Int, L:: Int)
    d = size(A,1)
    assert(d == size(A,2))
    Ws = [reshape(eye(Complex{Float64}, d), (d,d,1,1)) for k in 1:L]
    Ws[j] = reshape(A, (d,d,1,1))
    return mpo(Ws)
end

function kron(A :: MPO, B :: MPO)
    assert(A.L == B.L)
    L = A.L
    C = Array{Array{Complex{Float64}, 4}}(L)
    for j in 1:L
        AWj = A.W[j]
        BWj = B.W[j]
        @tensor CWj[s,r,sp,rp,al,bl,ar,br] := ( AWj[s,sp,al,ar]
                                              * BWj[r,rp,bl,br]
                                              )
        C[j] = reshape(CWj, (A.d      * B.d,
                             A.d      * B.d,
                             A.χ[j-1] * B.χ[j-1],
                             A.χ[j]   * B.χ[j],
                            ))
    end
    return mpo(C)
end

function matrix_element(φ :: MPS, A :: MPO, ψ :: MPS,)
    assert(φ.L == ψ.L)
    assert(φ.d == ψ.d)
    
    assert(A.d == φ.d)
    assert(A.L == φ.L)


    L = A.L
    
    ψW1  = squeeze(ψ.W[1],2)
    φW1c = conj(squeeze(φ.W[1],2))
    AW1  = A.W[1]
    
    @tensor C[c,l,a,f] := ( φW1c[s,f]
                        * AW1[s,sp,l,a]
                        * ψW1[sp,c]
                        )
    
    for j in 2:L
        
        ψWj = ψ.W[j]
        φWjc = conj(φ.W[j])
        AWj = A.W[j]
        
        #Is this the optimal contraction order?
        #(TensorOperations left-associates)
        @tensor D[c,l,a,f] :=  (    C[cl,l, al,fl]
                                * ψWj[sp,   cl, c]
                                * φWjc[s,   fl, f]
                                * AWj[s,sp, al, a]
                                )
        
        C = D
    end
    return squeeze(C)
end

function transpose(A :: MPO)
    W = Array{Array{Complex{Float64},4}}(A.L)
    for l in 1:A.L
        Wl = A.W[l]
        @tensor Wlp[s, sp, al, ar] := Wl[sp, s, al, ar]
        W[l] = Wlp
    end
    return mpo(W)
end
