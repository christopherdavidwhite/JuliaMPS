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
                    
