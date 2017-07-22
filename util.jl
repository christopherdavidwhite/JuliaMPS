import Base.squeeze

squeeze{T,N}(A :: Array{T,N}) = squeeze(A, (find(x -> x == 1, size(A))...))
