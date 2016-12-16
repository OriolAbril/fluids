function detM = fun5c(Ra)
global k
[M,~]=fun5a(k,Ra);
detM=abs(det(M));

end

