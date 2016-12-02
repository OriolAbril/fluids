function l= lmax(Re,k)
[~,eval,~] = ospp(k,Re);
ll=sort(real(eval),'descend');
l=ll(1);
end