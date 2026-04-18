function X = khatri_rao_b(C,U,L)

R = size(C,2);   

X = [];
index = 1;
for r = 1 : R
    ll = index + L(r) - 1;
    Ur = U(:, index:ll); 
    X = [X  kron(C(:,r),Ur)];
    index = ll+1;     
end
    
end