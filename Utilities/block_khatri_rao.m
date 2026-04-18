function X = block_khatri_rao(C,U)

R = size(C,2);    % number of blocks

X = [];
for r = 1 : R
    cr = C(:,r);
    Ur = U{1,r};
    X = [X  kron(C(:,r),Ur)];

end
    
end