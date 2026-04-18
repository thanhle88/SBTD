function AB=kron_onlineBTD(A,B)
    [I, J]=size(A);
    [K ,L]=size(B);
    A1 = reshape(A,[1 I 1 J]);
    B1 = reshape(B,[K 1 L 1]);
    AB =  reshape(bsxfun(@times,A1,B1),[I*K J*L]);
end