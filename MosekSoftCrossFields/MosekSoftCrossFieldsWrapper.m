% wrapping layer to cover the sparse matrix to triplet transformation
function x = MosekSoftCrossFieldsWrapper(A, b, D, a, psymbol, n, c)
    [Ai,Aj,As] = find(sparse(A));
    [Di,Dj,Ds] = find(sparse(D));
    A0 = size(A,1);
    A1 = size(A,2);
    D0 = size(D,1);
    D1 = size(D,2);

    x = MosekSoftCrossFields(Ai-1, Aj-1, As, A0, A1, full(b), Di-1, Dj-1, Ds, D0, D1, full(a), psymbol, n, full(c));
end
