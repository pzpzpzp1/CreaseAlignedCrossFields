% input As:n x n2 x m 
% output sparse bulkMat: n m x n2 m
function bulkMat = rank3tensor2blockdiag(As)

n = size(As,1);
n2 = size(As,2);
m = size(As,3);

A_n2m_n = reshape(permute(As,[1 3 2]),n2*m,n);
II = repmat([1:n*m]',n2,1);
Jt = reshape(repmat((0:(m-1))*n2,n,1),1,[]);
JJ = reshape((repmat(Jt,n2,1)+[1:n2]')',[],1)';
KK = A_n2m_n(:);
bulkMat = sparse(II,JJ,KK,n*m,n2*m);

end