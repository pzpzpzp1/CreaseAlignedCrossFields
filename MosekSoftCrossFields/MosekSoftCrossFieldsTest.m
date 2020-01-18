%% solves min_n |n|^p + c(n-n0)^2 in each element independently in parallel.
% verification of MosekSoftCrossFields compared to CVX solution. Requires CVX to run.
function MosekSoftCrossFieldsTest(A,b,D,a,p,n)
    addpath('..');
    maxNumCompThreads(96);

    if nargin == 0
%         nedges = randi(400)+2;
%         nfaces = randi(600)+2;
        nedges = randi(10)+2;
        nfaces = randi(10)+2;
        dimPerEdge = 9;
        D = rand(dimPerEdge * nedges, dimPerEdge * nfaces)-.5;
        A = rank3tensor2blockdiag(rand((dimPerEdge-2), dimPerEdge, nfaces)-.5);
        psymbol = randsample([-1 1 2 10:20:50],1); % -1 means inf
        a = rand(nedges,1);
        b = rand((dimPerEdge-2)*nfaces,1)-.5;
        n = randsample([0:.1:1],1);
        c = rand(nedges*dimPerEdge,1)-.5;
    end
    actualp = psymbol; if psymbol==-1; actualp = inf; end
    ap = a.^(1/actualp);

    %% compute with mex mosek 
    xmosek = MosekSoftCrossFieldsWrapper(sparse(A),b,sparse(D),a,psymbol,n,c);
    mosekE = norm( norms(reshape(D*xmosek+c,9,[])',2,2) .*ap,actualp);
    % norms(reshape(A*xmosek-b,dimPerEdge-2,[]),2,1) <= n
    
    %% compute with cvx
    cvx_begin
        cvx_solver mosek;
        cvx_precision best;
        variable xcvx(9*nfaces,1);
        variable ynorms(nedges,1);
        y = D*xcvx + c;
        pnormE = norm(ynorms.*ap,actualp);
        minimize pnormE
        subject to
            norms(reshape(A*xcvx-b,dimPerEdge-2,[]),2,1)<=n
            ynorms >= norms(reshape(y,9,[])',2,2); %nedges number of qcone constraints
    cvx_end
    cvxE = pnormE;
    
    fprintf('results: (xdiff: %f) (Ediff: %f) (p:%d n:%g) \n', norm(xcvx-xmosek), cvxE-mosekE, actualp, n);
    
end
