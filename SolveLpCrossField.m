% Inputs:
% X is the vertex list
% T is the triangle list
% mname is the mesh name. leave as '' if you don't know or don't care. This will be used to save the results into the Results folder unless it's ''.
% n is the softness on normal alignment. n=0 is hard alignment. n~.8 is no alignment
% p is the norm of the cross field. must be between 1 and infinity inclusive.
% Visualize is a toggle for whether or not to display the final cross field.
% isFixedTriangle is a boolean vector that is [size(T,1) x 1] and is true if you want to constrain the cross field on that triangle
% fixedTriangleFrames should be [3 x 3 x sum(isFixedTriangle)] and should be a rotation matrix (up to negation or orientation preserving permutation of columns). The column vectors make up the directions of the frame so one of them is recommended to (but isn't required to) be the normal of the triangle.
% Outputs:
% dirs1 is a [3x3xsize(T,1)] tensor where each 3x3 is a rotation matrix whose columns are vectors of the cross field or normals of the triangle.
% fname is the name of where the cross field is saved to file.
% data is the processed triangle mesh data.
function [dirs1, fname, data] = SolveLpCrossField(X, T, mname, n, p, Visualize, isFixedTriangle, fixedTriangleFrames)

    % verifications on pointwise dirichlet conditions
    if numel(isFixedTriangle)==0
        clear isFixedTriangle
        clear fixedTriangleFrames
    end
    if exist('isFixedTriangle','var')
        if~(exist('fixedTriangleFrames','var'))
            error('cant fix triangle frames without specifying the frame rotation values.');
        end
        assert(sum(isFixedTriangle) == size(fixedTriangleFrames,3));
        assert(size(isFixedTriangle,2)==1)
    end

    assert(p>=1);
    
    mosekError = false;
    shouldsave = numel(mname)~=0;
    if shouldsave
        if exist('isFixedTriangle','var')
            fname = sprintf('Results/mesh_%s_n_%g_p_%g_dc_%d.mat',mname,n,p,randi(1000));
        else
            fname = sprintf('Results/mesh_%s_n_%g_p_%g.mat',mname,n,p);
        end
    else
        fname = '';
    end
    
    if (exist(fname,'file')~=0 && nargin ~= 0)
        display('Found same result already computed in Results folder. Skipping.');
        fname = 'Skipped';
        dirs1 = [];
        data = [];
        return;
    end
    
    fiber = OctaMBO;
    data = getMeshData(X,T);
    
    % set up sph boundary
    vc = sqrt(5/12)*[0 0 0 0 0 0 0 0 1];
    vs = sqrt(5/12)*[1 0 0 0 0 0 0 0 0];
    vn = sqrt(7/12)*[0 0 0 0 1 0 0 0 0];
    a0 = vn + sin(pi/4)*vs + cos(pi/4)*vc;
    thresh = .665;
    
    SA = sum(data.triangleAreas);
    u0 = repmat(sparse([0;0;0;sqrt(7/12);0;0;0]),data.numTriangles,1);
    scaledu0 = u0; scaledu0(4:7:end)=sqrt(data.triangleAreas/SA)*sqrt(7/12);
    
    %% build wigner rots.
    D = OctaAlignMat(data.faceNormals);
    D(abs(D)<=1e-6)=0;
    N9 = rank3tensor2blockdiag(D);
    N2 = rank3tensor2blockdiag(D([1 9],:,:));
    N7 = rank3tensor2blockdiag(D(2:8,:,:));
    scaledN7 = rank3tensor2blockdiag(D(2:8,:,:).*permute(sqrt(data.triangleAreas/SA),[2 3 1]));
    
    prepocW = (1./data.edgeWeights);
    
    badFrameInds = (1:data.numTriangles)'; goodFrameInds=[]; 
    framesProj = zeros(data.numTriangles,9);
    if exist('isFixedTriangle','var')
        % account for dirichlet condition frames.
        badFrameInds = find(~isFixedTriangle);
        goodFrameInds = find(isFixedTriangle);
        framesProj(goodFrameInds,:) = Frames2Octa(fixedTriangleFrames)';
    end
    
    psymbol = p; if p==inf; psymbol = -1; end
    counter = 0; goodFramesPerMetaIter = 0;
    DM = kron(data.incidenceMatrix,speye(9));
    while numel(badFrameInds)>0 && counter < 10 && mosekError==false
        t1 = tic;
        counter = counter + 1;
        
        % fix the good/bad frames so they aren't variables anymore
        if numel(goodFrameInds~=0)
            goodFrameInds9Flat = unique(repmat((goodFrameInds-1)*9,1,9)+(1:9)); goodFrameInds9Flat = goodFrameInds9Flat(:);
        else; goodFrameInds9Flat = []; 
        end
        badFrameInds9Flat = unique(repmat((badFrameInds-1)*9,1,9)+(1:9)); badFrameInds9Flat = badFrameInds9Flat(:);

        
        AM = rank3tensor2blockdiag(D(2:8,:,badFrameInds));
        framesProjt = framesProj';
        cm = DM(:,goodFrameInds9Flat)*framesProjt(goodFrameInds9Flat); if numel(cm)==0; cm = zeros(data.numInteriorEdges*9,1); end;
        try 
            [xbad] = MosekSoftCrossFieldsWrapper(AM, u0(1:7*numel(badFrameInds)), DM(:,badFrameInds9Flat), prepocW, psymbol, n, cm);
        catch exception
            mosekError = true;
            if contains(exception.message,'Unknown') && exist('x','var')
                xbad = x(badFrameInds9Flat,:);
            else
                throw(exception);
            end
        end
        x(badFrameInds9Flat,:) = xbad;
        secondsElapsed=toc(t1);

        % store the result.
        framesUnprojRaw = full(reshape(x,9,[]))'; 
        framesUnproj = framesUnprojRaw./vecnorm(framesUnprojRaw,2,2);

        % Extremely rare frame projection failure. Random perturbation usually fixes it.
        randpert = zeros(9,numel(badFrameInds));
        while true
            try
                framesProj(badFrameInds,:) = fiber.proj(framesUnproj(badFrameInds,:)' + randpert)';
                break;
            catch 
                randpert = (rand(9,numel(badFrameInds))-.5)*1e-12;
                fprintf('Arff frame projection problem. Re-trying.\n');
            end
        end

        deviation = vecnorm(framesUnproj-framesProj,2,2); % higher deviation from oct variety means less good frame quality. Ambivalent to norm of frames.
        % threshold between good/bad frames
        goodFrameInds = find(deviation < thresh);
        badFrameInds = find(deviation >= thresh);

        if counter > 1 && numel(badFrameInds) == badFramesPerMetaIter(counter-1)
            break; % Guarantee termination.
        end

        dirs1 = Coeff2Frames(fiber.proj(framesProj')); 
        dirs2 = reshape([dirs1 -dirs1],3,[])';
        alignmentColors = vecnorm(full(reshape(N7*x - u0,7,[]))',2,2);

        goodFramesPerMetaIter(counter) = numel(goodFrameInds);
        badFramesPerMetaIter(counter) = numel(badFrameInds);
        frames{counter}.secondsElapsed = secondsElapsed;
        frames{counter}.frames = framesProj;
        frames{counter}.framesUnprojRaw = framesUnprojRaw;
        frames{counter}.goodFrameInds = goodFrameInds;
        frames{counter}.badFrameInds = badFrameInds;
        frames{counter}.alignmentColors = alignmentColors;
        if exist('cvx_status','var'); frames{counter}.cvx_status = cvx_status; end;
        if exist('details','var'); frames{counter}.details = details; end;
        frames{counter}.dirs1 = dirs1;
        frames{counter}.dirs2 = dirs2;
    end
    fprintf('Cross field completed. \n');
    lines = reshape(repmat(data.triangleBarycenters,1,6)', 3,[])';
    
    % visualize final cross field
    if Visualize
        dirs1 = Coeff2Frames(framesProj');
        dirs2 = reshape([dirs1 -dirs1],3,[])';
        figure; hold all; axis equal; rotate3d on;
        patch('Faces',data.triangles,'Vertices',data.vertices,'FaceColor','green');
        quiver3(lines(:,1),lines(:,2),lines(:,3),dirs2(:,1),dirs2(:,2),dirs2(:,3),'b');
    end
    
    if shouldsave
        toSave.frames = frames; % frame fields per iteration and which are good/bad
        toSave.mosekError = mosekError;
        toSave.mname = mname;
        toSave.metaIterationCounter = counter;
        toSave.goodFramesPerMetaIter = goodFramesPerMetaIter; % number of good/bad frames per meta iteration
        toSave.badFramesPerMetaIter = badFramesPerMetaIter;
        toSave.X = data.vertices;
        toSave.T = data.triangles;
        toSave.data = data;
        toSave.lines = lines;
        toSave.n = n;
        toSave.p = p;
        save(fname,'toSave');
    end
end







