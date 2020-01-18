% triangle meshes
function data = getMeshData(X,T)

data = [];
data.vertices = X;
data.triangles = double(T);
data.numVertices = size(X,1);
data.numTriangles = size(T,1);

normalf = cross( data.vertices(data.triangles(:,2),:)'-data.vertices(data.triangles(:,1),:)', ...
                 data.vertices(data.triangles(:,3),:)'-data.vertices(data.triangles(:,1),:)' );
d = sqrt( sum(normalf.^2,1) ); d(d<eps)=1;
data.faceNormals = (normalf ./ repmat( d, 3,1 ))';

%% extract surface mesh data.
rawedges = data.triangles(:,[1 2 2 3 3 1])';
edges = reshape(rawedges,2,[])';
data.edges = unique(sort(edges,2),'rows');
data.numEdges = size(data.edges,1);
T2V = sparse(repmat(1:data.numTriangles,3,1),data.triangles',ones(numel(data.triangles),1),data.numTriangles,data.numVertices);
V2E = sparse(repmat(1:data.numEdges,2,1),data.edges',ones(numel(data.edges),1),data.numEdges,data.numVertices)';
T2E = T2V*V2E; T2E = T2E == 2;
if any(sum(T2E)>2)
    nonManifoldEdges = find(sum(T2E)>2);
    error('Non Manifold Triangles Found!!!');
end
isBoundaryEdge = sum(T2E)==1;
[~, ia] = sort(isBoundaryEdge);
data.edges = data.edges(ia,:);
T2V = sparse(repmat(1:data.numTriangles,3,1),data.triangles',ones(numel(data.triangles),1),data.numTriangles,data.numVertices);
V2E = sparse(repmat(1:data.numEdges,2,1),data.edges',ones(numel(data.edges),1),data.numEdges,data.numVertices)';
T2E = T2V*V2E; T2E = T2E == 2;
isBoundaryEdge = sum(T2E)==1;
data.numInteriorEdges = sum(~isBoundaryEdge);
data.numBoundaryEdges = sum(isBoundaryEdge);
[ii1, jj1] = find(T2E');
[ii2, jj2] = find(T2E(:,1:data.numInteriorEdges));
data.triangles2edges = sort(reshape(ii1,3,[])',2);
data.edges2triangles = sort(reshape(ii2,2,[])',2); % data.numInteriorEdges x 2
data.isBoundaryEdge = isBoundaryEdge;

n=data.numTriangles;
adj = sparse(repmat((1:n)',3,1),data.triangles2edges(:),ones(3*n,1)); % tris x edges
[ii,jj] = find(adj(:,(data.numInteriorEdges+1):end));
data.edges2triangles = [data.edges2triangles; ii ii*0];

data.triangleBarycenters = (data.vertices(data.triangles(:,1),:)+data.vertices(data.triangles(:,2),:)+data.vertices(data.triangles(:,3),:))/3;

x1 = data.vertices(data.triangles(:,1),1);
y1 = data.vertices(data.triangles(:,1),2);
z1 = data.vertices(data.triangles(:,1),3);
x2 = data.vertices(data.triangles(:,2),1);
y2 = data.vertices(data.triangles(:,2),2);
z2 = data.vertices(data.triangles(:,2),3);
x3 = data.vertices(data.triangles(:,3),1);
y3 = data.vertices(data.triangles(:,3),2);
z3 = data.vertices(data.triangles(:,3),3);

A = sqrt((x1-x2).^2+(y1-y2).^2+(z1-z2).^2);
B = sqrt((x2-x3).^2+(y2-y3).^2+(z2-z3).^2);
C = sqrt((x3-x1).^2+(y3-y1).^2+(z3-z1).^2);
S = (A+B+C)/2;
data.triangleAreas = sqrt(S.*(S-A).*(S-B).*(S-C));
data.triangleEdgeBLength = B; % an arbitrary edge per triangle

data.edgeLengths = sqrt(sum((data.vertices(data.edges(:,1),:) - data.vertices(data.edges(:,2),:)).^2,2));
data.dualEdgeLengths = sqrt(sum((data.triangleBarycenters(data.edges2triangles(1:data.numInteriorEdges,1),:) - data.triangleBarycenters(data.edges2triangles(1:data.numInteriorEdges,2),:)).^2,2));
data.primalOverDualWeight = data.edgeLengths(1:data.numInteriorEdges)./data.dualEdgeLengths;

data.primalIncidence = sparse(repmat(1:data.numEdges,1,2), data.edges(:), [ones(data.numEdges,1); -ones(data.numEdges,1)], data.numEdges, data.numVertices);

% dual incidence. IExT
data.incidenceMatrix = sparse(1:data.numInteriorEdges, data.edges2triangles(1:data.numInteriorEdges,1), ones(data.numInteriorEdges,1), data.numInteriorEdges, data.numTriangles);
data.incidenceMatrix = data.incidenceMatrix - sparse(1:data.numInteriorEdges, data.edges2triangles(1:data.numInteriorEdges,2), ones(data.numInteriorEdges,1), data.numInteriorEdges, data.numTriangles);
data.dualGraphL = data.incidenceMatrix'*data.incidenceMatrix;

data.triangle2verts = T2V; 
data.vertNormals = data.triangle2verts'*(data.faceNormals.*data.triangleAreas);
data.vertNormals = data.vertNormals ./ vecnorm(data.vertNormals,2,2);


data.edgeWeights =  1./data.primalOverDualWeight;
assert(all(data.edgeWeights>=0));


%% compute triXtri2edge
% adj = sparse(repmat((1:n)',3,1),data.triangles2edges(:),ones(3*n,1)); % tris x edges
%adj_edgelabeled = sparse(repmat((1:n)',3,1),data.triangles2edges(:),data.triangles2edges(:)); % tris x edges
adj_trilabeled = sparse(repmat((1:n)',3,1),data.triangles2edges(:),repmat([1:data.numTriangles]',3,1)); % tris x edges
% compute tri x tri 2 tet
edgeXedge2tri = (adj'*adj_trilabeled);
edgeXedge2tri(1:data.numEdges+1:end)=0;
data.edgeXedge2tri = edgeXedge2tri;

% compute edges to triangles. one tri will be 0 for boundary edges
data.triXtri2edge = sparse(data.edges2triangles(1:data.numInteriorEdges,1),data.edges2triangles(1:data.numInteriorEdges,2),1:data.numInteriorEdges,data.numTriangles,data.numTriangles);
data.triXtri2edge = data.triXtri2edge + data.triXtri2edge';
data.triXedge2tri = sparse(data.edges2triangles(1:data.numInteriorEdges,1),1:data.numInteriorEdges,data.edges2triangles(1:data.numInteriorEdges,2),data.numTriangles,data.numEdges) +...
    sparse(data.edges2triangles(1:data.numInteriorEdges,2),1:data.numInteriorEdges,data.edges2triangles(1:data.numInteriorEdges,1),data.numTriangles,data.numEdges);



