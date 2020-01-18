%% build: This requires TBB, Mosek, gptoolbox, arff

% checks that we're in the right folder
currdir = pwd; 
assert(strcmp(currdir(end-23:end-1),'CreaseAlignedCrossField'))

% builds solver code
addpath('MosekSoftCrossFields');
cd('MosekSoftCrossFields');
mexbuild;
cd(currdir);

% check that gptoolbox is installed: https://github.com/alecjacobson/gptoolbox
% we only need it for the readOBJ function. 
assert(contains(path,'gptoolbox\mesh') || contains(path,'gptoolbox/mesh'));

% check that arff is installed: https://github.com/dpa1mer/arff
% this is included as a submodule already. 
% Follow the build instructions for arff.
assert(contains(path,'arff'));

%% run example
% fname = 'Meshes/Cyl3Inter_denser.obj';
fname = 'Meshes/notch5.obj';
% fname = 'Meshes/twistcube90.obj';

[~,mname,~] = fileparts(fname);
normalAlignment = 0;
pnorm = 2;
ShouldVisualize = true;
[X,T] = readOBJ(fname);
isFixedTriangle = [];
fixedTriangleFrames = [];
[dirs1, fname, data] = SolveLpCrossField(X, T, mname, normalAlignment, pnorm, ShouldVisualize, isFixedTriangle, fixedTriangleFrames);











