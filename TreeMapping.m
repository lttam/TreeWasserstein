function TX = TreeMapping(XX, WW, TM)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Tam Le
% RIKEN AIP
% October 24th, 2019
% tam.le@riken.jp
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% COMPUTE tree mapping vector for empirical measure (XX, WW) in tree metric TM
% TREE WASSERSTEIN among empirical measures (XX, WW) is equivalent to L1
% DISTANCE among tree mapping vector (TX)
% When we have TreeMapping --> Tree Wasserstein distance is a built-in function

% Main ideas:
% For each given support, we can travel over the tree to get the path of edges
% (had for parallel)
% --- PARALLELIZED ---
% Find the neareset neighbor from given supports to leaf-nodes in tree metric TM
% (by nearest neighbors with the option of using kd-tree)

%%%%%%%%
% INPUT:

% XX: cell of supports for N empirical measures
% Each element XX{ii} is a matrix of supports (N_i x dim) where N_i is the
% number of supports with dim dimensions.
%
% WW: corresponding weights for supports in XX
%
% TM: Tree metric (from BuildTreeMetric_HighDim)
% TM: Tree structure T
% Structure of TM
% ----TM.Vertex: set of vertices (4 components: parent, child, pos(ition), path)
% TM.nVertices : number of vertices
% TM.Vertex_ParentId
% TM.Vertex_ChildId
% TM.Vertex_Pos
% TM.Vertex_EdgeIdPath
%
% ----TM.Edge: set of Edges (3 components: lowNode, highNode, weight)
% TM.Edge_LowNode
% TM.Edge_HighNode
% TM.Edge_Weight
%
% TM.Level_sID: start ID (nodes) at each height level
% TM.Level_eID: ending ID (nodes) at each height level
%
% TM.LeavesIdArray: array of leaf IDs

%%%%%%%%
% OUTPUT:

% TX: tree mapping vector in R_m (m: number of edges in tree)

%%%%%%%%
% SUMMARY: Empirical measures
% XX: 1xm cell --> each cell: n_i x dim
% WW: 1xn cell --> each cell: n_i x 1 double

% number of empirical measures
N = length(XX);

% dimension of supports
dim = size(XX{1}, 2);

% index for each empirical measures (when gathering all supports together
% for building the tree)
nSupports = 0; %count the total number of supports in XX
sIDArray = zeros(N, 1); % starting index
eIDArray = zeros(N, 1); % ending index

% get start-end ID of each empirical measure
for ii = 1:N
    sIDArray(ii) = nSupports + 1; %starting index
    nSupports = nSupports + size(XX{ii}, 1);
    eIDArray(ii) = nSupports; %ending index    
end

% gathering ALL SUPPORTS & WEIGHTS
allXX = zeros(nSupports, dim);  
allWW = zeros(nSupports, 1);

for ii = 1:N
    allXX(sIDArray(ii):eIDArray(ii), :) = XX{ii};
    allWW(sIDArray(ii):eIDArray(ii)) = WW{ii};
end

% using MATLAB toolbox (adaptive between kdtree & exhaustive search)
% Euclidean distance metric
% Idx = knnsearch(X,Y) finds the nearest neighbor in X for each query point in Y 
% and returns the indices of the nearest neighbors in Idx, a column vector. 
% Idx has the same number of rows as Y.
% Input: N x dim
% Output: N x 1 (column vector)

% id of all leaves: TM.LeavesIDArray
allLeaves = TM.Vertex_Pos(TM.LeavesIDArray, :);
idLeaves = knnsearch(allLeaves, allXX);
allIdVertices = TM.LeavesIDArray(idLeaves);

TX = zeros(N, TM.nVertices - 1); % edge representation vector
for ii = 1:N
    tmpVector = zeros(1, TM.nVertices - 1);
    for jj = sIDArray(ii):eIDArray(ii)
        idEdges = TM.Vertex_EdgeIdPath{allIdVertices(jj)};
        tmpVector(idEdges) = tmpVector(idEdges) + allWW(jj); %empirical measures
    end
    TX(ii, :) = tmpVector .* TM.Edge_Weight'; %row vector
end

end




        
        
        
        
        
        
        
        
        
        
        
        
        
       