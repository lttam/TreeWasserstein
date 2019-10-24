function [TM, XX_VertexID] = BuildTreeMetric_HighDim_V2(XX, L, KC)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Tam Le
% RIKEN AIP
% October 24th, 2019
% tam.le@riken.jp
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%%%%%%%
% GOAL:
% Building tree metric for high-dimensional space
% Partition space by the farthest clustering of Gonzalez.

% Apply for discrete empirical measure \mu = \sum_i a_i \delta_{x_i}
% where weight a_i >=0, sum_i a_i = 1, and \delta_{\cdot} is a Dirac
% function.

% NOTE: heigh-level: from 0 --> L
%%%%%%%%
% INPUT:

% XX: cell of supports for N empirical measures
% Each element XX{ii} is a matrix of supports (N_i x dim) where N_i is the
% number of supports with dim dimensions.

% L: the predefined highest level of tree
% (as a stopping condition)

% KC: number of clusters (for the farthest clustering of Gonzalez)
% (when one partitions space by clustering)

%%%%%%%%
% OUTPUT:

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

% XX_VertexID: corresponding for XX in tree.

%%%%%%%%
% SUMMARY: Empirical measures
% XX: 1xm cell --> each cell: n_i x dim
% WW: 1xn cell --> each cell: n_i x 1 double

MAXNUM_TREE = KC^(L+1);

% maximum number of clusters for each height level
% At level \ell --> maximum number of clusters: KC^(\ell)
KCArray = KC.^(1:L)'; % L x 1

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

% gathering ALL SUPPORTS
allXX = zeros(nSupports, dim);  

for ii = 1:N
    allXX(sIDArray(ii):eIDArray(ii), :) = XX{ii};
end

% allXX (all supports) for clustering
nXX = nSupports;

% Node at parent level (centroids of clusters) --> for finding children nodes
% kcCenterPP: set of current parent nodes 
% ===<PP>: denotes for the Parent nodes (current leaf nodes)
% only the center for the first level (root node)
kcCenterPP = mean(allXX)'; % column vector (dim x 1)
numPP = 1; % number of Parent nodes at height level 0

% id of cluster for each empirical measure
% cluster ID starts from 0
% ===<ZZ: cluster ID for "supports" at parent level>
idZZPP = zeros(1, N); % all zeros (all supports in the same cluster id=0)

TM.nVertices = 0;
% Initialization (memory allocation)
TM.Vertex_ParentId = zeros(MAXNUM_TREE, 1);
TM.Vertex_ChildId = cell(MAXNUM_TREE, 1);
TM.Vertex_Pos = zeros(MAXNUM_TREE, dim);
TM.Vertex_EdgeIdPath = cell(MAXNUM_TREE, 1);
%--
TM.Edge_LowNode = zeros(MAXNUM_TREE, 1);
TM.Edge_HighNode = zeros(MAXNUM_TREE, 1);
TM.Edge_Weight = zeros(MAXNUM_TREE, 1);
%--
TM.Level_sID = zeros(L + 1, 1); % from level 0 --> level L
TM.Level_eID = zeros(L + 1, 1); % from level 0 --> level L

% add root into tree
TM.nVertices = 1;
% no parent --> id_parent = 0
TM.Vertex_Pos(1, :) = kcCenterPP; %column vector
TM.Vertex_EdgeIdPath{1} = []; %path is empty
% child nodes: ---> is processing inside the loop

TM.Level_sID(1) = 1;
TM.Level_eID(1) = 1;

% for each height level on the tree T
% (idLL ==> height level idLL + 1
for idLL = 1:L
    
    % HIERARCHICAL-STRUCTURE
    % initialize cluster_Id for child nodes
    % % ===<LL>: denotes for the children nodes (will be obtained by clustering)
    idZZLL = zeros(1, nXX); % cluster id for each support
    % all centroids at height level idLL in kcCenterLL (children nodes)
    % KCArray(idLL): number of clusters at level idLL
    kcCenterLL = zeros(dim, KCArray(idLL)); % having at most KCArray(idLL)
    nkcCenterLL = 0; %count the current real number of centroids at height level (idLL+1)
        
    % weight for each edge --- child node (kcCenterLL) & its corresponding
    % parent node (kcCenterPP)
    % wLL = zeros(1, KCArray(idLL)); %row vector
    
    % save/init info to TM
    TM.Level_sID(idLL+1) = TM.Level_eID(idLL) + 1;%assume >= 1
    TM.Level_eID(idLL+1) = TM.Level_eID(idLL); %will need to update for each idCCPP

    % consider for each cluster at parent level (ID cluster: 0 --> numPP-1)
    % then perform clustering to find children nodes for each parent node
    % idCCPP: ID cluster at parent nodes
    % ==<CC: cluster ID at children nodes>
    for idCCPP = 0:(numPP-1)
        
        % index of parent vertex
        idVertexPP = TM.Level_sID(idLL) + idCCPP;
        
        % First time for clustering or not?
        % idZZ_idCCPP: cluster ID --idCCPP-- for "supports" (idZZ) at parent level.
        if idLL == 1
            % Don't need to find ID for this cluster
            % only 1 cluster at height level 0
            idZZ_idCCPP = 1:(nXX); 
        else
            %extract 1 cluster from parent's clusters (Cluster ID: idCCPP)
            idZZ_idCCPP = find(idZZPP == idCCPP); 
        end
        
        % we have collected support data points on cluster ID --idCCPP--
        % --> perform clustering to obtain children nodes
        % If there is no support data points --> skip
% %         if ~isempty(idZZ_idCCPP)
        % ALSO SKIP if there is only 1 data point (using for small finite
        % number of supports)
        if length(idZZ_idCCPP) > 1
            
            % extract supports for the considered cluster ID from parent's clusters
            allZZ_idCCPP = allXX(idZZ_idCCPP, :); % supports_clusterID_idCCPP x dim
            
            % % function [K, rx, clusterIndex, clusterCenter, numPoints, clusterRadii] = 
            % %                    ...figtreeKCenterClustering(d, N, x, K)
            % % ===(THIRD-PARTY TOOLBOX: Gonzalez's farthest-point clustering algorithm)
            % % Input
            % %    * d --> data dimensionality.
            % %    * N --> number of source points.
            % %    * x --> d x N matrix of N source points in d dimensions 
            % %    * kMax --> maximum number of clusters.
            % % Ouput
            % %    * K --> actual number of clusters (less than kMax if duplicate pts exist)
            % %    * rx --> maximum radius of the clusters.
            % %    * clusterIndex --> vector of length N where the i th element is the 
            % %                cluster number to which the i th point belongs. 
            % %                ClusterIndex[i] varies between 0 to K-1.
            % %    * clusterCenters --> d x K matrix of K cluster centers 
            % %    * numPoints --> number of points in each cluster.
            % %    * clusterRadii --> radius of each cluster.
            
            % perform clustering to find children nodes (children's clusters)
            % rKCLL_idCCPP: ==<rKC(real #clusters)>
            % idZZLL_idCCPP: idZZ("supports")LL("child level") from "idCCPP" clusterID at parent level
            % kcCenterLL_idCCPP: kcCenter("centroids") ...
            [rKCLL_idCCPP, ~, idZZLL_idCCPP, kcCenterLL_idCCPP, ~, ~] = ...
                figtreeKCenterClustering(dim, length(idZZ_idCCPP), allZZ_idCCPP', KC);
            
            
            % Compute weight --> update to the global weight (wLL)!
            % at parent level, kcCenterPP(:, ii+1) is centroid (node) for cluster ID: ii 
            % rKCLL_idCCPP: number of children nodes
            ppMM = repmat(kcCenterPP(:, idCCPP + 1), 1, rKCLL_idCCPP);
            % Using Euclidean distance for edges
            wLL_idCCPP = sqrt(sum((kcCenterLL_idCCPP - ppMM).^2, 1));
            % wLL(nkcCenterLL + (1:rKCLL_idCCPP)) = wLL_idCCPP;
            
            
            setID_0 = find(wLL_idCCPP == 0);
            if length(setID_0) > 0 % exist 0-length

                kcCenterLL_idCCPP(:, setID_0) = []; % delete 0-length-edge clusters
                wLL_idCCPP(setID_0) = [];
                
                % RELABEL CLUSTER-ID
                % -1: 0-length-edge cluster
                % 0 --> rKCLL_idCCPP - 1
                % 0 --> rKCLL_idCCPP - length(setID_0) - 1
                               
                % set -1
                clusterID_ZeroLength = setID_0 - 1;
                allID_Zero = [];
                for iiCC_Zero = 1:length(clusterID_ZeroLength)
                    tmp = find(idZZLL_idCCPP == clusterID_ZeroLength(iiCC_Zero));
                    allID_Zero = [allID_Zero; tmp'];
                end
                idZZLL_idCCPP(allID_Zero) = -1;

                % RELABEL
                clusterID_NonZero = 0:(rKCLL_idCCPP-1);
                clusterID_NonZero(setID_0) = [];
                for iiCC_NonZero = 1:length(clusterID_NonZero)
                    if clusterID_NonZero(iiCC_NonZero) ~= (iiCC_NonZero - 1)
                        idZZLL_idCCPP(idZZLL_idCCPP == clusterID_NonZero(iiCC_NonZero)) = iiCC_NonZero - 1;
                    end
                end
                
                rKCLL_idCCPP = rKCLL_idCCPP - length(setID_0); % reduce #clusters
            end
            
            
            % update to global index (idZZLL): id of supports at child level
            % each cluster ID: idCCPP ---> have at most KC clusters
            
            idZZLL(idZZ_idCCPP) = nkcCenterLL + idZZLL_idCCPP;
            if length(setID_0) > 0 % exist 0-length
                idZZLL(idZZ_idCCPP(allID_Zero)) = -1;
            end

            
            if rKCLL_idCCPP > 0
                
                kcCenterLL(:, nkcCenterLL+(1:rKCLL_idCCPP)) = kcCenterLL_idCCPP;
            
                % update current real number of centroids (children nodes)
                nkcCenterLL = nkcCenterLL + rKCLL_idCCPP;

                % save to tree metric TM
                % rKCLL_idCCPP (number of new children nodes)
                % ID for new nodes (vertices): TM.Level_eID(idLL+1 +
                % [1:rKCLL_idCCPP]
                % kcCenterLL_idCCPP: centroids (node pos)

                % == Vertices
                TM.nVertices = TM.nVertices + rKCLL_idCCPP;
                idNewVertices = TM.Level_eID(idLL+1) + (1:rKCLL_idCCPP);
                TM.Vertex_ParentId(idNewVertices) = idVertexPP;
                TM.Vertex_ChildId{idVertexPP} = idNewVertices; % OK for parent node, later for children nodes
                TM.Vertex_Pos(idNewVertices, :) = kcCenterLL_idCCPP';
                % == Edges
                idNewEdges = idNewVertices - 1;
                TM.Edge_LowNode(idNewEdges) = idVertexPP;
                TM.Edge_HighNode(idNewEdges) = idNewVertices;
                TM.Edge_Weight(idNewEdges) = wLL_idCCPP;
                % path --vertex
                for ii = 1:rKCLL_idCCPP
                    % Path (Edge ID) from root to each node
                    TM.Vertex_EdgeIdPath{idNewVertices(ii)} = [TM.Vertex_EdgeIdPath{idVertexPP}, idNewEdges(ii)];
                end
                % == Level
                TM.Level_eID(idLL+1) = TM.Level_eID(idLL+1) + rKCLL_idCCPP;
            end
        end
    end
    
    % UPDATE clusterPP at parent level by information from children level
    % --> using for the next height level
    idZZPP = idZZLL; %cluster ID for each point! (values starting from 0)
    
    %update kcCenterPP (centroids) -- real centroids
    kcCenterPP = kcCenterLL(:, 1:nkcCenterLL);
    numPP = nkcCenterLL;

end

TM.LeavesIDArray = (TM.Level_sID(L+1):TM.Level_eID(L+1));

% TM.nVertices
% Initialization (memory allocation)
TM.Vertex_ParentId = TM.Vertex_ParentId(1:TM.nVertices); %zeros(MAXNUM_TREE, 1);
TM.Vertex_ChildId = TM.Vertex_ChildId(1:TM.nVertices); %cell(MAXNUM_TREE, 1);
TM.Vertex_Pos = TM.Vertex_Pos(1:TM.nVertices, :); %zeros(MAXNUM_TREE, dim);
TM.Vertex_EdgeIdPath = TM.Vertex_EdgeIdPath(1:TM.nVertices); %cell(MAXNUM_TREE, 1);
%--
TM.Edge_LowNode = TM.Edge_LowNode(1:(TM.nVertices - 1)); %zeros(MAXNUM_TREE, 1);
TM.Edge_HighNode = TM.Edge_HighNode(1:(TM.nVertices - 1)); %zeros(MAXNUM_TREE, 1);
TM.Edge_Weight = TM.Edge_Weight(1:(TM.nVertices - 1)); %zeros(MAXNUM_TREE, 1);

if nargout > 1
    XX_VertexID = cell(1, N);
    for ii = 1:N
        % for each empirical measure
        % sIDArray(ii) % starting index
        % eIDArray(ii) % ending index
        % idZZPP (starting from 0 ---> nkcCenterLL-1)
        % --> corresponding to TM.Level_sID(L+1) --> TM.Level_eID(L+1)
        % XX{ii}
        
        %corresponding set of idVertex
        XX_VertexID{ii} = TM.Level_sID(L+1) + idZZPP(sIDArray(ii):eIDArray(ii));
    end
end

end


















