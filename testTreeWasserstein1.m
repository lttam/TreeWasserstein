% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Tam Le
% RIKEN AIP
% October 24th, 2019
% tam.le@riken.jp
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clear all
clc

% compute tree metric from input data
% then compute tree-Wasserstein (TW) distance matrix for the same input data

load('Subset_1000.mat');

% parameter of tree metric
L = 5; % deepest level
KC = 4; % number of clusters for the farthest-point clustering

% building tree metric by the farthest-point clustering
disp('...Computing the tree metric from input data');
tic
[TM, TX] = BuildTreeMetric_HighDim_V2(XX, L, KC);
runTime = toc;
disp(['......running time: ' num2str(runTime)]);

disp('...Computing tree representation for input data');
tic
% mapping vector on tree
XX_TM = zeros(length(XX), length(TM.Edge_Weight));
for ii = 1:length(XX)
    % set of vertices of tree -- corresponding for probability XX
    XX_idVV = TX{ii};
    WW_idVV = WW{ii};
    for jj = 1:length(XX_idVV)
        XX_TM(ii, TM.Vertex_EdgeIdPath{XX_idVV(jj)}) = XX_TM(ii, TM.Vertex_EdgeIdPath{XX_idVV(jj)}) + WW_idVV(jj);
    end
end
% weighted mapping
XX_TMWW = XX_TM .* repmat(TM.Edge_Weight', length(XX), 1);
runTime = toc;
disp(['......running time: ' num2str(runTime)]);


disp('...Computing l1-distance for tree representation data');
tic
% compute TW distance matrix for XX
% L1 distance
DD_XX = zeros(length(XX), length(XX));
for ii = 1:(length(XX)-1)
    % L1 distances between ii^th id and (ii+1 : length(XX))^th ids
    tmp = sum(abs(repmat(XX_TMWW(ii, :), length(XX) - ii, 1) - XX_TMWW((ii+1):length(XX), :)), 2);
    DD_XX(ii, (ii+1):length(XX)) = tmp';
    DD_XX((ii+1):length(XX), ii) = tmp;
end
runTime = toc;
disp(['......running time: ' num2str(runTime)]);


disp('FINISH!');

