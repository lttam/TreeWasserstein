% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Tam Le
% RIKEN AIP
% October 24th, 2019
% tam.le@riken.jp
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clear all
clc

% compute tree metric from input data (1)
% then compute tree-Wasserstein (TW) distance matrix for a different input
% data (2)

% input data (1) --> for building tree metric
load('Subset_200.mat');

% parameter of tree metric
L = 5; % deepest level
KC = 4; % number of clusters for the farthest-point clustering

% building tree metric by the farthest-point clustering
disp('...Computing the tree metric from input data (1)');
tic
[TM, TX] = BuildTreeMetric_HighDim_V2(XX, L, KC);
runTime = toc;
disp(['......running time: ' num2str(runTime)]);

% input data (2) --> for computing tree-Wasserstein distance matrix
% using the tree metric built from input data (1)
load('Subset_1000.mat');

disp('...Computing tree representation for input data (2)');
tic
XX_TMWW = TreeMapping(XX, WW, TM);
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

