function [edgePot,edgeStruct]=CreateGridUGMModel(NumFils, NumCols, K, lambda)
%
%
% NumFils, NumCols: image dimension
% K: number of states
% lambda: smoothing factor
tic
% DEFINE MODEL STRUCTURE
nNodes = NumFils*NumCols;
adj = sparse(nNodes,nNodes);

% Add Down Edges
ind = 1:nNodes;
exclude = sub2ind([NumFils NumCols],repmat(NumFils,[1 NumCols]),1:NumCols); % No Down edge for last row
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+1)) = 1;

% Add Right Edges
ind = 1:nNodes;
exclude = sub2ind([NumFils NumCols],1:NumFils,repmat(NumCols,[1 NumFils])); % No right edge for last column
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+NumFils)) = 1;

% Add Up/Left Edges
adj = adj+adj';
edgeStruct = UGM_makeEdgeStruct(adj,K);

% DEFINE THE POTENTIALS
edgePot = zeros(K,K,edgeStruct.nEdges);

%Potts model
aux = ones(K)*exp(-lambda(2));
aux(logical(eye(K))) = exp(-lambda(1));
    
   
for e = 1:edgeStruct.nEdges
    edgePot(:,:,e) = aux;    
end

toc

end