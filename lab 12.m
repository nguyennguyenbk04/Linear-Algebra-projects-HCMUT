% Google PageRank algorithm on the example of random network

%% Load the network data
load("AdjMatrix.mat");
%% Check  ratio of non-zero elements
nnzAdjMatrix = (nnz(AdjMatrix) / numel(AdjMatrix)) * 100;
%% Dimensions of the matrix AdjMatrix
[m,n] = size(AdjMatrix);
%% Display a small amount of network
NumNetwork=500;
AdjMatrixSmall=AdjMatrix(1:NumNetwork,1:NumNetwork);

for j=1:NumNetwork
    coordinates(j,1)=NumNetwork*rand;
    coordinates(j,2)=NumNetwork*rand;
end;

gplot(AdjMatrixSmall,coordinates,'k-*');

%% Check the amount of links originating from each webpage
NumLinks=sum(AdjMatrixSmall,2);

%% Create a matrix of probabilities (Google matrix)
% Element (i,j) of the matrix shows the probability of moving from i-th
% page of the network to jth page. It is assumed that the user can follow
% any link on the page with a total probability of 85% (all hyperlinks are 
% equal), and jump (teleport) to any other page in the network with a total
% probability of 15% (again, all pages are equal).

alpha=0.15; 
GoogleMatrix=zeros(NumNetwork,NumNetwork);
for i=1:NumNetwork
    if NumLinks(i)~=0
       GoogleMatrix(i,:)=AdjMatrixSmall(i,:)./NumLinks(i);
    else
       GoogleMatrix(i,:)=1./NumNetwork;
    end;
end;

GoogleMatrix=(1-alpha)*GoogleMatrix+alpha*ones(NumNetwork,NumNetwork)./NumNetwork;

%% Check that all the rows in the GoogleMatrix matrix sum to 1
SumGoogleMatrix=sum(GoogleMatrix,2);

%% Finding an eigenvector corresponding to 1 (why is there sucj an eigenvector)?
w0=ones(1,NumNetwork)./NumNetwork;

w1=w0*GoogleMatrix;
w2=w1*GoogleMatrix;
w3=w2*GoogleMatrix;
w5=w0*(GoogleMatrix^5);
w10=w0*(GoogleMatrix^10);

% Check the difference between v30 and v20. Observe that it is pretty
% small
deltaw=w10-w5;
%% Compute the eigenvalues and the right eigenvectors 
[right_eigenvector, right_eigenvalue] = eig(GoogleMatrix);
% Explain the result
right_eigenvalue = diag(right_eigenvalue);
%%find eigenvalue 1
tolerance = 1e-10;
idx = find(abs(real(right_eigenvalue) - 1) < tolerance & abs(imag(right_eigenvalue)) < tolerance, 1);

% Check if an eigenvalue 1 was found
if isempty(idx)
    disp('No eigenvalue 1 found.');
else
   right_eigenvector = right_eigenvector(:, idx);
end
%% Compute the eigenvalues and the left eigenvectors
[left_eigenvector, left_eigenvalue] = eig(GoogleMatrix');
% Explain the result
left_eigenvalue = diag(left_eigenvalue);
%%find eigenvalue 1
tolerance = 1e-10;
idx = find(abs(real(left_eigenvalue) - 1) < tolerance & abs(imag(left_eigenvalue)) < tolerance, 1);

% Check if an eigenvalue 1 was found
if isempty(idx)
    disp('No eigenvalue 1 found.');
else
   left_eigenvector = left_eigenvector(:, idx);
end
%% Separate the eigenvector corresponding to the eigenvalue 1 and scale it
u1 = left_eigenvector;
u1=abs(u1)/norm(u1,1);
%% Select the maximum element and the corresponding element. 
%Which page is the most important in the network?
[MaxRank,PageMaxRank]=max(u1);
%% Check if it's the most popular (most linked to page):
MostLinks=sum(GoogleMatrix, 1);
[MaxLinks, PageMaxLinks]=max(MostLinks);



