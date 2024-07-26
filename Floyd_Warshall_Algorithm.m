% Example Graph
D = [
    0, 4, inf, 5, inf;
    inf, 0, 1, inf, 6;
    2, inf, 0, 3, inf;
    inf, inf, 1, 0, 2;
    1, inf, inf, 4, 0;
];
[d_result]= floydWarshall(D);
disp("The adjacency matrix that represent the shortest path between all the pairs of vertices: ")
disp(d_result);


function [D] = floydWarshall(D)
for k = 1:length(D)
  D = min(D,D(:,k) + D(k,:));
  disp("Intermediate node: ");
  disp(k);
  disp(D);
end
end