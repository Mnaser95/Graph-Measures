%% Clearing workspace
clear all
clc

%% User inputs
R_arr = load('');       % Import your adjacency matrix
threshold_level = 0.3;  % This is the edge weight threshold to remove all edges of lesser weight

%% Preprocessing
% Thresholding
G_before = graph(R_arr);                                     % Storing the original data before applying thresholding
weights_column = G_before.Edges.Weight;                      % This is the column containing the edge weights
edges_column = G_before.Edges.EndNodes;                      % This is the column containing the nodes connecting each edges
edges_to_remove = find(weights_column <= threshold_level);
G_after = rmedge(G_before, edges_to_remove);                 % Remove those edges

% Remove self loops
G_after = rmedge(G_after, 1:numnodes(G_after), 1:numnodes(G_after));

% Remove isolated noods
isolated_nodes = find(degree(G_after) == 0);
G_after = rmnode(G_after, isolated_nodes); 

% Generating unweighted matrix and graph
unweighted_matrix = adjacency(G_after);
G_unweighted=graph(unweighted_matrix);
num_nodes = numnodes(G_after);                               % This is the number of nodes after thresholding
num_edges = numedges( G_after );                             % This is the number of edges after thresholding

%% Feature extraction
% Centrality
centrality_eigen_values = centrality(G_unweighted, 'eigenvector'); 
centrality_rank_values = centrality(G_unweighted, 'pagerank'); 
centrality_betw_values = centrality(G_unweighted, 'betweenness'); 
centrality_clos_values = centrality(G_unweighted, 'closeness'); 
centrality_deg_values = centrality(G_unweighted, 'degree'); 
[max_edges_hub, hub_node] = max(centrality_deg_values);
%%%%%%%%%%%%%%
% Charecteristic Path Length (CPL)
shortest_path_distances = distances(G_unweighted);                % returns the shortes path between each pair of nodes
CPL = sum(shortest_path_distances(:))/(num_nodes *(num_nodes -1));
%%%%%%%%%%%%%%
% Global Efficiency (EGlob)
shortest_path_distances_inf_diagonal=shortest_path_distances;
shortest_path_distances_inf_diagonal(1:num_nodes+1:end) = inf;    % to isolate diagonal values
sum_inverse_shortest_paths = sum(1 ./ shortest_path_distances_inf_diagonal , 'all');
EGlob= sum_inverse_shortest_paths / (num_nodes * (num_nodes- 1));
%%%%%%%%%%%%%%
% Local Efficiency at node (Eloc_node) and Average Local Efficiency for graph(Eloc_graph)
Eloc_node = zeros(1,num_nodes);
for node = 1:num_nodes
    neighbor_values = neighbors(G_unweighted, node);
    k = numel(neighbor_values);
    if k < 2
        local_efficiency(node) = 0; % Handles nodes with insufficient neighbors
    else
        shortest_distances = zeros(1, k-1);
        for i = 1:k-1
            for j = i+1:k
                path = shortestpath(G_unweighted, neighbor_values(i), neighbor_values(j));
                shortest_distances(i) = shortest_distances(i) + length(path) - 1; 
            end
        end
        Eloc_node(node) = 1 / (sum(1./shortest_distances) / (k*(k-1)));
    end
end
Eloc_graph=sum(Eloc_node)/(num_nodes);
%%%%%%%%%%%%%%
% Clustering Coefficient (CC) and Average CC (CC_avg)            
CC = zeros(1, num_nodes);
for node = 1:num_nodes
    neighbors_list=neighbors(G_unweighted, node);
    subgraph_list = subgraph(G_unweighted, neighbors_list);
    num_edges_subgraph = numedges(subgraph_list);
    if num_edges_subgraph < 2
        CC(node) = 0;
    else
        possibleEdges = (num_edges_subgraph * (num_edges_subgraph - 1)) / 2;
        CC(node) = num_edges_subgraph / possibleEdges;
    end
end
CC_avg=mean(CC);
%%%%%%%%%%%%%%
% Network Density (ND) 
max_possible_edges = (num_nodes) * (num_nodes - 1)) / 2;
ND = num_edges / max_possible_edges;
%%%%%%%%%%%%%%
% Assortativity (ASSO)
ASSO = assortativity(unweighted_matrix, 0); % 0 means indirected graph
%%%%%%%%%%%%%%
% Modularity (Modu)
num_elements_1 = ceil(num_nodes / 2);
num_elements_2 = floor(num_nodes - num_elements_1);
random_arr = [ones(1, num_elements_1), 2 * ones(1, num_elements_2)]; % To be replaced
Modu= modularity(unweighted_matrix, random_arr);
%%%%%%%%%%%%%%
% Degree Correlations (DC)
DC = zeros(1, numnodes(G_unweighted));
for node = 1:num_nodes
    neighbors_list = neighbors(G_unweighted, node);
    neighbor_degrees = centrality_deg_values(neighbors_list);
    DC(node) = sum(neighbor_degrees) / (numel(neighbors_list) * centrality_deg_values(node));
end
%%%%%%%%%%%%%%
% Transitivity (Trans)
triangles_matrix = unweighted_matrix^3;
trace_triangles = trace(triangles_matrix);
num_connected_triples = sum(sum(triangles_matrix)) / 2; 
Trans = trace_triangles / (num_connected_triples * 3);
%%%%%%%%%%%%%%
% Small Worldness (SW)
C_random = (2 * num_edges) ./ (1:num_nodes * (1:num_nodes - 1));
C_norm = (CC_avg - C_random) / (1 - C_random);
diameter = max(max(shortest_path_distances));
PL_norm = CPL ./ (diameter  .* (diameter - 1));
SW=C_norm./PL_norm;
%%%%%%%%%%%%%%


%% Store the additional computed properties
graph_properties.G_after  =  G_after;
graph_properties.CentEigen= centrality_eigen_values;
graph_properties.CentRank = centrality_rank_values;
graph_properties.CentBetw = centrality_betw_values;
graph_properties.CentClos = centrality_clos_values;
graph_properties.CentDeg=   centrality_deg_values;
graph_properties.NumEdgeHub = max_edges_hub;
graph_properties.HubNode =  hub_node;
graph_properties.GlobEff =  EGlob;
graph_properties.LocEffNode =  Eloc_node;
graph_properties.LocEffGraph =  Eloc_graph;
graph_properties.ClustCoeff =  CC;
graph_properties.ClustCoeffAvg=  CC_avg;
graph_properties.NetDen =  ND;
graph_properties.Assor =  ASSO;
graph_properties.Modul=  Modu;
graph_properties.DegCorr =  DC;
graph_properties.Transi =  Trans;
graph_properties.SmaWor =  SW;

