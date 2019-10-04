% demo_MarovChain.m
clear all;

% Add your folder of 'complex_network_measures' 
% We need 'graph_generation.m' and 'cnm' toolbox. 
addpath('D:\Dropbox\My Code\ÀÓ¼±Èñ ±³¼ö´Ô ¼ö¾÷ ÄÚµå\complex_network_measures'); 

network_name = {'RE','SW','RA','SF','HY'};

% Generate artificial complex networks with sparsity 0.2
sparsity = 0.2;
[RE,SW,RA,SF,HY] = graph_generation(sparsity);
% number of nodes 
p = size(RE,1); 


%%% Markov Chain 
% The random walker moves 10 times from the node with the largest number of edges.
figure; 
% Color of nodes 
color_list = colormap(hot(120)); 
color_list = color_list(end-20:-1:1,:); 
for j = 1:5,
    if j == 1,
        A = RE;
    elseif j == 2,
        A = SW;
    elseif j == 3,
        A = RA;
    elseif j == 4,
        A = SF;
    else
        A = HY;
    end

    % Construct a graph for plot 
    [row,col] = find(A); % find the index of nodes of connected edges (row,col)  
    tind = find(row < col); % because A is a symmetric matrix 
    row = row(tind); col = col(tind); 
    G = graph(row,col); 

    % Estimate transition matrix 
    P = A./repmat(sum(A,1),[p 1]); 
    
    % Estimate the degree of nodes 
    degree = sum(A); 
    
    % Choose a node with the largest number of edges 
    [max_degree,ind_node] = max(degree); 
    
    % initial state 
    x = zeros(p,1); 
    x(ind_node) = 1; 

    for i = 1:10, 
        % Determine the size of nodes 
        node_size = max(ceil(x/max(x)*7),1);
        % Determine  the color of nodes 
        node_color = color_list(max(ceil(x/max(x)*100),1),:); 
        
        subplot(5,10,(j-1)*10+i), 
        h = plot(G,'Layout','force','MarkerSize',node_size,'NodeColor',node_color); 
        layout(h,'force','UseGravity',true,'Iterations',1000); 
        labelnode(h,1:p,'')
        if i == 5, 
            xlabel(network_name{j});
        end
        if j == 1, 
            title(['Step ' num2str(i-1)]);
        end
        set(gca,'FontSize',14); 
        
        % Estimate the next state 
        x = P*x; 
    end 
end

% The random walker moves 10 times from the node with the smallest number of edges.
figure; 
% Color of nodes 
color_list = colormap(hot(120)); 
color_list = color_list(end-20:-1:1,:); 
for j = 1:5,
    if j == 1,
        A = RE;
    elseif j == 2,
        A = SW;
    elseif j == 3,
        A = RA;
    elseif j == 4,
        A = SF;
    else
        A = HY;
    end

    % Construct a graph for plot
    [row,col] = find(A); % find the index of nodes of connected edges (row,col)  
    tind = find(row < col); % because A is a symmetric matrix 
    row = row(tind); col = col(tind); 
    G = graph(row,col); 

    % Estimate transition matrix 
    P = A./repmat(sum(A,1),[p 1]); 
    
    % Estimate the degree of nodes 
    degree = sum(A); 
    
    % Choose a node with the largest number of edges 
    [min_degree,ind_node] = min(degree); 
    
    % initial state 
    x = zeros(p,1); 
    x(ind_node) = 1; 

    for i = 1:10, 
        % Determine the size of nodes 
        node_size = max(ceil(x/max(x)*7),1);
        % Determine  the color of nodes 
        node_color = color_list(max(ceil(x/max(x)*100),1),:); 
        
        subplot(5,10,(j-1)*10+i), 
        h = plot(G,'Layout','force','MarkerSize',node_size,'NodeColor',node_color); 
        layout(h,'force','UseGravity',true,'Iterations',1000); 
        labelnode(h,1:p,'')
        if i == 5, 
            xlabel(network_name{j});
        end
        if j == 1, 
            title(['Step ' num2str(i-1)]);
        end
        set(gca,'FontSize',14); 
        
        % Estimate the next state 
        x = P*x; 
    end 
end


% Plot stationary distribution 
eps = 10^(-9); 

figure; 
% Color of nodes 
color_list = colormap(hot(120)); 
color_list = color_list(end-20:-1:1,:); 
for j = 1:5,
    if j == 1,
        A = RE;
    elseif j == 2,
        A = SW;
    elseif j == 3,
        A = RA;
    elseif j == 4,
        A = SF;
    else
        A = HY;
    end
    
    % Estimate transition matrix 
    P = A./repmat(sum(A,1),[p 1]); 
    
    % x = Px 
    [U,S,V] = svd(eye(p)-P); 
    ind_x = find(diag(S) < eps);
    % stationary distribution 
    x = abs(V(:,ind_x)); 
    
    % Construct a graph 
    [row,col] = find(A); % find the index of nodes of connected edges (row,col)  
    tind = find(row < col); % because A is a symmetric matrix 
    row = row(tind); col = col(tind); 
    G = graph(row,col); 

    % Determine the size of nodes 
    node_size = max(ceil(x/max(x)*7),1);
    % Determine  the color of nodes 
    node_color = color_list(max(ceil(x/max(x)*100),1),:); 
        
    subplot(1,5,j),
    h = plot(G,'Layout','force','MarkerSize',node_size,'NodeColor',node_color); 
    layout(h,'force','UseGravity',true,'Iterations',1000); 
    labelnode(h,1:p,'')
    xlabel(network_name{j});
    set(gca,'FontSize',14); 
end 
