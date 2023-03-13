function [neigh] = neighbours(c,n_agents)
neigh = zeros(n_agents);
for i = 1:n_agents
    for j = 1:n_agents
        idx = ismember(c{i}',c{j}', 'rows');
        if sum(idx)>0
            neigh(i,j) = 1;
        else
            neigh(i,j) = 0;
        end
    end
end
% c = 1:size(A, 2);
% d = c(idx); % is your answer