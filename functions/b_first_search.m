function [ visited, lvl ] = b_first_search( A, root)


visited = int16([]);
active = int16([]);
queue = int16([]);

visited = [visited; root];
queue = [queue; root];

while ~isempty(queue)
active = queue(1);
offsprings = find(A(:,active)~=0);
offsprings = offsprings(offsprings~=active);
if ~isempty(intersect(visited, offsprings))
[~, b] = ismember(visited, offsprings);
b(b==0)= [];
offsprings(b) = [];
clear b
end

visited = [visited; offsprings];
queue = [queue; offsprings];
queue = queue(queue~=active);


clear offsprings
end

levels = zeros(length(visited),1);

A(A == 2)=1;

D = graphallshortestpaths(sparse(A));
for i = 1:length(visited)

levels(i,1) = D(visited(i),root);
    
end

lvl(visited,1)=levels;

end

