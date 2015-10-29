%%%%%%%%%%          Calculating the shortest path     %%%%%%%%%%%%%%%%
function [path] = Shortest(cost, start, terminal)

global nextNode pre numNode INFINITY;
distance = ones(1, numNode) * INFINITY;
flag = zeros(1,numNode);
former = zeros(1,numNode);
for u=1:1:numNode
    distance(u) = cost(start, u);
    former(u) = -1;
    if (distance(u) < INFINITY - 1)
        former(u) = start;
    end
end

distance(start) = 0;
flag(start) = 2;
for j=2:1:numNode
    minimum = INFINITY;
    miniNode = start;
    for i=1:1:numNode
        if (flag(i) < 1) && (distance(i) < minimum)
            minimum = distance(i);
            miniNode = i;
        end
    end
    if (miniNode == start)
        break;
    end
    flag(miniNode) = 2;
    for i=1:1:numNode
        if (flag(i) < 1) && (cost(miniNode,i) < INFINITY - 1)
            if (distance(i) > distance(miniNode) + cost(miniNode, i))
                distance(i) = distance(miniNode) + cost(miniNode, i);
                former(i) = miniNode;
            end
        end
    end
end
path(1) = terminal;
i = terminal;
count = 1;
while (former(i) > 0)
    count = count + 1;
    i = former(i);
    %oppositePath(length(oppositePath)+1) = i;
    path(count) = i;
end

path = path(end:-1:1);
return;
end