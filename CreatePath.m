%%%%%%%%%%          Create feasible pathes     %%%%%%%%%%%%%%%%
function [pathSet] = CreatePath(cost, start, terminal, Quantity)
    global numNode;
    pathSet = zeros(5,numNode);
    tempCost = cost;
    tempPath = Shortest(cost, start, terminal);
    %%%%%  To be optimized!!!!   %%%%%
    for i = 1:1:length(tempPath)
        pathSet(1,i) = tempPath(i);
    end
    
    for index = 2:1:5
        i = 1;
        j = 2;
        len = length(tempPath);
        while (j <= len)
            u = tempPath(i);
            v = tempPath(j);
            if (u <= 0) || (v <= 0)
                disp 0;
                return;
            end
            if (u == v)
                disp 1;
                return;
            end
            tempCost(u,v) = tempCost(u,v) * 3;
            tempCost(v,u) = tempCost(v,u) * 3;
            i = i + 1;
            j = j + 1;
        end
        tempPath = Shortest(tempCost, start, terminal);
        %pathSet = tempPath;
        for k = 1:1:length(tempPath)
            pathSet(index,k) = tempPath(k);
        end
    end
    return;
end
