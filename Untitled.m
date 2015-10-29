%s=fileread('20h10r.json');
%s = fileread('hosts20.json');
%s
%res = parse_json(s);
%res
%res.v2

%%%%%%%%%%      environment setting   %%%%%%%%
flowNum = 5000;
INFINITY = 100000;
global nextNode pre numNode INFINITY

%%%%%%%%%%          get the graph   %%%%%%%%%%
numNode = 0;
numSwitch = 0;
numHost = 0;
fileIn = fopen('input_20_10.txt','r');
numNode = fscanf(fileIn, '%d', 1);
numHost = fscanf(fileIn, '%d', 1);
numSwitch = fscanf(fileIn, '%d', 1);

numLine = fscanf(fileIn, '%d', 1);

cost = ones(numNode, numNode) * INFINITY;
pre = zeros(1, numNode);
nextNode = zeros(numNode, numNode);
for line = 1:1:numLine
    v = fscanf(fileIn, '%d', 1);
    u = fscanf(fileIn, '%d', 1);
    if (cost(v,u) == INFINITY)
        cost(v,u) = 1;
        cost(u,v) = 1;
        pre(u) = pre(u) + 1;
        nextNode(u, pre(u)) = v;
        pre(v) = pre(v) + 1;
        nextNode(v, pre(v)) = u;
    end 
end

%%%%%%%%%%          generate the flows     %%%%%%%%%%%%%%%
flowStart = zeros(1, flowNum);
flowTerminal = zeros(1, flowNum);
for i = 1:1:flowNum
    randArray = randperm(numHost);
    u = randArray(1);
    v = randArray(2);
    flowStart(i) = u;
    flowTerminal(i) = v;
end

%%%%%%%%%%          create feasible pathes    %%%%%%%%%%%%

pathSet = CreatePath(cost, 10, 4, 3);
pathSet
%{
for u = 1:1:numNode
    for v = 1:1:numNode
        if u ~= v
            pathSet(u,v) = CreatePath(cost, u, v, 5);
        else
         %   pathSet(u,v) = [];
        end
        
    end
end
%}


    






