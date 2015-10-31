%s=fileread('20h10r.json');
%s = fileread('hosts20.json');
%s
%res = parse_json(s);
%res
%res.v2

%%%%%%%%%%      environment setting   %%%%%%%%
global nextNode pre numNode INFINITY
FLOWNUM = 5000;
INFINITY = 100000;
QUANTITY = 5;
TIMEMODIFY = 20;
TIMEINSERT = 20;
TIMELIMIT = 5000;


%%%%%%%%%%          get the graph   %%%%%%%%%%
numNode = 0;
numSwitch = 0;
numHost = 0;
fileIn = fopen('input_20_10.txt','r');
numNode = fscanf(fileIn, '%d', 1);
numHost = fscanf(fileIn, '%d', 1);
numSwitch = fscanf(fileIn, '%d', 1);

numLine = fscanf(fileIn, '%d', 1);

numRoad = 0;
cost = ones(numNode, numNode) * INFINITY;
pre = zeros(1, numNode);
nextNode = zeros(numNode, numNode);
capacity = zeros(numNode, numNode);
roadX = zeros(1, numNode * numNode);
roadY = zeros(1, numNode * numNode);
for line = 1:1:numLine
    v = fscanf(fileIn, '%d', 1);
    u = fscanf(fileIn, '%d', 1);
    if (cost(v,u) == INFINITY)
        numRoad = numRoad + 1;
        roadX(numRoad) = u;
        roadY(numRoad) = v;
        cost(v,u) = 1;
        cost(u,v) = 1;
        pre(u) = pre(u) + 1;
        nextNode(u, pre(u)) = v;
        pre(v) = pre(v) + 1;
        nextNode(v, pre(v)) = u;
        capacity(u,v) = randi([20*1500,30*1500]);
        capacity(v,u) = capacity(u,v);
    end 
end

%%%%%%%%%%          generate the flows     %%%%%%%%%%%%%%%
flowStart = zeros(1, FLOWNUM);
flowTerminal = zeros(1, FLOWNUM);
flowSize = zeros(1, FLOWNUM);
for i = 1:1:FLOWNUM
    randArray = randperm(numHost);
    u = randArray(1);
    v = randArray(2);
    flowStart(i) = u;
    flowTerminal(i) = v;
    flowSize(i) = randi([20,30]);
end

%%%%%%%%%%          create feasible pathes    %%%%%%%%%%%%

%pathSet1 = CreatePath(cost, 10, 4, 5);
%pathSet1

pathSet = zeros(numNode, numNode, QUANTITY, numNode);
for u = 1:1:numNode
    for v = 1:1:numNode
        if u ~= v
            pathSet(u,v,:,:) = CreatePath(cost, u, v, QUANTITY);
        else
         %   pathSet(u,v) = [];
        end
        
    end
end
inPath = zeros(numNode, numNode, QUANTITY, numRoad);
for u = 1:1:numNode
    for v = 1:1:numNode
        for p = 1:1:QUANTITY
            tempPath = pathSet(u, v, p);
            i = 1;
            j = 2;
            uu = tempPath(i);
            vv = tempPath(j);
            while (vv ~= 0)
                for tempRoad = 1:1:numRoad
                    uuu = roadX(tempRoad);
                    vvv = roadY(tempRoad);
                    if ((uu == uuu) && (vv == vvv))||((uu == vvv) && (vv == uuu))
                        inPath(u, v, p, tempRoad) = 1;
                    end
                end
                i = i + 1;
                j = j + 1;
                uu = tempPath(i);
                vv = tempPath(j);
            end
        end
    end
end



%%%%%%%%%%%   obtain the original(OSPF) state  %%%%%%%%%%%%%%
flowOSPF = zeros(1, FLOWNUM);
for i = 1:1:FLOWNUM
    flowOSPF(i) = randi(QUANTITY);
end

%%%%%%%%%%      obtain the update time on each switch (flow,path) %%%%%%%%
timeCost = zeros(numNode, FLOWNUM,  QUANTITY);
for i = 1:1:FLOWNUM
    for p = 1:1:QUANTITY
        u = flowStart(i);
        v = flowTerminal(i);
        w = pathSet(u, v, p, 2);
        j = 2;
        tempNode = pathSet(u,v,p,j);
        while (tempNode ~= 0)
            if (p ~= flowOSPF(i))
                if (tempNode == w)
                    timeCost(tempNode, i, p) = TIMEMODIFY;
                else 
                    timeCost(tempNode, i, p) = TIMEINSERT;
                end
            end
            j = j + 1;
            tempNode = pathSet(u,v,p,j);
        end
    end
end

%%%%%%%%%%      obtain the GRSU solution   %%%%%%%%%%%%%%%%%%
NN = QUANTITY * FLOWNUM + 1;

f = zeros(NN, 1);
f(end) = 1;

A = zeros(numNode + numRoad, NN);
for i = 1:1:numNode
    for j = 1:1:FLOWNUM
        for p = 1:1:QUANTITY
            A(i, (j-1)*QUANTITY + p) = timeCost(i, j, p);
        end
    end
end
for i = 1:1:numRoad
    for j = 1:1:FLOWNUM
        for p = 1:1:QUANTITY
            u = flowStart(j);
            v = flowTerminal(j);
            if (inPath(u,v,p,i) > 0)
                A(i, (j-1)*QUANTITY + p) = flowSize(j);
            end
        end
    end
    A(i, -1) = -Capacity(roadX(i), roadY(i));
end

b = zeros(numNode + numRoad, 1);
for i = 1:1:numNode
    b(i, 1) = TIMELIMIT;
end

Aeq = zeros(FLOWNUM, NN);
for i = 1:1:FLOWNUM
    for j = 1:1:QUANTITY
        Aeq(i,(i-1) * QUANTITTY + j) = 1;
    end
end

beq = ones(FLOWNUM, 1);

lb = zeros(NN, 1);
ub = [ones(NN-1, 1); 10];

options = optimset('Display' , 'off' , 'LargeScale' , 'off' , 'Simplex' , 'on');
result = linprog(f,A,b,Aeq,beq,lb,ub,[],options);

%%%%%%%%%%%%%%      obtain the solution of MCF    %%%%%%%%%%%%%%%%



    






