%s=fileread('20h10r.json');
%s = fileread('hosts20.json');
%s
%res = parse_json(s);
%res
%res.v2

format long g;
global nextNode pre numNode INFINITY
DELTA = 1E-8;
%%%%%%%%%%      environment setting   %%%%%%%%
%%%%%%%%%%%     version 1 (small)     %%%%%%%%
%{
FLOWNUM = 5000;
INFINITY = 10000000;
QUANTITY = 5;
TIMEMODIFY = 0.2;
TIMEINSERT = 0.2;
TIMELIMIT = 5;
%}

%%%%%%%%%%%     version 2 (big)     %%%%%%%%%%
FLOWNUM = 200000;
INFINITY = 10000000;
QUANTITY = 3;
TIMEMODIFY = 0.2;
TIMEINSERT = 0.2;
TIMELIMIT = 5;


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
    if (cost(v,u) >= INFINITY - DELTA)
        numRoad = numRoad + 1;
        roadX(numRoad) = u;
        roadY(numRoad) = v;
        cost(v,u) = 1;
        cost(u,v) = 1;
        pre(u) = pre(u) + 1;
        nextNode(u, pre(u)) = v;
        pre(v) = pre(v) + 1;
        nextNode(v, pre(v)) = u;
        capacity(u,v) = randi([2*15,3*15]);
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
    flowSize(i) = randi([20,30]) / 1000;
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
        if u ~= v
        for p = 1:1:QUANTITY
            i = 1;
            j = 2;
            uu = pathSet(u,v,p,i);
            vv = pathSet(u,v,p,j);
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
                uu = pathSet(u,v,p,i);
                vv = pathSet(u,v,p,j);
            end
        end
        end
    end
end



%%%%%%%%%%%   obtain the original(OSPF) state  %%%%%%%%%%%%%%
flowOSPF = zeros(1, FLOWNUM);
for i = 1:1:FLOWNUM
    flowOSPF(i) = randi(QUANTITY);
end
lamdaOSPF = 0;
for i = 1:1:numRoad
    lamdaTemp = 0;
    for j = 1:1:FLOWNUM
        p = flowOSPF(j);
        u = flowStart(j);
        v = flowTerminal(j);
        if inPath(u,v,p,i) - DELTA > 0
            lamdaTemp = lamdaTemp + flowSize(j);
        end
    end
    u = roadX(i);
    v = roadY(i);
    lamdaTemp = lamdaTemp / capacity(u,v);
    lamdaOSPF = max(lamdaOSPF, lamdaTemp);
end
lamdaOSPF

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
            if tempNode <= numHost 
                break;
            end
            if (p ~= flowOSPF(i))
                if (tempNode == w)
                    timeCost(tempNode, i, p) = TIMEMODIFY;
                else 
                    timeCost(tempNode, i, p) = TIMEINSERT;
                end
            end
            j = j + 1;
            tempNode = pathSet(u,v,p,j);
            if tempNode <= numHost 
                break;
            end
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
            if (inPath(u,v,p,i) - DELTA > 0)
                A(numNode + i, (j-1)*QUANTITY + p) = flowSize(j);
            end
        end
    end
    A(numNode + i, end) = -capacity(roadX(i), roadY(i));
end

b = zeros(numNode + numRoad, 1);
for i = 1:1:numNode
    b(i, 1) = TIMELIMIT;
end
b

%Aeq = zeros(FLOWNUM, NN);
Aeq = sparse(FLOWNUM, NN);
for i = 1:1:FLOWNUM
    for j = 1:1:QUANTITY
        Aeq(i,(i-1) * QUANTITY + j) = 1;
    end
end

beq = ones(FLOWNUM, 1);

lb = zeros(NN, 1);
ub = [ones(NN-1, 1); 10];
tic;
%options = optimset('Display' , 'off' , 'LargeScale' , 'off' , 'Simplex' , 'on');
options = optimset('Display' , 'final' , 'LargeScale' , 'on' , 'Simplex' , 'on');
resultGRSU = linprog(f,A,b,Aeq,beq,lb,ub,[],options);
%resultGRSU = linprog(f,A,b,Aeq,beq,lb,ub);
toc;


lamda = resultGRSU(end);
lamda

b0 = A * resultGRSU;
b0
beq0 = Aeq * resultGRSU;
beq0
%b0
%A(:,1:10)

A = zeros(numNode, NN);
for i = 1:1:numNode
    for j = 1:1:FLOWNUM
        for p = 1:1:QUANTITY
            A(i, (j-1)*QUANTITY + p) = timeCost(i, j, p);
        end
    end
end
timeDelay = A * resultGRSU;
timeDelay
delayCountGRSU = 0;
for i = 1:1:numNode
    if (timeDelay(i) > TIMELIMIT)
        delay0 = timeDelay(i);
       % delay0
       % TIMELIMIT
        delayCountGRSU = delayCountGRSU + 1;
    end
end
delayCountGRSU

%%%%%%%%%%%%%%      obtain the solution of MCF    %%%%%%%%%%%%%%%%


f = zeros(NN, 1);
f(end) = 1;

A = zeros(numRoad, NN);

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
    A(i, end) = -capacity(roadX(i), roadY(i));
end

b = zeros(numRoad, 1);

%Aeq = zeros(FLOWNUM, NN);
Aeq = sparse(FLOWNUM, NN);
for i = 1:1:FLOWNUM
    for j = 1:1:QUANTITY
        Aeq(i,(i-1) * QUANTITY + j) = 1;
    end
end

beq = ones(FLOWNUM, 1);

lb = zeros(NN, 1);
ub = [ones(NN-1, 1); 10];
tic;
options = optimset('Display' , 'off' , 'LargeScale' , 'on' , 'Simplex' , 'on');
resultMCF = linprog(f,A,b,Aeq,beq,lb,ub,[],options);
toc;
A = zeros(numNode, NN);
for i = 1:1:numNode
    for j = 1:1:FLOWNUM
        for p = 1:1:QUANTITY
            A(i, (j-1)*QUANTITY + p) = timeCost(i, j, p);
        end
    end
end
timeDelay = A * resultMCF(1:NN);
delayCountMCF = 0;
for i = 1:1:numNode
    if (timeDelay(i) > TIMELIMIT)
        delayCountMCF = delayCountMCF + 1;
    end
end
delayCountMCF
lamdaMCF = resultMCF(end);
lamdaMCF





