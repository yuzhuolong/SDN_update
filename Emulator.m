%s=fileread('20h10r.json');
%s = fileread('hosts20.json');
%s
%res = parse_json(s);
%res
%res.v2

%format long g;
function [lamdaOSPF, lamdaOPT, lamdaGRSU, lamdaMCF, maxTimeDelayMCF] = Emulator(FLOWNUM, QUANTITY, TIMELIMIT, INPUTFILE)
global nextNode pre numNode INFINITY
DELTA = 1E-8;
%%%%%%%%%%      environment setting   %%%%%%%%
%%%%%%%%%%%     version 1 (small)     %%%%%%%%

%FLOWNUM = 5000;
INFINITY = 10000000;
%QUANTITY = 5;
TIMEMODIFY = 0.02;
TIMEINSERT = 0.02;
%TIMELIMIT = 5;
%INPUTFILE = 'input_20_10.txt';


%%%%%%%%%%%     version 2 (big)     %%%%%%%%%%
%{
FLOWNUM = 50000;
INFINITY = 10000000;
QUANTITY = 3;
TIMEMODIFY = 0.02;
TIMEINSERT = 0.02;
TIMELIMIT = 5;
INPUTFILE = 'input_200_100.txt';
%}

%%%%%%%%%%          get the graph   %%%%%%%%%%
numNode = 0;
numSwitch = 0;
numHost = 0;
fileIn = fopen(INPUTFILE,'r');
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
    %flowSize(i) = randi([20,30]) / 1000;
end
BIGFLOWNUM = FLOWNUM / 5;
for i = 1:1:BIGFLOWNUM
    flowSize(i) = randi([20,30]) * 4 / 1000;
end
for i = BIGFLOWNUM+1:1:FLOWNUM
    flowSize(i) = randi([20,30]) / 1000 / 4;
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
%b

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


lamdaOPT = resultGRSU(end);
lamdaOPT


%%%%%%%%%%%%%%%         Greedy Part             %%%%%%%%%%%%%%%%%%%%%%%%
weight = zeros(numNode, numNode);
deployTime = zeros(numNode, FLOWNUM);
timeStack = zeros(1, numNode);

tempWeight = zeros(numNode, numNode);
tempDeployTime = zeros(numNode, FLOWNUM);
tempTimeStack = zeros(1, numNode);

z = zeros(1, FLOWNUM);
for i = 1:1:FLOWNUM
    p = flowOSPF(i);
    z(i) = 1 - resultGRSU((i-1) * QUANTITY + p);
end
[zTemp, zDescend] = sort(z, 'descend');
countUnchange = 0;
for i = 1:1:FLOWNUM
    j = zDescend(i);      %%% choose the flow j  
    %FLOWNUM
    tempY = resultGRSU((j-1)*QUANTITY + 1 : j*QUANTITY);
    [yTemp, yDescend] = sort(tempY, 'descend');
    
    tempWeight = weight;
    tempDeployTime = deployTime;
    tempTimeStack = timeStack; 
    flag = 1;
    %%%%%%%%%%    try to contain the best flow   %%%%%%%%%%%
    for k = 1:1:QUANTITY
        p = yDescend(k);    %%%% choose path NO.p
        %%%%%%%%%%%%%   try   %%%%%%%%%%%%%%%%%
        u = flowStart(j);
        v = flowTerminal(j);
        ii = 1;
        jj = 2;
        uu = pathSet(u,v,p,ii);
        vv = pathSet(u,v,p,jj);
        maxValue = 0;
        while (vv ~= 0)
            tempWeight(uu,vv) = tempWeight(uu,vv) + flowSize(j);
            tempWeight(vv,uu) = tempWeight(uu,vv);
            if (tempWeight(uu,vv) > capacity(uu,vv))
                flag = 0;
                break;
            end
            if ((jj > 2) && (vv >numHost))
                tempDeployTime(vv, j) = tempTimeStack(vv);
                tempTimeStack(vv) = tempTimeStack(vv) + TIMEINSERT;
                maxValue = max(maxValue, tempDeployTime(vv,j) + TIMEINSERT);
                if (tempTimeStack(vv) > TIMELIMIT)
                    flag = 0;
                    break;
                end
            end
            ii = ii + 1;
            jj = jj + 1;
            uu = pathSet(u,v,p,ii);
            vv = pathSet(u,v,p,jj);
        end
        if (flag == 0)
            continue;
        end
        uu = pathSet(u,v,p,2);
        tempDeployTime(uu, j) = maxValue;
        tempTimeStack(uu) = tempDeployTime(uu,j) + TIMEMODIFY;
        if (tempTimeStack(uu) > TIMELIMIT)
            flag = 0;
            continue;
        end
        weight = tempWeight;
        deployTime = tempDeployTime;
        timeStack = tempTimeStack;
        break;
    end
    if (flag == 0)
        countUnchange = countUnchange + 1;
        p = flowOSPF(j);
        u = flowStart(j);
        v = flowTerminal(j);
        for ii = 1:1:numRoad
            if (inPath(u,v,p,ii) - DELTA > 0)
                uu = roadX(ii);
                vv = roadY(ii);
                weight(uu,vv) = weight(uu,vv) + flowSize(j);
                weight(vv,uu) = weight(uu,vv);
            end
        end
    end
end
countUnchange
lamdaGRSU = 0;
for u = 1:1:numNode
    for v = 1:1:numNode
        lamdaGRSU = max(lamdaGRSU, weight(u,v) / capacity(u,v));
    end
end
lamdaGRSU
%b0 = A * resultGRSU;
%b0
%beq0 = Aeq * resultGRSU;
%beq0
%b0
%A(:,1:10)

%{
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
%}
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
maxTimeDelayMCF = 0;
for i = 1:1:numNode
    if (timeDelay(i) > TIMELIMIT)
        delayCountMCF = delayCountMCF + 1;
        maxTimeDelayMCF = max(maxTimeDelayMCF, timeDelay(i));
    end
end
delayCountMCF
maxTimeDelayMCF
lamdaMCF = resultMCF(end);
lamdaMCF





