%s=fileread('20h10r.json');
%s = fileread('hosts20.json');
%s
%res = parse_json(s);
%res
%res.v2

%format long g;
function [lamdaOSPF, lamdaOPT, lamdaGRSU, lamdaMCF, maxTimeDelayMCF] = Emulator(FLOWNUM, QUANTITY, TIMELIMIT, INPUTFILE)

global numNode numSwitch numHost nextNode pre capacity
global roadX roadY cost numRoad numLine
global flowStart flowTerminal flowSize INFINITY DELTA
global inPath pathSet flowOSPF

%%%%%%%%%%      environment setting   %%%%%%%%
%%%%%%%%%%%     version 1 (small)     %%%%%%%%

%FLOWNUM = 5000;

%QUANTITY = 5;
TIMEMODIFY = 0.015;
TIMEINSERT = 0.01;
%TIMELIMIT = 5;
%INPUTFILE = 'input_20_10.txt';
TIMESEG = 0.005;
NUMSEG = TIMELIMIT / TIMESEG;
NUMINSERTSEG = TIMEINSERT / TIMESEG;
NUMMODIFYSEG = TIMEMODIFY / TIMESEG;

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





%%%%%%%%%%%   obtain the original(OSPF) state  %%%%%%%%%%%%%%

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
ub = [ones(NN-1, 1); 100000];
tic;
%options = optimset('Display' , 'off' , 'LargeScale' , 'off' , 'Simplex' , 'on');
options = optimset('Display' , 'on' , 'LargeScale' , 'on' , 'Simplex' , 'on','maxiter',90000000);
[resultGRSU, fval, exitflag, output, lambdaFin] = linprog(f,A,b,Aeq,beq,lb,ub,[],options);
exitflag
%resultGRSU = linprog(f,A,b,Aeq,beq,lb,ub);
toc;


lamdaOPT = resultGRSU(end);
lamdaOPT


%%%%%%%%%%%%%%%         Greedy Part             %%%%%%%%%%%%%%%%%%%%%%%%
weight = zeros(numNode, numNode);
%deployTime = zeros(1, numNode);
%timeStack = zeros(1, numNode);
timeStack = zeros(numNode, NUMSEG);

tempWeight = zeros(numNode, numNode);
%tempDeployTime = zeros(1, numNode);
%tempTimeStack = zeros(1, numNode);
timeStack = zeros(numNode, NUMSEG);

%{
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
%}

z = zeros(1, FLOWNUM);
for i = 1:1:FLOWNUM
    p = flowOSPF(i);
    z(i) = 1 - resultGRSU((i-1) * QUANTITY + p);
end
[zTemp, zDescend] = sort(z, 'descend');
countUnchange = 0;
tempInsertSeg = zeros(1, NUMINSERTSEG);
tempModifySeg = zeros(1, NUMMODIFYSEG);
tic
for i = 1:1:FLOWNUM
    j = zDescend(i);      %%% choose the flow j  
    %FLOWNUM
    tempY = resultGRSU((j-1)*QUANTITY + 1 : j*QUANTITY);
    [yTemp, yDescend] = sort(tempY, 'descend');
    
    tempWeight = weight;
    %tempDeployTime = deployTime;
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
        maxValueSeg = 1;
        while (vv ~= 0)
            
            tempWeight(uu,vv) = tempWeight(uu,vv) + flowSize(j);
            tempWeight(vv,uu) = tempWeight(uu,vv);
            if (tempWeight(uu,vv) > capacity(uu,vv))
                flag = 0;
                break;
            end
            
            if ((jj > 2) && (vv >numHost))
                %%tempDeployTime(vv, j) = tempTimeStack(vv);
                %%tempTimeStack(vv) = tempTimeStack(vv) + TIMEINSERT;
                found = 0;
                probe = 1;
                while (probe <= NUMSEG - NUMINSERTSEG+1)
                    tempSum = sum(tempTimeStack(vv, probe:probe+NUMINSERTSEG-1));
                    if (tempSum - DELTA > 0)
                        probe = probe + tempSum;
                        continue;
                    end
                    found = 1;
                    break;
                end
                %{
                for probe = 1:1:(NUMSEG - NUMINSERTSEG)
                    %tempInsertSeg = tempTimeStack(vv, probe:probe+NUMINSERTSEG);
                    %if (sum(tempInsertSeg(1:NUMINSERTSEG)) - DELTA > 0)
                    if (sum(tempTimeStack(vv, probe:probe+NUMINSERTSEG)) - DELTA > 0)
                        continue;
                    end
                    found = 1;
                    break;
                end
                %}
                if (found == 0)
                    flag = 0;
                    break;
                end
                tempTimeStack(vv, probe:probe+NUMINSERTSEG) = 1;
               % tempDeployTime(vv) = (probe - 1) * TIMESEG;
                maxValueSeg = max(maxValueSeg, probe + NUMINSERTSEG);
                
                
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
        found = 0;
        
        probe = maxValueSeg;
        while (probe <= NUMSEG-NUMMODIFYSEG + 1)
            tempSum = sum(tempTimeStack(uu, probe:probe+NUMMODIFYSEG-1));
            if (tempSum - DELTA > 0)
                probe = probe + tempSum;
                continue;
            end
            found = 1;
            break;
        end
        %{
        for probe = maxValueSeg:1:NUMSEG - NUMMODIFYSEG
             %tempModifySeg = tempTimeStack(uu, probe:probe+NUMMODIFYSEG);
             %if (sum(tempModifySeg(1:NUMMODIFYSEG)) - DELTA > 0)
             if (sum(tempTimeStack(uu, probe:probe+NUMMODIFYSEG)) - DELTA > 0)
                continue;
             end
             found = 1;
             break;
        end
        %}
        if (found == 0)
            flag = 0; 
            continue;
        end
        tempTimeStack(uu, probe:probe+NUMMODIFYSEG) = 1;
        %tempDeployTime(uu) = (probe - 1) * TIMESEG;
        %tempDeployTime(uu, j) = maxValueSeg;
        %tempTimeStack(uu) = tempDeployTime(uu,j) + TIMEMODIFY;
        %if (tempTimeStack(uu) > TIMELIMIT)
        %    flag = 0;
        %    continue;
        %end
        weight = tempWeight;
       % deployTime = tempDeployTime;
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

toc
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
%{
for i = 1:1:FLOWNUM
    tempMax = 0;
    maxJ = 0;
    for j = 1:1:QUANTITY
        if (resultMCF((i-1)*QUANTITY+j) > tempMax)
            tempMax = resultMCF((i-1)*QUANTITY+j);
            maxJ = j;
        end
    end
    if maxJ == 0
        pause(5000);
        break;
    end
    flowMCF(i) = maxJ;
end
%}
for i = 1:1:FLOWNUM
    randNum = rand;
    j = 1;
    while (randNum > resultMCF((i-1) * QUANTITY + j)) && (j < QUANTITY)
        randNum = randNum - resultMCF((i-1) * QUANTITY + j);
        j = j + 1;
    end
    flowMCF(i) = j;
end
lamdaMCF = 0;
for i = 1:1:numRoad
    lamdaTemp = 0;
    for j = 1:1:FLOWNUM
        p = flowMCF(j);
        u = flowStart(j);
        v = flowTerminal(j);
        if inPath(u,v,p,i) - DELTA > 0
            lamdaTemp = lamdaTemp + flowSize(j);
        end
    end
    u = roadX(i);
    v = roadY(i);
    lamdaTemp = lamdaTemp / capacity(u,v);
    lamdaMCF = max(lamdaMCF, lamdaTemp);
end
lamdaMCF
lamdaMCF = resultMCF(end);
lamdaMCF

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






