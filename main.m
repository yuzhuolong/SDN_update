%s=fileread('20h10r.json');
%s = fileread('hosts20.json');
%s
%res = parse_json(s);
%res
%res.v2

format long g;

global numNode numSwitch numHost nextNode pre capacity
global roadX roadY cost numRoad numLine
global flowStart flowTerminal flowSize INFINITY DELTA
global inPath pathSet flowOSPF

QUANTITY = 5;
INPUTFILE = 'input_20_10.txt';
FLOWNUMMIN = 4000;
FLOWNUMLEN = 2000;
FLOWNUMMAX = 10000;
TIMEMIN = [0.5, 0.5,    1, 1];
TIMEMAX = [  2,   2,  2.5, 4];
TIMELEN = [0.5, 0.5,  0.5, 1];

%FLOWNUMMIN = 30000;
%FLOWNUMLEN = 10000;
%FLOWNUMMAX = 60000;
%TIMEMIN = [0.5,0.5,1,1];
%TIMEMAX = [2,2,2.5,4];
%TIMELEN = [0.5,0.5,0.5,1];
%%%%%%%%%%          get the graph   %%%%%%%%%%
INFINITY = 100000000;
DELTA = 1E-8;
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
        capacity(u,v) = randi([80,110]);
        capacity(v,u) = capacity(u,v);
    end 
end


%%%%%%%%%%          create feasible pathes    %%%%%%%%%%%%

pathSet = zeros(numNode, numNode, QUANTITY, numNode);
for u = 1:1:numHost
    for v = 1:1:numHost
        if u ~= v
            pathSet(u,v,:,:) = CreatePath(cost, u, v, QUANTITY);
        else
         %   pathSet(u,v) = [];
        end
        
    end
end

inPath = zeros(numNode, numNode, QUANTITY, numRoad);

for u = 1:1:numHost
    for v = 1:1:numHost
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
size(inPath)






picX = zeros(1, 4);
picY1 = zeros(1, 4);
picY2 = zeros(1, 4);
picY3 = zeros(1, 4);
picY4 = zeros(1, 4);


sheetNum = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];








num = 0;

for j = FLOWNUMMIN:FLOWNUMLEN:FLOWNUMMAX
    FLOWNUM = j;
    
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
    flowOSPF = zeros(1, FLOWNUM);
    for i = 1:1:FLOWNUM
        flowOSPF(i) = randi(QUANTITY);
    end

    num = num + 1;
    count = 0;
    for i = TIMEMIN(num):TIMELEN(num):TIMEMAX(num)
        count = count + 1;
        TIMELIMIT = i;
        [lamdaOSPF, lamdaOPT, lamdaGRSU, lamdaMCF, maxTimeDelayMCF] = Emulator(FLOWNUM, QUANTITY, TIMELIMIT, INPUTFILE);
        picX(count) = i;
        picY1(count) = lamdaOSPF;
        picY2(count) = lamdaOPT;
        picY3(count) = lamdaGRSU;
        picY4(count) = lamdaMCF;
        %picY5(count) = maxTimeDelayMCF;
    end
    
    picY5(num) = maxTimeDelayMCF;
    figure(num);
    plot(picX,picY1,'-square',picX,picY2,'-diamond',picX,picY3,'-o',picX,picY4,'-^','LineWidth',2);
    xlabel('T0(s)');
    ylabel('Lambda');
    legend('OSPF','OPT','GRSU','MCF',-1);
    title(['FlowNum=' ,num2str(j)]);
    saveas(gcf, num2str(num),'fig');
    matrix = [picX;picY1;picY2;picY3;picY4]';
    %colnames = {'num of flows(k)','OSPF','OPT','GRSU','MCF'};
    %t = uitable(matrix, colnames);
    xlswrite('result.xls',matrix,sheetNum(num));
end
picY5

