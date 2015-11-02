%s=fileread('20h10r.json');
%s = fileread('hosts20.json');
%s
%res = parse_json(s);
%res
%res.v2

format long g;

QUANTITY = 5;
picX = zeros(1, 5);
picY1 = zeros(1, 5);
picY2 = zeros(1, 5);
picY3 = zeros(1, 5);
picY4 = zeros(1, 5);


sheetNum = [1,2,3,4,5];
INPUTFILE = 'input_20_10.txt';
for i = 3:2:3
    TIMELIMIT = i;
    count = 0;
    for j = 2000:2000:10000
        count = count + 1;
        FLOWNUM = j;
        [lamdaOSPF, lamdaOPT, lamdaGRSU, lamdaMCF, maxTimeDelayMCF] = Emulator(FLOWNUM, QUANTITY, TIMELIMIT, INPUTFILE);
        picX(count) = j / 1000;
        picY1(count) = lamdaOSPF;
        picY2(count) = lamdaOPT;
        picY3(count) = lamdaGRSU;
        picY4(count) = lamdaMCF;
        picY5(count) = maxTimeDelayMCF;
    end
    num = (i-1) / 2;
    figure(num);
    plot(picX,picY1,'-square',picX,picY2,'-diamond',picX,picY3,'-o',picX,picY4,'-^','LineWidth',2);
    xlabel('Number of Flows(k)');
    ylabel('Lambda');
    legend('OSPF','OPT','GRSU','MCF',-1);
    title(['T0=' ,num2str(i),'s']);
    matrix = [picX;picY1;picY2;picY3;picY4]';
    %colnames = {'num of flows(k)','OSPF','OPT','GRSU','MCF'};
    %t = uitable(matrix, colnames);
    xlswrite('result.xls',matrix,sheetNum(num));
end
picY5

%{
sheetNum = [1,2,3,4,5];
INPUTFILE = 'input_200_100.txt';
for i = 5:2:11
    TIMELIMIT = i;
    count = 0;
    for j = 20000:20000:80000
        count = count + 1;
        FLOWNUM = j;
        [lamdaOSPF, lamdaOPT, lamdaGRSU, lamdaMCF, maxTimeDelayMCF] = Emulator(FLOWNUM, QUANTITY, TIMELIMIT, INPUTFILE);
        picX(count) = j / 10000;
        picY1(count) = lamdaOSPF;
        picY2(count) = lamdaOPT;
        picY3(count) = lamdaGRSU;
        picY4(count) = lamdaMCF;
        
    end
    num = (i-1) / 2;
    figure(num);
    plot(picX,picY1,'-square',picX,picY2,'-diamond',picX,picY3,'-o',picX,picY4,'-^','LineWidth',2);
    xlabel('Number of Flows(k)');
    ylabel('Lambda');
    legend('OSPF','OPT','GRSU','MCF',-1);
    title(['T0=' ,num2str(i),'s']);
    matrix = [picX;picY1;picY2;picY3;picY4]';
    %colnames = {'num of flows(k)','OSPF','OPT','GRSU','MCF'};
    %t = uitable(matrix, colnames);
    xlswrite('resultBig.xls',matrix,sheetNum(num));
end
%}
