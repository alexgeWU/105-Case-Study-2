clear;
close all;

load('COVID_STL.mat');
load('mockdata.mat');

deltaStart = datetime(2021,06,30,"Format",'yyyy-MM-dd');
deltaStart = dates(1,:) == deltaStart;
deltaStart = find(deltaStart);
deltaEnd = datetime(2021,10,27,"Format",'yyyy-MM-dd');
deltaEnd = dates(1,:) == deltaEnd;
deltaEnd = find(deltaEnd);

deltaDates = dates(1,deltaStart:deltaEnd);
deltaCases = cases_STL(1,deltaStart:deltaEnd);
deltaDeaths = deaths_STL(1,deltaStart:deltaEnd);

omicronStart = datetime(2021,10,27,"Format",'yyyy-MM-dd');
omicronStart = dates(1,:) == omicronStart;
omicronStart = find(omicronStart);
omicronEnd = datetime(2022,03,23,"Format",'yyyy-MM-dd');
omicronEnd = dates(1,:) == omicronEnd;
omicronEnd = find(omicronEnd);

omicronDates = dates(1,omicronStart:omicronEnd);
omicronCases = cases_STL(1,omicronStart:omicronEnd);
omicronDeaths = deaths_STL(1,omicronStart:omicronEnd);

figure();
subplot(3,1,1);
hold on; 
plot(deltaDates,deltaCases);
legend('show','cases','location','northwest');
subplot(3,1,2);
plot(deltaDates,deltaDeaths);
legend('show','deaths','location','northwest');
subplot(3,1,3);
plot(deltaDates,ones(1,18)*POP_STL-deltaCases-deltaDeaths);
legend('show','sus','location','northwest');

delta0 = [POP_STL-deltaCases(1)-deltaDeaths(1);
          deltaCases(1);
          0;
          deltaDeaths(1)];

StoI = .002;
ItoR = .0005;
ItoD = .00015;
RtoI = .002;
deltaA = [1-StoI 0           0      0;
          StoI   1-ItoR-ItoD RtoI   0;
          0      ItoR        1-RtoI 0;
          0      ItoD        0      1];

deltax = zeros(4,18);
deltax(:,1) = delta0;
for i = 2:18
    deltax(:,i) = deltaA * deltax(:,i-1);
end

figure();
subplot(3,1,1);
hold on; 
plot(deltaDates,deltax(2,:));
legend('show','cases','location','northwest');
subplot(3,1,2);
plot(deltaDates,deltax(4,:));
legend('show','deaths','location','northwest');
subplot(3,1,3);
plot(deltaDates,deltax(1,:));
legend('show','sus','location','northwest');

figure();
subplot(3,1,1);
hold on; 
plot(omicronDates,omicronCases);
legend('show','cases','location','northwest');
subplot(3,1,2);
plot(omicronDates,omicronDeaths);
legend('show','deaths','location','northwest');
subplot(3,1,3);
plot(omicronDates,ones(1,22)*POP_STL-omicronCases-omicronDeaths);
legend('show','sus','location','northwest');

omicron0 = [POP_STL-omicronCases(1)-omicronDeaths(1);
          omicronCases(1);
          0;
          omicronDeaths(1)];

StoI = .002;
ItoR = .0005;
ItoD = .00015;
RtoI = .002;
omicronA = [1-StoI 0           0      0;
            StoI   1-ItoR-ItoD RtoI   0;
            0      ItoR        1-RtoI 0;
            0      ItoD        0      1];

StoI = .0138;
ItoR = .0005;
ItoD = .00027;
RtoI = StoI;
omicronB = [1-StoI 0           0      0;
            StoI   1-ItoR-ItoD RtoI   0;
            0      ItoR        1-RtoI 0;
            0      ItoD        0      1];

StoI = .00148;
ItoR = .0005;
ItoD = .00018;
RtoI = StoI;
omicronC = [1-StoI 0           0      0;
            StoI   1-ItoR-ItoD RtoI   0;
            0      ItoR        1-RtoI 0;
            0      ItoD        0      1];

% split the omnicron into 3 sections with section 1 from omicronStart to
% omicronEnd1 then from End1 to End2 then from End2 to omicronEnd this is
% made by visually picking points where the slope changed darastically 
omicronEnd1 = datetime(2021,12,22,"Format",'yyyy-MM-dd');
omicronEnd2 = datetime(2022,02,02,"Format",'yyyy-MM-dd');
omicronEnd1 = dates(1,:) == omicronEnd1;
omicronEnd2 = dates(1,:) == omicronEnd2;
omicronEnd1 = find(omicronEnd1);
omicronEnd2 = find(omicronEnd2);

omicronx = zeros(4,22);
omicronx(:,1) = omicron0;
for i = 2:9
    omicronx(:,i) = omicronA * omicronx(:,i-1);
end

for i = 10:15
    omicronx(:,i) = omicronB * omicronx(:,i-1);
end

for i = 16:22
    omicronx(:,i) = omicronC * omicronx(:,i-1);
end

figure();
subplot(3,1,1);
hold on; 
plot(omicronDates,omicronx(2,:));
legend('show','cases','location','northwest');
subplot(3,1,2);
plot(omicronDates,omicronx(4,:));
legend('show','deaths','location','northwest');
subplot(3,1,3);
plot(omicronDates,omicronx(1,:));
legend('show','sus','location','northwest');

StoI = .0012;
ItoR = .0001;
ItoD = .00015;
RtoI = .002;
omicronA = [1-StoI 0           0      0;
            StoI   1-ItoR-ItoD RtoI   0;
            0      ItoR        1-RtoI 0;
            0      ItoD        0      1];

StoI = .00374;
ItoR = .0001;
ItoD = .00015;
RtoI = StoI;
omicronB = [1-StoI 0           0      0;
            StoI   1-ItoR-ItoD RtoI   0;
            0      ItoR        1-RtoI 0;
            0      ItoD        0      1];

StoI = .00143;
ItoR = .0005;
ItoD = .00018;
RtoI = StoI;
omicronC = [1-StoI 0           0      0;
            StoI   1-ItoR-ItoD RtoI   0;
            0      ItoR        1-RtoI 0;
            0      ItoD        0      1];

% split the omnicron into 3 sections with section 1 from omicronStart to
% omicronEnd1 then from End1 to End2 then from End2 to omicronEnd this is
% made by visually picking points where the slope changed darastically 
omicronEnd1 = datetime(2021,12,22,"Format",'yyyy-MM-dd');
omicronEnd2 = datetime(2022,02,02,"Format",'yyyy-MM-dd');
omicronEnd1 = dates(1,:) == omicronEnd1;
omicronEnd2 = dates(1,:) == omicronEnd2;
omicronEnd1 = find(omicronEnd1);
omicronEnd2 = find(omicronEnd2);

omicronx = zeros(4,22);
omicronx(:,1) = omicron0;
for i = 2:9
    omicronx(:,i) = omicronA * omicronx(:,i-1);
end

for i = 10:15
    omicronx(:,i) = omicronB * omicronx(:,i-1);
end

for i = 16:22
    omicronx(:,i) = omicronC * omicronx(:,i-1);
end

figure();
subplot(3,1,1);
hold on; 
plot(omicronDates,omicronx(2,:));
legend('show','cases','location','northwest');
subplot(3,1,2);
plot(omicronDates,omicronx(4,:));
legend('show','deaths','location','northwest');
subplot(3,1,3);
plot(omicronDates,omicronx(1,:));
legend('show','sus','location','northwest');

figure();
subplot(2,1,1);
plot(1:400,newInfections);
legend('show','infections','location','northwest');
subplot(2,1,2);
plot(1:400,cumulativeDeaths);
legend('show','deaths','location','northwest');
