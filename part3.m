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

% semi stable line till 120 then it increases steeply till 134 then it
% gradually evens out, split fmincon into three parts each using the
% optimized version and simulation from the one before, first from 0 to 120
% then from 120 to 134 then from 134 to 400.


% Define initial parameters and transitions
StoI1 = .002;
StoV1 = .01;
ItoR1 = .005;
ItoD1 = .001;
RtoI1 = .0002;
RtoV1 = .005;
VtoI1 = .000025;

StoI2 = .025;
StoV2 = .05;
ItoR2 = .06;
ItoD2 = .0024;
RtoI2 = .027;
RtoV2 = .06;

StoI3 = .021;
StoV3 = .046;
ItoR3 = .03;
ItoD3 = .000175;
RtoI3 = .03;
RtoV3 = .04;
vaxStart = 150;
initialGuess = [StoI1 StoV1 ItoR1 ItoD1 RtoI1 RtoV1 VtoI1 ... 
                StoI2 StoV2 ItoR2 ItoD2 RtoI2 RtoV2 ... 
                StoI3 StoV3 ItoR3 ItoD3 RtoI3 RtoV3 vaxStart];

fun = @(Aparams)errorNorm(Aparams, cumulativeDeaths, newInfections);

A = [1,1,0,0,0,0,0 ...
     0,0,0,0,0,0 ...
     0,0,0,0,0,0,0; ...

     0,0,0,0,0,0,0 ...
     1,1,0,0,0,0 ...
     0,0,0,0,0,0,0; ...

     0,0,0,0,0,0,0 ...
     0,0,0,0,0,0 ...
     1,1,0,0,0,0,0
     
     0,0,1,1,0,0,0 ...
     0,0,0,0,0,0 ...
     0,0,0,0,0,0,0; ...

     0,0,0,0,0,0,0 ...
     0,0,1,1,0,0 ...
     0,0,0,0,0,0,0; ...

     0,0,0,0,0,0,0 ...
     0,0,0,0,0,0 ...
     0,0,1,1,0,0,0
     
     0,0,0,0,1,1,0 ...
     0,0,0,0,0,0 ...
     0,0,0,0,0,0,0; ...

     0,0,0,0,0,0,0 ...
     0,0,0,0,1,1 ...
     0,0,0,0,0,0,0; ...

     0,0,0,0,0,0,0 ...
     0,0,0,0,0,0 ...
     0,0,0,0,1,1,0];
b = [1;1;1;1;1;1;1;1;1];
Aeq = [];
beq = [];
lb = [.004,0,0,0,0,0,0 ...
      0,0,0,0,0,0 ...
      0,0,0,0,0,0,1];
ub = [.006 ones(1,18) * .05 300];

function error = errorNorm(A, cumulativeDeaths, newInfections)
    [simulatedInfections, simulatedDeaths,~] = simulate(A);

    % RMS error
    error = (sqrt(mean((newInfections - simulatedInfections).^2))) + ...
            (sqrt(mean((cumulativeDeaths - simulatedDeaths).^2)));
end

% Simulation model
function [simulatedInfections, simulatedDeaths,simulatedVaccinated] = simulate(Aparams)
    StoI1 = Aparams(1);
    StoV1 = Aparams(2);
    ItoR1 = Aparams(3);
    ItoD1 = Aparams(4);
    RtoI1 = Aparams(5);
    RtoV1 = Aparams(6);
    VtoI1 = Aparams(7);
    StoI2 = Aparams(8);
    StoV2 = Aparams(9);
    ItoR2 = Aparams(10);
    ItoD2 = Aparams(11);
    RtoI2 = Aparams(12);
    RtoV2 = Aparams(13);
    StoI3 = Aparams(14);
    StoV3 = Aparams(15);
    ItoR3 = Aparams(16);
    ItoD3 = Aparams(17);
    RtoI3 = Aparams(18);
    RtoV3 = Aparams(19);
    vaxStart = Aparams(20);

    n_days = 400;
    x0 = [1; 0; 0; 0; 0];
    x = zeros(5, n_days);
    simulatedDeaths = zeros(1, n_days);
    simulatedInfections = zeros(1, n_days);
    
    A = [1-StoI1, 0, 0, 0, 0;
         StoI1, 1-ItoR1-ItoD1, RtoI1, 0, 0;
         0, ItoR1, 1-RtoI1, 0, 0;
         0, ItoD1, 0, 1, 0;
         0, 0, 0, 0, 1];
    tempx0 = x0;
    x0 = A * x0;
    x(:, 1) = x0;
    simulatedDeaths(1) = x0(4);
    tempA = zeros(5,5);
    tempA(2,1) = StoI1;
    tempA(2,3) = RtoI1;
    tempx0 = tempA * tempx0;
    simulatedInfections(1) = tempx0(2);

    for i = 2:120
        A = [1-StoI1, 0, 0, 0, 0;
        StoI1, 1-ItoR1-ItoD1, RtoI1, 0, 0;
        0, ItoR1, 1-RtoI1, 0, 0;
        0, ItoD1, 0, 1, 0;
        0, 0, 0, 0, 1];
        tempA = zeros(5,5);
        tempA(2,1) = StoI1;
        tempA(2,3) = RtoI1;
        if i >= vaxStart
            A(1,1) = A(1,1) - StoV1; A(5,1) = StoV1;
            A(3,3) = A(3,3) - RtoV1; A(5,3) = RtoV1;
        end
        if i >= 100
            A(2,5) = VtoI1;
            A(5,5) = A(5,5) - VtoI1;
            tempA(2,5) = VtoI1;
        end
        tempx0 = x0;
        x0 = A * x0;
        x(:, i) = x0;
        simulatedDeaths(i) = x0(4) + x(4,i-1);
        tempx0 = tempA * tempx0;
        simulatedInfections(i) = tempx0(2);
    end
    
    for i = 121:134
        A = [1-StoI2, 0, 0, 0, 0;
        StoI2, 1-ItoR2-ItoD2, RtoI2, 0, 0;
        0, ItoR2, 1-RtoI2, 0, 0;
        0, ItoD2, 0, 1, 0;
        0, 0, 0, 0, 1];
        tempA = zeros(5,5);
        tempA(2,1) = StoI2;
        tempA(2,3) = RtoI2;        
        if i >= vaxStart
            A(1,1) = A(1,1) - StoV2; A(5,1) = StoV2;
            A(3,3) = A(3,3) - RtoV2; A(5,3) = RtoV2;
        end
        if i >= 100
            A(2,5) = VtoI1;
            A(5,5) = A(5,5) - VtoI1;
            tempA(2,5) = VtoI1;
        end
        tempx0 = x0;
        x0 = A * x0;
        x(:, i) = x0;
        simulatedDeaths(i) = x0(4) + x(4,i-1);
        tempx0 = tempA * tempx0;
        simulatedInfections(i) = tempx0(2);
    end

    for i = 135:400
        A = [1-StoI3, 0, 0, 0, 0;
        StoI3, 1-ItoR3-ItoD3, RtoI3, 0, 0;
        0, ItoR3, 1-RtoI3, 0, 0;
        0, ItoD3, 0, 1, 0;
        0, 0, 0, 0, 1];
        tempA = zeros(5,5);
        tempA(2,1) = StoI3;
        tempA(2,3) = RtoI3;
        if i >= vaxStart
            A(1,1) = A(1,1) - StoV3; A(5,1) = StoV3;
            A(3,3) = A(3,3) - RtoV3; A(5,3) = RtoV3;
        end
        if i >= 100
            A(2,5) = VtoI1;
            A(5,5) = A(5,5) - VtoI1;
            tempA(2,5) = VtoI1;
        end
        tempx0 = x0;
        x0 = A * x0;
        x(:, i) = x0;
        simulatedDeaths(i) = x0(4) + x(4,i-1);
        tempx0 = tempA * tempx0;
        simulatedInfections(i) = tempx0(2);
    end
    simulatedVaccinated = x(5,:);
end

options = optimoptions('fmincon', 'Algorithm','interior-point', 'MaxIterations', 100000);
[optimalParams, fval] = fmincon(fun, initialGuess, A, b, Aeq, beq, lb, ub, [], options);

disp(['Optimal vaccination start day: ', num2str(optimalParams(20))]);
disp(['Optimal breakthrough infection rate: ', num2str(optimalParams(7))]);
disp(['Optimal Aparams: ', num2str(optimalParams(1:7))]);
disp(['Optimal Aparams: ', num2str(optimalParams(8:14))]);
disp(['Optimal Aparams: ', num2str(optimalParams(15:20))]);


% Simulate with optimal parameters for final output
% [simulated_infections, simulated_deaths] = simulate(initialGuess);
[simulated_infections, simulated_deaths,simulatedVaccinated] = simulate(optimalParams);

figure();
subplot(2,1,1);
plot(1:400, newInfections);
hold on;
plot(1:400, simulated_infections, 'r--');
xlabel('Days');
ylabel('New Infections');
title('newInfections: Observed vs Simulated');
legend('show','Observed New Infections','Simulated New Infections','location','northeast');
hold off;

subplot(2,1,2);
hold on;
plot(1:400, cumulativeDeaths);
plot(1:400, simulated_deaths, 'r--');
xlabel('Days');
ylabel('Cumulative Deaths');
title('cumulativeDeaths: Observed vs Simulated');
legend('show','Observed Cumulative Deaths','Simulated Cumulative Deaths','location','southeast');
hold off;

figure();
plot(1:400,simulatedVaccinated);
xlabel('Days');
ylabel('Fraction of Population Vaccinated');
title('Fraction of Population Vaccinated Over Time');
legend('show','Vaccinated','location','southeast');
hold off;
