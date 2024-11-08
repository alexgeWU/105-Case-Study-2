clear;
close all;

load('mockdata.mat');

% figure();
% subplot(2,1,1);
% plot(1:400,newInfections);
% legend('show','infections','location','northwest');
% subplot(2,1,2);
% plot(1:400,cumulativeDeaths);
% legend('show','deaths','location','northwest');

% semi stable line till 120 then it increases steeply till 134 then it
% gradually evens out, split fmincon into three parts each using the
% optimized version and simulation from the one before, first from 0 to 120
% then from 120 to 134 then from 134 to 400.


%% fmincon parameters
% specifies function used, its the errorNorm() function created later
fun = @(Aparams)errorNorm(Aparams, cumulativeDeaths, newInfections);

% Define initial guesses and gives then a name as a comment
initialGuess = [
    .0138;     % StoI1
    .048654;   % StoV1
    .0005;     % ItoR1
    .0027;     % ItoD1
    .0002      % ItoS1
    .01;       % RtoI1
    .05;       % RtoV1

    .00001;    % VtoI1
    .02238874; % StoI2
    150        % vaxStart
];

% linear inequality constraints such that each coloumn is at max 1
% and logically infection rates are lower for recovered people and even
% lower for vaccinated people
A = [1,1,0,0,0,0,0,0,0,0; ...  % StoI1 and StoV1 <= 1
     0,0,1,1,1,0,0,0,0,0; ...  % ItoR1 and ItoD1 <= 1
     0,0,0,0,0,1,1,0,0,0; ...  % RtoI1 and RtoV1 <= 1
     0,1,0,0,0,0,0,0,1,0; ...  % StoI2 and StoV1 <= 1
     0,0,-1,0,1,0,0,0,0,0; ... % ItoR1 >= ItoS1
     0,0,0,0,0,-1,0,1,0,0; ... % RtoI1 >= VtoI1
     -1,0,0,0,0,1,0,0,0,0]; ...% StoI1 >= RtoI1

b = [1;1;1;1;0;0;0];

% there are no linear equality constrains
Aeq = [];
beq = [];

% Lower and upper bounds for optimization
lb = [.001; % for the first 120 or so days the minimum value of the given 
            % newInfections is .0015 so the value we simulate shouldn't
            % vary that much on the lower bound
      0;    % vaccination rate shouldn't be negative
      0;    % infected to recovered rate shouldn't be negative
      0;    % infected to dead rate shouldn't be negative
      0;    % infected to susceptible shouldn't be negative
      0;    % recovered to infected rate shouldn't be negative
      0;    % recovered to vaccinated rate shouldn't be negative

      0;    % vaccinated to infected (breakout rate) shouldn't be negative
      0.01; % the second increased infection rate appears to be higher than
            % the first infection time period by a substantial margin,
            % between day 0 and 120 the max is .0078 so we set this bound
            % above that      
      130]; % vaccinations are assumed to start sometime around the when 
            % the cumulative deaths starts to even out this is the earliest
            % point where it starts to even out

            % realistically no change will be greater than 50% of any given
            % group so we set that as a very general upper bound.
ub = [.007; % the highest newInfections day before the jump at 120 days
            % is just below this

      .5;   
      .5;    
      .5;    
      .5;
      .5; 
      .001; % if someone recovers, they are less likely to get infected
            % compared to the susceptible group, so the ub for RtoI is set 
            % to the lb of StoI

      .001; % if someone is vaccinated, they are less likely to get infected
            % compated to the susceptible group, so the ub is the lb of
            % StoI
      .5;  % 
      155]; % 

%% Function to Optimize for fmincon
% Using the input parameters it calls the simulate function which simulates
% the 400 days with the given parameters, then it sets error to the root
% mean square error between simulated and given new infections plus 
% simulated and given cumulative deaths
function error = errorNorm(A, cumulativeDeaths, newInfections)
    [simulatedInfections,simulatedDeaths,~] = simulate(A);

    % RMS error
    error = (sqrt(mean((newInfections - simulatedInfections).^2))) + ...
            (sqrt(mean((cumulativeDeaths - simulatedDeaths).^2)));
end

%% Main Simulation Model
% Simulation model as a seperate function from error so that we can
% visualize the answers and the simulation after fmincon gives the best
% parameters by calling simulate ourselves.
function [simulatedInfections, simulatedDeaths,x] = simulate(Aparams)
    % gets parameters into their own variables
    StoI1 = Aparams(1);
    StoV  = Aparams(2);
    ItoR  = Aparams(3);
    ItoD  = Aparams(4);
    ItoS  = Aparams(5);
    RtoI  = Aparams(6);
    RtoV  = Aparams(7);
    VtoI  = Aparams(8);
    StoI2 = Aparams(9);
    vaxStart = Aparams(10);

    % sets up initial length of simulation, inital assumed to be all
    % susceptible, and initiates other variables to calculate.
    numberDays = 400;
    x0 = [1; 0; 0; 0; 0];
    x = zeros(5, numberDays);
    simulatedDeaths = zeros(1, numberDays);
    simulatedInfections = zeros(1, numberDays);
    
    % First day simulated seperately so that we can calculate cumulative
    % deaths from numberDays >= 2 buy adding to the day before
    
    % Constructs matrix A a coe fficient matrix with the given parameters
    A = [1-StoI1, ItoS, 0, 0, 0;
         StoI1, 1-ItoR-ItoD-ItoS, RtoI, 0, 0;
         0, ItoR, 1-RtoI, 0, 0;
         0, ItoD, 0, 1, 0;
         0, 0, 0, 0, 1];
    tempx0 = x0;
    x0 = A * x0;
    x(:, 1) = x0;
    simulatedDeaths(1) = x0(4);
    tempA = zeros(5,5);
    tempA(2,1) = StoI1;
    tempA(2,3) = RtoI;
    tempx0 = tempA * tempx0;
    simulatedInfections(1) = tempx0(2);

    % from day 2 till the day with higher infection rates simulate the
    % parameters
    for i = 2:round(vaxStart)
        A = [1-StoI1, ItoS, 0, 0, 0;
             StoI1, 1-ItoR-ItoD-ItoS, RtoI, 0, 0;
             0, ItoR, 1-RtoI, 0, 0;
             0, ItoD, 0, 1, 0;
             0, 0, 0, 0, 1];

        % Constructs another coefficient matrix that only includes the
        % changes from other variables to infected
        tempA = zeros(5,5);
        tempA(2,1) = StoI1;
        tempA(2,3) = RtoI;
       
        % if vaccinations have started, add them to the coefficient matrix
        if i >= vaxStart
            A(1,1) = A(1,1) - StoV; A(5,1) = StoV;
            A(3,3) = A(3,3) - RtoV; A(5,3) = RtoV;
        end

        % if breakout infections have started
        % add them to the coefficient matrixies
        if i >= 100
            A(2,5) = VtoI;
            A(5,5) = A(5,5) - VtoI;
            tempA(2,5) = VtoI;
        end

        tempx0 = x0;
        tempx0 = tempA * tempx0;
        simulatedInfections(i) = tempx0(2);

        x0 = A * x0;
        x(:, i) = x0;
        simulatedDeaths(i) = x0(4) + x(4,i-1);     
    end
    
    % Simulates the parameters with different StoI and RtoI
    for i = round(vaxStart)+1:400
        A = [1-StoI2, ItoS, 0, 0, 0;
        StoI2, 1-ItoR-ItoD-ItoS, RtoI, 0, 0;
        0, ItoR, 1-RtoI, 0, 0;
        0, ItoD, 0, 1, 0;
        0, 0, 0, 0, 1];

        % Constructs another coefficient matrix that only includes the
        % changes from other variables to infected
        tempA = zeros(5,5);
        tempA(2,1) = StoI2;
        tempA(2,3) = RtoI;  

        if i >= vaxStart
            A(1,1) = A(1,1) - StoV; A(5,1) = StoV;
            A(3,3) = A(3,3) - RtoV; A(5,3) = RtoV;
        end
        if i >= 100
            A(2,5) = VtoI;
            A(5,5) = A(5,5) - VtoI;
            tempA(2,5) = VtoI;
        end

        % calulates new infections without other changes to the infected
        % amount by using only variables that are "toI" using a tempx0
        % and the tempA created and maintained throughout the function.
        tempx0 = x0;
        tempx0 = tempA * tempx0;
        simulatedInfections(i) = tempx0(2);

        % calculates the new x at day i and stores the cumulative death by
        % adding the current x to the sum of xs before in i-1
        x0 = A * x0;
        x(:, i) = x0;
        simulatedDeaths(i) = x0(4) + x(4,i-1);
    end
end

%% Using fmincon and Results
% uses fmincon to determine the optimal x values within the constrains
optimalParams = fmincon(fun, initialGuess, A, b, Aeq, beq, lb, ub, []);

% displays the results
disp(['Optimal vaccination start day: ', num2str(optimalParams(10))]);
disp(['Optimal breakthrough infection rate (vaxbreak): ', num2str(optimalParams(8))]);
disp(['Optimal Aparams: ', num2str(optimalParams(1:7)')]);
disp(['Optimal Aparams: ', num2str(optimalParams(8:10)')]);

%%% Simulate with optimal parameters for final output

% used for barebones inital hand tuning
% [simulated_infections, simulated_deaths,x] = simulate(initialGuess);

% simulates the optimal parameters over the 400 days
[simulatedInfections, simulatedDeaths,x] = simulate(optimalParams);

vaxpop = x(5,:);
vaxbreak = optimalParams(10);

% Plots the results
figure();
subplot(2,1,1);
plot(1:400, newInfections);
hold on;
plot(1:400, simulatedInfections, 'r--');
xlabel('Days');
ylabel('New Infections');
title('newInfections: Observed vs Simulated');
legend('show','Observed New Infections','Simulated New Infections','location','northeast');
hold off;

subplot(2,1,2);
hold on;
plot(1:400, cumulativeDeaths);
plot(1:400, simulatedDeaths, 'r--');
xlabel('Days');
ylabel('Cumulative Deaths');
title('cumulativeDeaths: Observed vs Simulated');
legend('show','Observed Cumulative Deaths','Simulated Cumulative Deaths','location','southeast');
hold off;


figure();
plot(1:400,x(5,:));
legend('show','vaccinated');
title('vaxpop');
