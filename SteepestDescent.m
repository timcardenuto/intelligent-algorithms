%% Input and Output Membership Functions
clc
clear
x = [0:1:100];

cold = zeros(1,length(x));
cool = zeros(1,length(x));
warm = zeros(1,length(x));
slow = zeros(1,length(x));
medium = zeros(1,length(x));
fast = zeros(1,length(x));
for i=1:length(x)
    % X input MF's
    cold(i) = max([0 min([1 (30-x(i))/(30-0)])]);
    cool(i) = max([0 min([(x(i)-0)/(20-0) (70-x(i))/(70-20)])]);
    warm(i) = max([0 min([(x(i)-40)/(70-40) 1])]);
    
    % Z output MF's
    slow(i) = max([0 min([(x(i)-15)/(20-15) (25-x(i))/(25-20)])]);
    medium(i) = max([0 min([(x(i)-55)/(60-55) (65-x(i))/(65-60)])]);
    fast(i) = max([0 min([(x(i)-75)/(80-75) (85-x(i))/(85-80)])]);
end

subplot(2,1,1);
hold on
plot(x,cold);
plot(x,cool);
plot(x,warm);
xlabel('Universe X - Temperature (F)');
ylabel('MF Grade');
title('Input X Membership Functions');
legend('Cold', 'Cool', 'Warm');

subplot(2,1,2);
hold on
plot(x,slow);
plot(x,medium);
plot(x,fast);
xlabel('Universe Y - Fan Speed (RPM)');
ylabel('MF Grade');
title('Output Y Membership Functions');
legend('Slow', 'Medium', 'Fast')


%% Normal Product Inference Engine (PIE) with CA Defuzzification

% Current Z fuzzy set center points for defuzzification
% this is the parameter that needs to be tuned by the BLS and RLS
fast_c = 80;
medium_c = 60;
slow_c = 20;

% Crisp inputs are used as fuzzy singleton facts
xfacts = [0 10 30 40 50 70];
% Optimal output for given inputs (training data)
y_optimal = [150 128.5 100 100 72.7 50];

% "Fuzzify" the inputs at the given measurement values
% this output known as "firing strength" or "rule weight"
xfacts_cold_MFgrade = zeros(1,length(xfacts));
xfacts_cool_MFgrade = zeros(1,length(xfacts));
xfacts_warm_MFgrade = zeros(1,length(xfacts));
for i=1:length(xfacts)
    xfacts_cold_MFgrade(i) = max([0 min([1 (30-xfacts(i))/(30-0)])]);
    xfacts_cool_MFgrade(i) = max([0 min([(xfacts(i)-0)/(20-0) (70-xfacts(i))/(70-20)])]);
    xfacts_warm_MFgrade(i) = max([0 min([(xfacts(i)-40)/(70-40) 1])]);    
end

% Rule 1, if temp is cold then blower is fast
wR1 = xfacts_cold_MFgrade;      % only 1 input per rule
% R1_fast = wR1' * fast;        % would have needed this if not using CA
wR1_center = wR1 * fast_c;      % instead, using center point of output MF

% Rule 2, if temp is cool then blower is medium
wR2 = xfacts_cool_MFgrade;
% R2_medium = wR2' * medium;
wR2_center = wR2 * medium_c;

% Rule 3, if temp is warm then blower is slow
wR3 = xfacts_warm_MFgrade;
% R3_slow = wR3' * slow;
wR3_center = wR3 * slow_c;

sum_wR_center = wR1_center + wR2_center + wR3_center;
sum_wR = wR1 + wR2 + wR3;     % in this case = 1
% "Defuzzify" ouput using center average 
y_ca = sum_wR_center ./ sum_wR;  

figure()
hold on
plot(xfacts,y_ca)               % rule surface
plot(xfacts,y_optimal, '--k') 
xlabel('Temperature (F)');
ylabel('Fan Speed (RPM)');
title('Fuzzy Temperature Control Output');
legend('Current Profile', 'Optimal Profile')


%% Parallel Steepest Descent (SD)
% Estimate better y output membership function center points using SD

xfacts = [0 10 30 40 50 70];                % input training data
y_optimal = [150 128.5 100 100 72.7 50];    % output training data
thetahat_parallel = [80; 60; 20];                    % initial estimates

% Because I'm using a static step size, I chose an identical cost_delta to
% avoid overstepping. That's what happens here if you choose something like
% cost_delta = 1d-10
step_size = 1d-3;       % denoted as eta (n) or lamda in the notes
cost_delta = 1d-3;      % loop exit criteria
cost_value = 0;         % e(j)
prev_cost_value = -1d2; % e(j-1)

y_ca_defuzz = zeros(1,length(xfacts));  % defuzzified outputs using Center Average method
cost_gradient = zeros(3,length(xfacts));
save_cost_parallel = [];
save_y_defuzz_parallel = [];
% thetahat = zeros(3,length(xfacts));

% Loop while cost (error) is still changing by more than a chosen delta
% The assumption is that the cost difference decreases each iteration
% Iterations are called 'j' in notes, but we don't track them b/c we're
% not indexing anything based on them
while (abs(cost_value - prev_cost_value) > cost_delta)
    
    prev_cost_value = cost_value;   % must update prev_cost_value to converge!
    wR = zeros(3,length(xfacts));   % wR is a matrix of rule weights, rows are each Rule, columns are each input xfact   
    for i=1:length(xfacts)          % num data I/O, parallel SD uses all data at once
    
        xfact_cold_MFgrade = max([0 min([1 (30-xfacts(i))/(30-0)])]);
        xfact_cool_MFgrade = max([0 min([(xfacts(i)-0)/(20-0) (70-xfacts(i))/(70-20)])]);
        xfact_warm_MFgrade = max([0 min([(xfacts(i)-40)/(70-40) 1])]);    
        wR(:,i) = [xfact_cold_MFgrade xfact_cool_MFgrade xfact_warm_MFgrade]';
        
        % Calculate cost update, notes page 34
        % the cost function actually uses the "Defuzzified" ouput 
        y_ca_defuzz(i) = sum(wR(:,i) .* thetahat_parallel) ./ sum(wR(:,i));
        %cost_value = .5 * ((y_ca_defuzz(i) - y_optimal(i))^2);
        cost_value = y_ca_defuzz(i) - y_optimal(i);
        
        % Calculate gradiant, notes page 37
        % this is partial derivative of cost function
        cost_gradient(:,i) = cost_value * wR(:,i) / sum(wR(:,i));
        
        % Calculate new y output MF center
        % thetanew = theta - n*g
        % Must hard code which direction to go, positive/negative
        % if you choose wrong, the algorithm won't converge
        thetahat_parallel = thetahat_parallel - step_size * cost_gradient(:,i);
    end
    save_y_defuzz_parallel = [save_y_defuzz_parallel y_ca_defuzz'];
    save_cost_parallel = [save_cost_parallel cost_value];
end

thetahat_parallel    % the new estimated output centers
y_parallel = y_ca_defuzz


%% Serial Steepest Descent (SD)
% Estimate better y output membership function center points using SD

xfacts = [0 10 30 40 50 70];                % input training data
y_optimal = [150 128.5 100 100 72.7 50];    % output training data
thetahat_serial = [80; 60; 20];                    % initial estimates

step_size = 1d-2;       % denoted as eta (n) or lamda in the notes
cost_delta = 1d-10;     % loop exit criteria
cost_value = zeros(1,length(xfacts));       % e(j)
prev_cost_value = 1d5;                      % e(j-1)

wR = [];    % wR is a matrix of rule weights, rows are each Rule, columns are each input xfact
y_ca_defuzz = zeros(1,length(xfacts));  % defuzzified outputs using Center Average method
cost_gradient = zeros(3,length(xfacts));
% thetahat = zeros(3,length(xfacts));

save_costxy = [];
save_y_defuzz_all = [];
for i=1:length(xfacts)  % num data I/O, this is serial SD training loop
    
    xfact_cold_MFgrade = max([0 min([1 (30-xfacts(i))/(30-0)])]);
    xfact_cool_MFgrade = max([0 min([(xfacts(i)-0)/(20-0) (70-xfacts(i))/(70-20)])]);
    xfact_warm_MFgrade = max([0 min([(xfacts(i)-40)/(70-40) 1])]);    
    wR = [wR [xfact_cold_MFgrade xfact_cool_MFgrade xfact_warm_MFgrade]'];
    
    % Loop while cost (error) is still changing by more than a chosen delta
    % The assumption is that the cost difference decreases each iteration
    % Iterations are called 'j' in notes, but we don't track them b/c we're
    % not indexing anything based on them
    cost_value(i) = 0;
    prev_cost_value = -1d2;
    save_costx = [];
    save_y_defuzz = [];
    j = 0;
    while (abs(cost_value(i) - prev_cost_value) > cost_delta && j < 1d3)
        
        prev_cost_value = cost_value(i);    % must update prev_cost_value to converge!

        % Calculate cost update, notes page 34
        % the cost function actually uses the "Defuzzified" ouput 
        y_ca_defuzz(i) = sum(wR(:,i) .* thetahat_serial) ./ sum(wR(:,i));
        %cost_value(i) = .5 * ((y_ca_defuzz(i) - y_optimal(i))^2);
        cost_value(i) = y_ca_defuzz(i) - y_optimal(i);
        
        % Calculate gradiant, notes page 37
        % this is partial derivative of cost function
        cost_gradient(:,i) = cost_value(i) * wR(:,i) / sum(wR(:,i));
        
        % Calculate new y output MF center
        % thetanew = theta - n*g
        % Must hard code which direction to go, positive/negative
        % if you choose wrong, the algorithm won't converge
        thetahat_serial = thetahat_serial - step_size * cost_gradient(:,i);  
        
        save_costx = [save_costx cost_value(i)];
        save_y_defuzz = [save_y_defuzz y_ca_defuzz(i)];
        j = j + 1;  % iteration limit counter
    end
    save_costxy = [save_costxy; save_costx];
    save_y_defuzz_all = [save_y_defuzz_all; save_y_defuzz];
end

thetahat_serial    % the new estimated output centers
y_serial = y_ca_defuzz

%% Data Plots

figure()
hold on
scatter(xfacts,y_optimal, 'k') 
plot(xfacts,y_serial, 'r')
plot(xfacts,y_parallel, 'b')
xlabel('Temperature (F)');
ylabel('Fan Speed (RPM)');
title('Fuzzy Temperature Control Output');
legend('Optimal Profile', 'Serial Steepest Descent Tuned Profile', 'Parallel Steepest Descent Tuned Profile')

figure()
hold on
for i=1:length(xfacts)
    plot(save_costxy(i,:))
end
xlabel('Iteration');
ylabel('Cost');
title('Serial SD Cost VS Iteration');
legend('xfact 1', 'xfact 2', 'xfact 3', 'xfact 4', 'xfact 5', 'xfact 6')

figure()
hold on
for i=1:length(xfacts)
    plot(save_y_defuzz_all(i,:))
end
xlabel('Iteration');
ylabel('y Output');
title('Serial SD y Output VS Iteration');
legend('xfact 1', 'xfact 2', 'xfact 3', 'xfact 4', 'xfact 5', 'xfact 6')

figure()
hold on
plot(save_cost_parallel)
xlabel('Iteration');
ylabel('Cost');
title('Parallel SD Cost VS Iteration');
legend('all xfacts')

figure()
hold on
for i=1:length(xfacts)
    plot(save_y_defuzz_parallel(i,:))
end
xlabel('Iteration');
ylabel('y Output');
title('Parallel SD y Output VS Iteration');
legend('xfact 1', 'xfact 2', 'xfact 3', 'xfact 4', 'xfact 5', 'xfact 6')



%% try using new thetahats from serial SD - shows problem with serial

% Current Z fuzzy set center points for defuzzification
% this is the parameter that needs to be tuned by the BLS and RLS
fast_c = thetahat_serial(1);
medium_c = thetahat_serial(2);
slow_c = thetahat_serial(3);

% Crisp inputs are used as fuzzy singleton facts
xfacts = [0 10 30 40 50 70];
% Optimal output for given inputs (training data)
y_optimal = [150 128.5 100 100 72.7 50];

% "Fuzzify" the inputs at the given measurement values
% this output known as "firing strength" or "rule weight"
xfacts_cold_MFgrade = zeros(1,length(xfacts));
xfacts_cool_MFgrade = zeros(1,length(xfacts));
xfacts_warm_MFgrade = zeros(1,length(xfacts));
for i=1:length(xfacts)
    xfacts_cold_MFgrade(i) = max([0 min([1 (30-xfacts(i))/(30-0)])]);
    xfacts_cool_MFgrade(i) = max([0 min([(xfacts(i)-0)/(20-0) (70-xfacts(i))/(70-20)])]);
    xfacts_warm_MFgrade(i) = max([0 min([(xfacts(i)-40)/(70-40) 1])]);    
end

% Rule 1, if temp is cold then blower is fast
wR1 = xfacts_cold_MFgrade;      % only 1 input per rule
wR1_center = wR1 * fast_c;      % instead, using center point of output MF

% Rule 2, if temp is cool then blower is medium
wR2 = xfacts_cool_MFgrade;
wR2_center = wR2 * medium_c;

% Rule 3, if temp is warm then blower is slow
wR3 = xfacts_warm_MFgrade;
wR3_center = wR3 * slow_c;

sum_wR_center = wR1_center + wR2_center + wR3_center;
sum_wR = wR1 + wR2 + wR3;           % in this case = 1
y_ca = sum_wR_center ./ sum_wR;  % "Defuzzify" ouput using center average 

figure()
hold on
plot(xfacts,y_ca)               % rule surface
plot(xfacts,y_optimal, '--k') 
xlabel('Temperature (F)');
ylabel('Fan Speed (RPM)');
title('Fuzzy Temperature Control Output');
legend('Serial Tuned Profile', 'Optimal Profile')
