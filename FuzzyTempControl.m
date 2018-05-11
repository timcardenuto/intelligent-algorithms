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

% break;

%% Product Inference Engine (PIE) with CA Defuzzification

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

% break;

% put all rule weights into matrix so next section is more generic
wR = [wR1; wR2; wR3];


%% Batched Least Squares (BLS)
% assumes 1 input, 1 output

% Estimate better y output membership function center points using the
% training data 'y_optimal'
y_optimal = [150 128.5 100 100 72.7 50];

% matrix size = num rules (rows) x num data I/O (columns)
num_rules = length(wR(:,1));
xi = zeros(1,num_rules);
phi = [];

for i=1:length(xfacts)
    for j=1:num_rules
        xi(j) = wR(j,i)/sum(wR(:,i));
    end
    phi = [phi xi'];
end
phi = phi';

% new estimate for y output membership function centers using BLS
thetahat = inv(phi' * phi) * phi' * y_optimal'


%% Try again with optimal output y center points from BLS

% New fuzzy set center points for defuzzification, learned from BLS
fast_center = thetahat(1)
medium_center = thetahat(2)
slow_center = thetahat(3)

% re-using the same crisp inputs 'xfacts' as fuzzy singleton facts

% re-using the same "Fuzzified" inputs and same "rule weights"

% "Defuzzify" ouput using center average 
sum_wR_center = zeros(1,length(xfacts));
sum_wR = zeros(1,length(xfacts));
for i=1:num_rules
    sum_wR_center = sum_wR_center + wR(i,:) * thetahat(i);
    sum_wR = sum_wR + wR(i,:);
end
y_bls = sum_wR_center ./ sum_wR


%% Recursive Least Squares (RLS)

% initial estimates
% thetahat = zeros(num_rules,length(xfacts)+1);
% thetahat(:,1) = [80 60 20];
thetahat = [80; 60; 20];
P = 100*eye(num_rules);

xi = zeros(1,num_rules)';
y_rls = [];
wR = [];
for i=1:length(xfacts)  % num data I/O, this is RLS loop
    
    xfact_cold_MFgrade = max([0 min([1 (30-xfacts(i))/(30-0)])]);
    xfact_cool_MFgrade = max([0 min([(xfacts(i)-0)/(20-0) (70-xfacts(i))/(70-20)])]);
    xfact_warm_MFgrade = max([0 min([(xfacts(i)-40)/(70-40) 1])]);    
    wR = [wR [xfact_cold_MFgrade xfact_cool_MFgrade xfact_warm_MFgrade]'];
    
    % Regression vector update
    for j=1:num_rules
        xi(j) = wR(j,i)/sum(wR(:,i));
    end
        
    % Error matrix update
    P = P - (P*xi*xi'*P)/(1 + xi'*P*xi);
    
    % Parameter update - using training data
    thetahat = [thetahat (thetahat(:,i) + P*xi*(y_optimal(i) - xi'*thetahat(:,i)))];

    % "Defuzzify" ouput using center average 
    sum_wR_center = zeros(1,length(wR(1,:)));
    sum_wR = zeros(1,length(wR(1,:)));
    for j=1:num_rules
        sum_wR_center = sum_wR_center + wR(j,:) * thetahat(j,i+1);
        sum_wR = sum_wR + wR(j,:);
    end
    y_temp = sum_wR_center ./ sum_wR;
    y_rls = [y_rls y_temp(length(y_temp))]; 

end

figure()
hold on
plot(xfacts,y_bls, 'b')
plot(xfacts,y_rls, 'r')
plot(xfacts,y_optimal, '--k') 
xlabel('Temperature (F)');
ylabel('Fan Speed (RPM)');
title('Fuzzy Temperature Control Output');
legend('BLS Tuned Profile', 'RLS Tuned Profile', 'Optimal Profile')

break;

%% If you run the RLS for N number of input measurements to recalculate the
% output MF center points, you can batch Defuzzify the output...

% new estimate for y output membership function centers from RLS
thetahat = thetahat(:,length(thetahat));

% New fuzzy set center points for defuzzification, learned from BLS
fast_center = thetahat(1)
medium_center = thetahat(2)
slow_center = thetahat(3)

% re-using the same crisp inputs 'xfacts' as fuzzy singleton facts

% re-using the same "Fuzzified" inputs and same "rule weights"

% "Defuzzify" ouput using center average 
sum_wR_center = zeros(1,length(xfacts));
sum_wR = zeros(1,length(xfacts));
for i=1:num_rules
    sum_wR_center = sum_wR_center + wR(i,:) * thetahat(i);
    sum_wR = sum_wR + wR(i,:);
end
y_rls = sum_wR_center ./ sum_wR;  

figure()
hold on
plot(xfacts,y_bls, 'b')
plot(xfacts,y_rls, 'r')
plot(xfacts,y_optimal, '--k') 
xlabel('Temperature (F)');
ylabel('Fan Speed (RPM)');
title('Fuzzy Temperature Control Output');
legend('BLS Tuned Profile', 'RLS Tuned Profile', 'Optimal Profile')

