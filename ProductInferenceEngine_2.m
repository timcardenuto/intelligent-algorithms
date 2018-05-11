x = [-5:.1:5];
y = x;

% X input MF's
xsmall = 1./(1+exp(-(-3)*(x-0)));
xlarge = 1./(1+exp(-(3)*(x-0)));

subplot(3,1,1);
hold on
plot(x,xsmall);
plot(x,xlarge);
xlabel('Universe X');
ylabel('MF Grade');
title('Input X Membership Functions');
legend('Small', 'Large')

% Y input MF's
ysmall = 1./(1+exp(-(-0.9)*(y-0)));
ylarge = 1./(1+exp(-(0.9)*(y-0)));

subplot(3,1,2);
hold on
plot(y,ysmall);
plot(y,ylarge);
xlabel('Universe Y');
ylabel('MF Grade');
title('Input Y Membership Functions');
legend('Small', 'Large')

% Z output MF's
zneglarge = 1./(1+exp(-(-4)*(x-(-3.5))));
znegsmall = exp(-0.5*(((x-(-1.7))/1.5).^2));
zpossmall = exp(-0.5*(((x-(1.7))/1.5).^2));
zposlarge = 1./(1+exp(-(4)*(x-(3.5))));

subplot(3,1,3);
hold on
plot(x,zneglarge);
plot(x,znegsmall);
plot(x,zpossmall);
plot(x,zposlarge);
stem(-4,1,'b');         % center points, in same order
stem(-1.7,1,'b');
stem(1.7,1,'b');
stem(4,1,'b');
xlabel('Universe Z');
ylabel('MF Grade');
title('Output Z Membership Functions');
legend('NegLarge', 'NegSmall', 'PosSmall', 'PosLarge', 'Center Points')

%% Product Inference Engine (PIE) with CA Defuzzification

% Determine Z fuzzy set center points for defuzzification
% sigmoids are arbitrary, based on class slides
zneglarge_center = -4;  
zposlarge_center = 4;          
% gaussian center = crossover point z2 - z1 / 2
% determined from plotting and measuring, tops seem to be at the 'c' value
znegsmall_center = -1.7;
zpossmall_center = 1.7;

% Crisp inputs are used as fuzzy singleton facts
xfacts = [-5:.1:5];
yfacts = [-5:.1:5];

% Evaluate X and Y membership functions at these fact values
% this output known as "firing strength" ?
xfacts_xsmall_MFgrade = 1./(1+exp(-(-3)*(xfacts-0)));
xfacts_xlarge_MFgrade = 1./(1+exp(-(3)*(xfacts-0)));
yfacts_ysmall_MFgrade = 1./(1+exp(-(-0.9)*(yfacts-0)));
yfacts_ylarge_MFgrade = 1./(1+exp(-(0.9)*(yfacts-0)));

z_ca = [];

for i=1:length(yfacts),

    % Rule 1, if x is small and y is small then z is neglarge
    wR1 = xfacts_xsmall_MFgrade * yfacts_ysmall_MFgrade(i);    % this is the "rule weight" from all inputs (antecedents)?
    % R1_zneglarge = wR1' * zneglarge; % would have needed this if not using CA
    wR1_center = wR1 * zneglarge_center;

    % Rule 2, if x is small and y is large then z is negsmall
    wR2 = xfacts_xsmall_MFgrade * yfacts_ylarge_MFgrade(i);
    % R2_znegsmall = wR2' * znegsmall;
    wR2_center = wR2 * znegsmall_center;

    % Rule 3, if x is large and y is small then z is possmall
    wR3 = xfacts_xlarge_MFgrade * yfacts_ysmall_MFgrade(i);
    % R3_zpossmall = wR3' * zpossmall;
    wR3_center = wR3 * zpossmall_center;

    % Rule 4, if x is large and y is large then z is poslarge
    wR4 = xfacts_xlarge_MFgrade * yfacts_ylarge_MFgrade(i);
    % R4_zposlarge = wR4' * zposlarge;
    wR4_center = wR4 * zposlarge_center;
    
    sum_wR_center = wR1_center + wR2_center + wR3_center + wR4_center;
    sum_wR = wR1 + wR2 + wR3 + wR4;     % in this case = 1
    % ouput center average vector for 1 value of y
    z_ca = [z_ca; (sum_wR_center ./ sum_wR)];  
end

figure()
mesh(x,y,z_ca)  % rule surface
xlabel('Universe X');
ylabel('Universe Y');
zlabel('Universe Z');
title('Rule surface for simplified PIE with CA defuzzification');

break;


%% Plot everything

for i=1:length(xfacts),
    % Plot Rule 1
    figure();
    subplot(4,3,1);
    hold on
    plot(x,xsmall);
    stem(xfacts(i),xfacts_xsmall_MFgrade(i)); % Small MF grade for x 
    xlabel('Universe X');
    ylabel('MF Grade');
    title('if x is small');
    legend('Small', sprintf('x = %i', xfacts(i)))

    subplot(4,3,2);
    hold on
    plot(y,ysmall);
    stem(yfacts(i),yfacts_ysmall_MFgrade(i)); % Small MF grade for y 
    xlabel('Universe Y');
    ylabel('MF Grade');
    title('and y is small');
    legend('Small', sprintf('y = %i', yfacts(i)))

    subplot(4,3,3);
    hold on
    plot(x,zneglarge);
    area(x,R1_zneglarge(i,:));
    temp = ones(length(x)) * wR1(i);
    plot(x,temp, '--');    % Product of x and y MF grades = rule weight? 
    xlabel('Universe Z');
    ylabel('MF Grade');
    title('then z is neg-large');
    legend('NegLarge MF', 'R1 NegLarge Output', 'R1 weight')

    % Plot Rule 2
    subplot(4,3,4);
    hold on
    plot(x,xsmall);
    stem(xfacts(i),xfacts_xsmall_MFgrade(i)); % Small MF grade for x 
    xlabel('Universe X');
    ylabel('MF Grade');
    title('if x is small');
    legend('Small', sprintf('x = %i', xfacts(i)))

    subplot(4,3,5);
    hold on
    plot(y,ylarge);
    stem(yfacts(i),yfacts_ylarge_MFgrade(i)); % Small MF grade for y 
    xlabel('Universe Y');
    ylabel('MF Grade');
    title('and y is large');
    legend('Large', sprintf('y = %i', yfacts(i)))

    subplot(4,3,6);
    hold on
    plot(x,znegsmall);
    area(x,R2_znegsmall(i,:));
    temp = ones(length(x)) * wR2(i);
    plot(x,temp, '--');    % Product of x and y MF grades = rule weight? 
    xlabel('Universe Z');
    ylabel('MF Grade');
    title('then z is neg-small');
    legend('NegLarge MF', 'R2 NegLarge Output', 'R2 weight')

    % Plot Rule 3
    subplot(4,3,7);
    hold on
    plot(x,xlarge);
    stem(xfacts(i),xfacts_xlarge_MFgrade(i)); % Large MF grade for x 
    xlabel('Universe X');
    ylabel('MF Grade');
    title('if x is large');
    legend('Large', sprintf('x = %i', xfacts(i)))

    subplot(4,3,8);
    hold on
    plot(y,ysmall);
    stem(yfacts(i),yfacts_ysmall_MFgrade(i)); % Small MF grade for y 
    xlabel('Universe Y');
    ylabel('MF Grade');
    title('and y is small');
    legend('Small', sprintf('y = %i', yfacts(i)))

    subplot(4,3,9);
    hold on
    plot(x,zpossmall);
    area(x,R3_zpossmall(i,:));
    temp = ones(length(x)) * wR3(i);
    plot(x,temp, '--');    % Product of x and y MF grades = rule weight? 
    xlabel('Universe Z');
    ylabel('MF Grade');
    title('then z is pos-small');
    legend('pos-small MF', 'R3 pos-small output', 'R3 weight')

    % Plot Rule 4
    subplot(4,3,10);
    hold on
    plot(x,xlarge);
    stem(xfacts(i),xfacts_xlarge_MFgrade(i)); % Large MF grade for x 
    xlabel('Universe X');
    ylabel('MF Grade');
    title('if x is large');
    legend('Large', sprintf('x = %i', xfacts(i)))

    subplot(4,3,11);
    hold on
    plot(y,ylarge);
    stem(yfacts(i),yfacts_ylarge_MFgrade(i)); % Large MF grade for y 
    xlabel('Universe Y');
    ylabel('MF Grade');
    title('and y is large');
    legend('Large', sprintf('y = %i', yfacts(i)))

    subplot(4,3,12);
    hold on
    plot(x,zposlarge);
    area(x,R4_zposlarge(i,:));
    temp = ones(length(x)) * wR4(i);
    plot(x,temp, '--');    % Product of x and y MF grades = rule weight? 
    xlabel('Universe Z');
    ylabel('MF Grade');
    title('then z is pos-large');
    legend('pos-large MF', 'R4 pos-large output', 'R4 weight')
end

%%

% single point x,y example
figure()
hold on
area(x,R1_zneglarge(3,:));
temp = ones(length(x)) * wR1(3);
plot(x,temp, 'b--');
stem(zneglarge_center,wR1(3),'b');
area(x,R2_znegsmall(3,:));
temp = ones(length(x)) * wR2(3);
plot(x,temp, 'b--');
stem(znegsmall_center,wR2(3),'b');
area(x,R3_zpossmall(3,:));
temp = ones(length(x)) * wR3(3);
plot(x,temp, 'b--');
stem(zpossmall_center,wR3(3),'b');  
area(x,R4_zposlarge(3,:));
temp = ones(length(x)) * wR4(3);
plot(x,temp, 'b--');
stem(zposlarge_center,wR4(3),'b');

