c1 = [2,4,6];
c2 = [3,7,13];
sigma1 = [2,4,6];
sigma2 = [4,8,10];

x = [-10:.1:50];
mf = [];

for i = 1:length(c1)
    mf = [mf ts_gaussian(x,c1(i),sigma1(i),c2(i),sigma2(i))];
end

plot(x,mf)
xlabel('X = Universe of Discourse');
ylabel('Membership Grade');
title('2-sided Gaussian Membership Function');
legend('c1=2, c2=3, sigma1=2, sigma2=4', 'c1=4, c2=7, sigma1=4, sigma2=8', 'c1=6, c2=13, sigma1=6, sigma2=10')
%%
% you can find the index for the *closest* value to the crossover point .5, but it will rarely be precise
%a = abs(mf(:,1)-.5);
%[val idx] = min(a);
%mf(idx,1)

% mathematically...
% TECHNICALLY would need to only use this version of x1 (with -1* at the front) 
% where x is negative, then switch to what x2 looks like.
x1 = sqrt(log(.5)*(-2)) * sigma1 + c1;
if (x1 > c1)
    x1 = (-1)*sqrt(log(.5)*(-2)) * sigma1 + c1
end
x2 = sqrt(log(.5)*(-2)) * sigma2 + c2
bandwidth = x2 - x1