c1 = [6];
sigma1 = [6];
c2 = [-6];
sigma2 = [4];

x = [-20:1:20];

ua = exp(-0.5 * (((x-c1)/sigma1).^2));
ub = exp(-0.5 * (((x-c2)/sigma2).^2));

ua_cross_ub = zeros(length(x),length(x));

for i=1:length(x),
    for j=1:length(x),
       ua_cross_ub(i,j) = min([ua(i),ub(j)]);
    end
end

% hold on
% plot(x,ua)
% plot(x,ub)
mesh(x,x,ua_cross_ub);
xlabel('A Fuzzy Set');
ylabel('B Fuzzy Set');
zlabel('Cartesian Product AxB');
title('Tnorm min, Cartesian Product');
