x = [0:1:100];

mf = trapmf(x,[10 30 50 90]);

hold on
grid on
axis([ 0 100 0 1.25])
plot(x,mf)
xlabel('X = Universe of Discourse');
ylabel('Membership Grade');
title('Trapazoid Membership Function and Defuzzification');

% Center of Area (CoA) Defuzzification
% z_coa = defuzz(x,mf,'centroid') % using Fuzzy Logic Toolbox
z_coa = trapz(mf.*x)/trapz(mf); % trapz() 'estimates' integral of a trapazoid
stem(z_coa,1,'b')

% Mean of Max (MoM) Defuzzification
% z_mom = defuzz(x,mf,'mom') % using Fuzzy Logic Toolbox
index = find(mf==max(mf));  % all index(s) of mf that equal the max value
mean_index = mean(index);  % assuming the max area is contiguous...
z_mom = x(mean_index);     % ...can use index with input vector to find z
stem(z_mom,1,'r')

% Bisector of Area (BoA) Defuzzification
% z_boa = defuzz(x,mf,'bisector') % using Fuzzy Logic Toolbox
index=1;
while (logical(trapz(mf(1:index))~=trapz(mf(index:length(mf)))))
    index=index+1;
end
z_boa = x(index);
stem(z_boa,1,'g')

legend('Overall Output Membership Function', 'Center of Area (CoA) Defuzzification', 'Mean of Max (MoM) Defuzzification', 'Bisector of Area (BoA) Defuzzification')
