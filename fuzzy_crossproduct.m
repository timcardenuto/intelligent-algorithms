x=-20:1:20;
y=-20:1:20;

ux=exp((-x.^2)/30);
uy=exp((-y.^2)/90);

ux_cross_uy = zeros(length(x),length(y));

for i=1:length(x),
    for j=1:length(y),
       ux_cross_uy(i,j) = min([ux(i),uy(j)]);
    end
end

%hold on;
%plot(ux);
%plot(uy);
%break;
mesh(x,y,ux_cross_uy);
xlabel('x');
ylabel('y');
zlabel('Cartesian Product');