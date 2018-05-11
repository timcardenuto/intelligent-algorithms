x = [0:1:20];
y = x;

% triangle MF
ua1 = trimf(x,[5 10 15]);  
ub1 = trimf(y,[5 10 15]);  

% trapezoid MF
ua2 = trapmf(x,[3 8 12 17]);
ub2 = trapmf(y,[3 8 12 17]);

% min (mamdani)
mamdani1 = zeros(length(x),length(y));
for i=1:length(x),
    for j=1:length(y),
       mamdani1(i,j) = min([ua1(i),ub1(j)]);
    end
end
subplot(2,4,1);
mesh(x,y,mamdani1);
xlabel('A');
ylabel('B');
title('triangle min');

mamdani2 = zeros(length(x),length(y));
for i=1:length(x),
    for j=1:length(y),
       mamdani2(i,j) = min([ua2(i),ub2(j)]);
    end
end
subplot(2,4,5);
mesh(x,y,mamdani2);
xlabel('A');
ylabel('B');
title('trapezoid min');

% algebraic product
algebraic1 = zeros(length(x),length(y));
for i=1:length(x),
    for j=1:length(y),
       algebraic1(i,j) = ua1(i)*ub1(j);
    end
end
subplot(2,4,2);
mesh(x,y,algebraic1);
xlabel('A');
ylabel('B');
title('triangle algebraic product');

algebraic2 = zeros(length(x),length(y));
for i=1:length(x),
    for j=1:length(y),
       algebraic2(i,j) = ua2(i)*ub2(j);
    end
end
subplot(2,4,6);
mesh(x,y,algebraic2);
xlabel('A');
ylabel('B');
title('trapezoid algebraic product');

% bounded product
bounded1 = zeros(length(x),length(y));
for i=1:length(x),
    for j=1:length(y),
        bounded1(i,j) = max(0,(ua1(i)+ub1(j)-1));
    end
end
subplot(2,4,3);
mesh(x,y,bounded1);
xlabel('A');
ylabel('B');
title('triangle bounded product');

bounded2 = zeros(length(x),length(y));
for i=1:length(x),
    for j=1:length(y),
        bounded2(i,j) = max(0,(ua2(i)+ub2(j)-1));
    end
end
subplot(2,4,7);
mesh(x,y,bounded2);
xlabel('A');
ylabel('B');
title('trapezoid bounded product');

% drastic product
drastic1 = zeros(length(x),length(y));
for i=1:length(x),
    for j=1:length(y),
        if(ub1(j)==1)
            drastic1(i,j) = ua1(i);
        elseif(ua1(i)==1)
            drastic1(i,j) = ub1(j);
        else
            drastic1(i,j) = 0;
        end
    end
end
subplot(2,4,4);
mesh(x,y,drastic1);
xlabel('A');
ylabel('B');
title('triangle drastic product');

drastic2 = zeros(length(x),length(y));
for i=1:length(x),
    for j=1:length(y),
        if(ub2(j)==1)
            drastic2(i,j) = ua2(i);
        elseif(ua2(i)==1)
            drastic2(i,j) = ub2(j);
        else
            drastic2(i,j) = 0;
        end
    end
end
subplot(2,4,8);
mesh(x,y,drastic2);
xlabel('A');
ylabel('B');
title('trapezoid drastic product');
