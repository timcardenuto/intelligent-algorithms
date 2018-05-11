x = [0:1:20];
y = x;

% triangle MF
ua1 = trimf(x,[5 10 15]);  
ub1 = trimf(y,[5 10 15]);  

% trapezoid MF
ua2 = trapmf(x,[3 8 12 17]);
ub2 = trapmf(y,[3 8 12 17]);

% Zadeh's arithmetic rule
arithmetic1 = zeros(length(x),length(y));
for i=1:length(x),
    for j=1:length(y),
       arithmetic1(j,i) = min(1,(1-ua1(i)+ub1(j)));
    end
end
subplot(2,4,1);
mesh(x,y,arithmetic1);
xlabel('A');
ylabel('B');
title('Zadeh arithmetic rule - triangle');

arithmetic2 = zeros(length(x),length(y));
for i=1:length(x),
    for j=1:length(y),
       arithmetic2(j,i) = min(1,(1-ua2(i)+ub2(j)));
    end
end
subplot(2,4,5);
mesh(x,y,arithmetic2);
xlabel('A');
ylabel('B');
title('Zadeh arithmetic rule - trapezoid');

% Zadeh min-max rule
minmax1 = zeros(length(x),length(y));
for i=1:length(x),
    for j=1:length(y),
       minmax1(j,i) = max((1-ua1(i)),min(ua1(i),ub1(j)));
    end
end
subplot(2,4,2);
mesh(x,y,minmax1);
xlabel('A');
ylabel('B');
title('Zadeh min-max rule - triangle');

minmax2 = zeros(length(x),length(y));
for i=1:length(x),
    for j=1:length(y),
       minmax2(j,i) = max((1-ua2(i)),min(ua2(i),ub2(j)));
    end
end
subplot(2,4,6);
mesh(x,y,minmax2);
xlabel('A');
ylabel('B');
title('Zadeh min-max rule - trapezoid');

% Boolean implication
boolean1 = zeros(length(x),length(y));
for i=1:length(x),
    for j=1:length(y),
        boolean1(j,i) = max((1-ua1(i)),ub1(j));
    end
end
subplot(2,4,3);
mesh(x,y,boolean1);
xlabel('A');
ylabel('B');
title('Boolean implication - triangle');

boolean2 = zeros(length(x),length(y));
for i=1:length(x),
    for j=1:length(y),
        boolean2(j,i) = max((1-ua2(i)),ub2(j));
    end
end
subplot(2,4,7);
mesh(x,y,boolean2);
xlabel('A');
ylabel('B');
title('Boolean implication - trapezoid');

% Goguen's implication
goguen1 = zeros(length(x),length(y));
for i=1:length(x),
    for j=1:length(y),
        if(ua1(i)<=ub1(j))
            goguen1(j,i) = 1;
       else
            goguen1(j,i) = ub1(j)/ua1(i);
        end
    end
end
subplot(2,4,4);
mesh(x,y,goguen1);
xlabel('A');
ylabel('B');
title('Goguens implication - triangle');

goguen2 = zeros(length(x),length(y));
for i=1:length(x),
    for j=1:length(y),
        if(ua2(i)<=ub2(j))
            goguen2(j,i) = 1;
       else
            goguen2(j,i) = ub2(j)/ua2(i);
        end
    end
end
subplot(2,4,8);
mesh(x,y,goguen2);
xlabel('A');
ylabel('B');
title('Goguens implication - trapezoid');