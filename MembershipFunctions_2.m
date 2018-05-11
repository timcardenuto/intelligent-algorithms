x = [0:1:100]; %age

uyoung = exp(-(x/20).^2);
uold = exp(-((x-100)/30).^2);

hold on
plot(x,uyoung)
plot(x,uold)

% not very young and not very old
a = min((1-(uyoung.^2)),(1-(uold.^2)));
plot(x,a)

% very young or very old
b = max((uyoung.^2),(uold.^2));
plot(x,b)

xlabel('X - Universe of Discourse (Age)');
ylabel('Membership Grade');
title('Non-primary Calculations for Age Groups');
legend('young', 'old', 'not very young and not very old', 'very young or very old')
