function fitness = psofunc (x,y)
% single x,y value pair
% compute objective (fitness) function value
% + 7 at the end is really important, b/c the function has a large portion
% that is negative, if you don't raise the whole thing into the positive
% realm then you'll get negative fitness sums and the algorithm won't work
fitness = 3*((1-x)^2)*exp((-x^2)-((y+1)^2)) - 10*((x/5)-(x^3)-(y^5))*exp((-x^2)-(y^2)) - (1/3)*exp((-(x+1)^2)-(y^2));