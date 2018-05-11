clc
clear

rand('state',0);    % edit this to get different results every time

x=-5:0.313:5;       % the resolution here doesn't seem to make a difference
y=-5:0.313:5;       %  used 32 values to match the resolution of 'peaks'
z = zeros(2,length(x));
for i=1:length(x)
    for j=1:length(y)
        z(j,i) = (3*((1-x(i))^2)*exp((-x(i)^2)-((y(j)+1)^2)) - 10*((x(i)/5)-(x(i)^3)-(y(j)^5))*exp((-x(i)^2)-(y(j)^2)) - (1/3)*exp((-(x(i)+1)^2)-(y(j)^2)));
    end
end


% figure()
% hold on

%% Particle Swarm Optimization Algorithm

num_particles = [ 20 20 30 30 ];
fitness_delta = 1d-3;
phi1 = [1 2 1 2];
phi2 = [3 2.5 3 2.5];
vmax = [0.6 1 0.6 1];

for rr=1:4 % each simulation parameter setup

particle = [];
  
% initialize all particles
gbestfitness = 0;
for i=1:num_particles(rr);
    particle = [ particle struct('x',(rand()*10-5),'y',(rand()*10-5),'z',(rand()*10-5),'velx',0,'vely',0,'velz',0,'fitness',0,'pbestfitness',0,'pbestx',0,'pbesty',0,'pbestz',0) ];
    particle(i).fitness = 8 - psofunc(particle(i).x,particle(i).y);
    particle(i).pbestfitness = particle(i).fitness;
    particle(i).pbestx = particle(i).x;
    particle(i).pbesty = particle(i).y;
    particle(i).pbestz = particle(i).z;
    if (particle(i).fitness > gbestfitness)
        gbestfitness =  particle(i).pbestfitness;
        gbestx = particle(i).pbestx;
        gbesty = particle(i).pbesty;
        gbestz = particle(i).pbestz;
    end
end

% run optimization loop
last_gbestfitness = 1d5;
icnt = 1;
wmax = 1;
wmin = 0.4;
wsf = 0.999;
maxiter = 1d2;
w = wmax:-(wmax-wmin)/maxiter:wmin;
max_fitness = [];
average_fitness = [];
min_fitness = [];

figh = figure();
mesh(x,y,z);
hold on
xlabel('x');
ylabel('y');
zlabel('z');
title(['Fitness Function results for simulation ', num2str(rr)]);  

while ( abs(gbestfitness - last_gbestfitness) > fitness_delta )
    
    R = [rand() rand()];    % only updated in between runs, same random number applys to all particles?
    
    % re-draw every update
    clf(figh,'reset')
    figure(figh)
    mesh(x,y,z);
    hold on
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(['Fitness Function results for simulation ', num2str(rr)]);  

    fitness = [];
    for i=1:num_particles(rr) 
        
        % Update particle velocity
        particle(i).velx = w(icnt)*particle(i).velx + phi1(rr)*R(1)*(particle(i).pbestx - particle(i).x) +  phi2(rr)*R(2)*(gbestx - particle(i).x);
        particle(i).vely = w(icnt)*particle(i).vely + phi1(rr)*R(1)*(particle(i).pbesty - particle(i).y) +  phi2(rr)*R(2)*(gbesty - particle(i).y);
        particle(i).velz = w(icnt)*particle(i).velz + phi1(rr)*R(1)*(particle(i).pbestz - particle(i).z) +  phi2(rr)*R(2)*(gbestz - particle(i).z);

%         disp(['velocity update ', num2str(particle(i).velx), ' ', num2str(particle(i).vely)])
        
        % Constrain velocity update
        if (abs(particle(i).velx) > vmax(rr))
            particle(i).velx = vmax(rr)*(particle(i).velx/abs(particle(i).velx));    % etra term just forces the sign to be correct
        end
        if (abs(particle(i).vely) > vmax(rr))
            particle(i).vely = vmax(rr)*(particle(i).vely/abs(particle(i).vely));
        end
        if (abs(particle(i).velz) > vmax(rr))
            particle(i).velz = vmax(rr)*(particle(i).velz/abs(particle(i).velz));
        end
        
%         disp(['after constrain ', num2str(particle(i).velx), ' ', num2str(particle(i).vely)])
        
        % Update particle positions
        particle(i).x = particle(i).x + particle(i).velx;
        particle(i).y = particle(i).y + particle(i).vely;
        particle(i).z = particle(i).z + particle(i).velz;

        % plot particle
        plot3(particle(i).x,particle(i).y,8-particle(i).fitness,'*r')
        
        % Evaluate particle fitness
        %prevfitness = particle(i).fitness;
        particle(i).fitness = 8 - psofunc(particle(i).x,particle(i).y); % in example he subtracts from max possible value to turn into minimum function?
        
        % grab fitness statistics
        fitness = [ fitness particle(i).fitness ];
        
        % Check if current fitness is better than personal best (minimizing, does it go up or down???)
        if (particle(i).fitness > particle(i).pbestfitness) 
            particle(i).pbestfitness = particle(i).fitness;
            particle(i).pbestx = particle(i).x;
            particle(i).pbesty = particle(i).y;
            particle(i).pbestz = particle(i).z;
        end
        
        % Check if this particles best fitness is better than the global best fitness
        if (particle(i).pbestfitness > gbestfitness)
            last_gbestfitness = gbestfitness;   % store previous one for difference check
            gbestfitness = particle(i).pbestfitness;
            gbestx = particle(i).pbestx;
            gbesty = particle(i).pbesty;
            gbestz = particle(i).pbestz;
        end
    end
    
    % update inertia weigth
    icnt = icnt + 1;
    w(icnt) = max(wmin,w(icnt-1)*wsf);
    
    max_fitness = [ max_fitness max(fitness) ];
    average_fitness = [ average_fitness sum(fitness)/length(fitness) ];
    min_fitness = [ min_fitness min(fitness) ];
    
%     disp(['global best x,y,z   =  ',num2str(gbestx),'  ',num2str(gbesty),'  ',num2str(gbestz)])
%     disp(['global best fitness =  ',num2str(gbestfitness)])
%     disp(['original, 8-fitness =  ',num2str(8-gbestfitness)])
% 
    drawnow;
    pause(0.1);
end


%% Results
disp(['global best x,y,z   =  ',num2str(gbestx),'  ',num2str(gbesty),'  ',num2str(gbestz)])
disp(['global best fitness =  ',num2str(gbestfitness)])
disp(['original, 8-fitness =  ',num2str(8-gbestfitness)])
maximum = 8 - max(max_fitness)
average = 8 - sum(average_fitness)/length(average_fitness)
minimum = 8 - min(min_fitness)


%% Plot stuff

% clf
figure()
hold on
plot(8-max_fitness)
plot(8-average_fitness)
plot(8-min_fitness)
xlabel('Iteration');
ylabel('Fitness F(x)');
title(['Particle Swarm Fitness  -  population: ',num2str(num_particles(rr)),' max velocity: ',num2str(vmax(rr)),' phi1: ',num2str(phi1(rr)),' phi2: ',num2str(phi2(rr))]);
legend('max','average','min');

drawnow;
pause
end % simulation modifications, 4 different ones





