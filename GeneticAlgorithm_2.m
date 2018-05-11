clc
clear

rand('state',0);    % edit this to get different results every time

x=-3:0.188:3;       % the resolution here doesn't seem to make a difference
y=-3:0.188:3;       %  used 32 values to match the resolution of 'peaks'
z = zeros(2,length(x));
for i=1:length(x)
    for j=1:length(y)
        z(j,i) = 3*((1-x(i))^2)*exp((-x(i)^2)-((y(j)+1)^2)) - 10*((x(i)/5)-(x(i)^3)-(y(j)^5))*exp((-x(i)^2)-(y(j)^2)) - (1/3)*exp((-(x(i)+1)^2)-(y(j)^2)) +7;
    end
end

mesh(x,y,z)
xlabel('x');
ylabel('y');
zlabel('z');
title('Fitness Function to Maximize');


% figure()
% hold on

%% Genetic Algorithm
% for ss=1:20

max_generation = 25;
length_chromosome = 8;
population_size = 20;
pcross = 1.0;
pmutate = 0.01;

% set initial population & 
pop = zeros(population_size, length_chromosome);
x_decode = zeros(1, population_size);
y_decode = zeros(1, population_size);
x_mapped = zeros(1, population_size);
y_mapped = zeros(1, population_size);
fitness = zeros(1, population_size);

for ii=1:population_size,
    for jj=1:length_chromosome,
       pop(ii,jj) = (rand <= 0.5);  % generate chromosome bit
    end;
    
    % check that chromosome can be decomposed into substrings
    if (mod(length_chromosome,2)),
        disp('error, cant evenly derive substrings from chromosome')
        break;
    end;
    
    % decode chromosome substring x into decimal number
    x_substring = pop(ii,1:length_chromosome/2);
    accum = 0;
    powerof2 = 1;
    for jj=length(x_substring):-1:1,
        if ( x_substring(jj) ),
            accum = accum + powerof2;
        end
        powerof2 = 2*powerof2;
    end
    x_decode(ii) = accum;
    
    % decode chromosome substring y into decimal number
    y_substring = pop(ii,length_chromosome/2+1:length_chromosome);
    accum = 0;
    powerof2 = 1;
    for jj=length(y_substring):-1:1,
        if ( y_substring(jj) ),
            accum = accum + powerof2;
        end
        powerof2 = 2*powerof2;
    end
    y_decode(ii) = accum;
    
    % map decoded substring x to range [-3,3]
    x_mapped(ii) = ((3-(-3)) / ((2^4)-1)) * x_decode(ii) + (-3);
    if (x_mapped(ii) < -3 | x_mapped(ii) > 3)
        disp('error, x_mapped out of bounds')
    end
    
    % map decoded substring y to range [-3,3]
    y_mapped(ii) = ((3-(-3)) / ((2^4)-1)) * y_decode(ii) + (-3);
    if (y_mapped(ii) < -3 | y_mapped(ii) > 3)
        disp('error, y_mapped out of bounds')
    end
    
    % determine fitness from objective function
    fitness(ii) = objectivefunc(x_mapped(ii),y_mapped(ii));
end;

sum_fitness = sum(fitness);
average_fitness = sum_fitness/population_size;
max_fitness = max(fitness);
min_fitness = min(fitness);


% ga optimization loop
for kk=2:max_generation    
    
    % select mating pool M(i) from population for offspring reproduction
    parents = zeros(1, population_size);
    ii = 1;
	while (ii <= population_size),
        
        % select 1st parent
        jj = 0;
        partial_sum = 0;
        random_variable = rand * sum_fitness;
        while ((partial_sum <= random_variable) & (jj < population_size)),
            jj = jj + 1;
            partial_sum = partial_sum + fitness(jj);
        end;
        parents(ii) = jj;
       
        % select 2nd parent
        jj = 0;
        partial_sum = 0;
        random_variable = rand * sum_fitness;
        while ((partial_sum <= random_variable) & (jj < population_size)),
            jj = jj + 1;
            partial_sum = partial_sum + fitness(jj);
        end;
        parents(ii+1) = jj;
        
        % doing all these 2 children at a time inside population_size loop
        if (rand <= pcross)
            cross_point = randi([1,(length_chromosome-1)],1,1);
        else
            cross_point = length_chromosome;
        end
            
        % M(i) mutate - perturb the mated population stochastically
        for bb=1:cross_point   
            if (rand <= pmutate)
                child(ii,bb) = ~pop(parents(ii),bb);
            else
                child(ii,bb) = pop(parents(ii),bb);
            end
            if (rand <= pmutate)
                child(ii+1,bb) = ~pop(parents(ii+1),bb);
            else
                child(ii+1,bb) = pop(parents(ii+1),bb);
            end
        end
            
        % M(i) crossover - recombine genes of selected parents
        % TODO check this, his crossover includes a mutate (so there's two
        % mutates)...
        for jj=cross_point+1:length_chromosome
            if (rand <= pmutate)
                child(ii,jj) = ~pop(parents(ii+1),jj);
            else
                child(ii,jj) = pop(parents(ii+1),jj);
            end
            if (rand <= pmutate)
                child(ii+1,jj) = ~pop(parents(ii),jj);
            else
                child(ii+1,jj) = pop(parents(ii),jj);
            end
        end
        
        % decode child 1 chromosome x substring into decimal number
        x_substring(1,:) = child(ii,1:length_chromosome/2);
        accum = 0;
        powerof2 = 1;
        for jj=length(x_substring(1,:)):-1:1,
            if ( x_substring(1,jj) ),
                accum = accum + powerof2;
            end
            powerof2 = 2*powerof2;
        end
        x_decode(ii) = accum;

        % decode child 1 chromosome y substring into decimal number
        y_substring(1,:) = child(ii,length_chromosome/2+1:length_chromosome);
        accum = 0;
        powerof2 = 1;
        for jj=length(y_substring(1,:)):-1:1,
            if ( y_substring(1,jj) ),
                accum = accum + powerof2;
            end
            powerof2 = 2*powerof2;
        end
        y_decode(ii) = accum;
        
        % decode child 2 chromosome x substring into decimal number
        x_substring(2,:) = child(ii+1,1:length_chromosome/2);
        accum = 0;
        powerof2 = 1;
        for jj=length(x_substring(2,:)):-1:1,
            if ( x_substring(2,jj) ),
                accum = accum + powerof2;
            end
            powerof2 = 2*powerof2;
        end
        x_decode(ii+1) = accum;
        
        % decode child 2 chromosome y substring into decimal number
        y_substring(2,:) = child(ii+1,length_chromosome/2+1:length_chromosome);
        accum = 0;
        powerof2 = 1;
        for jj=length(y_substring(2,:)):-1:1,
            if ( y_substring(2,jj) ),
                accum = accum + powerof2;
            end
            powerof2 = 2*powerof2;
        end
        y_decode(ii+1) = accum;

        
        % map decoded substring x to range [-3,3]
        x_mapped(ii) = ((3-(-3)) / ((2^4)-1)) * x_decode(ii) + (-3);     % child 1
        if (x_mapped(ii) < -3 | x_mapped(ii) > 3)
            disp('error, x_mapped out of bounds')
        end
        x_mapped(ii+1) = ((3-(-3)) / ((2^4)-1)) * x_decode(ii+1) + (-3); % child 2
        if (x_mapped(ii+1) < -3 | x_mapped(ii+1) > 3)
            disp('error, x_mapped out of bounds')
        end
        
        % map decoded substring y to range [-3,3]
        y_mapped(ii) = ((3-(-3)) / ((2^4)-1)) * y_decode(ii) + (-3);     % child 1
        if (y_mapped(ii) < -3 | y_mapped(ii) > 3)
            disp('error, y_mapped out of bounds')
        end
        y_mapped(ii+1) = ((3-(-3)) / ((2^4)-1)) * y_decode(ii+1) + (-3); % child 2
        if (y_mapped(ii+1) < -3 | y_mapped(ii+1) > 3)
            disp('error, y_mapped out of bounds')
        end
    
        % determine fitness from objective function
        fitness(ii) = objectivefunc(x_mapped(ii),y_mapped(ii)); 
        fitness(ii+1) = objectivefunc(x_mapped(ii+1),y_mapped(ii+1));
        
        ii = ii + 2;    % create children in pairs
    end;

    pop = child;
    sum_fitness = sum(fitness);
    max_fitness = [max_fitness max(fitness)];
    average_fitness = [average_fitness sum_fitness/population_size];
    min_fitness = [min_fitness min(fitness)];
    
end

%% Results
maximum = max(max_fitness) - 7
average = sum(average_fitness)/length(average_fitness) - 7
minimum = min(min_fitness) - 7


%% Plot stuff

% clf
figure()
hold on
plot(max_fitness)
plot(average_fitness)
plot(min_fitness)
xlabel('Generation');
ylabel('Normalized Fitness F(x)');
title(['Genetic Algorithm Fitness Growth  -  pcross: ' num2str(pcross) '  pmutate: ' num2str(pmutate) '  popsize: ' num2str(population_size)]);
legend('max','average','min');

% drawnow
% pause
% end % simulation iterations ss





