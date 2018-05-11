clc
clear

rand('state',0);  % edit this to get different results every time

n = 2^30-1;
x = 0:1000:n;
% plot(x/n, (x/n).^2, 'r', x/n, (x/n).^10, 'b--', 'LineWidth', 2)
% xlabel('Normalized x (x/c)');
% ylabel('Normalized Fitness F(x)');
% title('Fitness Function to Maximize');
% legend('x^2', 'x^10');


%% Genetic Algorithm

max_generation = 25;

length_chromosome = 30;
population_size = 30;
pcross = 0.6;
pmutate = 0.0333;

% c = 2^length_chromosome - 1;
% Fx = (x/c).^10;

% determine encoding, how many bits
% num_bits = log(x(length(x)))/log(2); % = 30, hardcoded as length_chromosomes

% set initial population & 
pop = zeros(population_size, length_chromosome);
x_decode = zeros(1, population_size);
fitness = zeros(1, population_size);
for ii=1:population_size,
    for jj=1:length_chromosome,
       pop(ii,jj) = (rand <= 0.5);  % generate chromosome bit
    end;
    
    % decode chromosome into decimal number
    accum = 0;
    powerof2 = 1;
    for jj=length_chromosome:-1:1,
        if ( pop(ii,jj) ),
            accum = accum + powerof2;
        end
        powerof2 = 2*powerof2;
    end
    x_decode(ii) = accum;
    
    % determine fitness from objective function
    coeff = 2^length_chromosome - 1;
    fitness(ii) = (x_decode(ii)/coeff)^10;
end;


% evaluate fitness of initial population
% I'd have to loop again if I separated these steps
% Assumes chromosomes are encoded with LSB first
% x_decode = zeros(1, population_size);
% fitness = zeros(1, population_size);
% for ii=1:population_size,
%     powerof2 = 1;
%     for jj=1:length_chromosome,
%         if (pop(ii,jj)),
%            x_decode(ii) = x_decode(ii) + powerof2;    % decode chromosome bit
%         end;
%         powerof2 = 2*powerof2;
%     end;
%     % calculate fitness of decoded chromosome (run through Fx)
%     fitness(ii) = (x_decode(ii)/c).^10;
% end;

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
            cross_point = rnd(1,length_chromosome-1);
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
        
        % decode child 1 chromosome into decimal number
        accum = 0;
        powerof2 = 1;
        for jj=length_chromosome:-1:1,
            if ( child(ii,jj) ),
                accum = accum + powerof2;
            end
            powerof2 = 2*powerof2;
        end
        x_decode(ii) = accum;

        % determine fitness from objective function
        coeff = 2^length_chromosome - 1;
        fitness(ii) = (x_decode(ii)/coeff)^10;
        
         % decode child 2 chromosome into decimal number
        accum = 0;
        powerof2 = 1;
        for jj=length_chromosome:-1:1,
            if ( child(ii+1,jj) ),
                accum = accum + powerof2;
            end
            powerof2 = 2*powerof2;
        end
        x_decode(ii+1) = accum;

        % determine fitness from objective function
        coeff = 2^length_chromosome - 1;
        fitness(ii+1) = (x_decode(ii+1)/coeff)^10;
        
        ii = ii + 2;    % create children in pairs
    end;

    pop = child;
    
    
    
    
%     child = pop;
%     for ii=1:length(parents)-1,
% 
%         % mutate or crossover first ????
%         % M(i) mutate - perturb the mated population stochastically
%         for bb=1:length_chromosome   
%             if (rand <= pmutate)
%                 child(parents(ii),bb) = ~pop(parents(ii),bb);
%             end    
%             if (rand <= pmutate)
%             	child(parents(ii+1),bb) = ~pop(parents(ii+1),bb);
%             end
%         end
%         
%         %pop = child;
%         
%         % M(i) crossover - recombine genes of selected parents
%         if (rand <= pcross)
%             cross_point = randi([1 length_chromosome]);
%             for jj=cross_point:length_chromosome
%                 child(parents(ii),jj) = pop(parents(ii+1),jj); %right hand side should be child?
%                 child(parents(ii+1),jj) = pop(parents(ii),jj);
%             end
%         end
%         
%         ii = ii + 1;    % create children in pairs
%     end;


%     % update population
%     pop = child;
%     % evaluate new fitness 
%     x_decode = zeros(1, population_size);
%     for ii=1:population_size,
%         powerof2 = 1;
%         for jj=1:length_chromosome,
%            if (pop(ii,jj)),
%                x_decode(ii) = x_decode(ii) + powerof2;    % decode chromosome bit
%            end;
%            powerof2 = 2*powerof2;
%         end;
%         % calculate fitness of decoded chromosome (run through Fx)
%         fitness(ii) = (x_decode(ii)/c).^10;
%     end;

    sum_fitness = sum(fitness);
    max_fitness = [max_fitness max(fitness)];
    average_fitness = [average_fitness sum_fitness/population_size];
    min_fitness = [min_fitness min(fitness)];
    
end

%% Plot stuff

figure()
hold on
plot(max_fitness)
plot(min_fitness)
plot(average_fitness)
xlabel('Generation');
ylabel('Normalized Fitness F(x)');
title('Genetic Algorithm Fitness Growth, population 40');
legend('max', 'min', 'average');







