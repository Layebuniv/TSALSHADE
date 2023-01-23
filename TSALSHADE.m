%%%%%%%%%%%%%%%%%%%
% TSALSHADE: TSALSHADE: Improved LSHADE Algorithm with Tangent Search
% The code of LSHADE is used
%Author, inventor and programmer: Layeb Abdesslem
% In review paper
% e-Mail: abdesslem.layeb@univ-constantine2.dz
%%%%%%%%%%%%%%%%%%% 

% clc;
% clear all;

format long;

% cec 2022 problem parameters 
fhd=@cec22_test_func;
optimum=[300	400  600 800	900	1800	2000	2200	2300 2400	2600 2700 ];
dim = 20;
max_nfes = 10000 * dim;

val_2_reach = 10^(-16);
max_region = 100.0;
min_region = -100.0;
lu = [-100 * ones(1, dim); 100 * ones(1, dim)];
lb=-100;
ub=100;

for func = 1:12
  
 
  fprintf('\n-------------------------------------------------------\n')
  fprintf('Function = %d, Dimension size = %d\n', func, dim) 
 
  for run_id = 1 : 3
      rng('shuffle')
    %%  parameter settings for L-SHADE
    p_best_rate = 0.11;
    arc_rate = 1;
    memory_size = 5;
    pop_size = 18* dim; %15
   
    max_pop_size = pop_size;
    min_pop_size = 4;
     %%  parameter settings for TSA
     subpopr=1/6;  % subpoppoulation size
     tfm=0.01;     % tangent flight modification rate  
    %% Initialize the main population
    popold = repmat(lu(1, :), pop_size, 1) + rand(pop_size, dim) .* (repmat(lu(2, :) - lu(1, :), pop_size, 1));
    pop = popold; % the old population becomes the current population

    fitness = feval(fhd,pop',func);
    fitness = fitness';

    nfes = 0;
    bsf_fit_var = 1e+30;
    bsf_solution = zeros(1, dim);

    %%%%%%%%%%%%%%%%%%%%%%%% for out
    for i = 1 : pop_size
      nfes = nfes + 1;

      if fitness(i) < bsf_fit_var
	bsf_fit_var = fitness(i);

	bsf_solution = pop(i, :);
      end

	
      if nfes > max_nfes; break; end


    end
    %%%%%%%%%%%%%%%%%%%%%%%% for out

    memory_sf = 0.5 .* ones(memory_size, 1); %0.1
    memory_cr = 0.5 .* ones(memory_size, 1); %0.5 -0.8
    memory_pos = 1;

    archive.NP = arc_rate * pop_size; % the maximum size of the archive
    archive.pop = zeros(0, dim); % the solutions stored in te archive
    archive.funvalues = zeros(0, 1); % the function value of the archived solutions

    %% main loop
    while nfes < max_nfes
       
          if   (rand <=0.5 && nfes<=0.7*max_nfes)  || (rand <=1 && nfes<=0.3*max_nfes) %0.4

           for j=1:round(subpopr*pop_size)  
                X=popold(j,:);
                %%-----------------tangent flight operator----------------
               X=Tangent_Flight(X,bsf_solution,dim,nfes);
                    
               %%                   
                    
                    Xnew=popold(j,:);
                    id=randi(dim);
                    ind=find(rand(1,dim)<=tfm); %0.05
                    ind=[ind, id];
                    % end
                    Xnew(ind)=X(ind);
                    % -------------solution evaluation
                    %-------------Check boundries
                     
                     Xnew(Xnew>ub)=rand*(ub - lb) + lb;
                    Xnew(Xnew<lb)=rand*(ub - lb) + lb;
                    % Xnew= boundConstraint ( Xnew,X, lu);
                    % fitness= fun(Xnew);
                    fit=feval(fhd,Xnew',func) ;
                    if nfes > max_nfes; break; end
           
                    if  fit <fitness(j) 
                        popold(j,:) =Xnew ;
                       fitness(j) =fit;
                        if  fit< bsf_fit_var
                            bsf_solution=Xnew;
                            bsf_fit_var=fit   ;
                        end
                    end
                    nfes= nfes+1;
                    
                    if  nfes>max_nfes
                        break;
                    end
                    
                end
                
           end
  %% LSHADE operaor       
     if  rand<=1    && nfes>=0.5*max_nfes   %0.5  0.3
      pop = popold; % the old population becomes the current population
      [temp_fit, sorted_index] = sort(fitness, 'ascend');

      mem_rand_index = ceil(memory_size * rand(pop_size, 1));
      mu_sf = memory_sf(mem_rand_index);
      mu_cr = memory_cr(mem_rand_index);

      %% for generating crossover rate
      cr = normrnd(mu_cr, 0.2); %0.2
      term_pos = find(mu_cr == -1);
      cr(term_pos) = 0;
      cr = min(cr, 1);%0.8
      cr = max(cr, 0);

      %% for generating scaling factor
      sf = mu_sf + 0.05 * tan(pi * (rand(pop_size, 1) )); %0.01
      pos = find(sf <= 0);

      while ~ isempty(pos)
	sf(pos) = mu_sf(pos) + 0.05 * tan(pi * (rand(length(pos), 1) ));
	pos = find(sf <= 0);
      end

      sf = min(sf, 1); 
      
      r0 = [1 : pop_size];
      popAll = [pop; archive.pop];
      [r1, r2] = gnR1R2(pop_size, size(popAll, 1), r0);
      [r3, r4] = gnR1R2(pop_size, size(popAll, 1), r0);
      [r5, r6] = gnR1R2(pop_size, size(popAll, 1), r0);
      pNP = max(round(p_best_rate * pop_size), 2); %% choose at least two best solutions
      randindex = ceil(rand(1, pop_size) .* pNP); %% select from [1, 2, 3, ..., pNP]
      randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
      pbest = pop(sorted_index(randindex), :); %% randomly choose one of the top 100p% solutions

      %vi = pop + sf(:, ones(1, dim)) .* (pbest - pop + pop(r1, :) - popAll(r2, :));
      if rand<=0.2  && nfes>=0.7*max_nfes
          %% apply meanDE mutation--- Layeb, Abdesslem. Differential Evolution Algorithms with Novel Mutations, Adaptive Parameters and Weibull Flight Operator
       vi =(pop(r1, :) +popAll(r2, :))/2 + sf(:, ones(1, dim)) .* ((pbest +pop(r1, :))/2- 2*pop  + (pop(r1, :) +popAll(r2, :))/2);
     
      else
                   
       vi = pop +  sf(:, ones(1, dim)) .* (pbest - pop + pop(r1, :) - popAll(r2, :));
     
          end
      vi = boundConstraint(vi, pop, lu);

      mask = rand(pop_size, dim) > cr(:, ones(1, dim)); % mask is used to indicate which elements of ui comes from the parent
      rows = (1 : pop_size)'; cols = floor(rand(pop_size, 1) * dim)+1; % choose one position where the element of ui doesn't come from the parent
      jrand = sub2ind([pop_size dim], rows, cols); mask(jrand) = false;
      ui = vi; ui(mask) = pop(mask);

      children_fitness = feval(fhd, ui', func);
      children_fitness = children_fitness';

      %%%%%%%%%%%%%%%%%%%%%%%% for out
      for i = 1 : pop_size
	nfes = nfes + 1;

	if children_fitness(i) < bsf_fit_var
	  bsf_fit_var = children_fitness(i);
	  bsf_solution = ui(i, :);
	end

	 
	if nfes > max_nfes; break; end
      end
      %%%%%%%%%%%%%%%%%%%%%%%% for out

      dif = abs(fitness - children_fitness);

      %% I == 1: the parent is better; I == 2: the offspring is better
      I = (fitness > children_fitness);
      goodCR = cr(I == 1);  
      goodF = sf(I == 1);
      dif_val = dif(I == 1);

%      isempty(popold(I == 1, :))   
      archive = updateArchive(archive, popold(I == 1, :), fitness(I == 1));

      [fitness, I] = min([fitness, children_fitness], [], 2);
      
      popold = pop;
      popold(I == 2, :) = ui(I == 2, :);

      num_success_params = numel(goodCR);

      if num_success_params > 0 
	sum_dif = sum(dif_val);
	dif_val = dif_val / sum_dif;

	%% for updating the memory of scaling factor 
	memory_sf(memory_pos) = (dif_val' * (goodF .^ 2)) / (dif_val' * goodF);

	%% for updating the memory of crossover rate
	if max(goodCR) == 0 || memory_cr(memory_pos)  == -1
	  memory_cr(memory_pos)  = -1;
	else
	  memory_cr(memory_pos) = (dif_val' * (goodCR .^ 2)) / (dif_val' * goodCR);
	end

	memory_pos = memory_pos + 1;
	if memory_pos > memory_size;  memory_pos = 1; end
      end

      %% for resizing the population size
      plan_pop_size = round((((min_pop_size - max_pop_size) / max_nfes) * nfes) + max_pop_size);

      if pop_size > plan_pop_size %&&  nfes>=0.5*max_nfes
	reduction_ind_num = pop_size - plan_pop_size;
	if pop_size - reduction_ind_num <  min_pop_size; reduction_ind_num = pop_size - min_pop_size;end

	pop_size = pop_size - reduction_ind_num;
	for r = 1 : reduction_ind_num
	  [valBest indBest] = sort(fitness, 'ascend');
	  worst_ind = indBest(end);
	  popold(worst_ind,:) = [];
	  pop(worst_ind,:) = [];
	  fitness(worst_ind,:) = [];
	end
	  
	archive.NP = round(arc_rate * pop_size); 

	if size(archive.pop, 1) > archive.NP 
	  rndpos = randperm(size(archive.pop, 1));
	  rndpos = rndpos(1 : archive.NP);
	  archive.pop = archive.pop(rndpos, :);
	end
      end
    end
    end
    bsf_error_val = bsf_fit_var - optimum(func);
    if bsf_error_val < val_2_reach
        bsf_error_val = 0;
     end

    fprintf('%d th run, best-so-far error value = %1.8e\n', run_id , bsf_error_val)
%     outcome(func,run_id)= bsf_error_val;    
  end %% end 1 run


end %% end 1 function run
 

% for i=1:12
%     outcome(i,31)=(mean(outcome(i,1:30)));
%     outcome(i,32)=std(outcome(i,1:30));
%     outcome(i,33)=min(outcome(i,1:30));
%     outcome(i,34)=max(outcome(i,1:30));
% end
% 
% 
% 
% % save results to excell file
% writematrix(outcome,'TSAlshade2.xlsx','Sheet',1)