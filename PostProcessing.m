%load thetaout.txt;
%load Lambdaout.txt; 
%load Etaout.txt;   
% load gammaout.txt;
%load Workspace500.mat;
%load DesignMatrix.mat;

% matrici già salvate considerando burnin e thinning:
thetaout = readmatrix('Thetaout_seed_250.csv');  
Lambdaout = readmatrix('Lambdaout_seed_250.csv'); 
Etaout = readmatrix('Etaout_seed_250.csv');  

% queste dall'ultimo update del GS (con 10000 iterazioni e 1000 di burnin):
Thetatilde = readmatrix('thetatilde_seed_250.csv'); 
B = readmatrix('B_seed_250.csv');
thr1 = readmatrix('thresholds_seed_250.csv'); 

Y = readmatrix('Y.csv'); 

nrun = 10000; %500; 
burn = 1000; %20;
thin = 5;
its = burn:thin:nrun-1; 

q=5; 
num_righe = size(Thetatilde, 1);  % Numero di righe di Thetatilde
num_colonne = size(Thetatilde, 2); % Numero di colonne di The tatilde
THETA = zeros(num_righe, num_colonne); 
for i = 1:q
  index = find(bsxfun(@hypot, Thetatilde(: , 2*i-1), Thetatilde(:, 2*i)) >= thr1(:, i)); 
  if ~isempty(index)
    THETA(index, [2*i - 1 2*i]) = Thetatilde(index, [2*i - 1 2*i]);  
  end
end

% The following counts how many genes exhibit 24 hours periodicity
% (matrices from the last update of GS)
per24 = find(THETA(: , 1) == 0 & THETA(: , 2) == 0 & THETA(: , 3) == 0 & THETA(: , 4) == 0 ...
    & THETA(:, 5) == 0 & THETA(:, 6) == 0 & THETA(:, 7) == 0 & THETA(:, 8) == 0 ...
    & THETA(:, 9) ~= 0 & THETA(:, 10) ~= 0) ; length(per24) 

p = 500;                                    % Number of simulated protein profiles ( = number of genes)
T = 24;                                     % Number of observations per profile 
q = 5;                                      % Number of sin (cos) bases for a total of 2q bases -- no intercept
tij = [0:2:46]';                            % Time points at which the probes have been measured
m = max(tij) ;                              % Max hour of measurement
tij = tij / max(tij); 
cnt = 0; 
ib = sort(per24);
indicator = 0:500:4500;
indicatorLambda = 0:p:(12-1)*p; 
aperiodic = zeros(length(its)-1, T);


figure;
% Code to plot estimated curves
for h = 1:25
   i = ib(h); 
   cnt = cnt + 1; 
   subplot(5,5,cnt) 
   %thti = thetaout(its, indicator + i); 
   thti = thetaout(:, indicator + i); %considero tutte le righe perchè matrici già salvate con burnin+thinning
   % ghti = gammaout(its, indicator + i);
   for l = 1:length(its)-1
       %aperiodic(l, :) = reshape(Etaout(:, its(l)), T, 12)*Lambdaout(indicatorLambda + i, its(l)); 
       aperiodic(l, :) = reshape(Etaout(:, l), T, 12)*Lambdaout(indicatorLambda + i, l);
   end
   Eyi = (B*thti')' + aperiodic;        % REMOVED: + (C * ghti')'
   est = zeros(24,3); 
   %est = zeros(ni(i),3); 
   est(:,1) = mean(Eyi)'; 
   est(:,2:3) = prctile(Eyi,[2.5 97.5])';                                          % 95% pointwise credible interval
   plot(tij.*46, Y(:, i),'o', 'Color', 'k', 'MarkerSize',4)
   xlabel(['Time (hours)'])
   xlim([0 46])
   title(['Gene ' num2str(i)])
   line(tij.*46, Y(:, i),'LineStyle', '-','Color', 'k','LineWidth', 1.5 )
   line(tij.*46, est(:,1),'LineStyle', '-','Color', 'blue','LineWidth', 1.5)        % posterior mean estimate
   line(tij.*46, est(:,2),'LineStyle','--', 'Color', 'red', 'LineWidth',1)          % 99% pointwise intervals
   line(tij.*46, est(:,3),'LineStyle','--','Color', 'red','LineWidth',1)
end

% NO IDEA WHAT THIS DOES 
for h = 1:length(per24) 
    j = per24(h);
    %thti = thetaout(its, indicator + j);
    thti = thetaout(:, indicator + j);
    for l = 1:8
        countres(h, l) = (length(find(thti(:, l) == 0))/ size(thti,1)) * 100;
    end
    for l = 9:10
        countres(h,l) = (length(find(thti(:, l) ~= 0))/ size(thti,1)) * 100;
    end
end

% load DesignMatrix;
i = 400;
   %thti = thetaout(its, indicator + i);
   thti = thetaout(:, indicator + i);
   % ghti = gammaout(its, indicator + i);
   %thti = thetaout(its,((i-1)*(2*q+1) + 1):(i*(2*q + 1)));
   for l = 1:length(its)-1
       %aperiodic(l, :) = reshape(Etaout(:, its(l)), T, 15)*Lambdaout(indicatorLambda + i, its(l));
       aperiodic(l, :) = reshape(Etaout(:, l), T, 12)*Lambdaout(indicatorLambda + i, l);
   end
[mean(thti)' THETA(i,:)']
% [mean(ghti)' Gamma(i,:)']
%[mean(aperiodic)' (eta * Lambda(i,:)') abs(mean(aperiodic)') - abs(eta * Lambda(i,:)')]

% load trh1out.txt
% indthr1 = 0:p:(p-1)*q; 
% i = 50000; its = (burn / thin + 1): round(i / thin); 
% cnt = 0; 
% j = 373;
% for l = 1:length(indthr1)
%     cnt = cnt + 1;
%     subplot(3,2, cnt)
%     plot(trh1out(its, indthr1(l) + j));
%     hold on
%     plot([0 6000], [thr1(j, l) thr1(j, l)], 'Color', 'red', 'LineWidth', 1.5)
%     hold off
%     title(['Subejct ', num2str(j), ' thr_{theta} ', num2str(l)])
% end
% 
% % % load Factorout.txt
% % % prctile(Factorout(its), [2.5, 5, 97.5])
% 
% % % load sigmaout.txt
% % % mean(mean(sigmaout(its,:)))
% % % prctile(mean(sigmaout(its,:)), [2.5, 5, 97.5])
   
for h = 1:length(per24)
    j = per24(h);
    %thti = thetaout(its, indicator + j);
    thti = thetaout(:, indicator + j);
    %[mean(thti)' THETA(j,:)']
    % -- % Computing the probability that the protein is circardian % -- % 
    prob(h) = length(find(thti(: , 1) == 0 & thti(: , 2) == 0 & thti(: , 3) == 0 & thti(: , 4) == 0 ...
    & thti(:, 5) == 0 & thti(:, 6) == 0 & thti(:, 7) == 0 & thti(:, 8) == 0 ...
    & thti(:, 9) ~= 0 & thti(:, 10) ~= 0))/(length(its)-1) ; 
end

% calcolo di odd, preso da MCMC.m
odd = []; even = [];
for i = 1:q
  odd = [odd, 2*i - 1];                           % Index for sin bases
  even = [even, 2*i];                             % Index for cos bases
end  
dlmwrite('pcircardian.txt', [per24, countres(:, odd), prob'], 'delimiter', ' ', 'precision', 6);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% %    Computing the prob of a protein being circardian    % %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
indicator = 0:p:(2*q - 1)*p;
for h = 1:p
    %thti = thetaout(its, indicator + h); 
    thti = thetaout(:, indicator + h); 
    circardian(h) = length(find(thti(: , 1) == 0 & thti(: , 2) == 0 & thti(: , 3) == 0 & thti(: , 4) == 0 ...
    & thti(:, 5) == 0 & thti(:, 6) == 0 & thti(:, 7) == 0 & thti(:, 8) == 0 ...
    & thti(:, 9) ~= 0 & thti(:, 10) ~= 0))/(length(its)-1) ; 
end

ind = 1:p; ind = ind';
[A_C ordC ] = sort(circardian', 'descend');
circ_prot = ind(ordC);
beta = 1 - A_C;
list = [A_C 1 - A_C ordC];

% % % Choose a data-dependent kappa % %
% thr = 0.05; 
% kappa = .09;
% exp_fd = (list(:,2)'*(list(:,2) <= kappa));           % Expected number of false discoveries
% fdr = exp_fd / length(find(list(:, 2) <= kappa)) ;    % Expecte rate of false discoveries: needs to be <= thr;
% circa = ord(find(list(:, 2) <= kappa));               % Proteins identified as circardian
% 
% sensitivity = length(find(THETA(circa,1) == 0 & THETA(circa,2) == 0 & THETA(circa,3) == 0 & THETA(circa,4) == 0 & THETA(circa,5) == 0 & ...
%     THETA(circa,6) == 0 & THETA(circa,7) == 0 & THETA(circa,8) == 0 & THETA(circa,9) ~= 0 & THETA(circa,10) ~= 0)) / 25;
% 
% notcirc = setdiff(1:p, circa);              % Proteins not identified as circardian
% true_not_circ = p - 25;
% % Numb. of not circardian correctly identified as not circardian
% spec = (length(notcirc) - length(find(THETA(notcirc,1) == 0 & THETA(notcirc,2) == 0 & THETA(notcirc,3) == 0 & ...
%     THETA(notcirc,4) == 0 & THETA(notcirc,5) == 0 & THETA(notcirc,6) == 0 & THETA(notcirc,7) == 0 ...
%     & THETA(notcirc,8) == 0 & THETA(notcirc,9) ~= 0 & THETA(notcirc,10) ~= 0))) / true_not_circ;      
% 

% % ROC curve as function of different values of kappa % %
vect = 0:.005:1;
sensitivity=zeros(1,length(vect));
spec=zeros(1,length(vect));
for h = 1:length(vect)
 kappa = vect(h);
 circa = ordC(find(list(:, 2) <= kappa));     % Proteins identified as circardian
 sensitivity(h) = length(find(THETA(circa,1) == 0 & THETA(circa,2) == 0 & THETA(circa,3) == 0 & THETA(circa,4) == 0 & THETA(circa,5) == 0 & ...
    THETA(circa,6) == 0 & THETA(circa,7) == 0 & THETA(circa,8) == 0 & THETA(circa,9) ~= 0 & THETA(circa,10) ~= 0)) / 25;
 notcirc = setdiff(1:p, circa);             % Proteins not identified as circardian
 true_not_circ = p - 25;
 spec(h) = (length(notcirc) - length(find(THETA(notcirc,1) == 0 & THETA(notcirc,2) == 0 & THETA(notcirc,3) == 0 & ...
    THETA(notcirc,4) == 0 & THETA(notcirc,5) == 0 & THETA(notcirc,6) == 0 & THETA(notcirc,7) == 0 ...
    & THETA(notcirc,8) == 0 & THETA(notcirc,9) ~= 0 & THETA(notcirc,10) ~= 0))) / true_not_circ;       
end
figure;
plot(1 - spec, sensitivity)
xlim([0, 1])
xlabel(['False positive rate'])
ylabel(['True positive rate'])
title(['ROC curve'])
dlmwrite('ROCcircardianDepGibbs.txt', [(1-spec)' sensitivity'], 'delimiter', ' ', 'precision', 6);

% load ROCcircardianChainThree.txt
% plot(ROCcircardianChainThree(:, 1), ROCcircardianChainThree(:, 2), 'Color', 'red')
% hold on
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % %
% %  Plot histogram of circardian probabilities  % %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % %
foo = setdiff(1:p, per24);

figure;
subplot(1, 2, 1);
histogram(circardian(foo), 20, 'FaceColor', [0.4, 0.5, 0.6], 'EdgeColor', 'w');
title('Non-Circadian Genes', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Probability', 'FontSize', 10); 
ylabel('Frequency', 'FontSize', 10);  
xlim([0, 1]);                         
grid on;                              
set(gca, 'Box', 'off', 'FontSize', 10, 'LineWidth', 1.2); 

subplot(1, 2, 2);
histogram(circardian(per24), 20, 'FaceColor', [1, 0.85, 0.1], 'EdgeColor', 'w'); 
title('Circadian Genes', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Probability', 'FontSize', 10); 
ylabel('Frequency', 'FontSize', 10);  
xlim([0, 1]);                         
grid on;                              
set(gca, 'Box', 'off', 'FontSize', 10, 'LineWidth', 1.2); 

sgtitle('Circadian Probability Distributions', 'FontSize', 14, 'FontWeight', 'bold');

dlmwrite('circardian.txt', circardian, 'delimiter', ' ', 'precision', 6);
dlmwrite('per24.txt', per24, 'delimiter', ' ', 'precision', 6);

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % %
% %         Test for periodicity (not restricting to circardian)        % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for h = 1:p
    %thti = thetaout(its, indicator + h); 
    thti = thetaout(:, indicator + h);
    periodicity(h) = length(find((thti(: , 1) ~= 0 & thti(: , 3) == 0 & thti(: , 5) == 0 & thti(: , 7) == 0 & thti(: , 9) == 0) | ...
        (thti(: , 1) == 0 & thti(: , 3) ~= 0 & thti(: , 5) == 0 & thti(: , 7) == 0 & thti(: , 9) == 0) | ...
        (thti(: , 1) == 0 & thti(: , 3) == 0 & thti(: , 5) ~= 0 & thti(: , 7) == 0 & thti(: , 9) == 0) | ...
        (thti(: , 1) == 0 & thti(: , 3) == 0 & thti(: , 5) == 0 & thti(: , 7) ~= 0 & thti(: , 9) == 0) | ...
        (thti(: , 1) == 0 & thti(: , 3) == 0 & thti(: , 5) == 0 & thti(: , 7) == 0 & thti(: , 9) ~= 0))) / (length(its)-1);
end

ind = 1:p; ind = ind';
[A ord ] = sort(periodicity', 'descend'); % A=probabilità ordinate di essere un gene periodico, ord=indice originale corrispondente
% circ_prot = ind(ord);
beta = 1 - A;
list = [A 1 - A ord];


periodic_vars = zeros(length(per24) , q);
periodic_vars(1:length(per24), 5) = per24;
sensitivity=zeros(1,p);
for i = 1:p
    selected_vars = ord(1:i);
    TP(i) = length(intersect(selected_vars, periodic_vars));
    FP(i) = length(selected_vars)-TP(i);
    FN(i) = length(intersect(periodic_vars, setdiff(1:p, selected_vars)));
    TN(i) = length(intersect(setdiff(1:p,  periodic_vars), setdiff(1:p, selected_vars)));
    sensitivity(i) = TP(i) / (TP(i) + FN(i));
    fpr(i) = FP(i) / (FP(i) + TN(i));
    FDR(i) = FP(i) / (FP(i) + TP(i));
end

figure;
plot(fpr, sensitivity)
xlim([0, 1])
xlabel(['False positive rate'])
ylabel(['True positive rate'])
title(['ROC curve'])

figure;
plot(sensitivity, FDR) %plot(FDR, sensitivity)
xlim([0, 1])
ylabel(['False discovery rate']) %x
xlabel(['Sensitivity (power)'])  %y

dlmwrite('GibbsChain1.txt', [sensitivity' fpr' FDR'], 'delimiter', ' ', 'precision', 6);

% title(['ROC curve'])
% % periodicity_list = list(1:101, 3);
% % 
% % % For the proteins most likely to be periodic, compute the probability of
% % % simple periodicity
% % for j = 1:length(periodicity_list)
% %     h = periodicity_list(j);
% %     thti = thetaout(its, indicator + h); 
% %     simple_per(j, 1) = length(find((thti(: , 1) ~= 0 & thti(: , 3) == 0 & thti(: , 5) == 0 & thti(: , 7) == 0 & thti(: , 9) == 0))) / length(its);
% %     simple_per(j, 2) = length(find((thti(: , 1) == 0 & thti(: , 3) ~= 0 & thti(: , 5) == 0 & thti(: , 7) == 0 & thti(: , 9) == 0))) / length(its);
% %     simple_per(j, 3) = length(find((thti(: , 1) == 0 & thti(: , 3) == 0 & thti(: , 5) ~= 0 & thti(: , 7) == 0 & thti(: , 9) == 0))) / length(its);
% %     simple_per(j, 4) = length(find((thti(: , 1) == 0 & thti(: , 3) == 0 & thti(: , 5) == 0 & thti(: , 7) ~= 0 & thti(: , 9) == 0))) / length(its);
% %     simple_per(j, 5) = length(find((thti(: , 1) == 0 & thti(: , 3) == 0 & thti(: , 5) == 0 & thti(: , 7) == 0 & thti(: , 9) ~= 0))) / length(its);
% % end
% % 
% % for j = 1:length(periodicity_list)
% %     set(j, :) = [periodicity_list(j), find(simple_per(j, :) == max(simple_per(j,:)))];
% % end
% % per4_found = set(find(set(:, 2) == 1), 1);
% % per6_found = set(find(set(:, 2) == 2), 1);
% % per8_found = set(find(set(:, 2) == 3), 1);
% % per12_found = set(find(set(:, 2) == 4), 1);
% % per24_found = set(find(set(:, 2) == 5), 1);
% % 
% % [length(find(THETA(per4_found, 1) ~= 0 & THETA(per4_found, 3) == 0 & THETA(per4_found, 5) == 0 &...
% %     THETA(per4_found, 7) == 0 & THETA(per4_found, 9) == 0)), length(per4_found)]
% % 
% % [length(find(THETA(per6_found, 1) == 0 & THETA(per6_found, 3) ~= 0 & THETA(per6_found, 5) == 0 &...
% %     THETA(per6_found, 7) == 0 & THETA(per6_found, 9) == 0)), length(per6_found)]
% % 
% % [length(find(THETA(per8_found, 1) == 0 & THETA(per8_found, 3) == 0 & THETA(per8_found, 5) ~= 0 &...
% %     THETA(per8_found, 7) == 0 & THETA(per8_found, 9) == 0)), length(per8_found)]
% % 
% % [length(find(THETA(per12_found, 1) == 0 & THETA(per12_found, 3) == 0 & THETA(per12_found, 5) == 0 &...
% %     THETA(per12_found, 7) ~= 0 & THETA(per12_found, 9) == 0)), length(per12_found)]
% % 
% % [length(find(THETA(per24_found, 1) == 0 & THETA(per24_found, 3) == 0 & THETA(per24_found, 5) == 0 &...
% %     THETA(per24_found, 7) == 0 & THETA(per24_found, 9) ~= 0)), length(per24_found)]
% % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % Search by fixing a max to the FDR % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % thr = 0.05;                                  % Number of false discoveries we are willing to tolerate
% % kappa = .10;
% % 
% % exp_fd = (list(:,2)'*(list(:,2) <= kappa));  % Expected number of false discoveries
% % D = length(find(list(:,2) <= kappa))
% % fdr = exp_fd / D ;                          % Expected rate of false discoveries: needs to be <= thr;
% % periodicity_list = ord(find(list(:, 2) <= kappa));
% % 
% % % For the proteins most likely to be periodic, compute the probability of
% % % simple periodicity
% % for j = 1:length(periodicity_list)
% %     h = periodicity_list(j);
% %     thti = thetaout(its, indicator + h); 
% %     simple_per(j, 1) = length(find((thti(: , 1) ~= 0 & thti(: , 3) == 0 & thti(: , 5) == 0 & thti(: , 7) == 0 & thti(: , 9) == 0))) / length(its);
% %     simple_per(j, 2) = length(find((thti(: , 1) == 0 & thti(: , 3) ~= 0 & thti(: , 5) == 0 & thti(: , 7) == 0 & thti(: , 9) == 0))) / length(its);
% %     simple_per(j, 3) = length(find((thti(: , 1) == 0 & thti(: , 3) == 0 & thti(: , 5) ~= 0 & thti(: , 7) == 0 & thti(: , 9) == 0))) / length(its);
% %     simple_per(j, 4) = length(find((thti(: , 1) == 0 & thti(: , 3) == 0 & thti(: , 5) == 0 & thti(: , 7) ~= 0 & thti(: , 9) == 0))) / length(its);
% %     simple_per(j, 5) = length(find((thti(: , 1) == 0 & thti(: , 3) == 0 & thti(: , 5) == 0 & thti(: , 7) == 0 & thti(: , 9) ~= 0))) / length(its);
% % end
% % 
% % for j = 1:length(periodicity_list)
% %     set(j, :) = [periodicity_list(j), find(simple_per(j, :) == max(simple_per(j,:)))];
% % end
% % per4_found = set(find(set(:, 2) == 1), 1);
% % per6_found = set(find(set(:, 2) == 2), 1);
% % per8_found = set(find(set(:, 2) == 3), 1);
% % per12_found = set(find(set(:, 2) == 4), 1);
% % per24_found = set(find(set(:, 2) == 5), 1);
% % 
% % [length(find(THETA(per4_found, 1) ~= 0 & THETA(per4_found, 3) == 0 & THETA(per4_found, 5) == 0 &...
% %     THETA(per4_found, 7) == 0 & THETA(per4_found, 9) == 0)), length(per4_found)]
% % 
% % [length(find(THETA(per6_found, 1) == 0 & THETA(per6_found, 3) ~= 0 & THETA(per6_found, 5) == 0 &...
% %     THETA(per6_found, 7) == 0 & THETA(per6_found, 9) == 0)), length(per6_found)]
% % 
% % [length(find(THETA(per8_found, 1) == 0 & THETA(per8_found, 3) == 0 & THETA(per8_found, 5) ~= 0 &...
% %     THETA(per8_found, 7) == 0 & THETA(per8_found, 9) == 0)), length(per8_found)]
% % 
% % [length(find(THETA(per12_found, 1) == 0 & THETA(per12_found, 3) == 0 & THETA(per12_found, 5) == 0 &...
% %     THETA(per12_found, 7) ~= 0 & THETA(per12_found, 9) == 0)), length(per12_found)]
% % 
% % [length(find(THETA(per24_found, 1) == 0 & THETA(per24_found, 3) == 0 & THETA(per24_found, 5) == 0 &...
% %     THETA(per24_found, 7) == 0 & THETA(per24_found, 9) ~= 0)), length(per24_found)]
% % 


vect = 0:.0005:1; 
true = 112; % This number has to be verified
sensitivity=zeros(1,length(vect));
spec=zeros(1,length(vect));
for h = 1:length(vect)
 kappa = vect(h);
 circa = ord(find(list(:, 2) <= kappa));     % Proteins identified as periodic
 sensitivity(h) = length(find((THETA(circa , 1) ~= 0 & THETA(circa , 3) == 0 & THETA(circa , 5) == 0 & THETA(circa , 7) == 0 & THETA(circa , 9) == 0) | ...
        (THETA(circa , 1) == 0 & THETA(circa , 3) ~= 0 & THETA(circa , 5) == 0 & THETA(circa, 7) == 0 & THETA(circa , 9) == 0) | ...
        (THETA(circa , 1) == 0 & THETA(circa , 3) == 0 & THETA(circa , 5) ~= 0 & THETA(circa, 7) == 0 & THETA(circa , 9) == 0) | ...
        (THETA(circa , 1) == 0 & THETA(circa , 3) == 0 & THETA(circa , 5) == 0 & THETA(circa, 7) ~= 0 & THETA(circa , 9) == 0) | ...
        (THETA(circa , 1) == 0 & THETA(circa , 3) == 0 & THETA(circa , 5) == 0 & THETA(circa, 7) == 0 & THETA(circa , 9) ~= 0))) / true; 

 notcirc = setdiff(1:p, circa);             % Proteins not identified as circardian
 true_not_circ = p - true;
 spec(h) = (length(notcirc) - length(find((THETA(notcirc , 1) ~= 0 & THETA(notcirc , 3) == 0 & THETA(notcirc , 5) == 0 & THETA(notcirc , 7) == 0 & THETA(notcirc , 9) == 0) | ...
        (THETA(notcirc , 1) == 0 & THETA(notcirc , 3) ~= 0 & THETA(notcirc , 5) == 0 & THETA(notcirc, 7) == 0 & THETA(notcirc , 9) == 0) | ...
        (THETA(notcirc , 1) == 0 & THETA(notcirc , 3) == 0 & THETA(notcirc , 5) ~= 0 & THETA(notcirc, 7) == 0 & THETA(notcirc , 9) == 0) | ...
        (THETA(notcirc , 1) == 0 & THETA(notcirc , 3) == 0 & THETA(notcirc , 5) == 0 & THETA(notcirc, 7) ~= 0 & THETA(notcirc , 9) == 0) | ...
        (THETA(notcirc , 1) == 0 & THETA(notcirc , 3) == 0 & THETA(notcirc , 5) == 0 & THETA(notcirc, 7) == 0 & THETA(notcirc , 9) ~= 0)))) / true_not_circ;       
end
figure;
plot(1 - spec, sensitivity)
xlim([0, 1])
xlabel(['False positive rate'])
ylabel(['True positive rate'])
title(['ROC curve'])
dlmwrite('ROC.txt', [(1-spec)' sensitivity'], 'delimiter', ' ', 'precision', 6);

%%% PLOT TRAIETTORIE GENI IDENTIFICATI COME CIRCADIANI CON MAGGIORE
%%% PROBABILITA'
figure;
subplot(2,2,1)
plot(Y(:,ordC(1))) 
subplot(2,2,2)
plot(Y(:,ordC(2))) 
subplot(2,2,3)
plot(Y(:,ordC(3))) 
subplot(2,2,4)
plot(Y(:,ordC(4))) 

%%% PLOT TRAIETTORIE GENI CON BASSE PROBABILITA' DI ESSERE CIRCADIANI
figure;
subplot(2,2,1)
plot(Y(:,ordC(497))) 
subplot(2,2,2)
plot(Y(:,ordC(498))) 
subplot(2,2,3)
plot(Y(:,ordC(499))) 
subplot(2,2,4)
plot(Y(:,ordC(500))) 


writematrix(ordC, 'Indici_ordinati_geni_circadiani_synthetic.txt');
writematrix(A_C, 'probabilities_circadian_synthetic.txt');
