%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Circardian project (Synthetic data)           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This simulations generates data under the model and update
% all parameters but the THETA thresholds. 

clear all; clc;
rng(5000)                                   % Set seed number
p = 500;                                    % Number of simulated protein profiles ( = number of genes)
T = 24;                                     % Number of observations per profile 
q = 5;                                      % Number of sin (cos) bases for a total of 2q bases -- no intercept
tij = [0:2:46]';                            % Time points at which the probes have been measured
m = max(tij) ;                              % Max hour of measurement
tij = tij / max(tij);                       % Standardize for stability in analyzing data and ease in prior specification
tg = [0:0.1:46]'/ 46;                       % Fine grid of equally spaced measurement times for prediction

B = zeros(T, 2*q);                          % Design matrix for data in study
Bpred = zeros(length(tg), 2*q);             % Design matrix for estimation / prediction
lambda2 = [8, 12, 16, 24, 48];              % Range of periods for fixed basis functions (based on Fourier transforms for 2t)
periods = lambda2./2;                       % Periods on a unit-based increment
lambda = periods./46;                       % Periods for our time increments

for h = 1 : length(lambda)
    B(:, 2*h - 1) = sin(((2*pi)/(lambda(h))).*tij);
    Bpred(:, 2*h - 1) = sin(((2*pi)/(lambda(h))).*tg);
    B(:, 2*h) = cos(((2*pi)/(lambda(h))).*tij);
    Bpred(:, 2*h) = cos(((2*pi)/(lambda(h))).*tg);
end

% % -- % Generate true values for the model parameters % -- %
% -- % Generate matrix of factor loadings Lambda p \times k % -- %
% The number of non-zero elements in each column of Lambda  are chosen
% linearly between 2k and k +1 in a decreasing fashion. We randomly 
% allocate the location of the zeros in each column and simulate the 
% nonzero elements independently from a normal distribution with mean 0 
% and variance 9.
k = floor(log(p)*10);                                   % True number of factors
numeff = k + randperm(k);                               % Number of non-zero entries in each column of Lambda
numeff = numeff';                                       %G we transpose so all our vectors are columns
k = floor(log(q)*4);                                    %G arbitrary value, should be small (<10, >3)
Lambda = zeros(p,k);
% The number of non-zero entries varies between 63 and 114. Not many
% genes are correlated, but those that are are also strongly correlated
for h = 1:k
    temp = randsample(p, numeff(h)); 
    %G p = 500, numeff(h) is the number of non-null entries in that column e.g.
    %numeff(1) = 86 which means that in column 86 I will have 86 values !=0
    Lambda(temp,h) = normrnd(0,3,[numeff(h),1]);
    %G and then here I'm just filling the h-th column of lambda in 86
    %positions, sampling from a normal with mean 0 and sd 3 (as written in
    %paper (section 5.1))
end

%Omega = Lambda*Lambda';  

% -- % Generate low dimensional matrices W and Z % -- %
W = normrnd(0, 1, [2*q k]); % Z = normrnd(0, 1, [Ttilde k]);
Thetatilde = Lambda*W' + normrnd(0, 1, [p 2*q]);                   % Matrix of unshrunk coefficients
THETA = zeros(p, 2*q);                                             % Matrix of (fixed) basis functions coefficients
thr1 = unifrnd(0, 6, [p, q]);                                      % Matrix of thresholds for THETA (Ktheta = 6)
%G also here parameter of the uniform is arbitrary but must not be too big 

for i = 1:q
  index = find(bsxfun(@hypot, Thetatilde(: , 2*i-1), Thetatilde(:, 2*i)) >= thr1(:, i));
  %G here I'm saving the indexes of which values of thetatilde are above
  %the threshold I have defined above
  if length(index) ~= 0
  THETA(index, [2*i - 1 2*i]) = Thetatilde(index, [2*i - 1 2*i]); 
  %G I save in THETA only the values of Thetatilde that are above the
  %threshold I have saved, by using the indexes to extract from Thetatilde
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following counts how many genes exhibit 24/12/8/6/4 hours periodicit
% Can modify if we want more (or less): these numbers will change with rand
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%G so here we have built the matrix THETA s.t. genes with periodicity 24
%hrs are non-null in the last two columns, those with periodicity 12 in the
%penultimate two columns and so on and so forth. This makes it easy to
%extract but we're not fully sure why identifying them like these means
%that they're circadian with that period 

per24 = find(THETA(: , 1) == 0 & THETA(: , 2) == 0 & THETA(: , 3) == 0 & THETA(: , 4) == 0 ...
    & THETA(:, 5) == 0 & THETA(:, 6) == 0 & THETA(:, 7) == 0 & THETA(:, 8) == 0 ...
    & THETA(:, 9) ~= 0 & THETA(:, 10) ~= 0) ; length(per24)    % Period 24 -- 33 proteins
per12 = find(THETA(: , 1) == 0 & THETA(: , 2) == 0 & THETA(: , 3) == 0 & THETA(: , 4) == 0 ...
    & THETA(:, 5) == 0 & THETA(:, 6) == 0 & THETA(:, 7) ~= 0 & THETA(:, 8) ~= 0 ...
    & THETA(:, 9) == 0 & THETA(:, 10) == 0) ; length(per12)    % Period 12 -- 25 proteins
per8 = find(THETA(: , 1) == 0 & THETA(: , 2) == 0 & THETA(: , 3) == 0 & THETA(: , 4) == 0 ...
    & THETA(:, 5) ~= 0 & THETA(:, 6) ~= 0 & THETA(:, 7) == 0 & THETA(:, 8) == 0 ...
    & THETA(:, 9) == 0 & THETA(:, 10) == 0) ; length(per8)     % Period 8 -- 23 proteins
per6 = find(THETA(: , 1) == 0 & THETA(: , 2) == 0 & THETA(: , 3) ~= 0 & THETA(: , 4) ~= 0 ...
    & THETA(:, 5) == 0 & THETA(:, 6) == 0 & THETA(:, 7) == 0 & THETA(:, 8) == 0 ...
    & THETA(:, 9) == 0 & THETA(:, 10) == 0) ; length(per6)     % Period 6 -- 23 proteins
per4 = find(THETA(: , 1) ~= 0 & THETA(: , 2) ~= 0 & THETA(: , 3) == 0 & THETA(: , 4) == 0 ...
    & THETA(:, 5) == 0 & THETA(:, 6) == 0 & THETA(:, 7) == 0 & THETA(:, 8) == 0 ...
    & THETA(:, 9) == 0 & THETA(:, 10) == 0) ;  length(per4)    % Period 4 -- 13 proteins
periodic_vars = zeros(length(per24) , q);
periodic_vars(1:length(per4), 1) = per4;
periodic_vars(1:length(per6), 2) = per6;
periodic_vars(1:length(per8), 3) = per8;
periodic_vars(1:length(per12), 4) = per12;
periodic_vars(1:length(per24), 5) = per24;
%G here I save all the periodic genes I have identified, each column
%corresponds to a periodic interval so col1 = 4 hours, col2 = 6 hours,
%col3=8hours, col4 = 12 hours and col5 = 24 hours
%dlmwrite('periodicity.txt', periodic_vars, 'delimiter', ' ')
%G here I save all the genes with just a period of 24 hours 
%dlmwrite("circadian.txt", per24, 'delimiter', ' ')

% -- % Generate matrix of latent factors eta w dim T \times k % -- %
% The rows of eta correspond to the set of latent factors at time j's, thus
% the mouse-specific latent factors. We place a Gaussian prior on each row
% of eta
eta = zeros(T, k);
eta = normrnd(0, 1, [T, k]);

% Use sig_i^2 = 0.5 for every curve -- so cov becomes 0.5*eye(T) for each
% protein

Y = zeros(T, p);
for i = 1:p
    Y(:,i) = mvnrnd(B*THETA(i,:)' + eta*Lambda(i, :)', 0.5*eye(T)); % C*Gamma(i,:)' 
end

writematrix(Y, 'Y.csv')

%G this is the real generative step where we make the synthetic data using
%a normal centered in B*THETA(i,:)' + eta*Lambda(i, :)' and with variance
%0.5*I (this is done to introduce some noise, but the coefficient should
%not be too big, otherwise it becomes difficult to recognise periodicity)

% -- % Plot a few trajectories % -- %
ib = randsample(1:p, 6, 1);
cnt = 0;
for h = 1:25 % Visualise the trajectories of the first 25 circadian genes (lower if less than 25)
    i = per24(h);
    %i = ib(h);
    cnt = cnt + 1';
    subplot(5, 5, h)
    plot(tij.*46, Y(:, i),'o', 'Color', 'k', 'MarkerSize',4)
    line(tij.*46, Y(:, i), 'LineStyle', '-','Color', 'k','LineWidth', 1.5)
    xlim([0, 46])
    title(['True trajectory for protein ' num2str(i) ])
    xlabel(['Time (in hours)'])
%     if sum(Z(i, :)) == 1
%         line(tij.*46, B*THETA(:, i) + eta*Lambda(i, :)', 'Color', 'blue')
%         line(tij.*46, B*THETA(:, i), 'Color', 'red')
%         per = find(Z(i, :) == 1);
%         if per == 1
%         line([periods(1), periods(1)], [-10, 10]); line([periods(1)*2, periods(1)*2], [-10, 10]); line([periods(1)*3, periods(1)*3], [-10, 10])
%         end
%         if per == 2
%         line([periods(2), periods(2)], [-10, 10]); line([periods(2)*2, periods(2)*2], [-10, 10]); line([periods(2)*3, periods(2)*3], [-10, 10])
%         end
%         if per == 3
%         line([periods(3), periods(3)], [-10, 10]); line([periods(3)*2, periods(3)*2], [-10, 10]); line([periods(3)*3, periods(3)*3], [-10, 10])
%         end
%         if per == 4
%         line([periods(4), periods(4)], [-10, 10]); line([periods(4)*2, periods(4)*2], [-10, 10]); line([periods(4)*3, periods(4)*3], [-10, 10])
%         end
%         if per == 5
%         line([periods(5), periods(5)], [-10, 10]); line([periods(5)*2, periods(5)*2], [-10, 10]); 
%         end
%    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot one (random) circadian trajectory %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is to see how the different components of the model
% "modify" the trajectory
%figure()
%out = per24(end);

% First set of plots
%h1 = plot(tij.*46, B*THETA(out,:)', 'o', 'Color', 'red', 'MarkerSize',4);
%hold on
%line(tij.*46, B*THETA(out,:)', 'LineStyle', '-', 'Color', 'red', 'LineWidth', 1.5);

% Second set of plots
%h2 = plot(tij.*46, B*THETA(out,:)' + eta*Lambda(out,:)', 'o', 'Color', 'black', 'MarkerSize',4);
%line(tij.*46, B*THETA(out,:)' + eta*Lambda(out,:)', 'LineStyle', '-.', 'Color', 'black', 'LineWidth', 3);

% Third set of plots
%h3 = plot(tij.*46, B*THETA(out,:)' + eta*Lambda(out,:)' + error, 'o', 'Color', 'blue', 'MarkerSize',4);
%line(tij.*46, B*THETA(out,:)' + eta*Lambda(out,:)' + error, 'LineStyle', ':', 'Color', 'blue', 'LineWidth', 3);

%ylim([-8,8]); xlim([0, 46])
%xlabel(['Time (hours)'], 'FontSize', 24); ylabel(['Normalized Expression Level'], 'FontSize', 24)
%title('Protein 487', 'FontSize', 24)
%legend([h1, h2, h3], {'Fourier basis', 'Fourier basis + latent factor', 'Fourier basis + latent factor + error'}, 'FontSize', 24);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Visualise the basis functions %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Basis functions with 24 hours period %
plot(unique(tg).*46, Bpred(:, 10), 'Color', 'red')
hold on
plot(unique(tij).*46, B(: , 10), 'Color', 'blue')
line([24, 24], [-1, 1])
line([48, 48], [-1, 1])
plot(unique(tg).*46, Bpred(:, 9), 'Color', 'green')
hold on
plot(unique(tij).*46, B(: , 9), 'Color', 'black')
% 
% % Basis functions with 12 hours period %
% plot(unique(tg).*46, Bpred(:, 8), 'Color', 'red')
% hold on
% plot(unique(tij).*46, B(: , 8), 'Color', 'blue')
% line([12, 12], [-1, 1]); line([24, 24], [-1, 1]); line([36, 36], [-1, 1]); line([48, 48], [-1, 1])
% plot(unique(tg).*46, Bpred(:, 7), 'Color', 'green')
% hold on
% plot(unique(tij).*46, B(: , 7), 'Color', 'black')
% 
% % Basis functions with 8 hours period %
% plot(unique(tg).*46, Bpred(:, 6), 'Color', 'red')
% hold on
% plot(unique(tij).*46, B(: , 6), 'Color', 'blue')
% line([8,8], [-1, 1]); line([16, 16], [-1, 1]); line([24,24], [-1, 1]); line([32,32], [-1, 1])
% plot(unique(tg).*46, Bpred(:, 5), 'Color', 'green')
% hold on
% plot(unique(tij).*46, B(: , 5), 'Color', 'black')
% 
% % Basis functions with 6 hours period %
% plot(unique(tg).*46, Bpred(:, 4), 'Color', 'red')
% hold on
% plot(unique(tij).*46, B(: , 4), 'Color', 'blue')
% line([6,6], [-1, 1]); line([12, 12], [-1, 1]); line([18,18], [-1, 1]); line([24,24], [-1, 1])
% plot(unique(tg).*46, Bpred(:, 3), 'Color', 'green')
% hold on
% plot(unique(tij).*46, B(: , 3), 'Color', 'black')
% 
% % Basis functions with 4 hours period %
% plot(unique(tg).*46, Bpred(:, 2), 'Color', 'red')
% hold on
% plot(unique(tij).*46, B(: , 2), 'Color', 'blue')
% line([4,4], [-1, 1]); line([8, 8], [-1, 1]); line([12,12], [-1, 1]); line([16,16], [-1, 1]); line([20,20], [-1, 1])
% plot(unique(tg).*46, Bpred(:, 1), 'Color', 'green')
% hold on
% plot(unique(tij).*46, B(: , 1), 'Color', 'black')

save('Y')

% save(strcat('DesignMatrix'),'Y', 'tij', 'tg', 'p', 'q', 'T', 'B', 'eta', ...
%     'Lambda', 'THETA', 'Thetatilde', 'W', ...
%     'thr1', 'numeff', 'Bpred')
% 
% save(strcat('DesignMatrix'),'Y', 'tij', 'tg', 'p', 'q', 'T', 'Ttilde', 'B', 'C', 'eta', ...
%     'Lambda', 'THETA', 'Gamma', 'Thetatilde', 'Gammatilde','W', 'Z', ...
%     'thr1', 'thr2','numeff', 'Bpred')
