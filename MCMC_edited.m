%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Code for implementing circardian project with MGPSP % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main modifications wrt previous versions: update of the thresholds is 
% done via Gibbs step; update of latent parameter remains a MH step

clear; clc; 
warning off all;
tic;
load DesignMatrix;
% Set seed number
rng(5000)

% Need to start off normalizing the (real) data
 
p = size(Y, 2);                             % p = number of proteins
T = size(Y, 1);                             % T = number of time points
q = 5;                                      % Number of sin (cos) bases for a total of 2q bases -- no intercept
tij = [0:2:46]'/46;                         % Time points at which the probes have been measured
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

% % --- Define global constants --- % %
% ktr = 20; %G commented because it doesn't appear to be used anywhere else
rep = 1;
nrun = 10000;
burn = 1000;
thin = 5;
sp = (nrun - burn)/thin;                    % Number of posterior samples
% kinit = repmat(floor(log(p)*4),rep,1);    % Number of factors to start with (number of columns of Lambda)
k = 5;                                      % Number of factors to start with (for now)

b0 = 1; b1 = 0.0005;                        %G parameters for exponential - assumed to be arbitrary
epsilon = 1e-3;                             % Threshold limit
prop = 1.00;                                % Proportion of redundant elements within columns

% % --- Define hyperparameter values --- % %
as = 1; bs = 0.5;                           % Gamma hyperparameters for residual precision (true value res variance = 1 for every i)
df = 3;                                     % Gamma hyperparameters for t_{ij} %C for phi_{ij} (i.e. rho)?
ad1 = 2.1; bd1 = 1;                         % Gamma hyperparameters for delta_1
ad2 = 3.1; bd2 = 1;                         % gamma hyperparameters delta_h, h >= 2
adf = 1; bdf = 1;                           % Gamma hyperparameters for ad1 and ad2 or df

% % --- Initial values --- % %
sig = gamrnd(as, 1/bs, p, 1);                     % Residual variance (diagonal of sig^(-2)_i across proteins)
    
odd = []; even = [];
for i = 1:q
  odd = [odd, 2*i - 1];                           % Index for sin bases
  even = [even, 2*i];                             % Index for cos bases
end      

Lambda = zeros(p,k); % Factor loadings                          
eta = mvnrnd(zeros(T, k), eye(k));                % Latent factors (distrib. = mean 0 and identity cov matrix
W = mvnrnd(zeros(2*q, k), eye(k))'; % Low dim. matrix W
Thetatilde = Lambda*W;                            % Matrix of unshrunk coefficients
THETA = zeros(p, 2*q);                            % Matrix of (fixed) basis functions coefficients
Kappatheta = 5;                                   % Upper bound on the Unif prior on the thresholds (can modify this number as pleased)  
thr1 = unifrnd(0, Kappatheta, [p, q]);  % Matrix of thresholds for THETA

disp("Theta at the beginning:");
THETA

for i = 1:q
  index = find(bsxfun(@hypot, Thetatilde(: , 2*i-1), Thetatilde(:, 2*i)) >= thr1(:, i)); 
  index
  if length(index) ~= 0
  THETA(index, [2*i - 1 2*i]) = Thetatilde(index, [2*i - 1 2*i]);  
  end
end

disp("Theta at the end:");
THETA


phiih = gamrnd(df/2, 2/df, [p,k]);                         % Local shrinkage coefficients
delta = [gamrnd(ad1,bd1); gamrnd(ad2, bd2, [k-1,1])];      % Global shrinkage coefficients multilpliers
tauh = cumprod(delta);                                     % Global shrinkage coefficients
Plam = bsxfun(@times,phiih,tauh');                         % Precision of loadings rows %C i.e. 1/D
% This matrix stores the rejection probabilities 1 - mean(Accept) is the accept prob of the MH steps
% Acc1 = zeros(nrun, p); 
Acc2 = zeros(nrun, p);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % I will set and keep to ZERO all matrices related to the local bases     %
% % so they have no impact on the posterior update of the other parameters  %
% % and I will comment out below their updates %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- % Define output files % --- %
fidtheta = fopen('thetaout.txt', 'wt');
fidthetatilde = fopen('thetatildeout.txt', 'wt');
fidsigma = fopen('sigmaout.txt', 'wt');
fidthr1 = fopen('trh1out.txt', 'wt');

% --- % Initialize matrices of output files % --- %
thetaout = zeros(1,(2*q)*p);                % Thetatilde uses same string of theta
sigmaout = zeros(1,p);
thr1out = zeros(1, p*q);

type = '%3.3f\t';
stringa = '%3.3f\t';
for i = 2:length(thetaout) - 1
    stringa = strcat(stringa, ' ', type);
end
stringa = strcat(stringa, ' ' , '%3.3f\r\n');

stringa2 = '%3.3f\t';
for i = 2:length(sigmaout) - 1
    stringa2 = strcat(stringa2, ' ', type);
end
stringa2 = strcat(stringa2, ' ' , '%3.3f\r\n');

stringa4 = '%3.3f\t';
for i = 2:length(thr1out) - 1
    stringa4 = strcat(stringa4, ' ', type);
end
stringa4 = strcat(stringa4, ' ' , '%3.3f\r\n');

clear thetaout; clear sigmaout; clear thr1out; 

Lambdaout = zeros(p*15, nrun/thin);
Factorout = zeros(1, nrun/thin);
Kthetaout = zeros(1, nrun /thin);
Etaout = zeros(T*15, nrun/thin);
Wout = zeros(2*q*15, nrun/thin);

%  --- % Start Gibbs sampling % --- %
for i = 1:nrun
   
  % -- % Update error precisions % -- %
  Ytil = Y - B*THETA' - eta*Lambda'; % - C*Gamma'
  sig = gamrnd(as + 0.5 * T, 1./(bs+0.5*sum(Ytil.^2)))'; 
  %G alternative parametrization of Gamma
  
  % -- % Update eta % -- %        
  Lmsg = bsxfun(@times,Lambda,sig);
  Veta1 = eye(k) + Lmsg'*Lambda;                    % This is the cov matrix as in the paper
  Tchol = cholcov(Veta1); [Q,R] = qr(Tchol);
  S = inv(R); Veta = S*S';                          % Veta = inv(Veta1)
  Meta = ((Y - B*THETA')*Lmsg)*Veta;                % T x k REMOVED:  - C*Gamma'
  %G we already multiply by Veta since this parameter will be the mean of
  %etaj
  eta = Meta + normrnd(0,1,[T,k])*S';               % Update eta in a block
  %G faster and more stable way to create a normal with our desired
  %parameters

  % -- % Update of Lambda (rue & held) % -- %
  for h = 1:p
     Veta1 = sig(h)*eta'*eta + W*eye(2*q)*W' + diag(Plam(h,:));          % This is the cov matrix as in the paper; REMOVED: Z*eye(Ttilde)*Z'+
     % This should be V_lamda
     Tchol = cholcov(Veta1); [Q,R] = qr(Tchol);
     S = inv(R); Veta = S*S';                                             % Veta = inv(Veta1)
     Meta = (sig(h)*eta'*(Y(:,h) - B*THETA(h,:)') + ...                   % REMOVED: - C*Gamma(h,:)'
         W * eye(2*q) * Thetatilde(h,:)')' * Veta;    % 1 x k             % REMOVED:  + Z * eye(Ttilde) * Gammatilde(h,:)'  
     % This should be M_lamda
     Lambda(h,:) = Meta + normrnd(0,1,[1,k])*S'; 
  end

 % -- % Update phi_{ih}'s % -- %
  phiih = gamrnd(df/2 + 0.5, 1./(df/2 + bsxfun(@times,Lambda.^2,tauh')));
  %G check if 1/2 is missing or incorporated in bsxfun

  % -- % Update delta & tauh % -- %
  mat = bsxfun(@times,phiih,Lambda.^2);
  ad = ad1 + 0.5*p*k; bd = bd1 + 0.5*(1/delta(1))*sum(tauh.*sum(mat)');
  %G to remove the term zeta(h) (or better delta(h)) from the summation, we
  % multiply the total sum by 1/delta(h)
  delta(1) = gamrnd(ad,1/bd); %G delta = zeta in the paper
  tauh = cumprod(delta);

  for h = 2:k
    ad = ad2 + 0.5*p*(k-h+1); bd = bd2 + 0.5*(1/delta(h))*sum(tauh(h:end).*sum(mat(:,h:end))');
    delta(h) = gamrnd(ad,1/bd); tauh = cumprod(delta);
  end

 % -- % Update precision parameters % -- %
  Plam = bsxfun(@times,phiih,tauh');
  
 % -- % Update matrix W % -- %
 Veta1 = Lambda'*Lambda + eye(k);
 Tchol = cholcov(Veta1); [Q,R] = qr(Tchol);
 S = inv(R); Veta = S*S';
 for h = 1:2*q
   Meta = (Lambda' * Thetatilde(:, h))' * Veta; 
   W(: , h) = Meta + (normrnd(0,1,[1,k])*S');
 end
 
  % -- % Update of Thetatilde % -- %
 for h = 1:p
     Veta1 = sig(h) * B' * B + eye(2*q);            %G M_i in the paper
     Tchol = cholcov(Veta1); [Q,R] = qr(Tchol);
     S = inv(R); Veta = S*S';                       % Veta = inv(Veta1) 
     Ycent = Y(:, h) - eta*Lambda(h,:)';            % REMOVED: - C*Gamma(h,:)'
     Wtrans = W';
     Meta = (sig(h) * B' * Ycent + Wtrans * Lambda(h,:)')'*Veta;    %G m_i in the paper
     Thetatildeprop = Meta + normrnd(0,1, [1 2*q])*S';  % Generated Gammatilde_i from proposal distribution
     thetaprop = zeros(1, 2*q);
     index = find(bsxfun(@hypot, Thetatildeprop(odd)', Thetatildeprop(even)')' >= thr1(h, :)); 
     if length(index) ~= 0
      ind = sort([odd(index) even(index)]);
      thetaprop(ind) = Thetatildeprop(ind);
     end
     
     r = exp(-0.5*(Thetatildeprop' - Wtrans*Lambda(h,:)')'*(Thetatildeprop' - Wtrans*Lambda(h,:)') + ...
         0.5 * (Thetatilde(h,:)' - Wtrans*Lambda(h,:)')'*(Thetatilde(h,:)' - Wtrans*Lambda(h,:)') - ...
         0.5 * sig(h) * (Y(:, h) - B*thetaprop' - eta*Lambda(h,:)')' * (Y(:, h) - B*thetaprop' - eta*Lambda(h,:)') + ...
         0.5 * sig(h) * (Y(:, h) - B*THETA(h,:)' - eta*Lambda(h,:)')' * (Y(:, h) - B*THETA(h,:)' - eta*Lambda(h,:)') - ...
         0.5 * (Thetatilde(h,:)' - Meta')'*Veta1*(Thetatilde(h,:)' - Meta') + 0.5 * (Thetatildeprop' - Meta')'*Veta1*(Thetatildeprop' - Meta'));
     
     u = (rand(1) > min(1, r));  Acc2(i, h) = u;
     %G u = 0 --> accept, u = 1 --> reject
     thetacc = Thetatildeprop - (Thetatildeprop - Thetatilde(h,:)) * u;
     Thetatilde(h,:) = thetacc;
     THETA(h, :) = zeros(1, 2*q);
     index = find(bsxfun(@hypot, Thetatilde(h, odd)', Thetatilde(h,even)')' >= thr1(h, :)); 
     if length(index) ~= 0
      ind = sort([odd(index) even(index)]);
      THETA(h,ind) = Thetatilde(h,ind);
     end
 end
 
  % -- % Update of thresholds on THETA (thr1) % -- %
  for h = 1:q
        AA = reshape(B(:, setdiff(1: (2*q), 2*h-1:2*h)) * THETA(:, setdiff(1:2*q, 2*h-1:2*h))', T *p, 1);
        BB = AA + reshape(B(:, 2*h -1 : 2*h) * Thetatilde(:, 2*h -1 : 2*h)', T*p, 1);
        comp1 = bsxfun(@hypot, Thetatilde(:, 2*h - 1), Thetatilde(:, 2*h));
        Yprime = reshape(Y  - eta*Lambda', p*T, 1);         % REMOVED: - C*Gamma'
        nn = diag(reshape(Yprime - BB, T, p)' * reshape(Yprime - BB, T, p));
        Num = bsxfun(@times, exp(- 0.5 * bsxfun(@times, sig, nn)), comp1);
        dd = diag(reshape(Yprime - AA, T, p)' * reshape(Yprime - AA, T, p));
        Den = Num + bsxfun(@times, exp(- 0.5 * bsxfun(@times, sig, dd)), (Kappatheta - comp1));
        short2 = (Num ./ Den)';
        pc = find(bsxfun(@gt, comp1, Kappatheta));
         if length(pc) ~= 0
             thr1(pc, h) = unifrnd(0, Kappatheta, [length(pc), 1]);
         end
         npc = setdiff(1:p, pc);
         u = rand([1, length(npc)]);
         thr1(npc, h) = bsxfun(@le, u, short2(npc)).*unifrnd(0, comp1(npc)') + bsxfun(@gt, u, short2(npc)).*unifrnd(comp1(npc)', Kappatheta);
         THETA(:, 2*h-1 : 2*h) = zeros(p,2);
         index = find(comp1 >= thr1(: , h)); 
         if length(index) ~= 0
            THETA(index, 2*h - 1:2*h) = Thetatilde(index, 2*h-1:2*h);
         end
  end
  clear AA; clear BB;
   
 % -- % Make adaptations % -- %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % This part adapts the number of latent factors at the end of %
 %                      each run of the MCMC                   %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 prob = 1/exp(b0 + b1*i);                    % Probability of adapting
 uu = rand;
 lind = sum(abs(Lambda) < epsilon)/p;        % Proportion of elements in each column less than eps in magnitude
 vec = lind >=prop; num = sum(vec);          % Number of redundant columns

 if uu < prob
   if  i > 20 && num == 0 && all(lind < 0.995)
   k = k + 1;
   Lambda(:,k) = zeros(p,1);
   eta(:,k) = normrnd(0,1,[T,1]);  
   W(k, :) = normrnd(0,1,[1, 2*q]);
   phiih(:,k) = gamrnd(df/2,2/df,[p,1]);
   delta(k) = gamrnd(ad2,1/bd2);
   tauh = cumprod(delta);
   Plam = bsxfun(@times,phiih,tauh');
   elseif num > 0
   nonred = setdiff(1:k,find(vec));         % Non-redundant loadings columns
   k = max(k - num,1);
   Lambda = Lambda(:,nonred);
   W = W(nonred, :); % Z = Z(nonred, :);
   phiih = phiih(:,nonred);
   eta = eta(:,nonred);
   delta = delta(nonred);
   tauh = cumprod(delta);
   Plam = bsxfun(@times,phiih,tauh');
   end
 end

% --- % Save sampled values (after thinning) % --- %
if mod(i, thin) == 0 
   Lambdaout(1:p*k, i/thin)= reshape(Lambda, p*k,1);
   Etaout(1:T*k, i/thin) = reshape(eta, T*k,1);
   Factorout(1, i/thin) = k;
   Wout(1:(2*q*k), i/thin) = reshape(W, 2*q*k,1);
   thtout = THETA(:); fprintf(fidtheta, stringa, thtout); clear thtout;
   thtout = Thetatilde(:); fprintf(fidthetatilde, stringa, thtout); clear thtout;
   thtout = thr1(:); fprintf(fidthr1, stringa4, thtout); clear thtout;
   fprintf(fidsigma, stringa2, sig);  
   
end

[i, k]
end

fclose(fidtheta); fclose(fidsigma); % fclose(fidgamma);
fclose(fidthetatilde); % fclose(fidgammatilde); 
fclose(fidthr1); % fclose(fidthr2);

dlmwrite('Lambdaout.txt', Lambdaout, 'delimiter', '\t', 'precision', 6);
dlmwrite('Etaout.txt', Etaout, 'delimiter', '\t', 'precision', 6);
dlmwrite('Factorout.txt', Factorout, 'delimiter', '\t', 'precision', 6);
dlmwrite('Wout.txt', Wout, 'delimiter', '\t', 'precision', 6);
dlmwrite('Accept2.txt', Acc2, 'delimiter', '\t', 'precision', 6);

PostLambda = reshape(mean(Lambdaout'), p, 15);
PostEta = reshape(mean(Etaout'), T, 15);
dlmwrite('PostEta.txt', PostEta, 'delimiter', '\t', 'precision', 6);
dlmwrite('PostLambda.txt', PostLambda, 'delimiter', '\t', 'precision', 6);

toc;
% MAYBE NOT NECESSARY IF HAPPY WITH EVERYTHING
save('Workspace500.mat')

