%%% 2025-1-18

%%% We minimizes f(VL) to see if it converges to a point which it is 
%%% NOT a globally optimal point.
%%% When the objective value of the numerically convergence point
%%% fc >> 0, we see it is NOT a globally optimal point.

%%% randprob = true: use randomly generated data 'Lbar' and 'Lhat' (initial starting point for iteration).
%%% 'n' (the number of pi) and 'd' (the dimension of the ambinent space) can also be set. 
%%% randprob = false: load our picked counterexamples data with d=1, n=50 or d=2, n=100. 



%%% The sequence obtained by TR converges to a numerically second-order
%%% stationary point.

clear
profile clear
profile on

%seed = 100;
%rng(seed);
%rngsave = rng('shuffle');
%save('rngfile','rngsave');


%% Initializations
toler = 1e-5;
maxiter = 100;

randprob = true;   % use random data or load counterexamples data
flip_strategy = true;     %% can be set as 'true' only when randprob = true, d=1 or d=2
% when randprob = true, flip_strategy = true：use random data, and use
% 'flip' strategy to get a starting point. 
partial = false;  % use partial EDM data (distances) or not
verboseplot = false;  %plot or not
if randprob
	n = 200;  % number of pts
	d = 2;   % emb. dim
	ntests = 10;  % choose number of random problems to try
else
	ntests = 1;
	d = 2;   % choose 1 or 2
end
ncntexs = 0;  % count number of counterexamples
iter = 0;
iters = zeros(ntests,1);
nlngm = 0;
for ii = 1:ntests  % run ntests random problems
if randprob
% 	k = 3;
% 	n = 10;
%         [centers, Pbarprel] = constrwheel(k,d);  % points of the wheels
% 	baryc = sum(Pbarprel([2 4 7],:))/3;  % one of the opt baryc.
% 	Pbar = [Pbarprel;baryc];   %   n by d
% 	Phat = [Pbarprel;[-.4 .2]];   %   n by d

    
    if flip_strategy == 0
       Pbar = 10*randn(n,d);  %opt Pbar
       Pbar = Pbar-sum(Pbar)/n;
       Phat = 10*randn(n,d);  %starting point
       Phat = Phat-sum(Phat)/n;
       A= ([ones(1,n-1);-eye(n-1)]); 
       [V] = GS(A);
       Lbar = V'*Pbar;
       Lhat = V'*Phat;
    elseif flip_strategy == 1 && d==1%%%%%%%%%%%%%% 
       Pbar = randn(n,d);  %opt Pbar
       Pbar = Pbar - sum(Pbar)/n;
       [m, t] = max(abs(Pbar(:,1)));
       Pbar = Pbar - sum(Pbar)/n; 
       Phat = Pbar; Phat(1) = -Pbar(1);
       Phat = Phat - sum(Phat)/n;
       A= ([ones(1,n-1);-eye(n-1)]); 
       [V] = GS(A);
       Lbar = V'*Pbar;
       Lhat = V'*Phat;
    elseif flip_strategy == 1 %%%%%%%%%%%%%% 
       Lbar = randn(n-1,d);  %opt Pbar
       for jj=1:d
        Lbar(:,jj) = jj^2*Lbar(:,jj);
       end
       [m1, t1] = max(abs(Lbar(:,1)));
       %[m2, t2] = max(abs([Lbar(1:t1-1,1);0;Lbar(t1+1:n-1,1)]));
       Lhat = Lbar;
       Lhat(t1,:) = -Lbar(t1,:);
       %Lhat(t2,:) = -Lbar(t2,:);
       A = ([ones(1,n-1);-eye(n-1)]); 
       [V] = GS(A);
       Pbar = V*Lbar;
       Phat = V*Lhat;
    end   
	%Phat = randn(n,d);   % initial starting point for minimiz
else
	ntests = 1;
	if d == 1
        %%%%%%%%%%%%%%%%%%%% our known counterexample %%%%%%%%%%%%%%%%%%%%
            load ('Lbars');
	        n = 50;
            load ('Llngs');%% obtain llng, which is the lng minimizer
            A = ([ones(1,n-1);-eye(n-1)]); 
            [V] = GS(A);
            Pbar = V*Lbar;
            Phat = V*Llng+rand(n,d);%Phat is set as a perturbation from V*Llng

	elseif d == 2
      load ('lbars2');
      n=100;
      load ('llngs2');
      Lbar(1:n-1,1) = lbar2(1:n-1);
      Lbar(2:n-1,2) = lbar2(n:2*(n-1)-1);
      Llng = Lbar;
      Llng(6,:) = -Lbar(6,:);
      Llng(15,:) = -Lbar(15,:);
      A = ([ones(1,n-1);-eye(n-1)]); 
      [V] = GS(A);
      Pbar = V*Lbar;
      Phat = V*Llng;    
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% The adjacent matrix H %%%%%%%%%%%%%%%%%%%%%%%%%%% 
H = ones(n,n);
if partial
	%H(triu(ones(n),min(d+2,n))==1) = 0;  %Higher sparsity
    H(triu(ones(n),min(max(d+2,n/2),n))==1) = 0; %lower sparsity
    H = H+H';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lbar = V'*Pbar; Lhat = V'*Phat;
fprintf('\n End of the generation of Lbar Lhat Pbar Phat\n');
[Lc] = lngminTR(n,d,Lbar,Lhat,V,H,toler);
end