function [Lc] = lngminTR(n,d,Lbar,Lhat,V,H,toler)
en = ones(n,1);
%ed = ones(d,1);
Se = @(v)(en*v'+v*en');
K = @(B)(Se(diag(B))-2*B);%Lindenstrauss operator
%Ks = @(S)(2*(diag(sum(S))-S));%adjoint of K
Pbar = V * Lbar;
Phat = V * Lhat;
Bbar = Pbar*Pbar';
Dbar = K(Bbar);
%%% minimizing 1/2 norm(F_V(L) - Dbar)_F^2, F_V(L) = K_V(L) - Dbar
Lc = V'*Phat;
[fc,gc,Hc] = lngminFRobjgradHess(H,Lc,Dbar,d,V);
delta = max(4,norm(gc)/norm(Hc,'fro'));
delta2 = delta^2;
[uH,dH] = eig(Hc);
ddH=diag(dH);
lammin = min(ddH);
%dHlam = diag(dH);
%if lammin < 0
%	S = uH(:,1:length(Hc)-1)*diag(dHlam(2:end))*uH(:,1:length(Hc)-1)
%%%%%%%%%%%%%%%%%%%%%%%%  construct a convex reformulation of TRS %%%%%%%%%%%%%%%%%%%%%%%%
lambar = max(0,-lammin);
if lambar > 0    % neg eig for Hc
    NH = null(Hc - lammin*eye(length(Hc)));
    hardind = norm(NH'*gc); % not zero ==> easy case/move to bdry if hard case
    lamhat = lambar + hardind/delta;
    dS = sqrt(ddH + lamhat*ones(d*(n-1),1));
else
    dS = sqrt(ddH);
end
dS(dS < 1e-8) = 0;
dSpinv = zeros(d*(n-1),1);
dSpinv(dS > 1e-8) = 1./dS(dS > 1e-8);
S = uH*diag(dS)*uH';
gbar = uH*diag(dSpinv)*uH'*gc;
gbarlsq = lsqminnorm(S,gc);
if norm(gbarlsq) < norm(gbar)
    gbar = gbarlsq;
end
%%%%%now solve TRS  for Lc from given starting point
% cvx_begin quiet
% variable dl(d*(n-1),1);
%
% dual variable lamtrs
% minimize(norm(S*dl+gbar));
% subject to
%           %lamtrs : dl'*dl <= delta2;
%           lamtrs : norm(dl) <= delta;
% cvx_end
%%%%%%%%%%%%solve TRS by the (accelearated) projected gradient method %%%%%%%%%%%%%%%%%%%%%%
dl_p = zeros(d*(n-1),1);
dl = zeros(d*(n-1),1);
L=norm(S'*S, 2);
gradl = S'*(S*dl+gbar);
ngradl = norm(gradl);
ngc=norm(gc);
ii=0;
%%%%for solve TRS approximately
while ngradl>min(1e-4,ngc^2)*ngc &&  (gradl'*dl+ norm(gradl)*norm(dl) > min(1e-4,ngc^2)*ngc || abs(norm(dl)-delta)>toler)
    ii=ii+1;
    temp = dl+(ii-1)/(ii+1)*(dl-dl_p);%FISTA accelerated algorithm
    temp = temp-1/L*S'*(S*temp+gbar);
    normtemp = norm(temp);
    dl_p=dl;
    if normtemp > delta
        dl = temp/normtemp*delta;
    else
        dl = temp;
    end
    gradl = S'*(S*dl+gbar);
    ngradl = norm(gradl);
    if ii>1/toler
        fprintf('TRS not solve');
        break;
    end
end
lamtrs = norm(gradl)/delta;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if lammin < 0 && hardind < norm(gc)*1e-9 && (delta - norm(dl) ) > 1e-9*delta
    v = uH(:,1);  % eigenvec for lammin
    if -v'*dl >=0
        dlstep = ...
            (-v'*dl)+sqrt((v'*dl)^2-(norm(dl)^2-delta2));
    else
        dlstep = ...
            (-v'*dl)-sqrt((v'*dl)^2-(norm(dl)^2-delta2));
    end
    dl = dl + dlstep*v;
    fprintf('\nSTEP to BDRY: lammin<0=%g; cvx dual=%g\n', ...
        lammin,lamtrs)
    fprintf('%-8s','iter#');
    fprintf('%-15s','rel.stopcrit','fc','||gc||')
    fprintf('%-15s','delta','rtrs','lamminHessobj','lamtrs');
    fprintf('\n')
    if abs(norm(dl)-delta)>((1e-12)*delta)
        fprintf('error in step to bdry\n')
        keyboard
    end
end
lc = Lc(:);
lp = lc + dl;


fprintf('\n\n== NEW lngm search (prob size n=%i, emb dim d=%i) ====== \n', n,d)
fprintf('%-8s','iter#');
fprintf('%-15s','fc','||gc||','delta')
fprintf('%-15s','rtrs','lamminHessobj','lamtrs');
fprintf('\n')


iter = 0;Lp = zeros(n-1,d);
while norm(gc) > toler && fc>1 && iter < 1000
    %%%% change to make option ... only evaluate fc here
    iter = iter + 1;
    Lp = reshape(lp,n-1,d);
    [fp, ~, ~] = lngminFRobjgradHess(H,Lp,Dbar,d,V);
    deltaf = (fc-fp);
    deltaq =  - ((dl'*Hc*dl/2)+gc'*dl);  % our quad model does not have fc
    rtrs = deltaf/deltaq;
    if rtrs < 1/4
        delta = norm(dl)/4;
        %%%   dl is always active as we move to bdry till end
    elseif rtrs > 3/4 && (norm(dl) > delta/2)
        delta = 2*delta;
    end
    if rtrs >= 0
        lc = lp;
    end

    %%%% solve new TRS
    delta2 = delta^2;
    Lc = reshape(lc,n-1,d);
    [fc,gc,Hc]=lngminFRobjgradHess(H,Lc,Dbar,d,V);
    [uH,dH] = eig(Hc);
    ddH = diag(dH);
    lammin = min(ddH);
    %lammax = max(ddH);
    %lammean = mean(abs(ddH));
    lambar = max(0,-lammin);
    %%%%%%%%%%%%%%%%%%%%%%%%  construct a convex reformulation of TRS %%%%%%%%%%%%%%%%%%%%%%%%
    if lambar > 0    % neg eig for Hlc
        NH = null(Hc - dH(1,1)*eye(length(Hc)));
        hardind = norm(NH'*gc); % not zero ==> easy case/move to bdry if hard case
        lamhat = lambar + hardind/delta;
        dS = sqrt(ddH + lamhat*ones(d*(n-1),1));
    else
        dS = sqrt(ddH);
    end
    dS(dS < 1e-8) = 0;
    dSpinv = zeros(d*(n-1),1);
    dSpinv(dS > 1e-8) = 1./dS(dS > 1e-8);
    S = uH*diag(dS)*uH';
    gbar = uH*diag(dSpinv)*uH'*gc;

    %%%%%now solve TRS
    % cvx_begin sdp quiet
    % cvx_precision best
    % variable dl(d*(n-1),1);
    % dual variable lamtrs
    % minimize(norm(S*dl+gbar));
    % subject to
    %     lamtrs : norm(dl) <= delta;
    % cvx_end
    %%%%%% solve TRS by the (accelearated) projected gradient method %%%%%%%%%%%%%%%%%%%%%%
    dl_p = zeros(d*(n-1),1);
    dl = zeros(d*(n-1),1);
    L=normest(S'*S, 1e-6);
    gradl = S'*(S*dl+gbar);
    ii=0;
    ngradl=norm(gradl);
    ngc=norm(gc);
    while ngradl>min(1e-4,ngc^(1/2))*ngc && (gradl'*dl + norm(gradl)*norm(dl)> min(1e-4,ngc^(1/2))*ngc || abs(norm(dl)-delta)>toler)
        ii=ii+1;
        temp = dl+(ii-1)/(ii+1)*(dl-dl_p);
        temp = temp-1/L*S'*(S*temp+gbar);
        normtemp = norm(temp);
        dl_p = dl;
        if normtemp > delta
            dl = temp/normtemp*delta;
        else
            dl = temp;
        end
        gradl = S'*(S*dl+gbar);
        ngradl=norm(gradl);
        if ii>1/toler
            fprintf('\n PG fails\n');
            break;
        end
    end
    lamtrs = ngradl/delta;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%% if hard case hold, move the global minimizer of TRS to the boundry
    if lammin < 0 && hardind < norm(gc)*1e-6 && (delta - norm(dl)) > 1e-6*delta
        v = uH(:,1);  % eigenvec for lammin
        if -v'*dl >=0
            dlstep = ...
                (-v'*dl)+sqrt((v'*dl)^2-(norm(dl)^2-delta2));
        else
            dlstep = ...
                (-v'*dl)-sqrt((v'*dl)^2-(norm(dl)^2-delta2));
        end
        dl = dl + dlstep*v;
        fprintf('\nSTEP to BDRY: lammin<0=%g; cvx dual=%g\n', ...
            lammin,lamtrs)
        fprintf('%-8s','iter#');
        fprintf('%-15s','fc','||gc||')
        fprintf('%-15s','delta','rtrs','lamminHess','lamtrs');
        fprintf('\n')
        if abs(norm(dl)-delta) > (1e-12*delta)
            fprintf('error in step to bdry\n')
            keyboard
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 	%%% for a newton free step??????? temporary!!!!!
    %               %above done %  [fc,gc,Hc]=lngminFRobjgradHess(H,lc,Dbar,d,V);
    %  	      if norm(dl) < delta/2 || stopcrit < 1e2*toler
    % 		      %check for best search direction??????
    %  		      %fprintf('newton step???\n')
    %  	              dltemp = -Hc\gc;   %   Newton step??
    % 		      [R,flagc] = chol(Hc);
    %  		      if flagc == 0
    %  		      	      dltemp2 = -R\((R')\gc);
    % 		      else
    % 			      if min(eig(Hc)) < -n*eps(norm(Hc))
    %                       fprintf('roundoff error in TRS?\n')
    %  			              dltemp2 = -lsqminnorm(Hc,gc);
    %  			      end
    %  		      end
    %  		      [vals,indmin] = ...
    %  			    min([norm(Hc*dltemp+gc) norm(Hc*dltemp2+gc)  ...
    %  			              norm(Hc*dl+gc)]);
    %  			    %([norm(Hc*dltemp+gc) norm(Hc*dltemp2+gc)  ...
    % 			    %  norm(Hc*dl+gc) norm(dltemp-dl)]);
    %  		      if indmin == 1
    %  				      dl = dltemp;
    % 				      ndltemp = ndltemp + 1;
    %                       %fprintf('dl = dltemp\n');
    %  		      elseif indmin == 2
    %  				      dl = dltemp2;
    %  				      ndltemp2 = ndltemp2 + 1;
    %  				      %fprintf('dl = dltemp2 chol\n');
    %  		      end
    %           end
    lp = lc + dl;   % new lp

    fprintf('%-8i',iter);
    fprintf('%-15.2e',fc,ngc,delta,rtrs,lammin,lamtrs);
    fprintf('\n');
    if 0<abs(deltaf)/fp<1e-12
        fprintf('Rounding errors result in inaccurate calculations in deltaf');
        break;
    end
    %keyboard
end   % while loop
fprintf('\nAFTER While loop\n')
fprintf('%-8s','iter#');
fprintf('%-15s','fc','||gc||','delta')
fprintf('%-15s','rtrs','lamminHess','lamtrs');
fprintf('\n')
%%%%%%%%% check if a lngm is found
if norm(gc) < toler && fc>1
    fprintf('\n A lngm is found.\n');
    %nlngm = nlngm+1;
    %Bbar = Pbar*Pbar';
    %Dbar = K(Bbar);
    %Blngm = V*Lc*Lc'*V';
    %Dlngm = K(Blngm);
    if d==1
        figure(1)
        plot(Pbar);
        hold on;
        plot(V*Lc);
        hold off;
    end
    if d==2
        Plngm=V*Lc;
        figure(1)
        plot(Pbar(:,1),Pbar(:,2));
        hold on;
        plot(Plngm(:,1),Plngm(:,2));
        hold off;
    end
    keyboard;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if norm(gc)/(1+abs(fc)) < 1e2*toler && lammin > -1e2*toler && ...
% 	abs(fc) > 1e8*toler
% 	Lc = reshape(lc,n-1,d);
% 	Gc = V*Lc*Lc'*V'; Gc = (Gc+Gc')/2;
% 	Dc = K(Gc);
% 	fprintf('possible COUNTER EXAMPLE in lc,Lc \n')
% 	fprintf('normgc %g, relnorm gc %g; lammin %g rel lammin/lammean %g\n', ...
%       norm(gc),norm(gc)/(1+abs(fc)),lammin,lammin/(abs(lammean)+1))
%         ncntexs = 0;
%         if d == 1 && verboseplot
%           	   plot((V*lp)',zeros(1,n),'.')
% 	        title(['d=1; ',num2str(n),'  points on the line'])
%         end
% end

% iters(ii) = iter;
% %% for alternate Hessian calcs can use the following
% %[fc,gc,Hc,Hc1,Hc2,Jtilde] = lngminFRobjgradHess(H,lc,Dbar,d,V);
% Lc=reshape(lc,n-1,d);
% [fc,gc,Hc] = lngminFRobjgradHess(H,Lc,Dbar,d,V);
% [UHc,DHc] = eig(Hc);
% %%%%%%%%%%%%%test random orthog qrand matrix
% ldesc = UHc(:,1);   % direction of descent?? from eigenvector
% [qrand,~] = qr(randn(d));
% Lcq = Lc*qrand;
% lcq = Lcq(:);
% [fcq,gcq,Hcq] = lngminFRobjgradHess(H,Lc*qrand,Dbar,d,V);
% [UHcq,DHcq] = eig(Hcq);
% ts = linspace(-.01,.01,20);
% fcqs = zeros(length(ts),1);
% for ii = 1:length(ts)
% 	%lcts = lc+ts(ii)*ldesc;
% 	lcts = lc+ts(ii)*(lc-lcq);
%     Lcts=reshape(lcts,n-1,d);
% 	[fcqs(ii),gcq,Hcq] = lngminFRobjgradHess(H,Lcts,Dbar,d,V);
% end
% plot(fcqs)
% Pc = V*Lc;
% dP = V*reshape(ldesc,n-1,d);  % using eigenvector for 0 eigenvalue
% thirdderiv = trace(Ks(K((Pc*dP'+dP*Pc')))*(dP*dP'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fprintf('min eigs Hc1; Hc2: %g; %g\n',min(eig(Hc1)),min(eig(Hc2)))
fprintf('\n\n== lngm search (prob size n=%i, emb dim d=%i) ====== \n', n,d);
end
