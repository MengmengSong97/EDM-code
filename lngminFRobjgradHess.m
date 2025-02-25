function [gobj,gradgl,hessiangl] = lngminFRobjgradHess(H,L,Dbar,d,V)
%%% g(L)=f(VL), L: R(n-1)d , V=null(en): Rn(n-1)
%%% calculate the objective function value, gradient and Hessian matrix of g,
%%% regrading to L

%%%%%%%%%Intialize
n = length(Dbar);
l=vec(L);
en = ones(n,1);
tn = (n+1)*n/2;

%%% functions
Sy = @(X) (X+X')/2;      % Syms: Rnn --> Rnn
Se = @(v)(en*v'+v*en');
K = @(B)(Se(diag(B))-2*B);
Ks = @(S)(2*(diag(sum(S))-S));
M = @(P)(P*P'); % M(P) = PP'
F = @(P)(K(M(P)) - Dbar);


VL = V*L;
FVL = F(VL);
HFVL = H.*FVL;
KsHFVL = Ks(HFVL);
B = VL*VL';
D = Se(diag(B)) - 2*B;
gobj = norm(H.*(D - Dbar),'fro')^2/2;
%if nargout == 1    % never used?????
%gradgl = [];
%hessiangl = [];
%else

%%%%%%%%%%%%%Evaluation of the gradient of g(L) at l = vec(L)
%gradgl = vec( 2*V'*(Ks(H.*F(VL)))*VL ); % gradf: R(n-1)d --> R((n-1)*d)
temp = Ks(HFVL); % gradf: R(n-1)d --> R((n-1)*d)
temp = (temp+temp')/2;
gradgl = 2*vec( (V'*temp*V)*L ); % gradf: R(n-1)d --> R((n-1)*d)

%%%%%%%%%%%%Evaluation of the hessian of g(L), regarding to vec(L)%%%%%%%%%%%%%%
%%%% speed this up?  kron only once? external fn?????????
ked = (kron(eye(d),V))';
%JV = @(dL)(zsvec(K(Sy(VL*((V*dL)')))));   % alternate 4*JV'*JV == H1
%VCurvaturelH1 = @(dL)(4*ked*vec(((VL)'* ...
%	Sy(Ks(H.*K(Sy(VL*(V*dL)')))))'));
%VCurvaturelH2 = @(dL)(2*ked*vec(KsHFVL*V*dL));
%VCurvaturel = @(dL)( (4*ked*vec(((VL)'*Sy(Ks(H.*K(Sy(VL*(V*dL)')))))')) ...
%	          + ...
%	   (2*ked*vec(KsHFVL*V*dL))  );
% VCurvaturep: (Rnd,Rnd) --> Rnd
eynd=eye((n-1)*d);
hessiangl=zeros((n-1)*d);
%hessianglH1 = hessiangl;
%hessianglH2 = hessiangl;
%Jtilde = zeros(tn,(n-1)*d);
for i=1:(n-1)*d
    Ae=reshape(eynd(:,i),[n-1,d]);
    %    hessianglH1(:,i) = VCurvaturelH1(Ae);
    %    hessianglH2(:,i) = VCurvaturelH2(Ae);
    hessiangl(:,i) = (4*ked*vec(((VL)'*Sy(Ks(H.*K(Sy(VL*(V*Ae)')))))')) ...
        + ...
        (2*ked*vec(KsHFVL*V*Ae));   % skip for now?????
    %%following is less expensive for H1????
    %    Jtilde(:,i) = JV(Ae);   % alternate 4*Jtilde'*Jtilde === H1
end
%hessianglH1 = (hessianglH1+hessianglH1')/2;
%hessianglH2 = (hessianglH2+hessianglH2')/2;
%norm(hessianglH1+hessianglH2 - hessiangl,'fro')
%end  % nargout == 1
%hessiangl = 4*(Jtilde'*Jtilde) + hessianglH2;
hessiangl = (hessiangl+hessiangl')/2;
end
