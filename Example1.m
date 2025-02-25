%%% Data and Calculations in Example 1
%%% Error analysis are done
clc;clear;
n=50;d=1;
Pbar=[-0.60599420419914
   1.25253980362745
  -0.84238180586740
  -0.53579611085847
   0.19708931277192
  -1.79984782222111
  -0.41527865379273
   0.21915952071487
  -0.44909567056284
   0.31885876974828
   2.50030338660590
   1.49951591764291
   0.12255038262099
   0.16908428311885
  -1.11206775325418
   0.07998074823284
   0.03054981643026
  -0.76252052456168
  -0.65660027137029
  -0.82015747596304
   0.72087991038672
  -0.71053497218089
  -0.08107320462582
  -0.34225738304152
  -1.56521625701592
  -0.12539009887519
   0.41245452266308
  -0.70359762327541
   0.56572499282342
   0.36281513185751
   0.10360071656594
   0.31861446173879
   1.61101822550653
  -0.43395723595537
   1.05675299335403
  -0.14531333793808
   0.30177497095087
  -0.09375243537478
   0.60941253173112
   0.81204014650370
  -1.83584573330158
  -1.14955942388795
   0.48587665889069
   1.14619235440365
   0.00742509462435
   0.30446085136041
  -0.82672669605142
   0.43384263788938
  -0.05997259898909
   0.43041915039944];
Pc=-[0.43623325112033
  -1.22762020433635
   0.70010136135430
   0.35991415873168
  -0.28228560377420
   1.73081724658463
   0.23439486877264
  -0.29950524425661
   0.26875960498095
  -0.37830516857461
   1.99137198339376
  -1.47176630576290
  -0.22447090537459
  -0.26051580040020
   0.99896634828819
  -0.19153324524272
  -0.15316596774634
   0.61056805982301
   0.49218270988943
   0.67520088533096
  -0.72157601943584
   0.55233232235185
  -0.06496325807452
   0.16281245647958
   1.48497061031969
  -0.02888772296088
  -0.45429992108691
   0.54457520973120
  -0.58387966790691
  -0.41371997944492
  -0.20981187084355
  -0.37810961138561
  -1.58278065634917
   0.25328548668347
  -1.03696183043215
  -0.01239928459453
  -0.36466259224993
  -0.05471978005445
  -0.62201496705578
  -0.80516992708878
   1.76829199679283
   1.03991090257334
  -0.51554999274410
  -1.12363323905098
  -0.13511728276064
  -0.36680318321451
   0.68256330402195
  -0.47198785003490
  -0.08188748229023
  -0.46914820269597];
if abs(sum(Pbar))<1e-14 && abs(sum(Pc))<1e-14
    fprintf('\nThere are a centralized Pbar and a centralized Pc\n');
end
dis_diff=zeros(50,50);
for i=1:50
    for j=1:50
        dis_diff(i,j)=(Pc(i)-Pc(j))^2-(Pbar(i)-Pbar(j))^2;
    end
end
%Since Pc(i)-Pc(j) and Pbar(i)-Pbar(j) are in (-6,6) and are accurate,
%(Pc(i)-Pc(j))^2, (Pbar(i)-Pbar(j))^2 are accurate 
%for its shown digits, i.e., every one has at most 1e-14 error. 
%Thus, dis_diff(i,j) has at most 2e-14 error.
%We note that all dis_diff(i,j) are in (-20, 20).
%Thus, dis_diff(i,j)^2 has at most 2*20*2e-14+2e-14*2e-14<8e-13 error.
%Therefore, summation of dis_diff(i,j)^2 has at most 50*49*8e-13<1.96e-9 error.
f_Pc_lb=dis_diff(:)'*dis_diff(:)/2-1.96e-9;%lower bound of f_Pc
fprintf('\nA lower bound of f_L(Lc) is %g\n',f_Pc_lb);%It is known that f_L(Lc)=f(V*Lc)=f(Pc)
grad_Pc=zeros(50,1);
for i=1:50
for j=1:50
    grad_Pc(i)=grad_Pc(i)+4*dis_diff(i,j)*(Pc(i)-Pc(j));
end
end
%Since dis_diff(i,j) has at most 2e-14 error, and Pc(i,:)-Pc(j,:) are accurate in (-6, 6),
%4*dis_diff(i,j)*(Pc(i)-Pc(j)) has at most 4*6*2e-14<4.8e-13 error.
%Note that every absolute value of the element in
%4*dis_diff(i,j)*(Pc(i,:)-Pc(j,:)) is smaller than 172.
%The round-off error in 4*dis_diff(i,j)*(Pc(i)-Pc(j)) is no more than 1e-12.
%Thus, 4*dis_diff(i,j)*(Pc(i,:)-Pc(j,:)) has no more than 1e-12 error.
%Therefore, grad_Pc(i) has no more than 50*1e-12=5e-11 error.
%Since max(max(abs(grad_Pc))) is presented as 8.923981265063219e-08,
%norm(grad_Pc)^2 has no more than 50*(2*9e-8*5e-11+5e-11*5e-11)<4.6e-16 error.
norm_grads_ub=grad_Pc'*grad_Pc+4.6e-16;%an upper bound of \|norm(grad_Pc)\|^2
%norm(grad(f_L(Lc)))\le norm(grad(f(V*Lc)))=norm(grad(f(Pc)))
if 1e-7*1e-7>norm_grads_ub
fprintf('\nAn upper bound of norm(grad(f_L(Lc))) is %g\n',1e-7);
end
norm_relgrads_ub=norm_grads_ub/(1+f_Pc_lb);%an upper bound of \|norm(grad_Pc)\|^2
fprintf('\nAn upper bound of norm(grad(f_L(Lc))/(1+f_Pc_lb)) is %g\n',norm_relgrads_ub);
Hessian_Pc=zeros(50,50);
for i=1:50
for j=1:50
    if i~=j
        Hessian_Pc(i,j)=-4*dis_diff(i,j)-8*(Pc(i,1)-Pc(j,1))^2;
        Hessian_Pc(i,i)=Hessian_Pc(i,i)+4*dis_diff(i,j)+8*(Pc(i,1)-Pc(j,1))^2;
    end
end
end
%Since dis_diff(i,j) has at most 2e-14 error,
%all dis_diff(i,j) are in (-20, 20), and Pc(i)-Pc(j) are accurate in (-6, 6),
%every element of Hessian_Pc has at most 4*2e-14+8*1e-14<1.6e-13 error.
%Since max(max(abs(Hessian_Pc))) gives 1.420256399165363e+03,
%considering the round-off error of storing Hessian_Pc, Hessian_Pc has at most 1e-12 error.
A= ([ones(1,n-1);-eye(n-1)]); 
[V] = GS(A);
%Here V is not an accurate basis of null(en) as V(1,1)=1/sqrt(2) is an irrational number.
%However, every element of V has error no more than 1e-15.
%This can be verified by the following.
red_flag=0;
for i=1:49
x=V(i+1,i)-1e-15;
if (x/i)^2*i+x^2-1<0
    red_flag=1;
end
x=V(i+1,i)+1e-15;
if (x/i)^2*i+x^2-1>0
    red_flag=1;
end
end
if red_flag==0
    fprintf('\nEvery element of V has error no more than %g\n',1e-15);
end
Hessian_Lc=V'*Hessian_Pc*V;
eigHl=eig(Hessian_Lc);
%We first analyse the error in v'*Hessian_Pc*v, where v denotes any column of V.
%Since every element of Hessian_Pc has at most 1e-12 error, every element is no more than 1.5e+3, 
%and every element of v has at most 1e-15 error, Hessian_Pc(i,j)*v(j) has at most 
%1.5e+3*1e-15+1*1e-12+1e-15*1e-12<2.5e-12 error. Thus, every element of
%Hessian_Pc*v has at most 50*2.5e-12=5e-10 error.
%Since max(sum(abs(Hessian_Pc)))<3.6e+3, v'*Hessian_Pc*v has at most 50*(3.6e+3*1e-15+1*5e-10+1e-15*5e-10)w<2.6e-8 error. 
%Therefore,every element of Hessian_Lc has at most 2.6e-8 error.
%It is fair to think the numerical number min(eig(Hessian_Lc)) has no more than 0.5, 
%then, for all x\in\R^49 with the unit norm, we have x'*Hessian_Lc*x>-49*49*2.6e-8+min(eigHl)-0.5>210.
lambdamin_lb=-49*49*2.6e-8+min(eigHl)-0.5;
lambdamin_ub=49*49*2.6e-8+min(eigHl)+0.5;
fprintf('\nA lower bound of lambdamin(Hessian_Lc) is %g\n',lambdamin_lb);
fprintf('\nAn upper bound of lambdamin(Hessian_Lc) is %g\n',lambdamin_ub);

%Next we calculate gamma: the Lipschitz constant of Hessian
G=0;
for i=1:50
    for j=1:50
        G=G+abs(Pc(i)-Pc(j));
%Pc(i)-Pc(j) are accurate. 
    end
end
%Since G<2.2e+3, all (Pc(i)-Pc(j)) are accurate, G has round-off error at most 1e-11.
G_ub=G+1e-11;
fprintf('\nAn upper bound of sum_{i,j}norm(P(i)-P(j)) is %g\n',G_ub);
r=1e-3;
gamma=24*sqrt(2)*(G_ub+2*n*sqrt(n)*r);
%2*n*sqrt(n) has an error at most 1e-12, G_ub+2*n*sqrt(n) has an error at most 1e-11.
%gamma has an error at most 24*1e-14*2835+34*1e-11+1e-14*1e-11<1.1e-9
fprintf('\nAn upper bound of the Lipschitz of Hessian is %g\n',gamma+1.1e-9);



figure; 
scatter([1:10 12:n],[Pbar(1:10);Pbar(12:n)],'k','filled'); 
hold on; 
scatter(11,Pbar(11),'r','filled'); 
hold on; 
grid on; 
set(gca,'xtick',0:1:50);
hold on; 
scatter([1:10 12:n],[Pc(1:10);Pc(12:n)],'b'); 
xlim([0 50]);
ylim([-3 3]); 
hold on; 
scatter(11,Pc(11),'m'); 
hold on; 
legend('$\bar{p}_i, i\neq i_0$','$\bar{p}_{i_0}$','$\tilde{p}_i, i\neq i_0$','$\tilde{p}_{i_0}$','Interpreter','latex','FontSize',16);
xlabel('$i$','Interpreter','latex','FontSize',20);
ylabel('$p_i$','Interpreter','latex','FontSize',20);
hold off;