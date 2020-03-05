load System_Data_N.mat
%Satya 16D170026 Exam 24

global U_k

SetGraphics;

N_samples = 500;

n_st = length(Xs);
n_ip = length(Us);
n_op = length(Ys);
%initialising
xk = zeros(n_st, N_samples);
yk = zeros(n_op, N_samples);
uk = zeros(n_ip, N_samples);

Xk = zeros(n_st, N_samples);
Yk = zeros(n_op, N_samples);
Uk = zeros(n_ip, N_samples);
%perturbation
ip1 = 0.015*idinput(N_samples, 'rbs', [0 0.025]);
ip2 = 2*idinput(N_samples, 'rbs', [0 0.025]);

uk = [ip1 ip2]';
%noise
vk = mvnrnd (zeros(2, 1) , R_mat , N_samples);
vk = vk';
wk = mvnrnd(zeros(2 , 1), Q_mat , N_samples);
wk = wk';

%initial conditions
Xk(:,1) = Xs;
Yk(:,1) = C_mat*Xk(:,1) + vk(:,1);

xk(:,1) = Xk(:,1) - Xs;
yk(:,1) = Yk(:,1) - Ys;
uk(:,1) = Uk(:,1) - Us;
%time axis
kt = zeros(1,N_samples);
%for loop
for k = 1:1:N_samples-1
    kt(k) = k*samp_T;
    
    Uk(:,k) = Us + uk(:,k);
    
    U_k = Uk(:,k) + wk(:,k);
    
    [T, Xt] = ode45('dynamics', [0 samp_T], Xk(:,k));
    Xk(:,k+1) = Xt(end,:);
    
    Yk(:,k+1) = C_mat*Xk(:,k+1) + vk(:,k+1);
    yk(:,k+1) = Yk(:,k+1) - Ys;
      
    
end

k = N_samples;
kt(k) = k*samp_T;
Uk(:,k) = Us + uk(:,k);

%system identification

dat1 = iddata(Yk', Uk', samp_T);
model = pem(dat1);

A_mat = model.a;
B_mat = model.b;
D_mat = model.d;

cmod_lin = ss(A_mat, B_mat, C_mat, D_mat);
dmod = c2d(cmod, samp_T);


phy_est = dmod.a; % linearized phy
gama_est = dmod.b; % linearized gama

sys1 = ss(phy, gama, C_mat, D_mat, samp_T); % original system

ltiview({'step'},sys1);hold on; %step response of original system
ltiview({'step'},model); % step response of linearized system


