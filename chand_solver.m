%pkg load io
close all
% ============================
% PROBLEM SETUP
omega_bar = 1;
f = @(x) chand_fun(x,omega_bar,0);
df = @(x) chand_fun(x,omega_bar,1);
gsg_tol = 0.9; solver = 'newt';
iters = 0; maxiters = 100;
tol = 1e-8; err = tol + 1;
x = ones(1e3,1);
% ============================
plot_options = {'linewidth',1.5}; 
figure(1), hold on
% SOLVE 
pnm(f,df,x,'solver',solver,'plot',plot_options);
pnm(f,df,x,'solver',solver,'anderson',1,'plot',plot_options);
pnm(f,df,x,'solver',solver,'anderson',1,'safeguard',gsg_tol,'plot',[plot_options, '-.']);
legend(solver,['AA',solver,'(1)'],['\gammaAA',solver,'(1,0.9)']);
title(['Solve for CH Equation with ', solver]);
xlabel('iterations');
ylabel("||f(x)||");
set(gca,'fontsize',18);
