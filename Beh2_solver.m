%pkg load io
close all
% ============================
% PROBLEM SETUP
f = @(x) Beh2(x);
df = @(x) dBeh2(x);
gsg_tol = 0.9; solver = 'LM';
iters = 0; maxiters = 100;
tol = 1e-8; err = tol + 1;
x = [0.008;2.0];
% ============================
% SOLVE 
plot_options = {'linewidth',1.5}; 
figure(1), hold on
pnm(f,df,x,'solver',solver,'beh1','plot',plot_options);
pnm(f,df,x,'solver',solver,'beh1','anderson',1,'plot',plot_options);
pnm(f,df,x,'solver',solver,'beh1','anderson',5,'plot',plot_options);

pnm(f,df,x,'solver',solver,'beh1','anderson',1,'safeguard',gsg_tol,'plot',[plot_options,'-.']);

legend(solver,['AA',solver,'(1)'],['AA',solver,'(5)'],['\gammaAA',solver,'(1,0.9)']);
title(['Solve for Beh2 with ', solver]);
xlabel('iterations');
ylabel("||f'(x)^Tf(x)||");
set(gca,'fontsize',18);
