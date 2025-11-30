pkg load io
% ============================
% PROBLEM SETUP
f = @(x) Beh3(x);
df = @(x) dBeh3(x);
gsg_tol = 0.9; solver = 'LM';
iters = 0; maxiters = 100;
tol = 1e-8; err = tol + 1;
x = [pi;0.001];
% ============================
plot_options = {'linewidth',1.5}; 
figure(1), hold on
pnm(f,df,x,'solver',solver,'beh3','plot',plot_options);
pnm(f,df,x,'solver',solver,'beh3','anderson',1,'plot',plot_options);
pnm(f,df,x,'solver',solver,'beh3','anderson',1,'safeguard',gsg_tol,'plot',[plot_options, '-.']);
legend(solver,['AA',solver],['\gammaAA',solver]);
title(['Solve for Beh3 with ', solver]);
xlabel('iterations');
ylabel("||f'(x)^Tf(x)||");
set(gca,'fontsize',18);
