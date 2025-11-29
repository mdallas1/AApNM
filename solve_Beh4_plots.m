%pkg load io
close all
% ============================
% PROBLEM SETUP
f = @(x) Beh4(x);
df = @(x) dBeh4(x);
gsg_tol = 0.9; solver = 'LM';
iters = 0; maxiters = 100;
tol = 1e-8; err = tol + 1;
x = [0.01;0];
% ============================
plot_options = {'linewidth',1.5}; 
figure(1), hold on
pnm(f,df,x,'solver',solver,'beh4','plot',plot_options);
pnm(f,df,x,'solver',solver,'beh4','anderson',1,'plot',plot_options);
pnm(f,df,x,'solver',solver,'beh4','anderson',1,'safeguard',gsg_tol,'plot',[plot_options,'-.']);
pnm(f,df,x,'solver',solver,'beh4','anderson',5,'plot',plot_options);
legend(solver,['AA',solver],['\gammaAA',solver,'(1,0.9)'],...
	['AA',solver,'(5)']);
title(['Solve for Beh4 with ', solver, ' \mu = 5.0']);
xlabel('iterations');
ylabel("||f'(x)^Tf(x)||");
set(gca,'fontsize',18);

figure(2), hold on 
pnm(f,df,x,'solver',solver,'beh4_res','plot',plot_options);
pnm(f,df,x,'solver',solver,'beh4_res','anderson',2,'plot',plot_options);
pnm(f,df,x,'solver',solver,'beh4_res','anderson',3,'plot',plot_options);
pnm(f,df,x,'solver',solver,'beh4_res','anderson',4,'plot',plot_options);
pnm(f,df,x,'solver',solver,'beh4_res','anderson',5,'plot',plot_options);
legend(solver,['AA',solver,'(2)'],['AA',solver,'(3)'],...
	['AA',solver,'(4)'],['AA',solver,'(5)']);
title(['Solve for Beh4 with ', solver," and \mu = |f'(x)^Tf(x)|"]);
xlabel('iterations');
ylabel("||f'(x)^Tf(x)||");
set(gca,'fontsize',18);

