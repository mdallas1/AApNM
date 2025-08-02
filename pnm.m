function pnm(f,df,x0,varargin)
	%{
	-- pnm(f,df,x) solves f(x) = 0 
	using Anderson accelerated 
	perturbed Newton methods (pNM). 
	
	-- User may set various parameters 
		 described below. 
	 
	-- f and it's Jacobian df are 
		 function handles.
	-- x is the initial iterate.

	-- pnm(...,'solver',solver) 
		 solves f(x) = 0 using the 
		 user specified solver. The 
		 default solver is the 
		 Levenberg-Marquardt algorithm. 

	-- pnm(...,'anderson',m) applies 
		 Anderson acceleration to solver 
		 with algorithmic depth m.
	
	-- pnm(...,'safeguard',r) implements 
		 adaptive gamma safeguarding with 
		 tolerance r.

	-- pnm(...,'tol',tol) sets tolerance. 
		 Default tol = 1e-8. 

	-- pnm(...,'maxiters',maxiters) sets 
		 maximum iterations. Default is 100.

	-- pnm(...,'save','filename',<filename>) saves data to  
		 spreadsheet named <filename>. If no name is given, then 
		 a default one is used. 

	-- pnm(...,'plot',plot_options) plots 
		 residual history, plot_options is 
		 a cell array. 
	%}

	assert(isa(f,'function_handle'));
	assert(isa(df,'function_handle'));
	assert(isa(x0,'double'));

	[x_m,num_solves] = size(x0); 
	avg_iters = 0; avg_res = 0;

	p = inputParser(); 
	p.FunctionName = 'pnm' ; 

	default_tol = 1e-8;
	p.addParameter('tol',default_tol,@isdouble);

	default_maxiters = 100;
	p.addParameter('maxiters',default_maxiters,@isdouble);

	vld_solver = @(s) any(strcmp(s,{'newton','newt','lm','LM'}));
	default_solver = 'lm';
	p.addParameter('solver',default_solver,vld_solver);

	default_depth= 0;
	p.addParameter('anderson',default_depth,@isfloat);

	default_gsg_tol = 0; 
	p.addParameter('safeguard',default_gsg_tol,@isfloat);

	default_plot = {}; 
	p.addParameter('plot',default_plot,@(a) isa(a,'cell'));

	p.addSwitch('save');
	
	default_filename = ''; 
	p.addParameter('filename',default_filename,@ischar)

	p.addSwitch('inexact');

	p.parse(varargin{:}); 
	solver = p.Results.solver; m = p.Results.anderson; gsg_tol = p.Results.safeguard;
	tol = p.Results.tol; maxiters = p.Results.maxiters; 
	plot_options = p.Results.plot; 

	if any(strcmp('anderson',varargin));
		aa = 1;
	else 
		aa = 0;
	end 

	if any(strcmp('safeguard',varargin));
		gsg = 1;
	else 
		gsg = 0;
	end

	if any(strcmp('plot',varargin))
		plt = 1; 
	else 
		plt = 0; 
	end

	% - CREATE FILE TO SAVE DATA
	if p.Results.save
			if isempty(p.Results.filename);
				dat = datestr(now,31); 
				func_name = func2str(f);
				arr = {"solver",solver,'anderson',m,'safeguard',gsg_tol,'tol',tol,'maxiters',maxiters};
			if p.Results.inexact 
				arr{2} = ["inexact-",solver];
			end
			if aa
				filename = [func_name(6:end),'_',arr{2},'-AA','(',num2str(m),')','_',dat,'.csv'];
				if gsg 
					filename = [func_name(6:end),'_',arr{2},'-AA','(',num2str(m),')-',num2str(gsg_tol),'_',dat,'.csv'];
				end
			else 
				filename = [func_name(6:end),'_',arr{2},'_',dat,'.csv'];
			end
			cell2csv(filename,arr);
			end
		end

	% - SOLVE
	for k = 1:num_solves
		x = x0(:,k); 
		f = @(x) f(x); df = @(x) df(x);
		fx = f(x); dfx = df(x); 
		iters = 0; I = eye(length(x));
		res = norm(fx); 
		% ONLY FOR BEH PROBLEMS	
		%res = norm(dfx' * fx); disp("USING BEH RESIDUAL")
		
		damping = 1; 
		res_arr = []; x_arr = []; d_arr = []; xaa_arr = [];
		if res < tol 
			disp("Initial iterate satisfies tolerance.");
		end
		while (iters < maxiters && res > tol)
			res_arr = [res res_arr];
			if iters < 1
				x_arr = [x x_arr]; 
				if any(strcmp(solver,{'lm','LM'}))
					mu = 0.5*1e-8*norm(fx)^2; % LM parameter from KaYaFu03
					%mu = res; 
					if p.Results.inexact
						% inital forcing 0.5 follows Eisentat & Walker 
						forcing = 0.5;
						[d,~] = gmres((dfx'*dfx + mu * I),-dfx'*fx,20,min(forcing*res,forcing));
					else 
						d = - (dfx'*dfx + mu * I)\(dfx'*fx);
					end
				elseif any(strcmp(solver,{'newt','newton'}))
					if p.Results.inexact
						forcing = 0.5;
						[d,~] = gmres(dfx,-fx,20,min(forcing*res,forcing));
					else 
						d = -dfx\fx;
					end
				end
				x_arr = [x x_arr];
				d_arr = [d d_arr];
				x = x+d;
			else 
				if any(strcmp(solver,{'lm','LM'}))
					mu = min(mu,nfx^2);						
					%mu = res; 
					dfx_dot = dfx'*dfx;
					if p.Results.inexact
						% EISENSTAT & WALKER CHOICE 2
						forcing = min(0.9,0.5*(res_arr(:,1)/res_arr(:,2))^(1.5));
						[d,~] = gmres(dfx_dot + mu * I,-dfx'*fx,20,min(forcing*res_arr(:,1),forcing));
					else 
						d = - (dfx_dot + mu * I)\(dfx'*fx);
					end
				elseif any(strcmp(solver,{'newt','newton'}))
					if p.Results.inexact
						% EISENSTAT & WALKER CHOICE 2
						forcing = min(0.9,0.5*(res_arr(:,1)/res_arr(:,2))^(1.5));
						[d,~] = gmres(dfx,-fx,20,min(forcing*res_arr(:,1),forcing));
					else 
						d = -dfx\fx;
					end
				end
	
				% - ANDERSON	
				if aa 
					if gsg && res < 1e-1
						m = 1; 
					end
					mm = min(m,iters);
					x_arr = [x x_arr(:,1:mm)];
					d_arr = [d d_arr(:,1:mm)];
					%F = d_arr(:,1) - d_arr(:,2:end);
					%alpha = F \ d_arr(:,1);
					F = d_arr(:,1:end-1) - d_arr(:,2:end);
					gamma = F \ d_arr(:,1);
				end
	
				if gsg && (m < 2)
					[gamma,lambda,~] = gamma_safeguard(d_arr(:,1),d_arr(:,2),gamma,gsg_tol,1,eye(length(x)));
				end
	
				if aa 
						%x = x_arr(:,1) + damping*d_arr(:,1) - ( (x_arr(:,1) - x_arr(:,2:end)) + damping*F)*alpha;
						x = x_arr(:,1) + damping*d_arr(:,1) - ( (x_arr(:,1:end-1) - x_arr(:,2:end)) + damping*F)*gamma;	
				else 
					x = x+d;
				end
			end
			iters = iters+1;
			fx = f(x);
			nfx = norm(fx); 
			dfx = df(x);  
			res = nfx; 
			% ONLY FOR BEH PROBLEMS
			%res = norm(dfx' * fx); disp("USING BEH RESIDUAL"); 
			fprintf('%g,%g\n',iters,res)
		end
		res_arr = [res res_arr];
		x
		if p.Results.save 
			fid = fopen(filename,'a');
			dlmwrite(filename,[0:iters],'-append');
			dlmwrite(filename,flip(res_arr),'-append');
			fclose(fid);
		end
		avg_iters += iters;
		avg_res += res; 
	end
		avg_iters = avg_iters / num_solves ;
		avg_res = avg_res / num_solves;
		fid = fopen(filename,'a');
		dlmwrite(filename,[avg_iters avg_res],'-append');
		fclose(fid);

	if plt
		figure(2), hold on
		semilogy([0:iters],flip(res_arr),plot_options{:});
	end

end
			
