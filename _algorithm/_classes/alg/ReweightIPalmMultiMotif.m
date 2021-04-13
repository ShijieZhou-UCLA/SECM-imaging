classdef ReweightIPalmMultiMotif < Solver
% REWEIGHTIPALM Solve multiple iterates of IPalm, each with reweighted lda.
% For iterates k = 2,3,... the lda is determined by 
%                   lda^(k) = 1./(X^(k-1)+eps)
% where X^(k-1) is the solution of sparse map in last iterate.
% 
% obj = REWEIGHTIPALM(problem) Construct a solver for input problem. The 
% problem should contains properties {funcf, funch, proxf, gradh}.
%
% REWEIGHTIPALM methods:
%   SET_NITER_IPALM - Set number of iteration in each run of IPalm.
%
% For other methods/fields, see object Solver.
% 
% See also IPALM CALIBLASSO SOLVER

properties(Hidden)
    f_X    % function handle; null input function return sparse vector.
    prb    % Object CalibLasso; construct REWEIGHTIPALM with the problem.
    ipalm  % Object Ipalm; algorithm in each run. 
    niter_ipalm  % scalar; niter in each run of ipalm.
end

methods
    function obj = ReweightIPalmMultiMotif(problem)
    % REWEIGHTIPALM Construct Reweighted IPalm solver
        obj = obj@Solver(problem.vars);
        obj.f_X = @() obj.vars{1}.images; %%%Shijie
        obj.prb = problem;
        
        obj.plot_on = false; 
        obj.set_display(1);
        obj.set_maxiter(6);
        obj.set_niter_ipalm(30); %%%Shijie
    end
    
    function set_niter_ipalm(obj,niter_ipalm) 
    % obj.SET_NITER_IPALM(niter_ipalm) Set niter in each run of ipalm.
        obj.niter_ipalm = niter_ipalm; 
    end
end

methods (Access = protected)
    function iter(obj)
    % obj.ITER() Run one instance of IPALM
        X = obj.f_X();
        lda_array = cell(size(X));   %%%
%        Xsum = sum(cat(3,X{:}),3);   %%%
%        if norm(Xsum) ~= 0
        if obj.iiter > 1   %%%Shijie 
            for i = 1:size(X,1)
%                lda_array{i} = 1e-3*obj.objval./(X{i}+eps);   %%% Henry
%                lda_array{i} = 1e-5*obj.funchval./(X{i}+eps);  %%% 1tri_2
%                lda_array{i} = 1.225e-6*obj.funchval./(X{i}+eps);  %%% 8tri_2
%                lda_array{i} = 5.5e-6*obj.funchval./(X{i}+eps);  %%% 4tri_2
%                lda_array{i} = 2e-6*obj.funchval./(X{i}+eps);  %%% 4tri_1
                lda_array{i} = 3e-6*obj.funchval./(X{i}+eps);  %%% 4tri_1
            end
            obj.prb.set_lda_array(lda_array);
            obj.prb.set_init(obj.vars);
        end
        
        obj.ipalm = IPalm(obj.prb);
        obj.ipalm.set_display(5);
        obj.ipalm.set_maxiter(obj.niter_ipalm);
        figure();
        obj.ipalm.solve();
        close;
        obj.vars = obj.ipalm.vars;
    end
    
    function set_objval(obj); obj.objval = obj.ipalm.objval; end
    
    function set_funchval(obj); obj.funchval = obj.ipalm.funchval; end   %%%Shijie

    function stopping = stop_criterias(obj); stopping = false; end
end

end

