function qs = ik_ub(Ln,qn,xgn,lbound,ubound)

optFcn = @(q) finalJ_ub(q, Ln, xgn);

%optimopts = optimoptions('fmincon', 'Display', 'iter-detailed', 'SpecifyObjectiveGradient', true, 'HessianFcn', @(q, lambda) finalH_ub(qb,q,Ln,xgn));
optimopts = optimoptions('fmincon', 'Display', 'iter-detailed', 'SpecifyObjectiveGradient', true);

qs = fmincon(optFcn, qn, [], [], [], [], lbound, ubound, [], optimopts);
end

function [J, dJ] = finalJ_ub(q, Ln, xgn)
    J = J_ub_computable(q, Ln, xgn);
    dJ = dJ_ub_dq_computable(q, Ln, xgn).';
end

% function H_ub = finalH_ub(qb,q,Ln,xgn)
%     H_ub = ddJ_ub_ddq_computable(qb,q,Ln,xgn);
% end