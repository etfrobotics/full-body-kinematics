function qs = ik_ll(Ln,qn,xgn,lbound,ubound)

optFcn = @(q_l) finalJ_ll(q_l, Ln, xgn);

optimopts = optimoptions('fmincon', 'Display', 'iter-detailed', 'SpecifyObjectiveGradient', true, 'HessianFcn', @(q_l, lambda) finalH_ll(q_l,Ln,xgn));
% optimopts = optimoptions('fmincon', 'Display', 'iter-detailed', 'SpecifyObjectiveGradient', true);

qs = fmincon(optFcn, qn, [], [], [], [], lbound, ubound, [], optimopts);
end

function [J, dJ] = finalJ_ll(q_l, Ln, xgn)
    J = J_ll_computable(q_l, Ln, xgn);
    dJ = dJ_ll_dq_computable(q_l, Ln, xgn).';
end

function H_ll = finalH_ll(q_l,Ln,xgn)
    H_ll = ddJ_ll_ddq_computable(q_l,Ln,xgn);
end