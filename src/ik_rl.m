function qs = ik_rl(Ln,qn,xgn,lbound,ubound)

optFcn = @(q_l) finalJ_rl(q_l, Ln, xgn);

optimopts = optimoptions('fmincon', 'SpecifyObjectiveGradient', true, 'HessianFcn', @(q_l, lambda) finalH_rl(q_l,Ln,xgn),'MaxIterations',300);
%optimopts = optimoptions('fmincon', 'Display', 'iter-detailed', 'SpecifyObjectiveGradient', true, 'HessianFcn', @(q_l, lambda) finalH_rl(q_l,Ln,xgn),'MaxIterations',300);
% optimopts = optimoptions('fmincon', 'Display', 'iter-detailed', 'SpecifyObjectiveGradient', true);

qs = fmincon(optFcn, qn, [], [], [], [], lbound, ubound, [], optimopts);
end

function [J, dJ] = finalJ_rl(q_l, Ln, xgn)
    J = J_rl_computable(q_l, Ln, xgn);
    dJ = dJ_rl_dq_computable(q_l, Ln, xgn).';
end

function H_rl = finalH_rl(q_l,Ln,xgn)
    H_rl = ddJ_rl_ddq_computable(q_l,Ln,xgn);
end