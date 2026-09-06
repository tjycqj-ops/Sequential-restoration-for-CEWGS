function [H_lower,H_upper]=pump_head_projection_bounds(alpha,beta,q,omega_lower,omega_upper)
%PUMP_HEAD_PROJECTION_BOUNDS Project H=alpha*q*omega+beta*omega^2.
% Valid when dH/domega=alpha*q+2*beta*omega is nonnegative over the
% supplied operating interval, which holds for the current pump data.

if isnumeric(q)
    assert(all(alpha.*q+2*beta.*omega_lower>=-1e-12,'all'), ...
        'Pump head is not monotone in speed over the requested interval.');
end
H_lower=alpha.*q.*omega_lower+beta.*omega_lower.^2;
H_upper=alpha.*q.*omega_upper+beta.*omega_upper.^2;
end
