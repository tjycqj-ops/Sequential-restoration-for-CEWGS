function [ge_residual,le_residual]=pump_power_sca_dc_residuals( ...
    g,eff,P,H,q,H_ref,q_ref,q_scale,H_scale)
%PUMP_POWER_SCA_DC_RESIDUALS Two-sided convex restrictions of P=(g/eta)Hq.
% ge_residual<=0 is the SCA restriction for P>=physical power.
% le_residual<=0 is the SCA restriction for P<=physical power.
% Flow and head are normalized before the difference-of-convex split.

x=q./q_scale;
y=H./H_scale;
x_ref=q_ref./q_scale;
y_ref=H_ref./H_scale;
power_scale=g.*(1./eff).*q_scale.*H_scale;
z=P./power_scale;

minus_ref=x_ref-y_ref;
plus_ref=x_ref+y_ref;
minus_square_linear=0.25*minus_ref.^2+0.5*minus_ref.*((x-x_ref)-(y-y_ref));
plus_square_linear=0.25*plus_ref.^2+0.5*plus_ref.*((x-x_ref)+(y-y_ref));

% 0.25*(x+y)^2-0.25*(x-y)^2-z <= 0  <=>  z>=x*y.
ge_residual=0.25*(x+y).^2-minus_square_linear-z;
% z+0.25*(x-y)^2-0.25*(x+y)^2 <= 0  <=>  z<=x*y.
le_residual=z+0.25*(x-y).^2-plus_square_linear;
end
