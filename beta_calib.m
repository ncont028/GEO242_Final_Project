function [x] = beta_calib(in1)
%alpha_calib - calibration for the alpha parameter
% Uses least-square to find the best fit
% Tr_in : recurrence intervals [yr] (column vector)
% M0_in : seismic moments [dyne-cm] (column vector)
% alpha : alpha parameter
%%-------------------------------------------------------------------------
% Chen, K. H., Nadeau, R. M., and Rau, R.-J. (2007), 
% Towards a universal rule on the recurrence interval scaling of repeating 
% earthquakes? Geophys. Res. Lett., 34, L16308, doi:10.1029/2007GL030554.
%%-------------------------------------------------------------------------
% solving using linear squares regression
% Ax = b
% From Nadeau (1998) eqn. : log(Tr) = beta_T*log(M0) + alpha_T --> y=mx+b
% A matrix : log(M0) values
% x : beta_T & constant
% b : log(Tr)
% If the average rate  of seismic  slip is indpendent of M0:
%       beta_T = beta_d  (Nadeau 1998, eqn. 19)
% If this is not the case, then eqn. 15 is an approximation.

Tr_in = in1(:,1);
M0_in = in1(:,2);

A = [log10(M0_in),ones(length(M0_in),1)];
b = log10(Tr_in);
x = A\b;
% beta=x(1);
% alpha_T = x(2);

end