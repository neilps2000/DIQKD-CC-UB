function [id] = maximalViolation(d)

% This funciton gives the value of the maximal violation of the CGLMP 
% expression using the maximally entangled state in dimension d, following
% the paper "Bell inequalities for arbitrarily high dimensional systems" by
% Collins et al.

% Inputs:
% - d: dimension of the system

% Outputs:
% - id: aximal violation of the CGLMP expression

q = @(c) 1/(2*d^3*sin(pi*(c+0.25)/d)^2);

id=0;

for k=0:floor(d/2)-1
    id = id + (1-2*k/(d-1))*(q(k)-q(-(k+1)));
end

id = 4*d*id;

end

