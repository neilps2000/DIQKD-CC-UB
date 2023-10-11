function [pa] = paTerm(margiA_ideal,locWeight)

%Calculate PA term of the key rate for measurement settings x,y and 
% dimension d, following Appendix C.2.1 of the Lukanowski et al. "Upper 
% bounds on key rates in device-independent quantum key distribution based 
% on convex-combination attacks" paper

% Inputs:
% - margiA_ideal: Alice's marginal probabilities corresponding to the 
% ideal noise-free non-local state
% - locWeight: local weight of the convex decomposition of the observed
% probabilities

% Outputs:
% - pa: the PA term of the upper bound on the key rate 

pa = (1-locWeight)*shannonEntropy(margiA_ideal);

end

