function [ec] = ecTerm(margiB,joint)

% Calculate EC-term following Appendix C.1.1 of the Lukanowski et al. "Upper 
% bounds on key rates in device-independent quantum key distribution based 
% on convex-combination attacks" paper

% Inputs:
% - margiB: Bob's marginal probabilities
% - joint: matrix containing Alice and Bob's joint probability distribution

% Outputs:
% - ec: EC-term

% Calculate conditional probabilities p_A|B(a|b)
conditionalProbs = joint.*(1./margiB);
    
% Calculate EC term (equation 68 of Lukanowski et al.)
ec = margiB*shannonEntropy(conditionalProbs');

end

