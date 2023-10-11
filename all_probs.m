function [joint,margiA,margiB] = all_probs(d,alphas,betas,state)
% Given a set of settings {A_1,A_2,...,A_m} and {B_1,B_2,...,B_m} for Alice
% and Bob, this function returns the joint probabilities for each pair of 
% settings, and the marginal probabilities for Alice and Bob

% Inputs: 
% - d: dimension of the system, i.e. number of possible outcomes
% - alphas: angles defining Alice's measurments, 
% i.e.[\alpha_1,...,\alpha_m]
% - betas: angles defining Bob's measurments, [\beta_1,...,\beta_m]

% Outputs:
% - joint: dictionary whose keys are Alice and Bob's measurement settings. 
% If we want to know the joint distribution we get if Alice performs
% measurement A_x and Bob performs measurement B_y, we have to write 
% [joint,margiA,margiB] = all_probs(d,alphas,betas,state). Then, 
% joint("xy") is a d x d matrix whose element [a,b] is p(a,b|x,y).
% For this to work, we hava to write prob_xy = prob("xy") and then
% prob_xy(a,b). Notice that "xy" is a string.
% - margiA,margiB: Alice and Bob's marginal pribabilities

mA=length(alphas); % number of measurements (Alice)
mB=length(betas); % number of measurements (Bob)

% Define the dictionary and matrices to save the probabilities
joint = containers.Map;
margiA = zeros(d,mA);
margiB = zeros(d,mB);

% Go through all possible measurement settings and compute the
% corresponding probability
for x=1:1:mA
    for y=1:1:mB
        joint(string(x)+string(y)) = ...
            single_prob_for_settings(d,state,alphas(x),betas(y));
    end
end
for x=1:mA
    margiA(:,x) = sum(joint(string(x)+string(1)),2);
end
for y=1:mB
    margiB(:,y) = sum(joint(string(1)+string(y)),1);
end

end

