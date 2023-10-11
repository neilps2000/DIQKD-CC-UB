function [val] = bellValue(probs,mA,mB,d)

% This function calculates the value of the Collins-Gisin Bell expression 
% for a given probability vector (Equation 6 from the paper "Bell 
% inequalities for arbitrarily high dimensional systems" b y Collins 
% et al.)

% Inputs:
% - probs: probability vector for which we want to calculate the value of
% the Bell expression
% - mA,mB: number of measurement settings for Alice and Bob
% - d: dimension of the system, i.e. number of possible outcomes

% Outputs:
% - val: value of the Bell expression

% Define an anonymous function to get the joint probabilities for setting 
% x and y
probs_xy = @(x,y) probs_for_settings(probs,mA,mB,d,x,y);
% Initialise an anonymous function to get the sum of joint probabilities 
% for outcomes differing by k
probs_k = @(probs,k) 0; 
for j = 0:d-1 % loop over all possible outcomes
    % update probs_k by adding the joint probability for outcome j and 
    % j+k (mod d)
    probs_k = @ (probs,k) probs_k (probs,k) + probs(j+1,mod(j+k,d)+1); 
end

% get the joint probabilities for setting i (Alice) and j (Bob)
probs_11 = probs_xy(1,1);
probs_12 = probs_xy(1,2);
probs_21 = probs_xy(2,1);
probs_22 = probs_xy(2,2);

%Implement Bell expression from Collins-Gisin equation 6
val=0; % initialise an empty variable for val
for k = 0:(floor(d/2)-1 )
    val = val + (1-2*k/(d-1))*(probs_k(probs_11,k) ...
        + probs_k(probs_21,-k-1) ...
        + probs_k(probs_22,k) ...
        + probs_k(probs_12,-k) ...
        - probs_k(probs_11,-k-1) ...
        - probs_k(probs_21,k) ...
        - probs_k(probs_22,-k-1) ...
        - probs_k(probs_12,k+1)); 
end

end

