function [probs] = single_prob_for_settings(d,state,alpha,beta)
% This function computes the ideal (noiseless) probability settings for 
% a given state and measurement angles

% Inputs: 
% - d: dimension of the system, i.e. the number of outcomes for each 
% measurement setting
% - state: matrix representing the ideal quantum state used
% - alpha: angle defining Alice's measurement
% - beta: angle defining Bob's measurement

% Outputs:
% - probs: matrix containing the probabilities of each pair of outcomes for
% Alice and Bob

% Compute the probability for each of the possible outcomes
probs = zeros(d,d); % initialise an empty matrix
for a = 0:(d-1) % loop over all possible outcomes for Alice
    for b = 0:(d-1) % loop over all possible outcomes for Bob
        % Obtain the quantum operators
        operator = AB_proj_operators(d,a,b,alpha,beta);  
        % Compute the mean value
        mel = trace(state*operator);
        probs(a+1,b+1) = mel;
    end
end
probs=real(probs);
end