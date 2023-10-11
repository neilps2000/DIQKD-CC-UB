function [H,error] = shannonEntropy(probs)

% This function calculates the Shannon entropy of a set of probabilities 

% Inputs: 
%  - probs: a matrix of probabilities, each row represents a discrete 
% probability distribution 

% Outputs: 
% - H: a vector of Shannon entropy values, one for each row of probs 
% - error: a flag indicating whether there was an error in the input

error = 0; % initialise error flag to zero
n = size(probs,1); % get the number of rows in probs
H = zeros(n,1); % initialise H vector to zeros
%tol = 1e-8; % set a numerical tolerance

for i=1:n % loop over each row of probs
    % if abs(sum(probs(i,:)) - 1) > tol % check if the probabilities sum to 1
    %     error = 1;
    %     disp("Something went wrong: probabilities don't sum to 1");
    %     disp(probs(i,:));
    %     break;
    % end
    prob = probs(i,:); % get the ith row of probs
    prob = prob(prob>0); % remove zero probabilities
    H(i) = - prob*(log2(prob))'; % calculate the Shannon entropy
end

