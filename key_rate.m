function [kr,pa,ec] = key_rate(probsMeasured,probNL,locWeights,mA,mB,d,x,y)

% This function calculates an upper bound on the key rate given the
% probability distribution measured by Alice and Bob, the ideal noise-free
% probability distribution, the number of measurements mA and mB Alice and 
% Bob can perform, the dimension, and the key settings x and y used by 
% Alice and Bob

%Inputs:
% - probsMeasured: matrix whose rows correspond to the each of the observed
% probability vectors
% - probNL: probability vector of the noise-free non-local state
% - mA,mB: number of measurements available to Alice and Bob
% - d: dimension of the system
% - x,y: key settings used by Alice and Bob

% Outputs:
% - kr: upper bound on the key rate
% - pa: PA-term
% - ec: EC-term

n=size(probsMeasured,1); % number of observed probability vectors

pa = zeros(1,n); % initialise a zero vector
ec = zeros(1,n); % initialise a zero vector
for i=1:n
    % calculate the PA and EC terms for each of the probability vectors in 
    % probsMeasured
    [joint,~,margiB] = probs_for_settings(probsMeasured(i,:),mA,mB,d,x,y);
    [~,margiA_ideal] = probs_for_settings(probNL,mA,mB,d,x,y);
    pa(i) = paTerm(margiA_ideal,locWeights(i));
    ec(i) = ecTerm(margiB,joint);
end
% Calculate the key rates by subtracting the EC term from the PA term
kr=pa-ec;
end

