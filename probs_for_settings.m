function [joint,margiA,margiB] = probs_for_settings(probs,mA,mB,d,x,y)

% Returns marginal probabilities for Alice and Bob, as well as joint
% probability matrix for given measurement settings, given a row of
% probabilities probs

% Inputs:
% - probs: a row vector of probabilities for all possible outcomes 
% - mA: the number of measurement settings for Alice 
% - mB: the number of measurement settings for Bob 
% - d: the dimension of the system 
% - x: the measurement setting chosen by Alice 
% - y: the measurement setting chosen by Bob

% Outputs:
% - joint: matrix containing joint probabilities for each pair of outcomes
% for Alice and Bob
% - margiA: Alice's marginal probbilities
% - margiB: Bob's marginal probabilities

margiA = zeros(1,d); % initialise an empty matrix
margiB = zeros(1,d); % initialise an empty matrix

% Extract Alice's marginal probabilities for measurement setting x
margiA(1:d-1) = probs(1+(x-1)*(d-1):x*(d-1));
margiA(d) = 1-sum(margiA(1:d-1));

% Extract Bob's marginal probabilities for measurement setting y
margiB(1:d-1) = probs(mA*(d-1)+1+(y-1)*(d-1):mA*(d-1)+y*(d-1));
margiB(d) = 1-sum(margiB(1:d-1));

% Reshape joint probabilities into matrix form
joint = reshape(probs((mA+mB)*(d-1)+1 : end),mA*(d-1),mB*(d-1));

% Remove the elements of the matrix not corresponding to measurement 
% settings
joint = joint((x-1)*(d-1)+1:x*(d-1),:);
joint = joint(:,(y-1)*(d-1)+1:y*(d-1));
% Complete the joint probability matrix with missing elements corresponding
% to the last outcome, which is omitted in probs
joint = [joint, margiA(1:d-1)' - sum(joint,2)];
joint = cat( 1, joint, cat(2,margiB(1:d-1) - ...
    sum(joint(:,1:end-1),1),margiB(d) - sum(joint(:,end))) );

end

