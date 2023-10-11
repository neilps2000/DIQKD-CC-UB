function [] = plot_CC_UB(probNL,probL,mA,mB,d,x,y,Vrange)

% This function plots the key rate for a specified range of visibilities,
% given an ideal noise-free probability vector

% Inputs:
% - probNL: noise-free non-local probability distribution
% - probL: local deterministic strategies used for the decomposition
% - mA,mB: number of measurements Alice/Bob
% - d: dimension of the system
% - x,y: measurement settings for Alice/bob
% - Vrange: vector of visibilities for which the key rate is calculated

% Calculate the observed probabilities for different visibilities
probsMeasured = Vrange'*probNL+(1-Vrange')*[ones(1,(mA+mB)*(d-1))*1/d,...
    ones(1,mA*mB*(d-1)^2)*1/(d^2)];

% Find the local weight corresponding to each of the observed probabilities
% by calling the maxLocWeights function
[locWeights,error] = maxLocWeights(probNL,probL,probsMeasured,mA,mB,d);

if error == 1 % check for errors
    disp("Something went wrong: maxLocWeights");
else % calculate and plot the upper bound on the key rate as a function 
    % of the visibility
    kr = key_rate(probsMeasured,probNL,locWeights,mA,mB,d,x,y);
    plot(Vrange, kr);
    title('UB on key rate ('+string(mA)+string(mB)+string(d)+string(d)...
        +'-protocol)')
    xlabel('V');
    ylabel('r');
end

end

