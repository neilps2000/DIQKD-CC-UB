function [critV,locWeight,out] = critical_visibility(probIdeal,probNL,...
    probL,mA,mB,d,x,y,Vrange,tol)


% This function computes the critical visibility

% Inputs: 
% - probNL: a vector containing the ideal noise-free probability
% distribution
% - probL: a matrix with each row corresponding to a local deterministic 
% strategy 
% - mA: number of measurement settings for Alice 
% - mB: number of measurement settings for Bob 
% - d: dimension of the system, i.e. the number of outcomes for each 
% measurement setting 
% - x,y: key settings for Alice and Bob
% - Vrange: a vector containing two elements corresponding to the initial 
% lower and upper bounds for the critical visibility 
% - tol: tolerance for the key rate function 

% Outputs: 
% - critV: critical visibility 
% - locWeight: local weight in the convex decomposition of the observed
% probability vector
% - out: a flag indicating whether the computation was successful (1) or 
% not (2)

out = 0; % initialize the output flag
iter=0; % initialize the iteration counter
maxIter = 100; % set the maximum number of iterations

% Compute the observed probability vectors for the lower and upper bounds 
% of Vrange
probsMeasured = Vrange'*probIdeal+(1-Vrange')*[ones(1,(mA+mB)*(d-1))*1/d,...
    ones(1,mA*mB*(d-1)^2)*1/(d^2)];

% Call the maxLocWeights function to find the local weights for each 
% observed probability vector
[locWeight,error] = maxLocWeights(probNL,probL,probsMeasured,mA,mB,d);
if error == 1 % check if there was an error in maxLocWeights
    disp("Something went wrong: maxLocWeights");
else % if no error, proceed with the computation
    % Call the key_rate function to find the key rates for each observed 
    % probability vector
    kr = key_rate(probsMeasured,probNL,locWeight,mA,mB,d,x,y);
    while out==0
        % Compute the mean of Vrange as the current estimate of critical 
        % visibility
        critV = mean(Vrange);
        % Compute the observed probability vector for the current estimate 
        % of critical visibility
        probsMeasured = critV*probIdeal+(1-critV)*...
            [ones(1,(mA+mB)*(d-1))*1/d,ones(1,mA*mB*(d-1)^2)*1/(d^2)];
        % Call the maxLocWeights function to find the local weight for the 
        % current observed probability vector
        [locWeight,error] = maxLocWeights(probNL,probL,probsMeasured,...
            mA,mB,d);
        if error == 1 % check if there was an error in maxLocWeights
            disp("Something went wrong: maxLocWeights");
        else % if no error, proceed with the computation
            % Call the key_rate function to find the key rate for the 
            % current observed probability vector
            kr_new = key_rate(probsMeasured,probNL,locWeight,mA,mB,d,x,y);
            % Check if the absolute value of key rate is less than the
            % tolerance
            if abs(kr_new) < tol
                    out = 1; % set out to 1 to indicate success
            else % if not, update Vrange and kr accordingly
                if kr(1)*kr_new <= 0 % check if there is a sign change 
                    % between kr(1) and kr_new
                    kr(2) = kr_new; % update kr(2) to kr_new
                    Vrange(2) = critV; % update Vrange(2) to critV
                    critV = mean(Vrange); % update critV to mean of Vrange
                elseif kr(2)*kr_new <= 0 % check if there is a sign change 
                    % between kr(2) and kr_new
                    kr(1) = kr_new; % update kr(1) to kr_new
                    Vrange(1) = critV; % update Vrange(1) to critV
                    critV = mean(Vrange); % update critV to mean of Vrange
                else % if neither case holds, something went wrong
                    disp("Something went wrong: both key rates have the same sign")
                    out = 2; % set out to 2 to indicate failure
                end
            end
        end
        
        iter=iter+1; % increment the iteration counter
        if iter>=maxIter % check if the maximum number of iterations is 
            % reached
            disp("Something went wrong: maximum iterations exceeded");
            out = 2; % set out to 2 to indicate failure
        end

    end

end

end

