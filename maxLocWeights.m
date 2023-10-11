function [q,error] = maxLocWeights(probNL,probL,probsMeas,mA,mB,d)

% This function computes the maximum local weight q of the 
% decomposition of a set of probability vectors into a combination of a set 
% of local deterministic strategies and a given non-local strategy

% Inputs:
% - probNL: a well chosen vector of non-local probabilities
% - probL: the set of local deterministic strategies
% - probsMeas: the set of measured probabilities
% - mA,mB: number of measurements Alice/Bob
% - d: dimension of the system 

% Outputs:
% - q: maximum local weight of the decomposition
% - error: flag used for error management

error = 0; % initialise a variable for error management
q=zeros(1,size(probsMeas,1)); % initialise the local weights to zero

% Construct a matrix M by concatenating the matrix probL containing the 
% local deterministic with probNL
% This matrix represents the linear constraints for the 
% optimisation problem (i.e the distributed probability has to be a 
% linear combination of these probabilities)
M=[probL;probNL];
s = size(M,1); % number of rows in M

x=sdpvar(s,1); % define optimisation variable (sx1 column vector)
% corresponding to the probability of distributing each local/non-local
% probability vector

% We define the matrix A in order to sum the probabilities of 
% distributing a local probability vector (all but the last element of x)
A=eye(s);
A(s,s)=0;

obj=sum(A*x); % define the objective function of the optimisation problem

% Restrictions on the linear optimisation problem
F1 = x>=0; %The probabilities must be positive
F2 = sum(x)==1; %The probabilities must sum to one
F0 = F1+F2; 

for i=1:size(probsMeas,1) % loop over each of the observed probabilities
    % First, check if probsMeas(i,:) is local
    % If it is, the local weight is 1
    if bellValue(probsMeas(i,:),mA,mB,d) <= 2
        q(i) = 1;
    else
        F3 = (M'*x)==probsMeas(i,:)'; % the distributed probability must 
        % coincide with the measured probability
    
        F = F0+F3; % add this constraint to the previous ones
        
        % Define the optimisation problem using the function optproblem
        problem = optproblem(F,obj,...
            sdpsettings('solver','mosek','verbose',0));
        
        % Solve the problem and obtain diagnostics
        diagnostics = maximize(problem);
        if diagnostics.problem == 0
            % disp('MaxLocWeights: Solver thinks it is feasible');
            q(i) = value(obj); % set the local weight equal to the value of
            % the objective function
        elseif diagnostics.problem == 1
            disp('MaxLocWeights: Solver thinks it is infeasible')
            error = 1;
        else
            disp('MaxLocWeights: Something else happened')
            error = 1;
        end
    end
end
end

