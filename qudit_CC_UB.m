% The code implements the CGLMP Bell inequality for a two-qudit system
% It involves measuring two observables per party, each with d possible 
% outcomes
% The code aims to find the measurement settings maximising the violation 
% of the inequality by the maximally entangled state

clear % clear the workspace

format long

% Initialise the parameters
d=3; % dimension of the system
m=2; % number of measurements used for the Bell inequality
mA=m; % number of measurements performed by Alice
mB=m+1; % number of measurements performed by Bob

% We define the CGLMP Bell operator in terms of generic projective 
% operators
% The projective operators are defined by the function AB_proj_operators, 
% which takes the dimension, the outcome indices, and the measurement 
% angles as inputs

% First, define operators op_xy_k(x,y,k,XX) as the sum of the 
% AB_proj_operators for measurement angles XX and settings x,y, and outcome 
% indices differing by k
op_xy_k = @(x,y,k,XX) 0; % initialise an empty matrix
for j = 0:d-1
    % XX is a 2 x m matrix of measurement angles for Alice and Bob
    op_xy_k = @ (x,y,k,XX) op_xy_k (x,y,k,XX) ...
    + AB_proj_operators(d,j,mod(j+k,d),XX(1,x),XX(2,y)); 
end
% The Bell operator is a linear combination of projective operators for 
% different measurements and outcomes
bellOperator=@(XX) 0; % initialise an empty matrix
for k = 0:floor(d/2)-1
    bellOperator = @(XX) bellOperator(XX) - (1-2*k/(d-1))*(...
        (op_xy_k(1,1,k,XX) ...
        + op_xy_k(2,1,-k-1,XX) ...
        + op_xy_k(2,2,k,XX) ...
        + op_xy_k(1,2,-k,XX)) ...
        - op_xy_k(1,1,-k-1,XX) ...
        - op_xy_k(2,1,k,XX) ...
        - op_xy_k(2,2,-k-1,XX) ...
        - op_xy_k(1,2,k+1,XX) ...
        ); 
end

measAngles = [[0,-1/2];[1/4,-1/4]]; % set measurement angles

% Create the density matrix for the maximally entangled state
% The state is given by |Psi> = (1/sqrt(d)) sum_j |jj>
maxEnt = zeros(d); % initialise an empty matrix
for i=1:d % loop over all possible outcomes
    maxEnt(i,i) = 1/sqrt(d); % assign the diagonal elements to 1/sqrt(d)
end
maxEnt = reshape(maxEnt,1,[]); % reshape the matrix into the ket column |Psi>
maxEnt = maxEnt'*maxEnt; % calculate the density matrix |Psi><Psi|

[eigenvectors,eigenvalues] = eig(bellOperator(measAngles));
% Find the minimum eigenvalue of a matrix eigenvalues and its corresponding 
% index
[M,I] = min( min(eigenvalues,[],'linear'));
% Compute the density matrix of the state maximally violating the Bell
% inequality
maxVio = eigenvectors(:,I)*eigenvectors(:,I)';

% Find the optimal measurement angle meas_angle for Bob to maximise the 
% correlation with Alice's measurement setting 1
% This is equivalent to minimising the EC-term, calculated with the 
% function ecTerm, which takes the marginal and joint probabilities as 
% inputs
% The probabilities are given by the function single_prob_for_settings, 
% which takes the dimension, the state, and the measurement angles as 
% inputs
% We use a numerical optimisation method fminunc to find the optimal 
% angle for Bob
measProbs = @(meas_angle) single_prob_for_settings(d,maxEnt,...
    measAngles(1,1),meas_angle);
ec = @(meas_angle) ecTerm(sum(measProbs(meas_angle),2)',...
    measProbs(meas_angle));
option = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-6,...
    'MaxIter', 3000) ;
[optMeasAngle, maxCorrelation] = fminunc(ec,rand,option);

% Define the measurement angles for Alice and Bob
alphas = measAngles(1,:); % Alice's angles are the first row of optAngles
betas = [measAngles(2,:),optMeasAngle]; % Bob's angles are the second row of
% optAngles plus the optimal angle

% Find the noiseless probabilities for all possible settings and outcomes
% The probabilities are given by the function all_probs, 
% which takes the dimension, the angles, and the state as inputs
% The function returns the joint probabilities for each pair of settings,
% and the marginal probabilities for Alice and Bob
[joint_MV,margiA_MV,margiB_MV] = all_probs(d,alphas,betas,maxVio);
[joint_ME,margiA_ME,margiB_ME] = all_probs(d,alphas,betas,maxEnt);

% Define the probability vector, which eliminates redundancy
% The probability vector is a one-dimensional array that contains only the 
% marginal probabilities for outcomes 1 to d-1 for each setting for Alice 
% and Bob, and the joint probabilities for outcomes 1 to d-1 for each pair 
% of settings
% The probabilities for outcome d can be obtained by subtracting the sum of 
% the other outcomes from 1
probMV = zeros(1,(mA+mB)*(d-1)+mA*mB*(d-1)^2); % initialise empty array
allJointProbs = zeros((d-1)*mA,(d-1)*mB); % initialise an empty matrix
for x=1:mA % loop over all possible settings for Alice
    for y=1:mB % loop over all possible settings for Bob
        joint_xy = joint_MV(string(x)+string(y)); % get the joint probability 
        % matrix for setting x and y
        % Assign the submatrix of allJointProbs corresponding to
        % measurement settings x and y
        allJointProbs((x-1)*(d-1)+1:x*(d-1),(y-1)*(d-1)+1:y*(d-1)) = ...
            joint_xy(1:d-1,1:d-1);
    end
end
% Assign the flattened matrix of allJointProbs to the last part of 
% probNL
probMV((mA+mB)*(d-1)+1:end) = reshape(allJointProbs,1,[]);

% Assign the marginal probabilities to the corresponding part of probNL
probMV( 1:mA*(d-1) ) = reshape(margiA_MV(1:d-1,:), 1 ,[]);

% Assign the marginal probabilities to the corresponding part of probIdeal
probMV(mA*(d-1)+1:(mA+mB)*(d-1)) = reshape(margiB_MV(1:d-1,:), 1 ,[]);


probME = zeros(1,(mA+mB)*(d-1)+mA*mB*(d-1)^2); % initialise empty array
allJointProbs = zeros((d-1)*mA,(d-1)*mB); % initialise an empty matrix
for x=1:mA % loop over all possible settings for Alice
    for y=1:mB % loop over all possible settings for Bob
        joint_xy = joint_ME(string(x)+string(y)); % get the joint probability 
        % matrix for setting x and y
        % Assign the submatrix of allJointProbs corresponding to
        % measurement settings x and y
        allJointProbs((x-1)*(d-1)+1:x*(d-1),(y-1)*(d-1)+1:y*(d-1)) = ...
            joint_xy(1:d-1,1:d-1);
    end
end
% Assign the flattened matrix of allJointProbs to the last part of 
% probNL
probME((mA+mB)*(d-1)+1:end) = reshape(allJointProbs,1,[]);

% Assign the marginal probabilities to the corresponding part of probNL
probME( 1:mA*(d-1) ) = reshape(margiA_ME(1:d-1,:), 1 ,[]);

% Assign the marginal probabilities to the corresponding part of probIdeal
probME(mA*(d-1)+1:(mA+mB)*(d-1)) = reshape(margiB_ME(1:d-1,:), 1 ,[]);

% Find the local strategies that achieve the local bound
% The local strategies are given by the function probLoc, which takes the 
% number of measurements and the dimension as inputs
probL = probLoc(mA,mB,d);

% Calculate the critical visibility
% The visibility is a measure of how much noise is added to the state
% The critical visibility is the minimum visibility required achieve a
% positive bound on the key rate
Vrange = [0.75 0.825]; % set a range of visibilities to search
tol = 1e-8; % set a tolerance for the error
critV = critical_visibility(probMV,probMV,probL,mA,mB,d,1,mB,Vrange,tol);
disp(critV); % display the result

% Plot the upper bound on the key rate
Vrange = linspace(0.7,1,50); % define the range of visibilities to plot
plot_CC_UB(probMV,probMV,probL,mA,mB,d,1,mB,Vrange)
