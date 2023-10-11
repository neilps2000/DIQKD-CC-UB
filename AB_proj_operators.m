function [operator,A_op,B_op] = ...
    AB_proj_operators(d,a,b,alpha,beta)

% This function gives the projective operators for each outcome, which
% define the measurements. That is, Alice's measurement A_i is defined by
% the projective operators {\Pi_1,...,\Pi_d} (similarly for Bob)

% Inputs:
% - d: dimension of Alice and Bob's Hilbert spaces
% - a: the projector we want to calculate for Alice i.e. the a in \Pi_a
% - b: the same for Bob
% - alpha: angle defining Alice's measurement
% - beta: angle defining Bob's measurements

% Outputs:
% operator: matrix corresponding to the operator acting on the whole system
% A_op,B_op: matrices corresponding to the operators acting on Alice and
% Bob's parts of the system seperately

% Define the eigenstates
eigenstate_A = 0; % initialise an empty vector
eigenstate_B = 0; % initialise an empty vector
% Compute the eigenstates
for q = 1:d 
    % Define "Fock" component
    state = zeros(d,1); % initialise an empty vector
    state(q,1) = 1; % assign the q-th element to 1
    % add the weighted phase factor to eigenstate_A
    eigenstate_A = eigenstate_A + ...
        sqrt(1/d)*exp((2*pi*1j/d)*(q-1)*(a+alpha))*state;
    % add the weighted phase factor to eigenstate_B
    eigenstate_B = eigenstate_B + ...
        sqrt(1/d)*exp(-(2*pi*1j/d)*(q-1)*(b+beta))*state;
end
% Compute the projectors
A_op = eigenstate_A*eigenstate_A'; % use outer product to get A_op
B_op = eigenstate_B*eigenstate_B'; % use outer product to get B_op

% Compute the total projector
operator = kron(A_op,B_op); % use Kronecker product to get the total 
% operator

end