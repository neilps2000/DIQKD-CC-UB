function [X] = canonicalBasis(v,n)

% This function creates a matrix with all the possible combinations of
% n elements from v

% Inputs:
% - v: vector containing all the elements that we want to use
% - n: number of elements from v to put in each combination

% Outputs:
% - X: matrix with rows corresponding to all possibl combination of n
% elements from v

b=size(v,2);

for i=1:n
    X(:,i)=reshape(repmat(v,[b^(n-i),b^(n-1)/(b^(n-i))]),[],1);
end