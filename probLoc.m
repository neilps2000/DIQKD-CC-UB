function [dd] =probLoc(mA,mB,d)

% This function computes the local strategies that will be used for the
% decomposition
% It then calculates the value of the Bell expression for each strategy 
% using the function bellValue and then discards the strategies that are 
% don't achieve the local bound, which is 2 for the CGLMP inequality

% Inputs:
% - mA,mB: number of measurement settings for Alice and Bob
% - d: dimension of the system, i.e. number of possible outcomes

% Outputs:
% - dd: matrix in which each row represents a deterministic strategy
% [As1o1,As1o2,...As1o(nboutput-1),..., Asnbsettingso1,As1o2,...,
% AsnbSettingso(nbOutput-1),Bs1o1,Bs1o2,...,Bs1o(n-1),...,As1o1*Bs1o1,
% As1o1*Bs1o2,...,As1o1*BsnbSettingsonbOutput,As1o2*Bs101,...,...] 
% where Asioj is the probability of Alice obtaining the jth output for the 
% ith settings.

% We use the function canonicalBasis to obtain all possible outcomes for
% Alice
A=canonicalBasis([0,1],d-1);
% We remove any unormalised probabilities by setting rows summing to more 
% than one to be empty
A(sum(A,2)>1,:)=[];

% We repeat this matrix to account all of Alice's measurement settings
for t=1:mA
    AA=repmat(A,1,size(A,1)^(mA-1)/(size(A,1)^(t-1)));
    AA=reshape(AA',d-1,[])';
    AA=repmat(AA,(size(A,1))^(t-1),1);
    AAA(:,:,t)=AA;
end
% fffA contains all possible deterministic strategies for Alice
fffA = reshape(AAA,size(AAA,1),[],1);

% We go through the same process for Bob
for t=1:mB
    BB=repmat(A,1,size(A,1)^(mB-1)/(size(A,1)^(t-1)));
    BB=reshape(BB',d-1,[])';
    BB=repmat(BB,(size(A,1))^(t-1),1);
    BBB(:,:,t)=BB;
end
% fffB contains all all possible deterministic strategies for Bob
fffB = reshape(BBB,size(BBB,1),[],1);

% We create all the joint strategies for Alice and Bob
gg=repmat(fffA,1,size(fffB,1));
gg=reshape(gg',mA*(d-1),[])';
l=repmat(fffB,size(fffA,1),1);
M=[gg,l];

% We calculate the associated joint probabilities
dd=zeros(size(M,1),(mA+mB)*(d-1)+mA*mB*(d-1)^2);
dd(:,1:(mA+mB)*(d-1))=M;

for i=1:size(dd,1)
    dd(i,(mA+mB)*(d-1)+1:size(dd,2))=reshape(   (  M(i,1:mA*(d-1))'...
        * M(i,mA*(d-1)+1:(mA+mB)*(d-1))  ),1,[]  );
end

% We calculate value of Bell expression for each local strategy
bell = zeros(size(dd,1),1);
for i=1:size(dd,1)
    bell(i) = bellValue(dd(i,:),mA,mB,d);
end
% We discard strategies that don't achieve the local bound
tol = 1e-8; % Set a tolerance for numerical error
dd = dd(bell > 2-tol,:); % For the CGLMP expression, the local bound is 2
end

