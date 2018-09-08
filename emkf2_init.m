function [Q,P,pi,F,H] = emkf2_init(v,q)

% STATE TRANSITION MATRIX
F=[1 1 0 0;
   0 1 0 0;
   0 0 1 1;
   0 0 0 1];

% MEASUREMENT MATRIX
H=[1 0 0 0;
   0 0 1 0];

nObjects = length(q);
Q = zeros(2,2,nObjects);
P = zeros(4,4,nObjects);
pi = zeros(nObjects,1);

for i = 1:nObjects
   % SYSTEM NOISE COVARIANCE MATRIX
   Q(:,:,i) = eye(2)*q(i);
   
   % INITIAL STATE COVARIANCE MATRICES
   P(:,:,i) = eye(4)*v(i);
   
   % WEIGHTS
   pi(i) = 1/nObjects;
end