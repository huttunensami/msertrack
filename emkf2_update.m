function [xp,Pp,pi,w] = emkf2_update(xp,Pp,H,z,pi,R,weightThreshold) %#codegen

%% Get weights
[mu,C,w] = emkf2_getWeights(xp,Pp,H,z,pi,R);
%% CONSTANTS
a = 0.95;

nObjects = size(xp,2);

if ~nObjects
   %       x = xp;
   %       P = Pp;
   %       pi = [];
   %       w = [];
   return;
end

%% Kalman Filter Measurement Update ("Correct")
for i = 1:nObjects
   %% Remove negligible measurements
   ind = w(:,i)>weightThreshold;
   weights = w(ind,i);
   zw = z(:,ind);
   
   n = size(zw,2);
   HM = kron(ones(n,1),H);
   
   %% UPDATE FILTER
   RM = kron(inv(diag(weights(:)+1e-8)),R(:,:,i));
   [ xp(:,i), Pp(:,:,i), ~] = KalmanCorrect( xp(:,i), Pp(:,:,i), HM, zw, RM );
   
   %% Update Weights
   pi(i) = a*pi(i)+(1-a)*mean(w(:,i));
end
end