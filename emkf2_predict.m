function [x,P] = emkf2_predict(x,F,P,Q) %#codegen

Ga=[1/2 1 0 0;0 0 1/2 1]';
fading = 0;

nObjects = size(x,2);

%xp = zeros(size(x,1),nObjects);
%Pp = zeros(size(P,1),size(P,2),nObjects);

% if nObjects == 0
%    return;
% end

%%  Kalman Filter Time Update ("Predict")
for i = 1:nObjects
   [x(:,i), P(:,:,i)] = KalmanPredict(F,x(:,i),P(:,:,i),Ga,Q(:,:,i),fading);
end

% tic
% P2=squeeze(num2cell(P,[1 2]));
% 
% x2=num2cell(x,1);
% x2=reshape(x2,nObjects,1);
% 
% Ga2=repmat(Ga,[1 1 nObjects]);
% Ga2=num2cell(Ga2,[1 2]);
% Ga2=squeeze(Ga2);
% 
% Q2=squeeze(num2cell(Q,[1 2]));
% 
% F2=repmat(F,[1 1 nObjects]);
% F2=squeeze(num2cell(F2,[1 2]));
% 
% tmp=num2cell(repmat(0,nObjects,1));
% 
% [xo Po]=cellfun(@KalmanPredict,F2,x2,P2,Ga2,Q2,tmp,'UniformOutput',false);
% toc
% 
% pause;
end