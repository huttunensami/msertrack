function [mu,C,w] = emkf2_getWeights(xp,Pp,H,z,pi,R) %#codegen
 
coder.extrinsic('clear')

   nObjects = size(xp,2);
   
   if ~nObjects
      mu = [];
      C = [];
      w = [];
      return;
   end
   
   nMeasurements = size(z,2);
   mu = zeros(size(z,1),nObjects);
   C = zeros(size(mu,1),size(mu,1),nObjects);
   %w = zeros(nMeasurements,nObjects);
   p = zeros(nMeasurements,nObjects);
       
%% Calculate Weights
   parfor i = 1:nObjects
      mu(:,i) = H*xp(:,i);
      C(:,:,i) = H*Pp(:,:,i)*H' + R(:,:,i);
      %obj = gmdistribution(mu(:,i)',C(:,:,i));
      %p(:,i) = pdf(obj,z');
      p(:,i) = mdgauss(z,mu(:,i),C(:,:,i));
   end
   
   %p = p./max(p(:));
   
   pp = p*pi;
   
%    for j=1:nMeasurements
%       w(j,:) = pi'.*p(j,:)./(pp(j)+0.0000001);
%    end
   
   pps = repmat(pp,1,nObjects)+0.0000001;
   clear pp;
   pis = repmat(pi',nMeasurements,1);
   pis = pis.*p;
   clear p;
   %w = pis.*p./(pps);
   w = pis./pps;
   
   %w = bsxfun(@rdivide,w,max(w,[],1));
   