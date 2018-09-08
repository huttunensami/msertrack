function [ x, P, K ] = KalmanCorrect( xp, Pp, H, z, R ) %#codegen
%Kalman Filter Measurement Update ("Correct")
%
% disp('Kalman filter: ');
% tic
%% Kalman Gain
K = Pp*H'/(H*Pp*H' + R);
      
%% a posteriori state estimate
x = xp + K*(z(:)-H*xp);
      
%% a posteriori error covariance estimate 
P = (eye(length(x))-K*H)*Pp;

% toc
% 
% pause(0.5)
% disp('Information filter: ');
% tic
% %% Information filter
% Yp=inv(Pp);
% yp=inv(Pp)*xp;
% I=H'*inv(R)*H;
% i=H'*inv(R)*z(:);
% Y=Yp+I;
% y=yp+i;
% P2=inv(Y);
% x2=P2*y;
% toc