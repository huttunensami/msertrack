function [xp Pp] = KalmanPredict(Fi,x,P,Ga,Q,fadingFactor)
%Kalman Filter Time Update ("Predict")
%
%% a priori state estimate
    xp=Fi*x;
        
%% a priori error covariance estimate
    Pp=Fi*P*Fi'*exp(fadingFactor)+Ga*Q*Ga';