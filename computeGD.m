function Ang_Step = computeGD(A, E, Rj, hj, dAB, gammaA, gammaJ, channelParam, alphN, Angle,n, nUAV, typeA )
    % n:    time step
    % tht:  current theta value
    
    %dTht = (5e-3)*Angle;      % Infinitesimal approximation for derivative
    dTht = (5e-2)*Angle;      % Infinitesimal approximation for derivative
    
    UAVs_FD     = setNewPos_N(nUAV, Angle+dTht, hj, Rj, typeA);
    UAVs_F      = setNewPos_N(nUAV, Angle, hj, Rj, typeA);
    
%     WSC_FD = computeWSC_NSMRT_NUAV(A, E, UAVs_FD, dAB, gammaA, gammaJ, channelParam );
%     WSC_F  = computeWSC_NSMRT_NUAV(A, E, UAVs_F, dAB, gammaA, gammaJ, channelParam );
    WSC_FD = computeWSC_ZF_NUAV(A, E, UAVs_FD, dAB, gammaA, gammaJ, channelParam );
    WSC_F  = computeWSC_ZF_NUAV(A, E, UAVs_F, dAB, gammaA, gammaJ, channelParam );
    
    
    WSC_DEV = real((WSC_FD-WSC_F)/dTht);
    
    if alphN==0
        alphN = 1/n;
    end
    %alphN = 1/n;
    
    Ang_Step    = Angle - alphN*WSC_DEV;
    
    Ang_Step = mod(Ang_Step,2*pi);              % Angle 
    
end