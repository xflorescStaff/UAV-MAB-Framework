function [WSC_Max, Ang_Max, DeltaComp_Max] = optimalWSC_NOP(A, E, Rj, hj, dAB, gammaA, gammaJ, channelParam, thetaV, nUAV, typeA )
    WSC_Max = 0;
    Ang_Max = 0;
    DeltaComp_Max = 0;
    for ang=thetaV
        UAVs = setNewPos_N(nUAV, ang, hj, Rj, typeA);
        [WSC_Dum, DeltaComp] = computeWSC_NOP_NUAV(A, E, UAVs, dAB, gammaA, gammaJ, channelParam );
%         scatter(UAVs(:,1), UAVs(:,2)); xlim([-Rj, Rj])
%         fprintf('ang %.2f\t WSC %.2f\n',ang,WSC_Dum)
        if WSC_Dum > WSC_Max
            WSC_Max = WSC_Dum;
            Ang_Max = ang;
            DeltaComp_Max = DeltaComp;
        end 
    end 
%     fprintf('end\n')
end