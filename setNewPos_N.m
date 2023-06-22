function UAVs = setNewPos_N(nUAV, Ang, h, Rj, typeA)
    
    UAVs = zeros(nUAV, 3);
    
    switch typeA
        case 1
            % Angle of first UAV that preserves symmetry
            phi = pi - ((nUAV-1)/2)*Ang;                    
            for iU = 1:nUAV
                UAVs(iU,1) = Rj*cos(phi + (iU-1)*Ang);
                UAVs(iU,2) = Rj*sin(phi + (iU-1)*Ang);
                UAVs(iU,3) = h;
            end
        case 2
            % For the ZF Pattern
            phi = pi - ((nUAV-1))*Ang;
            for iU = 1:nUAV
                UAVs(iU,1) = Rj*cos(phi + (nUAV-iU)*Ang);
                UAVs(iU,2) = Rj*sin(phi + (nUAV-iU)*Ang);
                UAVs(iU,3) = h;
            end
    end
end