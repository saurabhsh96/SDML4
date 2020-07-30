%Routine to find the propogation constant
function kRho = findpropSubStrate(k0, er, h, kgRho, zeta0, mode)
     
    %delK
    dK = k0./500;
    
    %Derivative
    dKg = Den_Substrate(k0, er, h, kgRho, zeta0, mode);
    dKgDelP = Den_Substrate(k0, er, h, (kgRho + dK/2), zeta0, mode);
    dKgDelM = Den_Substrate(k0, er, h, (kgRho - dK/2), zeta0, mode);
    dKgDel = (dKgDelP - dKgDelM)./dK;
    
    %kRho
    kRho = kgRho - dKg./dKgDel;
end