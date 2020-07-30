%Routine to find the propogation constant
function kRho = findpropSuperStrate(k0, er, h, hs, kgRho, zeta0, mode)
     
    %delK
    dK = k0./500;
    
    %Derivative
    dKg = Den_SuperStrate(k0, er, h, hs, kgRho, zeta0, mode);
    dKgDelP = Den_SuperStrate(k0, er, h, hs, (kgRho + dK/2), zeta0, mode);
    dKgDelM = Den_SuperStrate(k0, er, h, hs, (kgRho - dK/2), zeta0, mode);
    dKgDel = (dKgDelP - dKgDelM)./dK;
    
    %kRho
    kRho = kgRho - dKg./dKgDel;
end