%Plotting and max function kRho
function [Dlargest, kRReq] = kRhoReq(k0, er, h, hs, kRho, zeta0, mode, num)
    Dlargest = Den_SuperStrate(k0, er, h, hs, kRho, zeta0, mode);
    figure();
    surface(real(kRho)./k0, imag(kRho)./k0, db(abs(1./Dlargest)./max(max(abs(1./Dlargest)))), 'linestyle','none');
    colorbar;
    title('1/Denominator function to identify singularities ' + mode + num2str(num));
    xlabel('Re(k_\rho)/k_0');
    ylabel('Im(k_\rho)/k_0');
    [row, col] = find(ismember((1./Dlargest), max(1./Dlargest(:))));
    kRReq = kRho(row, col);
end