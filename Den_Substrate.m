%Ground slab denominator function
function D = Den_Substrate(k0, er, h, kRho, zeta0, mode)
    
    %In Slab
    ks = sqrt(er).*k0;
    zetaS = zeta0./sqrt(er);
    
    %kZ %Top Riemann sheet
    kz0 = 1j*sqrt(-((k0^2)-(kRho.^2)));
    kzs = 1j*sqrt(-((ks^2)-(kRho.^2)));

    %Defining according to the mode
    if(mode == "TE")
        Z0 = (zeta0.*k0)./kz0;
        Zs = (zetaS.*ks)./kzs;
        %ZL = Zs.*(Z0 + 1j.*Zs.*tan(kzs.*hs))./(Zs + 1j.*Z0.*tan(kzs.*hs));
        Zup = Zs;
        Zdown = 1j*Z0.*tan(kz0.*h);
    else 
        Z0 = (zeta0.*kz0)./k0;
        Zs = (zetaS.*kzs)./ks;
%         ZL = Zs.*(Z0 + 1j.*Zs.*tan(kzs.*hs))./(Zs + 1j.*Z0.*tan(kzs.*hs));
        Zup = Zs;
        Zdown = 1j*Z0.*tan(kz0.*h);
    end 
    
    %Finding Dispersion equation
    D = Zup + Zdown;
end