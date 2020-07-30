%Analytical calculations
function [kroTE1, kroTM1, kroTM0] =  krhoAnalyticalSubStrate(k0,lam,h,er)   
%     Assumptions
%     tan(k0h) = kz0h - npi
%     hd = lam./4
%     kro = 0
%     zl = zeta0./er
%     D = zeta0./er + jz0(kz0h + npi)
%     kz0 = jzeta0./(er*h*z0) + npi./h)
%     zeta0 = 120*pi;
    hbar = h./lam;

    kz0TE1 = (k0./(pi.*sqrt(er).*(2*hbar).^2)).*((2*pi.*hbar.*sqrt(er) + 1j)...
        ./(1 + (1./(((2*pi.*hbar).^2).*er))));
%     Z0TE1 = (zeta0.*k0)./(kz0TE1);

    kz0TM0 = sqrt((1j.*k0)./(sqrt(er).*h));
%     Z0TM0 = (zeta0.*kz0TM0)./k0;

    kz0TM1 = (k0./(4.*hbar)).*(1 + sqrt(1 + 8.*1j.*hbar./(pi.*sqrt(er))));
%     Z0TM1 = (zeta0.*kz0TM1)./k0;

%     DTE1 = zeta/er + 1j.*Z0TE1.*(kz0TE1.*h - pi);
%     DTM0 = zeta/er + 1j.*Z0TM0.*(kz0TM0.*h - pi);
%     DTM1 = zeta/er + 1j.*Z0TM1.*(kz0TM1.*h - pi);
    
    kroTE1 = sqrt(k0.^2 - kz0TE1.^2);
    kroTM1 = sqrt(k0.^2 - kz0TM1.^2);
    kroTM0 = sqrt(k0.^2 - kz0TM0.^2);
end