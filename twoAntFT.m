%Fuction for current FT of 2 antennas
function CurrFT = twoAntFT(k0, kl, kw, L, W, J, d)
    FT = CurrentFT(k0, kl, kw, L, W, J);
    CurrFT(1,:,:) = squeeze(FT(1,:,:)).*2.*cos(kw.*d./2);
    CurrFT(2,:,:) = squeeze(FT(2,:,:)).*2.*cos(kw.*d./2);
    CurrFT(3,:,:) = squeeze(FT(3,:,:)).*2.*cos(kw.*d./2);
end