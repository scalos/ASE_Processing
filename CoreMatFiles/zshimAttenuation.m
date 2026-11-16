function F = zshimAttenuation(Gz, tDeph, sliceThk)
% zshimAttenuation  Through-slice dephasing attenuation for a rect slice profile.
%   Gz       : array of gradients (uT/m)
%   tDeph    : dephasing time(s) (s); scalar or broadcastable against Gz
%   sliceThk : slice thickness (m)
%   F = |sinc(gamma*dz*Gz*t/2)| with the UNNORMALIZED sinc sin(x)/x (not MATLAB's sinc).
    gamma = 2.67502e8;                         % rad/s/T
    x = gamma*sliceThk*(Gz*1e-6).*abs(tDeph)/2;
    F = abs(sin(x)./x);
    F(x == 0) = 1;
    F(isnan(Gz)) = NaN;
end
