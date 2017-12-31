% this function is used for calculate shear viscosity according to
% a rheological model proposed by Morris & Boulay(1999)
% By Rui Luo 2017/12/27

function [normalViscosityCartesian,normalViscosityPolar] = getNormalViscosity(concentrationCartesian,concentrationPolar)
    Kn = 0.75;
    phiMax = 0.585;
    concentration = concentrationCartesian(:,3);
    normalViscosity = Kn*(concentration/phiMax).^2./(1-concentration/phiMax).^2;
    normalViscosityCartesian = concentrationCartesian;
    normalViscosityCartesian(:,3) = normalViscosity;
    normalViscosityPolar = concentrationPolar;
    normalViscosityPolar(:,3) = normalViscosity;
end