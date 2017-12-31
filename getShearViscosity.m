% this function is used for calculate shear viscosity according to
% a rheological model proposed by Morris & Boulay(1999)
% By Rui Luo 2017/12/27

function [shearViscosityCartesian,shearViscosityPolar] = getShearViscosity(concentrationCartesian,concentrationPolar)
    Ks = 0.1;
    phiMax = 0.585;
    concentration = concentrationCartesian(:,3);
    shearViscosity = 1+2.5*concentration./(1-concentration/phiMax)...
        +Ks*(concentration/phiMax).^2./(1-concentration/phiMax).^2;
    shearViscosityCartesian = concentrationCartesian;
    shearViscosityCartesian(:,3) = shearViscosity;
    shearViscosityPolar = concentrationPolar;
    shearViscosityPolar(:,3) = shearViscosity;
end