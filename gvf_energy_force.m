function [ Eext, Fext ] = gvf_energy_force( I, Options )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Transform the Image into an External Energy Image
Eext = ExternalForceImage2D(I, Options.Wline, Options.Wedge, Options.Wterm, Options.Sigma1);

% Make the external force (flow) field.
Fx = ImageDerivatives2D(Eext, Options.Sigma2,'x');
Fy = ImageDerivatives2D(Eext, Options.Sigma2,'y');
Fext(:,:,1)= - Fx*2*Options.Sigma2^2;
Fext(:,:,2)= - Fy*2*Options.Sigma2^2;

% Do Gradient vector flow, optimalization
GVF = true;
if GVF
    Fext = GVFOptimizeImageForces2D(Fext, Options.Mu, Options.GIterations, Options.Sigma3);
end

end

