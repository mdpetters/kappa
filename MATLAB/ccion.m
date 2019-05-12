function [Vbs] = ccion(npl, nmi, nup, num, A, V, Vts, s)
%+
% NAME:
%	ccion
% PURPOSE:
%       computes the volume fraction of surfactant in the bulk of a
%       droplet. Equation is based on Eq. (8) in Petters and
%       Kreidenweis (2013) which built from Eq. 11 in Raatikainen and
%       Laaksonen (2012)
%
%       Petters, M. D. and Kreidenweis, S. M.: A single parameter
%       representation of hygroscopic growth and cloud condensation
%       nucleus activity Part 3: Including surfactant partitioning,
%       Atmos. Chem. Phys., 12, 22687-22712, 2012. 
%
%       Raatikainen, T. and Laaksonen, A.: A simplified treatment of
%       surfactant effects on cloud drop activation, Geosci. Model
%       Dev., 4, 107-116, doi:10.5194/gmd-4-107-2011
%	
% CALLING SEQUENCE:
%       Vbs = ccion(npl, nmi, nup, num, A, V, Vts, s)  
%
% INPUT:
%	See publications for def of npl, nmi, nup, num
%       A = surface area of droplet
%       V = volume of droplet
%       Vts = dry volume of droplet
%       s = structure with surfactant properties, defined in
%       sc_sfc.m
%
% OUTPUT:
%	volume of surfactant in the bulk
%
% DEPENDENCIES: 
%       cuberoot.pro 
%
% EXAMPLE:
%	see sc_sfc.pro
%
% REVISION HISTORY:
%	Markus Petters, 2012
%       Translated to MATLAB, 2015
%-
  k1 = nmi / num + npl / nup;              
  k2 = nup / num * nmi + num / nup * npl;  
  
  nts = s.e * Vts / s.alpha;

  a0 = k2 * nts * s.beta * V; 
  a1 = k2 * nts + (s.nu * nts - k2) * s.beta * V - k1 * A * s.Gmax;
  a2 = s.nu * nts - k2 - s.nu * s.beta * V - A * s.Gmax;
  a3 = - s.nu; 
  
  f = cuberoot([a0,a1,a2,a3]);              
  nb = f(f > 0); 
  Vbs = nb(1) * s.alpha;
end
