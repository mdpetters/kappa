function [Sm, sigm] = sc_sft(s, k2, Dd, counter)
%+
% NAME:
%	sc_sft.pro
% PURPOSE:
%       computes the critical supersaturation of a dry particle that
%       is composed of a solute and a surfactant. It implements the
%       equations from Petters and Kreidenweis (2013) which is base
%       on the analytical solution presented by Raatikainen and
%       Laaksonen (2012) 
%
%       Petters, M. D. and Kreidenweis, S. M.: A single parameter
%       representation of hygroscopic growth and cloud condensation
%       nucleus activity Part 3: Including surfactant partitioning,
%       Atmos. Chem. Phys., 13, 22687-22712, 2012. 
%
%       Raatikainen, T. and Laaksonen, A.: A simplified treatment of
%       surfactant effects on cloud drop activation, Geosci. Model
%       Dev., 4, 107-116, doi:10.5194/gmd-4-107-2011
%	
% CALLING SEQUENCE:
%       sc = sc_sfc(s, k2, Dd)   
%
% INPUT:
%       s = structure with surfactant properties, 
%       s = {R:R, T:T, alpha:alpha, nu:nu, Gmax:Gmax, beta:beta,
%            k:ksft, sigma0:sigma0, cmc:cmc, A:A, e:1.0d}
%         R: ideal gas constant
%         T: temperature
%         alpha: molar volume of surfactant
%         nu: number of dissociable ions (surfactant)
%         Gmax: max surface excess (surfactant)
%         beta: inverse activity coefficient (surfactant)
%         ksft: chemical kappa of surfactant 
%         sigma0: surface tension of pure water at T
%         cmc: critical micelle concentration
%         A = 8.69251d-6  : A - parameter [K3 m3 J-1]
%         e: volume fraction of surfactant
%         Dd = Particle dry diameter
%         k2 = kappa of second solute
%
%OUTPUT:
%	volume of surfactant in the bulk
%
%DEPENDENCIES: 
%      ccion.m
%      cuberoot.m 
%
%NOTES:
%      Keyword /COUNTER
%      the common counter ion properties are hard coded in this
%      routine. For a system different than SDS/NaCl this must be
%      modified. This note can be disregarded when modeling systems
%      with no common counter ions
%
%EXAMPLE:
%	see example.pro
%
%REVISION HISTORY:
%	Markus Petters, 2012
%       Translated to MATLAB from IDL, 2015
%-  
    Dmax = 20.0 * Dd;
    Sm = 1.0;  sigm = 0.0;
    D = Dd * 1.01;
    while D < Dmax
        D = D * 1.01;                  % do this diameter
        A = pi * D^2.0;                % surface2area of the droplet
        V = 1.0 / 6.0 * pi * D^3.0;    % volume of droplet
        Vts = 1.0 / 6.0 * pi * Dd^3.0; % volume of dry particle
        
        %% compute amount in bulk fraction
        try foo = counter;
            a_NaCl = (58.43 / 2.165) * 1e-6;         % common ion for NaCl
            npl = (1.0 - s.e) * Vts / a_NaCl; 
            Vbs = ccion(npl, 0, 1, 1, A, V, Vts, s); % change for others common 
        catch err
            g = s.e * Vts / s.alpha - s.beta * V - A * s.Gmax / s.nu;
            Vbs = s.alpha / (2.0) * (g + sqrt(g^2.0 + 4.0 * s.e * ...
                                              Vts * s.beta * V / ...
                                              s.alpha)); 
        end
                    
        %% compute surface tension
        sigma = s.sigma0 - s.R*s.T*s.Gmax*log(1.0 + Vbs/(s.alpha*V*s.beta)); 
        if sigma < s.cmc; sigma = s.cmc; end 
            
        %% compute mixed kappa and partitioning fraction
        xi = Vbs / (s.e * Vts);
        k = s.e * xi * s.k + (1 - s.e)*k2;
        aw = (D^3.0 - Dd^3.0)/(D^3.0 - Dd^3.0 * (1-k));
        Sx = aw * exp(s.A * sigma / (s.T * D));
        if Sx > Sm; Sm = Sx; end 
        if sigma > sigm; sigm = sigma; end
    end
end
