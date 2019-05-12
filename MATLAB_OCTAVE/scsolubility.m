function [smax] = scsolubility(Dd, Ci, ei, ks, T, sigma)
% +
% PURPOSE: finds the critical supersaturation from a dry diameter and
%          a set of kappa values, volume fractions, and solubilities. 
%
%
% AUTHOR: Markus Petters 
%         Department of Marine Earth and Atmospheric Sciences
%         NC State University (mdpetter@ncsu.edu)
%
%
% COMMENTS: This code is published in Petters and Kreidenweis,
%           ACPD, 2008. 
%
%-

A = 8.69251d-6*sigma/T;
g = 1.01;
Smax = 0;
while g < 20
    xi = Ci.*(g.^3.0 - 1.0)./ei;
    xi(find(xi > 1)) = 1;
    k = sum(ks .* (ei.*xi));
    xw = ((Dd.*g).^3.0 - Dd.^3.0)./((Dd.*g).^3.0 - Dd.^3.0.*(1.0-k));
    S = xw.*exp(A./(Dd.*g));
    if S > Smax; Smax = S; end;
    g = g*1.0005;
end

smax = (Smax-1)*100;
end
