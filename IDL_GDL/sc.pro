; +
;
; PROGRAM: function Sc
;
;
; PURPOSE: finds the critical supersaturation from a dry diameter and
;          a set of kappa values, volume fractions, and solubilities. 
;
;
; AUTHOR: Markus Petters (petters@atmos.colostate.edu)
;         Department of Atmospheric Science
;         Colorado State University, Fort Collins, CO
;
;
; COMMENTS: This code is published in Petters and Kreidenweis,
;           ACPD, 2008. Default values for T are 298.15K and sigma =
;           0.072 J m-2 
;
;
; EXAMPLE: This code can be executed directly from the IDL Python
;          interpreter ;
;
; MODIFICATION HISTORY:
;   
; Copyright (C) 2008, Markus Petters, Department of Atmospheric
; Science, Colorado State Univsersity
;
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.
;
;-


function sc, Dd, g, Ci, ei, ks, T = T, sigma = sigma
  if not keyword_set(T) then T = 298.15d
  if not keyword_set(sigma) then sigma = 0.072d
  A = 8.69251d-6*sigma/T
  xi = Ci*(g^3.0d - 1.0d)/ei
  i = where(xi gt 1) & if i[0] ne -1 then xi[i] = 1
  k = total(ks * (ei*xi))
  xw = ((Dd*g)^3.0 - Dd^3.0)/((Dd*g)^3.0 - Dd^3.0*(1.0-k))
  S = xw*exp(A/(Dd*g))
  return, (f = (g gt 20) ? S : max([S, sc(Dd,g*1.01,Ci,ei,ks,T=T,sigma=sigma)]))
end

print, sc(100d-9, 1+1d-11, [1d15, 0.1], [0.5, 0.5], [0.6, 0.2])
print, sc(100d-9, 1+1d-11, [1d15], [1], [0.6], T=273.15)

end
