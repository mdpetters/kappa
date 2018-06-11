; +
;
; PROGRAM: function findKappa
;
;
; PURPOSE: finds the kappa from a dry diameter and critical
;          supersaturation data pair. 
;
;
; AUTHOR: Markus Petters (petters@atmos.colostate.edu)
;         Department of Atmospheric Science
;         Colorado State University, Fort Collins, CO
;
;
; COMMENTS: This code is based on the work in Petters and Kreidenweis,
;           ACP, 2007. Default values for T are 298.15K and sigma =
;           0.072 J m-2, This code requires the routines "defined.pro"
;           and "init" in the same subdirectory
;
;           Note that this code also works for an array of D and s values 
;-


function findKappa, D, s, T = T, sv = sv
if not defined(sv) then sv = 0.072d
if not defined(T) then T = 298.15d

k = fltarr((n=n_elements(D))) 

for i = 0, n-1 do begin
    k0 = 0.2 & repeat begin
        Sk0 = S[i]-findSc(k0+0.01, D[i], T = T, sigma = sv)
        Sk1 = S[i]-findSc(k0, D[i], T = T, sigma = sv)
        Sk2 = S[i]-findSc(k0-0.01, D[i], T = T, sigma = sv)
        dSdk = (Sk2-Sk0)/0.02d
        k0 = k0 + Sk1/dSdk
    endrep until Sk1^2 le 1d-20
    k[i] = k0
endfor

return, k
end


