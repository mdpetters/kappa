; +
;
; PROGRAM: function findSc
;
;
; PURPOSE: finds the critical supersaturation from a dry diameter and
;          a kappa value. 
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
;-


function findSc, k, Ds, T = T, sigma = sigma
init
if not defined(T) then T = 298.15d
if not defined(sigma) then sigma = 0.072d
A = 4*!Mv*sigma/(!R*T*!rhow)

Dt = Ds & Sold = 0.0d & Snew = 0.1d
while Snew le 0.9 do begin
    Dt *= 1.1 & top = (Dt^3.0 - Ds^3.0) & bot = Dt^3 - Ds^3.0 * (1.0-k)
    aw = top/bot & Sold = Snew & Snew = aw*exp(A/Dt)
endwhile

while Snew gt Sold do begin
    Dt *= 1.0005 & top = (Dt^3.0 - Ds^3.0) & bot = Dt^3 - Ds^3.0 * (1.0-k)
    aw = top/bot & Sold = Snew & Snew = aw*exp(A/Dt)
endwhile

if k eq 0 then Snew = exp(A/Ds)
return,  (Snew-1)*100
end

