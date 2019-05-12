; +
;
; PROGRAM: function numgen
;
;
; PURPOSE: To generate n floating point numbers between [min, max]
;
;
; AUTHOR: Markus Petters (petters@uwyo.edu)
;         Department of Atmospheric Science
;         University of Wyoming
;
;
; COMMENTS: 
; Mar-28-2006: Added the keyword log this calculates the logarithmic
;              geometric series bewteen min and max. THis series is
;              useful for geometrically stepped grids in log space,
;              i.e. dN/dlogD = const. 
;-


function numgen, n, min = min, max = max, log = log

if not defined(log) then begin
    if not keyword_set(min) then min = 0
    if not keyword_set(max) then max = n
    
    
return, dindgen(n)*(max - min)/(n - 1) + min
endif else begin
    x = dblarr(n) & x[0] = min
    step = alog10(max/min)/(n-1)
    for i = 1, n-1 do x[i] = x[i-1]*10^step
    return, x
endelse
end
