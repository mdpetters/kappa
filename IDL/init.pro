;-
;
; PROGRAM: init.pro
;
;
; PURPOSE: create my system variables
;
;
; AUTHOR: Markus Petters (petters@uwyo.edu)
;         Department of Atmospheric Science
;         University of Wyoming
;
; COMMENTS:
;
;-  

pro init
common initializiation_block, system_variables_are_initialized

if not defined(system_variables_are_initialized) then begin
    defsysv, '!true', 1
    defsysv, '!false', 0
    defsysv, '!fptr', 1
    defsysv, '!fptr1', 2
    defsysv, '!g', 9.81
    defsysv, '!cpd', 1005.2d
    defsysv, '!cpv', 1900.0d
    defsysv, '!cpi', 2120.0d
    defsysv, '!cpw', 4218.0d
    defsysv, '!R', 8.314d
    defsysv, '!Rd', 287.05d
    defsysv, '!Rv', 461.15d
    defsysv, '!ep', 0.622d
    defsysv, '!ew0', 610.7d
    defsysv, '!Md', 28.964*1d-3
    defsysv, '!Mv', 18.015*1d-3
    defsysv, '!lv0', 2.501d6
    defsysv, '!lf0', 0.334d6
    defsysv, '!ls0', 2.834d6
    defsysv, '!rhow', 997.1d
    defsysv, '!NA', 6.023d23
    defsysv, '!T0', 273.15d
    defsysv, '!P0', 101325.0d
    defsysv, '!k', 1.380658d-23           
    defsysv, '!e', 1.60217733d-19         
endif    

system_variables_are_initialized = 1
end
