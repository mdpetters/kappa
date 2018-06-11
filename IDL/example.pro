init

; define the constants and setup the problem 
R = 8.314d                   ; universal gas constant
T = 298.15d                  ; temperature
alpha = (288.33/1.176)*1d-6  ; molar volume of compound A (SDS)
Gmax = 13.9*1d-3/(R*T)       ; Szyskowski parameter SDS
beta = 9.5e-1                ; Szyskowski parameter SDS
nu = 2.0d                    ; dissociation constant (SDS)  
cmc = 0.03d                  ; critical micelle concentration
ksft = 0.134                 ; kappa chem SDS from Ruehl et al. 
sigma0 = 0.072d              ; surface tension at standard state [J m-2]
A = 8.69251d-6               ; A - parameter [K3 m3 J-1]
e = 0.9d                     ; volume fraction of surfactant

;; surfactant properties in a structure
s = {R:R, T:T, alpha:alpha, nu:nu, Gmax:Gmax, beta:beta, k:ksft, $
     sigma0:sigma0, cmc:cmc, A:A, e:e}
Dd = 40d-9            ; dry diameter
kNaCl = 1.28d         ; kappa NaCl (here ideal for comparison)

result1 = sc_sft(s, kNaCl, Dd)
result2 = sc_sft(s, kNaCl, Dd, /counter)

print, 'Without common counter ion'
print, 'sc (%)    sigma (J m-2)'
print, result1.sc, result1.tsig

print, ''
print, 'With common counter ion'
print, 'sc (%)    sigma (J m-2)'
print, result2.sc, result2.tsig
end
