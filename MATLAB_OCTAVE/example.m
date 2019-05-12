
%% define common constants and setup the problem
R = 8.314;                   % universal gas constant
T = 298.15;                  % temperature
sigma0 = 0.072;              % surface tension at standard state [J m-2]
A = 8.69251e-6;              % A - parameter [K3 m3 J-1]
e = 0.9;                     % volume fraction of surfactant
alphaw = (18.0/0.997)*1e-6;  % molar volume of water

%% Initialize the calculation SDS
alpha = (288.33/1.176)*1e-6; % molar volume of compound A 
Gmax = 13.9*1e-3/(R*T);      % Szyskowski parameter SDS
beta = 9.5e-1;               % Szyskowski parameter 
nu = 2.0;                    % dissociation constant   
cmc = 0.03;                  % critical micelle concentration
ksft = 0.134;                % kappa chem  
sds = struct('R', R, 'T', T, 'alpha', alpha, 'nu', nu, 'cmc', cmc, ...
           'Gmax', Gmax, 'beta', beta, 'k', ksft, 'sigma0', sigma0, ...
           'A', A, 'e', e) ; 

Dd = 40e-9;     % dry diameter
kNaCl = 1.28;   % kappa NaCl (here ideal for comparison)

[Sc, sigm] = sc_sft(sds, kNaCl, Dd);
sc = (Sc-1)*100;
disp('Without common counter ion')
disp('sc (%)    sigma (J m-2)')
fprintf("%.4f   %.4f", sc, sigm)

[Sc, sigm] = sc_sft(sds, kNaCl, Dd, true);
sc = (Sc-1)*100;
fprintf('\n')
disp('With common counter ion')
disp("sc (%)    sigma (J m-2)")
fprintf("%.4f   %.4f", sc, sigm)
fprintf('\n')

eps = 0.01:0.001:1;
for i = 1:numel(eps)
    sds.e = eps(i);
    [Sc, sigm] = sc_sft(sds, kNaCl, Dd, true);
    sc(i) = (Sc-1)*100;
end

% Compare to Raatikainen and Laaksonen, Figure 2
% https://www.geosci-model-dev.net/4/107/2011/gmd-4-107-2011.pdf
plot(eps*100, sc)
xlabel('Surfactant dry mass fraction (%)')
ylabel('Critical supersaturation (%)')
print -dpdf example.pdf
