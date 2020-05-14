function [] = MertonJumpEuro(CallPut, AssetP, Strike, RiskFree, Time, Vol, Jumps, Gamma, MaxIter)
% Requires BlackScholesEuro.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the Merton (1976) Jump diffusion model for European Call/Put Option Values based
% on the following inputs:
% CallPut           =       Call = 1, Put = 0
% AssetP            =       Underlying Asset Price
% Strike            =       Strike Price of Option
% RiskFree          =       Risk Free rate of interest
% Time              =       Time to Maturity
% Vol               =       Volatility of the Underlying
% Jumps             =       Number of Jumps per Year
% Gamma             =       Percent of total volatility explained by jumps
% MaxIter           =       Max number of iterations to be used
% Please note that the use of this code is not restricted in anyway.
% However, referencing the author of the code would be appreciated.
% To run this program, simply use the function defined in the 1st line.
% http://www.global-derivatives.com 
% info@global-derivatives.com
% Kevin Cheng (Nov 2003)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Delta = sqrt((Gamma * (Vol ^ 2)) / Jumps);
A = sqrt(Vol ^ 2 - Jumps * Delta ^ 2);
Value = 0;

for i = 0:MaxIter
    VV = sqrt(A ^ 2 + Delta ^ 2 * (i / Time));
    Value = Value + (exp(-Jumps * Time) * (Jumps * Time) ^ i / factorial(i)) * BlackScholesEuro(CallPut, AssetP, Strike, RiskFree, Time, VV);
end



