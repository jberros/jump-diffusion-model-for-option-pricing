function [] = BinomialEuro(CallPut, AssetP, Strike, RiskFree, Div, Time, Vol, nSteps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the Cox, Ross & Rubinstein (1979) Binomial Tree for European Call/Put Option Values based
% on the following inputs:
% CallPut           =       Call = 1, Put = 0
% AssetP            =       Underlying Asset Price
% Strike            =       Strike Price of Option
% RiskFree          =       Risk Free rate of interest
% Div               =       Dividend Yield of Underlying
% Time              =       Time to Maturity
% Vol               =       Volatility of the Underlying
% nSteps            =       Number of Time Steps for Binomial Tree to take
% Please note that the use of this code is not restricted in anyway.
% However, referencing the author of the code would be appreciated.
% To run this program, simply use the function defined in the 1st line.
% http://www.global-derivatives.com 
% info@global-derivatives.com
% Kevin Cheng (Nov 2003)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = Time / nSteps;

if CallPut
    b = 1;
end
if ~CallPut
    b = -1;
end

RR = exp(RiskFree * dt);
Up = exp(Vol * sqrt(dt));
Down = 1 / Up;
P_up = (exp((RiskFree - Div) * dt) - Down) / (Up - Down);
P_down = 1 - P_up;
Df = exp(-RiskFree * dt);

for i = 0:nSteps
    state = i + 1;
    St = AssetP * Up ^ i * Down ^ (nSteps - i);
    Value(state) = max(0, b * (St - Strike));
end

for tt = nSteps - 1 : -1 : 0
    for i = 0:tt
        state = i + 1;
        Value(state) = (P_up * Value(state + 1) + P_down * Value(state)) * Df;
    end
end

Binomial = Value(1)