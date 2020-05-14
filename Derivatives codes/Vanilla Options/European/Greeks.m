function [] = Greeks(CallPut, AssetP, Strike, RiskFree, DividendY, Time, Volatility)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the Black-Scholes-Merton European Call/Put Greeks based
% on the following inputs:
% CallPut           =       For a call, the input is "1". Puts, use "0".
% AssetP            =       Underlying Asset Price
% Strike            =       Strike Price of Option
% RiskFree          =       Risk Free rate of interest
% DividendY         =       Dividend Yield of Underlying
% Time              =       Time to Maturity
% Volatility        =       Volatility of the Underlying
% Please note that the use of this code is not restricted in anyway.
% However, referencing to the global derivatives website would be appreciated.
% To run this program, simply use the function defined in the 1st line.
% http://www.global-derivatives.com 
% info@global-derivatives.com
% Kevin Cheng (Nov 2003)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = Volatility * sqrt(Time);                                
df = RiskFree -DividendY + 0.5 * Volatility ^ 2;                        % Computes the drift term
d1 = (log( AssetP / Strike ) + df * Time ) / dt;             % Calculates the d1 term used in Black-Scholes
d2 = d1 - dt;                                                % Calculates the d2 term used in Black-Scholes

% The cumulative normal distribution functions for use in computing calls
nd1 = normcdf(d1);
nd2 = normcdf(d2);
% The cumulative normal distribution functions for use in computing puts
nnd1 = normcdf(-d1);
nnd2 = normcdf(-d2);
nn1 = (1 / sqrt(2 * pi) * exp(-0.5 * d1 ^ 2));

% Computes call greeks
if CallPut
    Delta = nd1 * exp(-DividendY * Time)
    Gamma =(nn1 * exp(-DividendY * Time)) / (AssetP * dt)
    Theta = -((AssetP * nn1 * exp(-DividendY * Time) * Volatility) / (2 * sqrt(Time))) + (-DividendY * AssetP * nd1 * exp(-DividendY * Time)) - (RiskFree * Strike * exp(-RiskFree * Time) * nd2)
    Vega = AssetP * sqrt(Time) * nn1 * exp(-DividendY * Time)
    Rho1 = Strike * Time * exp(-RiskFree * Time) * nd2
    Rho2 = -AssetP * exp(-DividendY * Time) * Time * nd1
end
% Computes put greeks
if ~CallPut
    Delta = nd1 * exp(-DividendY * Time) - 1
    Gamma =(nn1 * exp(-DividendY * Time)) / (AssetP * dt)
    Theta = -((AssetP * nn1 * Volatility * exp(-DividendY * Time)) / (2 * sqrt(Time))) - (DividendY * AssetP * nnd1 * exp(-DividendY * Time)) + (RiskFree * Strike * exp(-RiskFree * Time) * nnd2)
    Vega = AssetP * sqrt(Time) * nn1 * exp(-DividendY * Time)
    Rho1 = -Strike * Time * exp(-RiskFree * Time) * nnd2
    Rho2 = AssetP * exp(-DividendY * Time) * Time * nnd1
end