function [] = BlackScholesMertonEuro(AssetP, Strike, RiskFree, DividendY, Time, Volatility)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the Black-Scholes-Merton European Call/Put Option Values based
% on the following inputs:
% AssetP            =       Underlying Asset Price
% Strike            =       Strike Price of Option
% RiskFree          =       Risk Free rate of interest
% DividendY         =       Dividend Yield of Underlying
% Time              =       Time to Maturity
% Volatility        =       Volatility of the Underlying
% Please note that the use of this code is not restricted in anyway.
% However, referencing the author of the code would be appreciated.
% To run this program, simply use the function defined in the 1st line.
% http://www.global-derivatives.com 
% info@global-derivatives.com
% Kevin Cheng (Nov 2003)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = Volatility * sqrt(Time);                                
df = RiskFree - DividendY + 0.5 * Volatility ^ 2;            % Computes the drift term
d1 = (log( AssetP / Strike ) + df * Time ) / dt;             % Calculates the d1 term used in Black-Scholes
d2 = d1 - dt;                                                % Calculates the d2 term used in Black-Scholes

% The cumulative normal distribution functions for use in computing calls
nd1 = normcdf(d1);
nd2 = normcdf(d2);
% The cumulative normal distribution functions for use in computing puts
nnd1 = normcdf(-d1);
nnd2 = normcdf(-d2);

% Computes call price
CallPrice = AssetP * exp(-DividendY * Time) * nd1 - Strike * exp(-RiskFree * Time) * nd2
% Computes put price
PutPrice = Strike * exp(-RiskFree * Time) * nnd2 - AssetP * exp(-DividendY * Time) * nnd1
