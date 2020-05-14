function volatility = JDimpv(S0, X, r, T, a, b, lambda, value)
%**************************************************************************
%Acknowledgement. 
%This function was inspired by the blsimpv function of the financial toolbox
%**************************************************************************
%JDimpv Log-Uniform Jump-Diffusion implied volatility.
%   Compute the implied volatility of an underlying asset from the market 
%   value of European call using a Log-Uniform Jump-Diffusion model.
%
%   volatility = JDimpv(S0, X, r, T, a, b, lambda, value)
%
% Inputs: 
%   S0      - Current price of the underlying asset.
%
%   X       - Strike (i.e., exercise) price of the option.
%
%   r       - Annualized continuously compounded risk-free rate of return over
%           the life of the option, expressed as a positive decimal number.
%
%   T       - Time to expiration of the option, expressed in years.
%
%   a       - Upper bound of the uniform jump-amplitude mark density  
%
%   b       - Lower bound of the uniform jump-amplitude mark density
%
%   lambda  - Annualized jump rate 
%
%   value   - Price (i.e., value) of a European option from which the implied
%     volatility of the underlying asset is derived.
%     *** This input is allowed to be a Matrix ***
%
% Output:
%   volatility - Implied volatility of the underlying asset derived from 
%     European option prices, expressed as a decimal number. If no solution
%     can be found, a NaN (i.e., Not-a-Number) is returned.
%
%   Example:
%      Consider a European call option trading at $7 with an exercise price
%      of $102 that expires in 3 months. Assume the underlying stock is 
%      trading at $98, the the risk-free rate is 2.5% per annum, the 
%      jump-amplitude mark density is on (a,b)=(-0.030,0.028) and the 
%      annualized jump rate is lambda be 60;
%      So (S0=98;X=102;r=0.025;T=0.25;a=-0.030;b=0.028;lambda=60; value=7;)
%      
%   volatility = JDimpv(S0, X, r, T, a, b, lambda, value)
% 
%       volatility =
% 
%           0.4112  % i.e 41.12% per annum.
%Note: The result is random from about the third decimal place due to the 
%Monte-Carlo simulation necessary to exectute JDprice            
%
%**************************************************************************
% Rodolphe Sitter - MSFM Student - The University of Chicago
% March 2009
%**************************************************************************

[nRows, nCols] = size(value);
volatility = nan(nRows * nCols, 1);

options = optimset('fzero');
options = optimset(options, 'TolX', 1e-6, 'Display', 'off');

for i = 1:length(volatility)
    if (S0 > 0) && (X > 0) && (T > 0)

    try
          volatility(i) = fzero(@objfcn, [0 9], options, ...
                                 S0, X, r, T, a, b, lambda, value(i));

                catch
          volatility(i) = NaN;
    end
    end
end
volatility = reshape(volatility, nRows, nCols);

% * * Jump Diffusion Implied Volatility Objective Function * * * * *
function delta = objfcn(volatility, S0, X, r, T, a, b,lambda, value)

callValue = JDprice(S0, X, r, T, volatility, a, b, lambda);
delta = value - callValue;


