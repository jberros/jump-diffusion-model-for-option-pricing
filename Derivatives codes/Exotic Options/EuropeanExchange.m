================ Matlab File ======================

function price = margrabeEuroExchageOptions(  S1, S2, sigma1, sigma2, D1, D2, T, r, rho) % European exchange options Margrabe model %  price = margrabeEuroExchageOptions(  S1, S2, sigma1, sigma2, D1, D2, T, r, rho) % % Input Parameters:
% ================
%   S1     = Sacrificed asset price
%   sigma1 = Sacrificed asset volatility
%   D1     = Sacrificed asset dividend yield
%   S2     = Exchanged asset price
%   sigma2 = Exchanged asset volatility
%   D2     = Exchanged asset dividend yield
%   r      = Risk free rate
%   rho    = Correlation coefficient between the two assets
%   T      = Time to maturity
%
% Output Parameters:
% =================
%   price = the european exchange price
%
%
% Example:
% =======
% >> price = euroExchageOptions(  130, 105, 0.2, 0.2, 0.06, 0.04, 1, 0.1, 0.5) %  price =  23.50829459339569
% Code adapted by Sione @ sionep@xtra.co.nz
% http://www.global-derivatives.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if(S1<0)
  error('S1 - must be positive.');
elseif(S2<0)
  error('S2 - must be positive.');
elseif(sigma1<=0)
  error('sigma1 must be greater than 0.');
elseif(sigma2<=0)
  error('sigma2 must be greater than 0.');
elseif(T<=0)
  error('T must be a positive number.');
elseif(r<=0)
  error('r must be a positive number.');
elseif(D1<=0)
  error('D1 must be a positive number.');
elseif(D2<=0)
  error('D2 must be a positive number.');
elseif(rho<0)
  error('rho must be at least zero.');
end

sigmaTot = sqrt( sigma1^2 + sigma2^2 - 2*rho*sigma1*sigma2);
d1 = (log(S1/S2) + (D2 - D1 + 0.5*sigmaTot^2)*T)/ (sigmaTot*sqrt(T));
d2 = d1 - sigmaTot*sqrt(T);

price = S1*exp(-D1*T)*normcdf(d1) - S2*exp(-D2*T)*normcdf(d2);

=============== End Matlab File =================

