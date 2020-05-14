=========== Matlab Code ============================

function price = fltStrikeLookback(callPut, S, S_max_min, T, D, r, sigma) % Floating Strike Lookback Options % price = fltStrikeLookback(callPut, S, S_max_min, T, D, r, sigma) % % Input Parameters:
% ================
%   callPut   = 1 (call option) or 0 (put option)
%   S         = current asset price 
%   S_max_min = if callPut is 1, then S_max_min is the minimum asset price
%               else S_max_min is the maximum asset price
%   r         = risk free rate
%   D         = dividend yield
%   sigma     = volatility
%   T         = time to maturity
%
% Output Parameter:
% =================
%   price = call option price if callPut is 1 or put option price if callPut
%           is 0
%
%
% Example from Matlab Prompt:
% ==========================
%  >> call = fltStrikeLookback(1, 100, 98, 0.501, 0.03, 0.1, 0.3) %
%     call = 17.18954547977881
%
%  >> put = fltStrikeLookback(0, 100, 98, 0.501, 0.03, 0.1, 0.3) %
%     put = 14.04526346904223
% Code adapted by Sione @ sionep@xtra.co.nz
% http://www.global-derivatives.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


if(callPut~=0)
  if(callPut~=1)
   error('callPut - must be a logical value of either 1 (call option) or 
0 (put option).');
  end
elseif(S_max_min<0)
  error('S_max_min - must be positive.');
elseif(S<0)
  error('S - must be positive.'); 
elseif(T<=0)
  error('T must be a positive number.');
elseif(r<=0)
  error('r must be a positive number.');
elseif(D<=0)
  error('D must be a positive number.');
elseif(sigma<=0)
  error('sigma must be a positive number.');
elseif(r==D)
  error('D must not be equal to  r.');
end

price      = NaN;
term1      = NaN;
term2      = NaN;
brackCoeff = S*exp(-r*T)*(sigma*sigma)/(2*(r - D));
brack      = NaN;
brack2     = NaN;
temp       = NaN;
dt         = sigma*sqrt(T);


if(callPut)
   a1 = ( log(S/S_max_min) + (r - D + 0.5*sigma*sigma)*T ) / dt;
   a2 = a1 - dt;
   term1 = S*exp(-D*T)*normcdf(a1);
   term2 = S_max_min*exp(-r*T)*normcdf(a2);  
   temp = -a1 + 2*(r - D)*sqrt(T)/sigma;
   brack = (S/S_max_min)^(-2*(r - D)/(sigma*sigma))*normcdf(temp);
   brack2 = exp((r-D)*T)*normcdf(-a1);
   price = term1 - term2 + brackCoeff*(brack - brack2);
else
   b1 = (log(S_max_min/S) + (D + r + 0.5*sigma*sigma)*T)/dt;
   b2 = b1 - dt;
   term1 = S_max_min*exp(-r*T)*normcdf(-b2);
   term2 = S*exp(-D*T)*normcdf(-b1);  
   temp = b1 - 2*(r - D)*sqrt(T)/sigma;
   brack = (S/S_max_min)^(2*(r - D)/(sigma*sigma))*normcdf(temp);
   brack2 = exp((r-D)*T)*normcdf(b1);
   price = term1 - term2 + brackCoeff*(-brack + brack2);
end

=========== End Matlab Code ===============