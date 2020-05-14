function [JDCallPrice,std_err]=JDprice(S0, X, r, T, vol, a, b, lambda)
%   JDprice Log-Uniform Jump-Diffusion price.
%   Compute European call option price using a Log-Uniform Jump-Diffusion
%   model.
%
%       [JDCallPrice,std_err]=JDprice(S0, X, r, T, vol, a, b,lambda)
%
%*************************************************************************
% ACKNOWLEDGMENTS:
% Thanks to Zongwu Zhu and Floyd B. Hanson for their paper
% "A Monte-Carlo Option-Pricing Algorithm for Log-Uniform".
% This function uses the Monte-Carlo method with antithetic
% and control variates techniques (abbreviated as AOCV) described in
% this paper. For more detail about the algorithm please refer to  
% Zongwu Zhu and Floyd B. Hanson paper.
%
% This function requires the blsprice function of the financial toolbox
%
%*************************************************************************  
%
%   Inputs:
%   S0   	- Current price of the underlying asset.
%
%   X       - Strike (i.e., exercise) price of the option.
%
%   r       - Annualized continuously compounded risk-free rate of return
%                 over the life of the option, expressed as a positive decimal
%                 number.
%
%   T       - Time to expiration of the option, expressed in years.
%
%   vol     - Annualized asset price volatility (i.e., annualized standard
%                 deviation of the continuously compounded asset return),
%                 expressed as a positive decimal number.
%  
%
%   a       - Upper bound of the uniform jump-amplitude mark density  
%
%   b       - Lower bound of the uniform jump-amplitude mark density
%
%  In a real market, the ratio a/b will be close to 1 and b+a will be very 
%  small. The parameters a,b and lambda should be estimated by MLE
%
%   lambda   - Annualized jump rate 
%
%   Fixed Parameter: m        - number of Monte-Carlo simulations
%   m=100;  % This function output faces stability problems for m>=500
%
%   Outputs:
%   JDCallPrice        - Price (i.e., value) of a European call option.
%
%   std_err            - Standard deviation of the error due to the
%                          Monte-Carlo simulation 
%                          std_err=std(sample)/sqrt(length(sample))
%
%**************************************************************************
%
%   Example:
%      Consider European stock options with an exercise price of $102 that
%      expire in 3 months. Assume the underlying stock is trading at $98, 
%      and has a volatility of 45% per annum, the the risk-free rate 
%      is 2.5% per annum. Let the jump-amplitude mark density be on (a,b)
%      (a,b)=(-0.030,0.028) and the annualized jump rate lambda be 60;
%      So S0=98;X=102;r=0.025;T=0.25;vol=0.45;a=-0.030;b=0.028; lambda=60;
%      Using this data, and doing m=100 Monte-Carlo simulations we get:
%
%      [JDCallPrice,std_err]=JDpricer(S0, X, r, T, vol
%      a, b, lambda)
%
%       JDCallPrice =
% 
%           7.7286  % in dollars
% 
%       std_err =
% 
%           0.0554
%
%       If we compare to the Black-Scholes price:
%       BlackScholesCallPrice=blsprice(S0, X, r, T, vol)
%       BlackScholes_Price =
%           7.3465
%
% So that JDCallPrice > BlackScholesCallPrice, which is a general result
% whish can be proven by Jensen’s inequality
%
%**************************************************************************
% Rodolphe Sitter - MSFM Student - The University of Chicago
% March 2009
%**************************************************************************

m=100; %Number of Monte-Carlo simulations 
if lambda==0
    JDCallPrice=blsprice(S0,X,r,T,vol);
    std_err=0;
else
alpha=10^(-6);
Jav=(exp(b)-exp(a))/(b-a)-1;
K=poissinv(1-alpha,lambda*T);
p0=poisspdf(0,lambda*T);
scaling=poisscdf(K,lambda*T);
BS_k0=BS(S0*exp(-lambda*Jav*T),0,X,T,r,r,vol); % k=0, no jumps
p=poisspdf(1:K,lambda*T);
U=rand(K,m);
Sk=zeros(K,m);
Sk_a=zeros(K,m);
S0_poiss=zeros(K,m);
S0_poiss_a=zeros(K,m);
BS_k=zeros(K,m);
BS_k_a=zeros(K,m);
for k=1:K;     
    Sk(k,:)=k*a+(b-a)*sum(U(1:k,:));
    Sk_a(k,:)=(a+b)*k-Sk(k,:);
    S0_poiss(k,:)=S0*exp(Sk(k,:)-lambda*Jav*T);
    S0_poiss_a(k,:)=S0*exp(Sk_a(k,:)-lambda*Jav*T);
    BS_k(k,:)=BS(S0_poiss(k,:),0,X,T,r,r,vol);
    BS_k_a(k,:)=BS(S0_poiss_a(k,:),0,X,T,r,r,vol);
end  
sample=(p0*BS_k0+p*BS_k)/scaling;
%     JDCallPrice=mean(sample);    % without using the AOCV method
%     std_err=std(sample)/sqrt(m); % without using the AOCV method
sample_a=(p0*BS_k0+p*BS_k_a)/scaling;
x=.5*(sample+sample_a);
y=.5*(p*exp(Sk)+p*exp(Sk_a));
VARy=.5*(exp(lambda*T*Jav)-2*exp(2*lambda*T*Jav)+exp(lambda*T*(exp(a+b)-1)));
beta=m/(m-1)*(mean(x.*y)-mean(x)*mean(y))/VARy;
Z=x-beta*(y-exp(lambda*T*Jav));
y2=y.*(2*mean(y)-y);
bias=1/(m-1)*(mean(x.*y2)-mean(x)*mean(y2))/VARy;
JDCallPrice=mean(Z)-bias;
std_err=std(Z)/sqrt(m);
end

