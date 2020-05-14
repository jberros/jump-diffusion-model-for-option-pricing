function Call=BS(S0,t,K,T,Rgrow,Rdisc,sigma)

F=S0.*exp(Rgrow.*T);
d1=log(F./K)./(sigma.*sqrt(T-t))+sigma.*sqrt(T)/2;
d2=log(F./K)./(sigma.*sqrt(T-t))-sigma.*sqrt(T)/2;
Call = exp(-Rdisc.*T).*(F.*normcdf(d1)-K.*normcdf(d2));
% Put=Call+K.*exp(-Rdisc.*T)-S0;

