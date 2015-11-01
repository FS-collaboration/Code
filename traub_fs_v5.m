%NOTES.
%  From [Cunningham PNAS 2004 SI].
%  x.  No persistent Na
%  x.  No h-current
%  x.  Updated fast Na activation (done)
%  x.  Updated fast Na inactivation (done)
%  x.  Updated KDR (done)

function [V,t,mNaF,hNaF,mKDR,mCaH,kV,mKM,ic] = traub_fs_v5(T, I0, gLNa, gLK, gL, gNaF, gKDR, gCaH, gKM, gKv3, EK, C, sigma, ic)

  dt = 0.005; %time step
 
  V = zeros(T,1);
  mNaF = zeros(T,1);
  hNaF = zeros(T,1);
  mKDR = zeros(T,1);
  mCaH = zeros(T,1);
  mKM = zeros(T,1);
  kV = zeros(T,1);

%intial conditions
  if isstruct(ic)
      V(1) = ic.V;
      mNaF(1) = ic.mNaF;
      hNaF(1) = ic.hNaF;
      mKDR(1) = ic.KDR;
      mCaH(1) = ic.mCaH;
      mKM(1)  = ic.mKM;
      kV(1)   = ic.kV;
  end
      
  for i=1:T-1
      
      noise   = sqrt(dt).*randn().*sigma;
      
      V(i+1) = V(i) + dt.*C*( -I0(i)  ...
          - gL*(70.0 + V(i)) + ...
          - gNaF.*mNaF(i).^3.0.*hNaF(i).*(-50.0 + V(i)) ...    ;NaF current  Rev 50.
          - (  gKDR.*mKDR(i).^4.0 ...                          ;KDR current  Rev -100.
             + gKM.*  mKM(i) ...
             + gKv3*   kV(i)).*(-EK(i) + V(i)) ...
          - gCaH.*mCaH(i).^2.0.*(-125 + V(i))) ...
          + noise;
      
      mNaF(i+1) = mNaF(i) + dt*(alpha_mNaF(V(i))*(1-mNaF(i))  - beta_mNaF(V(i))*mNaF(i));
      hNaF(i+1) = hNaF(i) + dt*(alpha_hNaF(V(i))*(1-hNaF(i))  - beta_hNaF(V(i))*hNaF(i));
      mKDR(i+1) = mKDR(i) + dt*(alpha_mKDR(V(i))*(1-mKDR(i))  - beta_mKDR(V(i))*mKDR(i));
      mCaH(i+1) = mCaH(i) + dt*((1.6.*(1 - mCaH(i)))./(1 + exp(-0.072.*(-5 + V(i)))) - (0.02.*mCaH(i).*(8.9 + V(i)))./(-1 + exp((8.9 + V(i))/5.)));
      kV(i+1)   = kV(i)   + dt*(alpha_KV(  V(i))*(1-  kV(i))  - beta_KV(  V(i))*kV(i));
      mKM(i+1)  = mKM(i)  + dt*(alpha_mKM( V(i))*(1- mKM(i))  - beta_mKM( V(i))*mKM(i));

  end
  

  ic = {};
  t = dt*(1:T);
  ic.V = V(end);
  ic.mNaF = mNaF(end);
  ic.hNaF = hNaF(end);
  ic.KDR = mKDR(end);
  ic.mCaH = mCaH(end);
  ic.mKM = mKM(end);
  ic.kV = kV(end);

end
    
function res = alpha_mNaF(V)    % NaF activation [Cunningham SI 2004].
  minf = 1.0 ./ (1.0 + exp( (-V-38.0)/10 ));
  if V < -30
      taum = 0.0125 + 0.1525 * exp((V+30.0)/10);
  else
      taum = 0.02   + 0.145  * exp((V-30.0)/10);
  end
  res = minf./taum;
end

function res = beta_mNaF(V)    % NaF activation  [Cunningham SI 2004].
  minf = 1.0 ./ (1.0 + exp( (-V-38.0)/10 ));
  if V < -30
      taum = 0.0125 + 0.1525 * exp((V+30.0)/10);
  else
      taum = 0.02   + 0.145  * exp((V-30.0)/10);
  end
  res = (1.0 - minf)./taum;
end

function res = alpha_hNaF(V)    % NaF inactivation [Cunningham SI 2004].
  hinf = 1.0 ./ (1.0 + exp( (V+58.3)/6.7 ));
  tauh = 0.225 + 1.125 / (1.0+exp((V+37.0)/15.0));
  res = hinf./tauh;
end

function res = beta_hNaF(V)    % NaF inctivation  [Cunningham SI 2004].
  hinf = 1.0 ./ (1.0 + exp( (V+58.3)/6.7 ));
  tauh = 0.225 + 1.125 / (1.0+exp((V+37.0)/15.0));
  res = (1.0 - hinf)./tauh;
end

function res = alpha_mKDR(V)  % KDR activation [Cunningham SI 2004].
  minf = 1.0 ./ (1.0 + exp( (-V-27.0)/11.5 ));
  if V < -10
      taum = 0.25 + 4.35 * exp((V+10)/10);
  else
      taum = 0.25 + 4.35 * exp((-V-10)/10);
  end
  res = minf./taum;
end

function res = beta_mKDR(V)  % KDR activation [Cunningham SI 2004].
  minf = 1.0 ./ (1.0 + exp( (-V-27.0)/11.5 ));
  if V < -10
      taum = 0.25 + 4.35 * exp((V+10)/10);
  else
      taum = 0.25 + 4.35 * exp((-V-10)/10);
  end
  res = (1.0 - minf)./taum;
end

function aK = alpha_KV(V)
a = 0.0189324;
b = -4.18371;
c = 6.42606;
aK = a*(-((V)+b)) / (exp(-((V)+b)/c)-1);
end

function bK = beta_KV(V)
d = 0.015857;
e = 25.4834;
bK = d*exp(-(V)/e);
end

function alpha = alpha_mKM(V)  % 	  ;M-current forward rate function [Traub, 2003].
  alpha = 0.02./(1.0 + exp((-20 - V)/5.));
end

function beta = beta_mKM(V)  % 	  ;M-current backward rate function [Traub, 2003].
  beta = 0.01.*exp((-43 - V)/18.);
end