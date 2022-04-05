function p = p2R(pv)
   % Probability density function for albedo of NEOs
   % from Wright et. al. 2016
   fd = 0.253;
   d  = 0.03;
   b  = 0.168; 
   p = fd*pv.*exp(-pv.^2./(2*d^2))./(d^2) +(1-fd).*pv.*exp(-pv.^2./(2*b^2))./(b^2);
end
