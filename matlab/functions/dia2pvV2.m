function [pv, dpvdd] = dia2pvV2(d, H)
   % Convert the diameter to albedo, given the absolute magnitude H
   pv    = (1329*1000*10.^(-H./5)./d).^2;
   dpvdd = -2*(1329*1000*10.^(-H./5)).^2./d.^3;
end