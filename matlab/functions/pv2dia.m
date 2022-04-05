function d = pv2dia(pv, H)
   % Convert the albedo in diameter, for a given value of the absolute magnitude H

   d = 1000*1329*10.^(-H/5)./sqrt(pv);
end
