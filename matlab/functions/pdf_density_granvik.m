function [x, pdf] = pdf_density_granvik(p)
   % Construct the probability density function for the density of a given NEO. 
   % In input, the vector of the probabilities of coming from each route is given.
   % The vector of the probabilities is ordered as follows:
   % 1 : nu6
   % 2 : 3:1
   % 3 : 5:5
   % 4 : Hungaria
   % 5 : Phocaea
   % 6 : 2:1
   % 7 : JFC

   % Interval of densities considered
   x = 400:5:7500;
   pdf = zeros(1, length(x));
   rhoC = 1200;
   rhoS = 2720;
   rhoX = 2350;
   sigmaC = 300;
   sigmaS = 540;
   sigmaX = 520;
   muC = log(rhoC);
   xxC = (1+sqrt(1+4*sigmaC^2/exp(2*muC)))/2;
   ssC = sqrt(log(xxC));
   muS = log(rhoS);
   xxS = (1+sqrt(1+4*sigmaS^2/exp(2*muS)))/2;
   ssS = sqrt(log(xxS));
   muX = log(rhoX);
   xxX = (1+sqrt(1+4*sigmaX^2/exp(2*muX)))/2;
   ssX = sqrt(log(xxX));
   pdfC = lognpdf(x,muC,ssC);
   pdfS = lognpdf(x,muS,ssS);
   pdfX = lognpdf(x,muX,ssX);
   for i=1:7
      switch i 
         % USE THE ALBEDO DISTRIBUTION BY Morbidelli et. al. 2020
         case 1 % nu6, inner belt 
            fC = 0.12;
            fS = 0.558;
            fX = 0.322;
         case 2 % 3:1, inner belt
            fC = 0.144;
            fS = 0.782;
            fX = 0.074;
         case 3 % 5:2, mid belt
            fC = 0.294;
            fS = 0.557;
            fX = 0.149;
         case 4 % Hungaria
            fC = 0.021;
            fS = 0.113;
            fX = 0.866;
         case 5 % Phocaea, inner belt
            fC = 0.501;
            fS = 0.452;
            fX = 0.047;
         case 6 % 2:1, outer belt
            fC = 0.399;
            fS = 0.2;
            fX = 0.461;
         case 7 % Jupiter Family Comets, outer belt
            fC = 1;
            fS = 0;
            fX = 0;
      end
      ff = p(i);
      pdfaux  = fC*pdfC + fS*pdfS + fX*pdfX;
      pdf = pdf + pdfaux*ff;
   end
   % Normalize with the integral
   pdf = pdf/trapz(x, pdf);
end
