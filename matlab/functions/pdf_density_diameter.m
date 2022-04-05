function [sample_rho, sample_D] = pdf_density_diameter(H_mean, H_std, proute, nsample)
   sample_rho = zeros(nsample,1);
   sample_D   = zeros(nsample,1);

   % Define the densities of C-, S-, X-complex and compute the mean and the standard
   % deviations of the corresponding lognormal distributions
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


   [pv, pdf_pv] = pdf_albedo_total(proute);

   j = 1;
   while(j<=nsample)
      % Generate a value of the magnitude, using a Gaussian distribution
      H = normrnd(H_mean, H_std);

      % Generate a value of the albedo pV
      pV_route = pdfrnd(pv, pdf_pv, 1);

      % Generate the value of the density
      D = pv2dia(pV_route, H);
      if(D == Inf)
         % If it is infinite, use another pV, H
         continue
      end
      
      % Generate the value of the density
      if( pV_route <= 0.1 )
         sample_rho(j) = lognrnd(muC, ssC);
      elseif( pV_route > 0.1 && pV_route <= 0.3 )
         sample_rho(j) = lognrnd(muS, ssS);
      else
         sample_rho(j) = lognrnd(muX, ssX);
      end

      sample_D(j) = D;
      j = j+1;
   end
   
end
