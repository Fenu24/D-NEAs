function [x, pdf] = pdf_cosgamma()
   % Probability density function for \cos\gamma from Tardioli et. al. 2017

   % Disretize [-1, 1]
   x = linspace(-1,1,2000);
   % Define the weights
   a =  1.12; 
   b = -0.32; 
   c =  0.13;
   pdf = a*x.^2 + b*x + c;
end
