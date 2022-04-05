function [x, pdf] = pdf_gamma()
   % Probability density function of the obliquity gamma, as given in Tardioli et al. 2017

   gamma = linspace(0, 180, 1000);
   a =  1.12; 
   b = -0.32; 
   c =  0.13;
   x = gamma;
   pdf = (a*cos(gamma*pi/180).^2+b*cos(gamma*pi/180)+c).*(sin(gamma*pi/180))*pi/180;
end
