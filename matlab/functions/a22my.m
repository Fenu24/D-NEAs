function yarko = a22my(A2, a, e)
   % Convert the A2 parameter from the JPL website into an acceleration
   % in AU/My
   
   % Sun gravitational constant
   GM   = 1.32712440018d20;
   % Conversion factors
   au2m = 149597870700.d0;
   s2d = 1/86400;
   % Compute the mean motion
   nstd = sqrt(GM/(a*au2m)^3)/s2d;
   % Convert the value of A2:
   p = a*(1-e^2);
   d = 2;
   d2y  = 1/365.25;
   y2my = 1/10^6; 
   d2my = d2y*y2my;
   yarko = 2*A2*(1-e^2)/p^d/nstd/d2my;
   %2 A2 (1-e^2)*(1/(a(1-e^2)))^2/n
end
