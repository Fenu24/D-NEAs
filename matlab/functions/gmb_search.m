function [p, s] = gmb_search(a, e, i, H, model)
   [n,m] = size(model);
   aeiH = model(:, 1:4);
   aux  = [a*ones(n,1), e*ones(n,1), i*ones(n,1), H*ones(n,1)]; 
   diff = aeiH-aux; 
   V = vecnorm(diff');
   [val, index] = min(V);
   % Here the order is 
   % nu6     5/2     2/1     HUN     3/1     PHO     JFC 
   prob = model(index,7:13);
   sigma = model(index,14:20);
   % Change the order to
   % nu6 3:1 5:2 Hungaria Phocaea 2:1 JFC
   p(1) = prob(1);
   p(2) = prob(5);
   p(3) = prob(2);
   p(4) = prob(4);
   p(5) = prob(6);
   p(6) = prob(3);
   p(7) = prob(7);
   
   s(1) = sigma(1);
   s(2) = sigma(5);
   s(3) = sigma(2); 
   s(4) = sigma(4); 
   s(5) = sigma(6); 
   s(6) = sigma(3); 
   s(7) = sigma(7); 

end
