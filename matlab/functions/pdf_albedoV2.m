function pd = pdf_albedoV2(pv, proute)
   % Albedo density distribution by Morbidelli etal 2020
   % pv: input albedo
   % p : probability of coming from each route 
   % Invexes for the probability vector are
   % 1 : nu6
   % 2 : 3:1
   % 3 : 5:5
   % 4 : Hungaria
   % 5 : Phocaea
   % 6 : 2:1
   % 7 : JFC

   % Define the albedo probabilities for each route. The first index
   % define the albedo category, while the second one the route
   % 1 : nu6
   % 2 : 3:1
   % 3 : 5:5
   % 4 : Hungaria
   % 5 : Phocaea
   % 6 : 2:1
   % 7 : JFC
   
   p_albedo = [0.120, 0.144, 0.294, 0.021, 0.501, 0.399, 1.0;
               0.558, 0.782, 0.557, 0.113, 0.452, 0.200, 0;
               0.322, 0.074, 0.149, 0.866, 0.047, 0.461, 0];
           
           
    integral = 0.104526;
    fact = (1./2.6).^((pv-0.3)./0.1)./integral;
    pd1 = p_albedo(1,:)*proute'/0.1;
    pd2 = p_albedo(2,:)*proute'/0.2;
    pd3 = p_albedo(3,:)*proute';
    
    pd = zeros(1,length(pv));
    
    id1 = (pv<=0.1);
    id2 = (pv>0.1 & pv<=0.3);
    id3 = (pv>3);
    
    pd(id1) = pd1;
    pd(id2) = pd2;
    pd(id3) = fact(id3)*pd3;
end