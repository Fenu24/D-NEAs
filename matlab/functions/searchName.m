function index = searchName(Namelist, name)
   [n,~] = size(Namelist);

   index = [];
   for i=1:n
      if(strcmp(Namelist{i},name))
         index = [index i];
      end
   end   

end
