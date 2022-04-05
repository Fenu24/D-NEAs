function res = isUnique(Namelist, name)
   index = searchName(Namelist, name);
   if(length(index)==1)
      res = true;   
   else
      res = false;
   end
end
