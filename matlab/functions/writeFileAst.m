function writeFileAst(name,a,e)
    folder = strcat('input/', name);
    filename = strcat('../',folder,'/gamma_est_mc.nml');
    line_a  = strcat('        semiaxm =', 32, num2str(a),',\n'); 
    line_e  = strcat('            ecc =', 32, num2str(e),',\n'); 
    line_fn = strcat('      filename  =', 32, '''', name, '.txt'',\n');
    fid = fopen(filename, 'w');
    fprintf(fid, '&asteroid \n');
    fprintf(fid, '           C    = 680.d0,\n');
    fprintf(fid, ' thermalCondMin = 0.0000001d0,\n');
    fprintf(fid, ' thermalCondMax = 100.0d0,\n');
    fprintf(fid, line_a); 
    fprintf(fid, line_e);
    fprintf(fid, '       absCoeff = 1.0d0,\n');
    fprintf(fid, '        emissiv = 0.984d0,\n');
    if(e<0.2)
       fprintf(fid, '         method = 1,\n');
    else
       fprintf(fid, '         method = 2,\n');
    end
    fprintf(fid, line_fn);
    fprintf(fid, '     comparison = 0,\n');
    fprintf(fid, '           expo = 0.0,\n');
    fprintf(fid, '/\n');
    fclose(fid);
end
