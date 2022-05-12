function gen_distribV2(D, D_sigma, H, H_sigma, dadt_mean, dadt_sigma,...
    P, P_sigma, p, folder, n)
   % gen_distrib(D, D_sigma, H, H_sigma, dadt_mean, dadt_sigma, P, P_sigma, p, folder, n)
   
   
   % Optional parameter
    if ~exist('n','var')
        n = 1.e4;
    end

   path = strcat('../',folder);
   
   % Sample diameter:
   if ~isnan(D_sigma)
        % If the diameter is measured, use a Gaussian distribution for the diameter, 
        % and the old version of the density distribution
        sample_diam  = normrnd(D, D_sigma, n,1);
        [rho, pdf_rho] = pdf_density_granvik(p);
        sample_rho   = pdfrnd(rho, pdf_rho, n);
   else
        % If the diameter is not known, the samples for diameter and density have to be
        % constructed using the same albedo distribution
        H_sigma = 0.3;
        [sample_rho, sample_diam] = pdf_density_diameter(H, H_sigma, p, 50000);
   end
   
   % Distribution gamma
   [gam, pdf_gam] = pdf_gamma();

   % Generate samples using the inverse transform sampling
   sample_gamma = pdfrnd(gam, pdf_gam, n);
   sample_dadt  = normrnd(dadt_mean, dadt_sigma, n, 1);
   if isnan(P)
       sample_P = generatePeriodDistribution(H, n);
   else
       sample_P = normrnd(P, P_sigma, n, 1);
   end
   
   %{
   figure(1)
   t = tiledlayout(1,2); % Requires R2019b or later
   nexttile
   plot(rho, pdf_rho, 'LineWidth', 2.5)
   xlabel('density, $\rho$ (kg/m$^3$)', 'interpreter', 'latex')
   ylabel('pdf', 'interpreter','latex')
   xlim([800 6000])
   set(gca, 'FontSize', 20) 
   nexttile
   plot(diam, pdf_D, 'LineWidth', 2.5)
   xlabel('diameter, D (m)', 'interpreter', 'latex')
   xlim([0 max(diam)/3])
   %xticks([0 20 40 60 80 100 120 140 160]);
   set(gca, 'FontSize', 20) 
   % Set the distances between the figures
   t.Padding = 'compact';
   t.TileSpacing = 'compact';
   set(gcf, 'Position', [662 403 1259 448])

   cmd = strcat('print -depsc', 32, path, '/pdf_input.eps');
   eval(cmd);
    %}
   
   figure(2)
   t = tiledlayout(1,5); % Requires R2019b or later
   nexttile
   hist(sample_diam, 40)
   xlabel('D (m)', 'interpreter', 'latex')
   ylabel('count')
   nexttile
   hist(sample_rho, 40)
   xlabel('$\rho$ (kg/m$^3$)', 'interpreter', 'latex')
   nexttile
   hist(sample_gamma, 40)
   xlabel('$\gamma$ (deg)', 'interpreter', 'latex')
   nexttile
   hist(sample_dadt, 40)
   xlabel('da/dt (AU/My)', 'interpreter', 'latex')
   nexttile
   hist(sample_P(sample_P>0), 40)
   xlabel('P (h)', 'interpreter', 'latex')

   cmd = strcat('print -depsc', 32, path, '/sample_input.eps');
   eval(cmd);

   % Set the distances between the figures
   t.Padding = 'compact';
   t.TileSpacing = 'compact';
   set(gcf, 'Position', [2 529 1914 426]);

   writematrix(sample_diam,  strcat(path,'/diam_mc.txt'));
   writematrix(sample_rho,   strcat(path,'/rho_mc.txt'));
   writematrix(sample_gamma, strcat(path,'/gamma_mc.txt'));
   writematrix(sample_dadt,  strcat(path,'/dadt_mc.txt'));
   writematrix(sample_P(sample_P>0),  strcat(path,'/period_mc.txt'));
    
end
