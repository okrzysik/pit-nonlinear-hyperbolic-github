% Plot fluxes and their derivatives for the various PDEs.

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%%%%%%%%%%%%%%%%
save_fig = ~true;
%%%%%%%%%%%%%%%

umax = 1.5;
umin = -umax;


u = linspace(umin, umax, 1e3);

%pde_ids = {'burgers', 'traffic', 'buckley-leverett'};
pde_ids = {'burgers', 'buckley-leverett'};

mycols     = {'r', 'b', [0, 0.5, 0], [0.75, 0, 0.95], [0.9290, 0.6940, 0.1250], [0.3010 0.7450 0.9330]};
mymarkers  = {'x' ,'>', '*', 'p', '<', 's', '+'};
linestyles = {'-'; '-.'}; 

for pde_idx = 1:numel(pde_ids)
   
    pde_id = pde_ids{pde_idx};
    
    % Create PDE object.
    if strcmp(pde_id, 'burgers')
        my_cons_law = burgers(pde_pa);

    elseif strcmp(pde_id, 'buckley-leverett')
        my_cons_law = buckley_leverett(pde_pa);

    end
    
    f      = my_cons_law.flux(u);
    fprime = my_cons_law.flux_jacobian(u);
    
    nameval = {'Color', mycols{pde_idx}, 'LineWidth', 3, 'LineStyle', linestyles{pde_idx}};
    
    figure(1)
    plot(u, f, nameval{:}, 'DisplayName', sprintf('%s', pde_ids{pde_idx}))
    hold on
    
    
    figure(2)
    plot(u, fprime, nameval{:}, 'DisplayName', sprintf('%s', pde_ids{pde_idx}))
    hold on
    plot(u, abs(fprime), nameval{:}, 'LineWidth', 1, 'HandleVisibility', 'Off')
end

figure(1)
title('$f(u)$')
xlabel('$u$')
lh = legend();
lh.set('Location', 'Best')
axis tight

if save_fig; figure_saver(gcf, './figures/pde_flux_example'); end

figure(2)
title('$f''(u)$')
xlabel('$u$')
axis tight
if save_fig; figure_saver(gcf, './figures/pde_flux_deriv_example'); end


% Save figure
function figure_saver(fig, fig_name)

    %set(gca, 'LooseInset', get(gca,'TightInset'))
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    saveas(gcf, strcat(fig_name, '.pdf'));
    %saveas(gcf, strcat(fig_name, '.eps'), 'epsc');
    
end