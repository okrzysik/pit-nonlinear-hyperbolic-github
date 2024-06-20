% Make a plot of material parameters c0 and Z0. They are expected to be in
% disc_pa under the fields c_cc and Z_cc.
%
function fh = acoustics_plot_material_parameters(figno, mesh_pa, disc_pa, invisible_col_bar)

    if nargin == 3
        invisible_col_bar = ~true;
    end

    fh = figure(figno);
    plot(mesh_pa.x_centers, disc_pa.Z_cc, '-k',  'DisplayName', '$Z_0$', 'LineWidth', 2)
    hold on
    plot(mesh_pa.x_centers, disc_pa.c_cc, '--', 'Color', [0.75, 0, 0.95], 'DisplayName', '$c_0$', 'LineWidth', 2)
    xlabel('$x$')
    axis tight
    lh = legend();
    lh.set('Location', 'best')
    title('material properties')
    data = [disc_pa.Z_cc(:); disc_pa.c_cc(:)];
    data_min = min(data);
    data_max = max(data);
    height   = data_max - data_min;
    tol = 0.05;
    if data_min ~= data_max
        ylim([data_min - tol*height, data_max + tol*height])
    end
    xmin     = mesh_pa.xmin;
    xmax     = mesh_pa.xmax;
    xlim([xmin, xmax])
    xticks([xmin:(xmax-xmin)/10:xmax])
    xticks([xmin:2*(xmax-xmin)/10:xmax])
    
    if invisible_col_bar
        ylabel('$t$', 'Color', 0.998*[1 1 1])
        cb = colorbar(); scale = 0.998; cb.Color = scale*[1 1 1]; myColorMap = scale*ones(1, 3); colormap(gcf, myColorMap); 
    end
    
    
% This is a check to make sure I worked out the formula consistently 
% for the paper since the way I implement it is more complicated...    
%     % For periodic media problem. if we have n layers, need to map each layer 
%     % to an integer domain, so the 1st layer is [0,1], the second is [1,2], the 
%     % third is [2,3], etc. Then, say the impedance is one constant on the
%     % even intervals [k, k+1] where k is even, and it's the other constant
%     % in the odd inervals [k, k+1] where k is odd. So, if x \in [0,1] I do
%     % x <- n_layers * x, then if mod( floor(x), 2 ) is 0 I'm in
%     % an even layer, and if it's 1 then I'm in an odd layer
%     n_layers = 16;
%     x = mesh_pa.x_centers;
%     Z = zeros(size(x));
%     EvenI = find( mod(floor(n_layers * x), 2) == 0 ); % Even layers
%     OddI  = find( mod(floor(n_layers * x), 2) == 1 ); % Odd layers
%     Z(EvenI) = 1;
%     Z(OddI) = 2;
%     hold on
%     plot(x, Z, '-ro')
end