function [] = Problem_1()

    %%%%%%
    % Solves the linear convection-diffusion equation using the Fourier pseudo-spectral
    % method.
    %
    % Ryan Skinner, November 2015
    %%%
    
    Set_Default_Plot_Properties();
    
    for IC_str = {'a','b'}
    for n = [5,6]
    
    % Spatial domain.
    nn = 2^n + 1;
    Ly = 2.0;
    y = linspace(0, Ly, nn)';
    
    % Temporal domain.
    dt = 0.001;
    if strcmp(IC_str,'a')
        t_final = 2;
    else
        t_final = 2;
    end
    t = [0:dt:t_final];
    
    % Physical parameters.
    Re = 1;
    Pr = 25;
    alpha = 1 / (Re * Pr);
    v = sin(pi*y);

    % Initialize solution.
    T = nan(nn,length(t));
    if strcmp(IC_str,'a')
        T(:,1) = cos(2*pi*y) .* sin(pi*y);
    else
        T(:,1) = cos(2*pi*y);
    end

    %%%
    % Solve problem numerically.
    %%%
    
    for t_n = 1:(length(t)-1)
        
        T_n = T(:,t_n);
        
        dTdy   = find_dfdn(  T_n',nn,Ly)';
        d2Tdy2 = find_d2fdn2(T_n',nn,Ly)';
        
        % Update solution.
        T(:,t_n+1) = T_n + dt * (alpha * d2Tdy2 - v .* dTdy);

    end

    %%%
    % Process results.
    %%%
    
    % Time evolution of solution for a single Courant number.
    n_plot = 31;
    cmap = winter(n_plot);
    step_numbers = round(linspace(1,length(t),n_plot));
    hf = figure();
    set(hf,'Position',[100,500,900,300]);
    hold on;
    plot_handles = [];
    for t_n = 1:length(step_numbers)
        tmp = sprintf('t = %.2f', t(step_numbers(t_n)));
        hp = plot(y, T(:,step_numbers(t_n)), 'DisplayName', tmp, 'Color', cmap(t_n,:));
        if(mod(t_n-1,5) == 0)
            plot_handles(end+1) = hp;
        end
    end
    title(sprintf('IC = (%s), n = %.0f',IC_str{1},n));
    xlabel('x');
    ylabel('u');
    ylim([-1,1]);
    xlim([0,Ly]);
    hleg = legend(plot_handles);
    set(hleg,'Location','eastoutside');
    
    end
    end
    
    disp('Done.');
    
end













