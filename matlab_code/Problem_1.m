function [] = Problem_1()

    %%%%%%
    % Solves the linear convection-diffusion equation using the Fourier pseudo-spectral
    % method.
    %
    % Ryan Skinner, November 2015
    %%%
    
    Set_Default_Plot_Properties();
    
    % Spatial domain.
    n = 5;
    nn = 2^n + 1;
    Ly = 2.0;
    y = linspace(0, Ly, nn)';
    
    % Temporal domain.
    dt = 0.001;
    tfinal = 1;
    t = [0:dt:tfinal];
    
    % Physical parameters.
    Re = 1;
    Pr = 25;
    alpha = 1 / (Re * Pr);
    v = sin(pi*y);

    % Initialize solution.
    T = nan(nn,length(t));
%     T(:,1) = cos(2*pi*y) .* sin(pi*y); IC_str = '(a)';
    T(:,1) = cos(2*pi*y); IC_str = '(b)';

    %%%
    % Solve problem numerically.
    %%%
    
    for time_n = 1:(length(t)-1)
        
        T_n = T(:,time_n);
        
        dTdy   = find_dfdn(  T_n',nn,Ly)';
        d2Tdy2 = find_d2fdn2(T_n',nn,Ly)';
        
        % Update solution.
        T(:,time_n+1) = T_n + dt * (alpha * d2Tdy2 - v .* dTdy);

    end

    %%%
    % Process results.
    %%%
    
%     surf(y,t,T');
%     xlabel('y');
%     ylabel('t');
    
    % Time evolution of solution for a single Courant number.
    n_plot = 21;
    cmap = winter(n_plot);
    step_numbers = round(linspace(1,length(t),n_plot));
    hf = figure();
    set(hf,'Position',[100,500,900,300]);
    hold on;
    plot_handles = [];
    for time_n = 1:length(step_numbers)
        tmp = sprintf('t = %.2f', t(step_numbers(time_n)));
        hp = plot(y, T(:,step_numbers(time_n)), 'DisplayName', tmp, 'Color', cmap(time_n,:));
        if(mod(time_n-1,4) == 0)
            plot_handles(end+1) = hp;
        end
    end
    title(sprintf('IC = %s, n = %.0f',IC_str,n));
    xlabel('x');
    ylabel('u');
    ylim([-1,1]);
    xlim([0,Ly]);
    hleg = legend(plot_handles);
    set(hleg,'Location','eastoutside');
    
    disp('Done.');
    return
    
end













