function [] = Problem_2()

    %%%%%%
    % Solves the linear convection-diffusion equation using the FTCS explicit method.
    %
    % Ryan Skinner, November 2015
    %%%
    
    Set_Default_Plot_Properties();
    
    cases = {{'a',5},{'a',6},{'b',5},{'b',6}};
    for case_i = 1:4
            
    IC_str = cases{case_i}{1};
         n = cases{case_i}{2};
    
    % Spatial domain.
    nn = 2^n + 1;
    Ly = 2.0;
    y = linspace(0, Ly, nn)';
    
    % Temporal domain.
    dt = 0.001;
    if strcmp(IC_str,'a')
        t_final = 1;
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
    
    n_plot = 31;
    cmap = jet(n_plot);
    step_numbers = round(linspace(1,length(t),n_plot));
    hf = figure(case_i);
    set(hf,'Position',[100,500,900,300]);
    hold on;
    plot_handles = [];
    for t_n = 1:length(step_numbers)
        tmp = sprintf('t = %.2f', t(step_numbers(t_n)));
        if t_n == 1
            hp = plot(y, T(:,step_numbers(t_n)), 'k-.', 'LineWidth', 3, 'DisplayName', tmp);
        else
            hp = plot(y, T(:,step_numbers(t_n)), 'DisplayName', tmp, 'Color', cmap(t_n,:));
        end
        if(mod(t_n-1,5) == 0)
            plot_handles(end+1) = hp;
        end
    end
    title(sprintf('IC = (%s), n = %.0f',IC_str,n));
    xlabel('y');
    ylabel('T');
    ylim([-1,1]);
    xlim([0,Ly]);
    hleg = legend(plot_handles);
    set(hleg,'Location','eastoutside');
    
    end
    
    disp('Done.');
    
end













