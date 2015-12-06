function [] = Problem_2()

    %%%%%%
    % Solves the linear convection-diffusion equation using the FTCS explicit method.
    %
    % Ryan Skinner, December 2015
    %%%
    
    Set_Default_Plot_Properties();
    
    % Display results from Problem 1 at only a few time steps.
    n_plot = 11;
    T_prob1 = Problem_1(n_plot);
    
    cases = {{'a',5},{'a',6},{'b',5},{'b',6}};
    for case_i = 1:4
            
    IC_str = cases{case_i}{1};
         n = cases{case_i}{2};
    
    % Spatial domain.
    nn = 2^n + 1;
    Ly = 2.0;
    y = linspace(0, Ly, nn)';
    dy = y(2) - y(1);
    
    % Temporal domain.
    dt = 0.001;
    if strcmp(IC_str,'a')
        t_final = 4;
    else
        t_final = 4;
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
        
        Tn = T(:,t_n);
        
        for i = 1:nn
            if i == 1 ; Tnim1=Tn(nn); else Tnim1=Tn(i-1); end % Periodic BCs
            if i == nn; Tnip1=Tn( 1); else Tnip1=Tn(i+1); end % Periodic BCs
            T(i,t_n+1) = Tn(i) + dt * ( alpha * (Tnip1 - 2*Tn(i) + Tnim1) / (dy^2) ...
                                       - v(i) * (Tnip1           - Tnim1) / (2*dy) );
        end

    end

    %%%
    % Process results.
    %%%
    
    step_numbers = round(linspace(1,length(t),n_plot));
    figure(case_i);
    hold on;
    for t_n = 2:length(step_numbers)
        plot(y, T(:,step_numbers(t_n)), 'k--');
    end
    
    cmap = jet(n_plot);
    hf = figure(length(cases) + case_i);
    set(hf,'Position',[100,500,900,300]);
    hold on;
    plot_handles = [];
    for t_n = 2:length(step_numbers)
        tmp = sprintf('t = %.2f', t(step_numbers(t_n)));
        method_err = (T_prob1{case_i}(:,step_numbers(t_n)) - T(:,step_numbers(t_n)));
        hp = plot(y(2:end-1), method_err(2:end-1), 'DisplayName', tmp, 'Color', cmap(t_n,:));
        if(mod(t_n-1,5) == 0 || n_plot <= 11)
            plot_handles(end+1) = hp;
        end
    end
    title(sprintf('IC = (%s), n = %.0f',IC_str,n));
    xlabel('y');
    ylabel('T_{Fourier} - T_{FTCS}');
    ylim([-0.05,0.05]);
    xlim([0,Ly]);
    hleg = legend(plot_handles);
    set(hleg,'Location','eastoutside');
    
    end
    
    disp('Done.');
    
end













