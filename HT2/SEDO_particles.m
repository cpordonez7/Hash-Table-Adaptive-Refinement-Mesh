function [y] = SEDO_particles(t,y0,BETAC,MUC)
    %PDE functions Lagranguianas particles
    % y = zeros(size(y0));
    xc = 0.0;
    yc = 0.0;

    Double_number_initial_condit = length(y0);
    Number_initial_condit = Double_number_initial_condit/2;

    % d = zeros(Num_cond_ini,1);
    y = zeros(size(y0));

    d = sqrt((y0(1:2:Double_number_initial_condit)-xc).^2+(y0(2:2:Double_number_initial_condit)-yc).^2);
    
    vel_part_ = vel_part(d,t,BETAC,MUC);   %column vector

    for iter=1:Number_initial_condit
        y(2*iter-1) = (y0(2*iter)-yc)*vel_part_(iter)/d(iter);
        y(2*iter) = -(y0(2*iter-1)-xc)*vel_part_(iter)/d(iter);
    end
end

function vel_part = vel_part(d,t,BETAC,MIUC)

    vel_part = (BETAC./d).*(1-exp(-(d.^2)/(MIUC*t)));
end
