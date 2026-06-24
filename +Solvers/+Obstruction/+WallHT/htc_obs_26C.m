function [htc_w,htc_fs] = htc_obs_26C(W,H,L,th_w,th_ww,W_w,th_fs,th_wfs,power,G,e_w,e_ww,x_in,z,stdev_f,y_plus_vs)
% this function approximates the htc's of an annular flow in a rectangular channel
% flowing past a cylindrical obstacle on one side of the channel

%% inputs
    % test section
        % W = [m] test section flow channel width (also width of fto glass heater)
        % H = [m] test section flow channel height
        % L = [m] length of heater (length of active test section)
    % wake geometry (caused by obstruction)
        %th_w = [m] thickness of wake base film (determined from experiment or correlations)
        %th_ww = [m] thickness of waves in wake (determined from experiment or correlations)
        %W_obs = [m] width of obstuction
        %W_w = [m] width of wake
    % free stream geometry
        %th_fs = [m] thickness of free stream base film
        %th_wfs = [m] thickness of free stream waves
    % heat flux input to heater
        %power = [W] power input to one heater (total power in the case of only one heater)
    % input parameters based on experimental flow conditions
        %G = [kg/m^s-s] fluid mass flux
        %e_w = [-] fraction of film mass found in wake (this is a guess value for now, the team is working very hard to correlate this parameter)
        %e_ww = [-] fraction of film mass found in waves (guess value, will be correlated later)
        %x_in = [-] inlet quality of fluid
        %z = [m] location of interest
        %stdev_f = [m] standard deviation of film thickness
        %y_plus_vs = [-] thickness of viscous sublayer in wall units (it is common practice to use a value of around 5)
    
% caluated parameters based on geometry and heat inputs
    per = 2*W+2*H; %[m] wetted perimeter
    A_c = W*H; %[m^2] cross-sectional area
    D_h = 2*W*H/(W+H); %[m] hydraulic diameter of test section
    per_w = W_w; %[m] wetted perimeter of wake (wake on one side of test section channel)
    A_c_w = W_w*th_w; %[m] cross-sectional area of wake
    W_fs = W-W_w; %[m] width of free stream on one side
    per_fs = W_fs+W+2*H; %[m] wetted perimeter of free stream
    A_c_fs = W_fs*th_fs+W*th_fs+2*H*th_fs; %[m] cross-sectional area of free stream
    power_t = 2*power; %[W] power input total (half this in the case of only one heater)
    q_flux = power/(L*W); %[W/m^2] heat flux into heater
    
%Fluid Properties
    % These properties all depend solely on the saturation temperature of the
    % fluid at the desired pressure and thus, they only need to be changed if
    % the saturation temperature is changed. If this is the case, these
    % properties will need to be determined from the appropriate property table
    % or from EES.
    fluid = 'R245fa'; %fluid used in flow loop
    T_sat = 26; %[C] saturation temperature of fluid at desired pressure
    P_sat = 153335; %[Pa] saturation pressure of fluid
    mu_l = 0.0004007; %[kg/m-s] dynamic viscosity of liquid phase of fluid
    mu_v = 0.00001019; %[kg/m-s] dynamic viscosity of vapor phase of fluid
    rho_l = 1336; %[kg/m^3] density of liquid phase
    rho_v = 8.829; %[kg/m^3] density of vapor phase
    Pr_l = 6.563; %[-] prandtl number of liquid phase
    Pr_v = 0.6956; %[-] prandtl number of vapor phase
    k_l = 0.08085; %[W/m-K] thermal conductivity of liquid phase
    k_v = 0.01403; %[W/m-K] thermal conductivity of vapor phase
    sigma_l = 0.0135; %[N/m] surface tension of liquid phase
    dh_fg = 189747; %[J/kg] enthalpy of vaporization of fluid
    nu_l = mu_l/rho_l; %[m^2/s] kinematic viscosity of liquid phase
    nu_v = mu_v/rho_v; %[m^2/s] kinematic viscosity of vapor phase
    v_fg = (1/rho_v)-(1/rho_l); %[m^3/kg] specific volume difference between liquid and vapor phase
    
% calculated parameters based on flow condition inputs
    enth_data = table2array(readtable('enthalpy_quality_26Csat_R245fa.csv')); %enthalpy data versus quality (table made from EES)
    h_in = enth_data(find(abs(x_in - enth_data(:,1))==min(abs(x_in - enth_data(:,1)))),2); %[J/kg] inlet enthalpy: depends on inlet quality and saturation pressure
    x = x_in+q_flux*(2*W*z)/(G*A_c*dh_fg); %[-] quality at point of interest
    m_dot_t = G*A_c; %[kg/s] total mass flow rate
    h_out = q_flux*(2*W*z)/(G*A_c)+h_in; %[J/kg] enthalpy at location of interest
    m_dot_l = m_dot_t*(1-x); %[kg/s] mass flow rate of liquid phase
    m_dot_v = m_dot_t-m_dot_l; %[kg/s] mass flow rate of vapor phase
    u_v = m_dot_v/(rho_v*A_c); %[m/s] average velcotiy of vapor phase
    Re_v = rho_v*u_v*D_h/mu_v; %[-] reynolds number of vapor
    
% Droplet Entrainment Fraction
    % calculated using correlation from Cioncolini and Thome (2011). this
    % correlation requires (1) an iterative method to solve or (2) a
    % predictor-corrector method as outlined in by Cioncolini and Thome. The
    % predictor-corrector method is much easier to implement and is thus used
    % here. This correlation captures about 6 out 10 data points within plus or
    % minus 30% of the true measured value of entrainment.
    % predictor step
        We_p = rho_v*u_v^2*D_h/sigma_l; %[-] weber number predictor (based on true density of vapor)
        e_p = (1+279.6*We_p^(-0.8395))^(-2.209); %[-] entrained droplet fraction predictor
    % corrector step
        rho_c = (x+e_p*(1-x))/((x/rho_v)+(e_p*(1-x)/rho_l)); %[kg/m^3] vapor core density corrector (droplets plus vapor)
        We_c = rho_c*u_v^2*D_h/sigma_l; %[-] weber number corrected
        e_d = (1+279.6*We_c^(-0.8395))^(-2.209); %[-] entrained droplet fraction corrected
        m_dot_d = e_d*m_dot_l; %[kg/s] mass flow rate of liquid phase entrained in droplets
        m_dot_f = m_dot_l-m_dot_d; %[kg/s] mass flow rate of liquid phase in film
    
% Calculation of Momentum and Thermal Boundary Layers
    % friction factor calculation
        f = 0.005*(1+300*(stdev_f/D_h)); %[-] friction factor for shear stress and film-wall interface
    % wake region
        m_dot_w = e_w*m_dot_f; %[kg/s] mass flow rate of film found in wake
        m_dot_ww = e_ww*m_dot_w; %[kg/s] mass flow rate of wake found in waves
        m_dot_wbf = m_dot_w-m_dot_ww; %[kg/s] mass flow rate of wake found in base film
        u_wbf = m_dot_wbf/(rho_l*A_c_w); %[m/s] velocity of base film in wake
        tau_w = 0.5*rho_v*(u_v-u_wbf)^2*f; %[Pa] shear stress at vapor-liquid interface in wake
        u_star_w = sqrt(tau_w/rho_l); %[-] normalized wake velocity
        dmw = y_plus_vs*nu_l/u_star_w; %[m] momentum boundary layer thickness in wake
        dtw = dmw/(Pr_l^(1/3)); %[m] thermal boundary layer thickness in wake
    % free stream region
        m_dot_fs = m_dot_f-m_dot_w; %[kg/s] mass flow rate of film found in free stream
        m_dot_fsw = e_ww*m_dot_fs; %[kg/s] mass flow rate of free stream found in waves
        m_dot_fsb = m_dot_fs-m_dot_fsw; %[kg/s] mass flow rate of free stream found in base film
        u_fsb = m_dot_fsb/(rho_l*A_c_fs); %[m/s] velocity of base film in free stream
        tau_fs = 0.5*rho_v*(u_v-u_fsb)^2*f; %[Pa] shear stress at vapor-liquid interface in free stream
        u_star_fs = sqrt(tau_fs/rho_l); %[-] normalized wake velocity
        dmf = y_plus_vs*nu_l/u_star_fs; %[m] momentum boundary layer thickness in wake
        dtf = dmf/(Pr_l^(1/3)); %[m] thermal boundary layer thickness in wake
    
% Critical Film Thickness
    % correlation taken from Chun et al (2003)
    G_w = m_dot_w/A_c_w; %[kg/m^2-s] mass flux in wake
    G_fs = m_dot_fs/A_c_fs; %[kg/m^2-s] mass flux in free stream
    th_crit_w = (q_flux/(dh_fg*G_w))^0.35*v_fg*(mu_l)^2*10^(8.8*(mu_v/mu_l)^0.617)/sigma_l; %[m] minimum stable film thickness in wake
    th_crit_fs = (q_flux/(dh_fg*G_fs))^0.35*v_fg*(mu_l)^2*10^(8.8*(mu_v/mu_l)^0.617)/sigma_l; %[m] minimum stable film thickness in free stream
    
% Calculation of Heat Transfer Coefficients
    % single phase vapor htc correlation 
        htc_dry = (k_v/D_h)*0.023*Re_v^0.8*Pr_v^0.4; %[W/m^2-K] dry wake htc (Dittus-Boelter)
        f_v = (0.79*log(Re_v)-1.64)^(-2); %[-] friction factor for turbulent vapor
        htc_dry_G = k_v/D_h*(f_v/8)*(Re_v-1000)*Pr_v/(1+12.7*(f_v/8)^(1/2)*(Pr_v^(2/3)-1)); %[W/m^2-K] dry wake htc (Gnielinski)
    % wake region
        htc_w_thin = k_l/th_w; %[W/m^2-K] htc of wake when wake thickness is smaller than thermal boundary layer thickness
        htc_w_thick = k_l/dtw; %[W/m^2-K] htc of wake when wake thickness is larger than thermal boundary layer thickness
        if th_w > dtw & th_w > th_crit_w %if statment determines htc based on film thickness condition
            htc_w = htc_w_thick;
        elseif th_w < dtw & th_w > th_crit_w
            htc_w = htc_w_thin;
        else
            htc_w = htc_dry;
        end
        htc_w; %[W/m^2-K] modelled htc of wake based on film thickness
    % free stream region
        htc_fs_thin = k_l/th_fs; %[W/m^2-K] htc of wake when wake thickness is smaller than thermal boundary layer thickness
        htc_fs_thick = k_l/dtf; %[W/m^2-K] htc of wake when wake thickness is larger than thermal boundary layer thickness
        if th_fs > dtf & th_fs > th_crit_fs %if statment determines htc based on film thickness condition
            htc_fs = htc_fs_thick;
        elseif th_fs < dtw & th_fs > th_crit_fs
            htc_fs = htc_fs_thin;
        else
            htc_fs = htc_dry;
        end
        htc_fs; %[W/m^2-K] modelled htc of wake 
    
end
    
    