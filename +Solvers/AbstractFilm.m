classdef (Abstract) AbstractFilm < Solvers.AbstractField
    %ABSTRACTFILM Base class for liquid film models in three- and four-field solvers
    %
    % This abstract class defines shared methods and properties for liquid film modeling.
    % It includes calculations for film thickness, entrainment, shear forces, and velocity models.
    % Specific entrainment models (e.g., Okawa, Govan) and friction models are implemented.

    properties (SetAccess={?Solvers.AbstractField,?Solvers.AbstractSolver})

        DZ           (1,1) double  {mustBeNumeric}                         =0      % Axial step size [m]
        inputSet                   {isa(inputSet,'Inputs.InputSet')}               % :class:`Inputs.InputSet` object containing geometry, model, and boundary conditions
        fluid                      {isa(fluid,'Inputs.FluidProperties')}           % :class:`Inputs.FluidProperties` object containing fluid thermophysical properties
        mix          (1,1)         {isa(mix, 'Solvers.Mixture.Mixture')}   = NaN   % :class:`Solvers.Mixture.Mixture` object of the mixture solver
    end


    methods

        function absfilm = AbstractFilm(inputSet, fluid)
            %AbstractFilm Constructor for AbstractFilm class
            %
            % Initializes :class:`Inputs.InputSet` and :class:`Inputs.FluidProperties`

            if nargin > 0
                % Store inputSet as object property
                absfilm.inputSet = inputSet;
                absfilm.fluid  = fluid;
            end

            % Overload copyable properties
            %mix.flowProperties = {'W','U','H'};
        end

        function wl = WL(absfilm,zIdx)
            %WL Film mass flow rate per unit perimeter [kg/s/m]

            if nargin < 2, zIdx = (1:absfilm(1).NZ).'; end

            perim  = absfilm.inputSet.geometry.PERIM;                      % [m] Perimeter

            wl = absfilm.W(zIdx,:)./perim;
        end

        function thick = THICK(absfilm,zIdx)
            %THICK Film thickness [m]

            if nargin < 2, zIdx = (1:absfilm(1).NZ).'; end

            rhof  = absfilm.fluid.RHOF;                                    % [kg/m^3] Saturated liquid density

            thick = absfilm.WL(zIdx)./absfilm.U(zIdx,:)./rhof;
        end

        function thickHH = THICKHH(absfilm, zIdx)
            %THICKHH Film thickness according to Henstock and Hanratty 1976
            if nargin < 2, zIdx = (1:absfilm(1).NZ).'; end
            RE = max(absfilm.RE(zIdx), 1E-6);                              % Film Reynolds Number
            tau_i = abs(absfilm.FWALL(zIdx));                              % [Pa] Film wall shear stress
            tau_w = abs(absfilm.FVAPOR(zIdx));                             % [Pa] Film interfacial shear stress
            rhof  = absfilm.fluid.RHOF;                                    % [kg/m^3] Saturated liquid density
            mu_f = absfilm.fluid.MUF;                                      % [Pa*s] Saturated liquid viscosity
            nu_f = mu_f / rhof;                                            % [m^2/s] Saturated liquid kinematic viscosity

            tau_c = (1/3) * tau_i + (2/3) * tau_w;                         % Representative fluid shear stress
            u_l_star = sqrt(tau_c/rhof);                                   % Film friction velocity
            delta_plus =((0.707 * RE.^0.5).^2.5 + (0.037*RE.^0.9).^2.5).^0.4;% Film thickness in inner coordinates according to H&H
            thickHH = delta_plus*nu_f./u_l_star;                            % Film thickness 
        end

        function re = RE(absfilm,zIdx)
            %RE Film Reynolds number [-]

            if nargin < 2, zIdx = (1:absfilm(1).NZ).'; end

            muf   = absfilm.fluid.MUF;                                     % [kg/m^3] Saturated liquid viscosity

            re = abs(4.*absfilm.WL(zIdx)./muf);
        end

        function ment = MENT(absfilm,zIdx)
            %MENT Film entrainment mass flux [kg/m^2/s]
            % 
            % Based on selected :attr:`Inputs.Model.ENTRAINMENT` model.

            if nargin < 2, zIdx = (1:absfilm(1).NZ).'; end

            model = absfilm.inputSet.model;

            Wf = absfilm.W(zIdx,:);
            negfilm = find(Wf<0);
            Wf = abs(Wf);

            switch model.ENTRAINMENT
                case InputEnums.ENTRAINMENT.NONE
                    % Suppress film entrainment
                    ment = zeros(length(zIdx),1);

                case InputEnums.ENTRAINMENT.GOVAN
                    % Govan & Hewitt film entrainment model
                    vapor = absfilm.mix.vapor;
                    rhof  = absfilm.fluid.RHOF;                            % [kg/m^3] Saturated liquid density
                    rhog  = absfilm.fluid.RHOG;                            % [kg/m^3] Saturated vapor density
                    muf   = absfilm.fluid.MUF;                             % [kg/m^3] Saturated liquid viscosity
                    mug   = absfilm.fluid.MUG;                             % [kg/m^3] Saturated vapor viscosity
                    sig   = absfilm.fluid.SIGMA;                           % [N/m] Surface tension
                    hdiam = absfilm.inputSet.geometry.HDIAM;               % [m] Hydraulic diameter
                    area  = absfilm.inputSet.geometry.AREA;                % [m^2] Coolant area
                    perim = absfilm.inputSet.geometry.PERIM;               % [m^2] Coolant area

                    k = 5.75e-5; n1 = 0.316; n2 = 0.632;                   % Model constants
                    Wfc = muf.*exp(5.8504+0.4249*mug/muf*sqrt(rhof/rhog)).*perim./4; % [kg/s] Critical film flow rate
                    ment = k*((Wf./perim-Wfc./perim).^2*16/(rhof*sig*hdiam)).^n1.*(rhof/rhog)^n2.*vapor.W(zIdx)./area; % [kg/m^2/s] Entrainment mass flux
                    ment(Wf<=Wfc) = 0;                                     % Set to 0 below critical film flowrate

                case InputEnums.ENTRAINMENT.OKAWA2003
                    % Okawa et al. 2003 film entrainment model
                    coefs = [320 0.111 4.79E-4 1];
                    ment  = absfilm.OKAWAMENT(zIdx, coefs);

                case InputEnums.ENTRAINMENT.OKAWA2004
                    % Okawa et al. 2004 film entrainment model
                    coefs = [320 0 0.0310 2.3 0.0675 1.2 0.2950 0.5];
                    ment  = absfilm.OKAWAMENT(zIdx, coefs);

                case InputEnums.ENTRAINMENT.OKAWA2004MOD
                    % Modified Okawa et al. (2004) from Adamsson and Le Corre (2011)
                    coefs = [320 0 0.0310 2.3 0.0675 1.2];
                    ment  = absfilm.OKAWAMENT(zIdx, coefs);

                case InputEnums.ENTRAINMENT.OKAWAMFVAL
                    % Modified Okawa et al. (2004) with Rodarte (2015) data for MFVAL test facility
                    coefs = [320 0 0.0310 2.3 0.0387 0.39];
                    ment  = absfilm.OKAWAMENT(zIdx, coefs);
                
                case InputEnums.ENTRAINMENT.OKAWACT
                    % Modified Okawa et al. (2004) based on fit from Ciancolini and
                    % Thome Correlation (2012)
                    coefs = [320 0 0.0310 2.3 0.0386 0.908];
                    ment  = absfilm.OKAWAMENT(zIdx, coefs);

                case InputEnums.ENTRAINMENT.OKAWARD
                    % Modified Okawa et al. (2004) based on fit for
                    % refrigerants from Rodarte (2015)
                    coefs = [320 0 0.0310 2.3 0.0319 0.8496];
                    ment  = absfilm.OKAWAMENT(zIdx, coefs);
                    
                case InputEnums.ENTRAINMENT.OKAWAGEN
                    % Generic Okawa model
                    coefs = model.OKAWACOEFS;
                    ment  = absfilm.OKAWAMENT(zIdx, coefs);
            end

            ment(negfilm)=-ment(negfilm);
            ment = -absfilm.mix.AFDISTR(0,ment,zIdx);                      % [kg/m^2/s] Entrainment mass flux, in annular flow region only
        end

        function Mtot = MTOT(absfilm,drop,zIdx)
            %MTOT Total film mass transfer [kg/m^2/s]

            if nargin < 3, zIdx = (1:absfilm(1).NZ).'; end

            Mtot  = absfilm.MEVAP(zIdx,:)+absfilm.MENT(zIdx)+drop.MDEP(zIdx);
        end

        function Cw = CW(absfilm,zIdx)
            %CW Film wall friction factor [-]
            %
            % Based on selected :attr:`Inputs.Model.THINFILMFRIC` model.

            if nargin < 2, zIdx = (1:absfilm(1).NZ).'; end

            model = absfilm.inputSet.model;
            C = 0.005;                                                     % [-]  Wall Friction constant
            Re_transition = model.RETRANSITION;                            % [-]  Reynolds number at which transition to turbulence is assumed
            f_w_lam = model.FWLAM;                                         % Factor in numerator in laminar wall friction calculation
            RE = max(absfilm.RE(zIdx), 1E-6);                              % [-] Film Reynolds number
            

            switch model.THINFILMFRIC
                case InputEnums.THINFILMFRIC.TURBULENT
                    %
                    Cw = absfilm.CW_TURB_CALC(zIdx, C);                    % Call private method

                case InputEnums.THINFILMFRIC.LAMINAR
                    %
                    Cw = absfilm.CW_LAM_CALC(zIdx, C);                     % Call private method

                case InputEnums.THINFILMFRIC.TRANSITION
                    if RE < Re_transition
                        Cw = f_w_lam./RE;                                  % Laminar wall friction factor
                    else
                        Cw = (3.6 * log10(6.9./RE)).^(-2);                  % Colebrook equation for turbulent wall friction factor
                    end 
                case InputEnums.THINFILMFRIC.TRACE
                    f_lam = f_w_lam./RE;                                   %Laminar wall friction factor for pipe flow
                    f_turb = (3.6 * log10(6.9./RE)).^(-2);                 %Turbulent friction factor for a smooth pipe according to Haalands approximation of the Colebrook equation:
                                                                           %S.E. Haaland, "Simple and Explicit Formulas for the Friction Factor in Turbulent Pipe Flow," J. Fluids Eng., 105, 89-90, 1983.
                    
                    mask = RE > 100;                                       %For RE<100 the TRACE model diverges from CW=16/RE so use CW=16/RE directly for small RE
                    Cw = f_lam;
                    Cw(mask) = (f_lam(mask).^3 + f_turb(mask).^3).^(1/3); %TRACE model for annular flow friction factor. 
                   
            end

            Cw  = absfilm.mix.AFDISTR(C,Cw,zIdx);
        end

        function Fwall = FWALL(absfilm,zIdx)
            %FWALL Film wall shear stress [N/m^2]

            if nargin < 2, zIdx = (1:absfilm(1).NZ).'; end

            Fwall  = -0.5.*absfilm.CW(zIdx).*absfilm.fluid.RHOF.*absfilm.U(zIdx,:).^2;
        end

        function Uwall = UWALL(absfilm,zIdx)
            %UWALL Film wall friction velocity [m/s]

            if nargin < 2, zIdx = (1:absfilm(1).NZ).'; end

            rho_ls = absfilm.fluid.RHOF;                                   % [kg/m^3] Saturated liquid density

            Uwall  = (-absfilm.FWALL(zIdx)./rho_ls).^0.5;
        end

        function thick = YPLUS2THICK(absfilm, yplus, zIdx)
            %YPLUS2THICK Thickness based on wall unit value [m]
            %Uses film velocity as friction velocity

            if nargin < 3, zIdx = (1:absfilm(1).NZ).'; end

            rho_ls = absfilm.fluid.RHOF;                                   % Saturated liquid mass density
            mu_ls = absfilm.fluid.MUF;                                     % Saturated liquid viscosity
            nu_ls = mu_ls./rho_ls;                                         % Saturated liquid kinematic viscosity

            thick = yplus./absfilm.UWALL(zIdx).*nu_ls;                     % Converted thickness [m]
        end

        function thick = YPLUSTHICKTAUW(absfilm, yplus, zIdx)
            %YPLUSTHICKTAUW Thickness based on wall unit value [m]
            % Uses the wall shear stress based friction velocity

            if nargin < 3, zIdx = (1:absfilm(1).NZ).'; end

            rho_ls = absfilm.fluid.RHOF;                                   % Saturated liquid mass density
            mu_ls = absfilm.fluid.MUF;                                     % Saturated liquid viscosity
            nu_ls = mu_ls./rho_ls;                                         % Saturated liquid kinematic viscosity
            tau_w = abs(absfilm.FWALL(zIdx));                                   % Wall sheat stress
            u_star = (tau_w/rho_ls).^(0.5);                                % Friction velocity

            thick = yplus./u_star.*nu_ls;                     % Converted thickness [m]
        end

        function Cv = CV(absfilm,zIdx)
            %CV Film/vapor interfacial friction factor [-]
            %
            % Based on selected :attr:`Inputs.Model.VAPORFRIC` model.

            if nargin < 2, zIdx = (1:absfilm(1).NZ).'; end

            model = absfilm.inputSet.model;
            nwall = absfilm.inputSet.geometry.NWALL;                       % Number of walls
            %C = model.VAPORFRICCST;                                        % [-] Friction constant
            Re_vap = absfilm.mix.vapor.RE(zIdx);                            %Vapor Reynolds Number
            C = (3.6 * log10(6.9./Re_vap)).^(-2);
            

            switch model.VAPORFRIC
                case InputEnums.VAPORFRIC.CONSTANT
                    %
                    Cv = repmat(C,length(zIdx),nwall);                     % [-]

                case InputEnums.VAPORFRIC.WALLIS
                    %
                    vf = absfilm.mix.vapor.VF(zIdx);                       % [-]
                    Cv = C.*(1+75.*(1-vf));                                % [-]
                    Cv = repmat(Cv,1,nwall);                               % [-]

                case InputEnums.VAPORFRIC.WALLISTHICK
                    %
                    thick = abs(absfilm.THICK(zIdx));                      % [m] Film thickness
                    Cv = absfilm.CV_WALLISTHICK_CALC(thick,C);             % Call private method

                case InputEnums.VAPORFRIC.SMOOTH
                    Re_vap = absfilm.mix.vapor.RE(zIdx);                   %Vapor Reynolds Number
                    Cv = (3.6 * log10(6.9./Re_vap)).^(-2);                 %Turbulent friction factor for a smooth pipe according to Haalands approximation of the Colebrook equation:
                                                                           %S.E. Haaland, "Simple and Explicit Formulas for the Friction Factor in Turbulent Pipe Flow," J. Fluids Eng., 105, 89-90, 1983.
            end

            Cv  = absfilm.mix.AFDISTR(0,Cv,zIdx);                          % [-]
        end

        function Fvapor = FVAPOR(absfilm,zIdx)
            %FVAPOR Vapor shear stress on film [N/m^2]

            if nargin < 2, zIdx = (1:absfilm(1).NZ).'; end

            UVAP = absfilm.mix.vapor.U(zIdx);                              % [m/s] Vapor velocity

            Fvapor = 0.5.*absfilm.CV(zIdx).*absfilm.fluid.RHOG.*(UVAP-absfilm.U(zIdx,:)).^2; % [N/m^2]
        end

        function Fbuoy = FBUOY(absfilm,zIdx)
            %FBUOY Film buoyancy force [N/m^2]

            if nargin < 2, zIdx = (1:absfilm(1).NZ).'; end

            thick = abs(absfilm.THICK(zIdx));                              % [m] Film thickness
            DPDZ = absfilm.mix.DP.Tot(zIdx)/absfilm.DZ;                    % [Pa/m] Pressure gradient

            Fbuoy = thick.*(DPDZ);                                         % [N/m^2]

            Fbuoy = absfilm.mix.AFDISTR(0,Fbuoy,zIdx);
        end

        function Fgrav = FGRAV(absfilm,zIdx)
            %FGRAV Gravitational force on film [N/m^2]

            if nargin < 2, zIdx = (1:absfilm(1).NZ).'; end

            model = absfilm.inputSet.model;
            thick = abs(absfilm.THICK(zIdx));                              % [m] Film thickness

            Fgrav = -thick.*(model.G*cos(model.ANGLE*pi/180)*absfilm.fluid.RHOF); % [N/m^2]

            Fgrav = absfilm.mix.AFDISTR(0,Fgrav,zIdx);
        end

        function Fdep = FDEP(absfilm,drop,zIdx)
            %FDEP Drop deposition shear force [N/m^2]

            if nargin < 3, zIdx = (1:absfilm(1).NZ).'; end

            dep  = drop.MDEP(zIdx);                                        % [kg/m^2/s] Drop deposition mass flux

            Fdep = (drop.U(zIdx)-absfilm.U(zIdx,:)).*dep;                  % [N/m^2]

            Fdep = absfilm.mix.AFDISTR(0,Fdep,zIdx);
        end

        function Ftot = FTOT(absfilm,drop,zIdx)
            %FTOT Total film forces per unit wall area [N/m^2]

            if nargin < 3, zIdx = (1:absfilm(1).NZ).'; end

            Ftot  = absfilm.FWALL(zIdx)+absfilm.FVAPOR(zIdx)+absfilm.FBUOY(zIdx)+absfilm.FGRAV(zIdx)+absfilm.FDEP(drop,zIdx);  % [N/m^2]
        end

        function Ualgebr = UALGEBR(absfilm,zIdx)
            %UALGEBR Film velocity using simple algebraic model [m/s]
            %
            % From :cite:t:`ADAMSSON2014316` (eq. 21 or 23).

            if nargin < 2, zIdx = (1:absfilm(1).NZ).'; end

            TAUW  = absfilm.mix.TAUW(zIdx);                                % [Pa] Wall shear stress
            Cw = absfilm.CW(zIdx);

            Ualgebr = sqrt((2*TAUW/absfilm.fluid.RHOF)./Cw);               % [m/s]

            Ualgebr  = absfilm.mix.AFDISTR(absfilm.mix.liquid.U(zIdx),Ualgebr,zIdx);
        end

        function Uequil = UEQUILS(absfilm,zIdx)
            %UEQUILS Film velocity assuming equilibrium between wall and vapor
            % shear  [m/s]

            if nargin < 2, zIdx = (1:absfilm(1).NZ).'; end

            iter(1).U = absfilm.U(zIdx,:);
            iter(1).Ftot = absfilm.FVAPOR(zIdx)+absfilm.FWALL(zIdx);

            iter(2).U = iter(1).U+0.1;
            absfilm.U(zIdx,:)=iter(2).U;
            iter(2).Ftot = absfilm.FVAPOR(zIdx)+absfilm.FWALL(zIdx);

            eps = 1.0;
            for k = 3:100
                Uiter = iter(k-2).U-iter(k-2).Ftot.*(iter(k-1).U-iter(k-2).U)./(iter(k-1).Ftot-iter(k-2).Ftot);
                iter(k).U = (1-eps).*iter(k-1).U+eps.*Uiter;
                absfilm.U(zIdx,:)=iter(k).U;
                iter(k).Ftot = absfilm.FVAPOR(zIdx)+absfilm.FWALL(zIdx);
                err = max(abs(iter(k).Ftot),[],'all');
                if err<1E-3, break; end
            end
            if err > 1E-3
                fprintf('%s UEQUILS model : not converged -> err=%0.4f\n',class(absfilm), err);
            end

            Uequil = absfilm.mix.AFDISTR(absfilm.mix.liquid.U(zIdx),absfilm.U(zIdx,:),zIdx);
        end

        function Uequil = UEQUIL(absfilm,drop,zIdx)
            %UEQUIL Film velocity using full force equilibrium model [m/s]

            if nargin < 3, zIdx = (1:absfilm(1).NZ).'; end

            iter(1).U = absfilm.U(zIdx,:);
            iter(1).Ftot = absfilm.FTOT(drop,zIdx);

            iter(2).U = iter(1).U+0.1;
            absfilm.U(zIdx,:)=iter(2).U;
            iter(2).Ftot = absfilm.FTOT(drop,zIdx);

            eps = 1.0;
            for k = 3:100
                Uiter = iter(k-2).U-iter(k-2).Ftot.*(iter(k-1).U-iter(k-2).U)./(iter(k-1).Ftot-iter(k-2).Ftot);
                iter(k).U = (1-eps).*iter(k-1).U+eps.*Uiter;
                absfilm.U(zIdx,:)=iter(k).U;
                iter(k).Ftot = absfilm.FTOT(drop,zIdx);
                err = max(abs(iter(k).Ftot),[],'all');
                if err<1E-3, break; end
            end
            if err > 1E-3
                fprintf('%s UEQUIL model : not converged -> err=%0.4f\n',class(absfilm), err);
            end

            Uequil = absfilm.mix.AFDISTR(absfilm.mix.liquid.U(zIdx),absfilm.U(zIdx,:),zIdx);
        end

        function X = X(absfilm,zIdx)
            %X Vapor quality in the near-wall region [-]
            %
            % Vapor quality in the near-wall region, considering only the
            % film for the liquid phase. The near-wall region is defined in
            % the mixture solver.
            %
            % Inputs:
            %
            % - drop  — :class:`Solvers.ThreeField.Drop` object
            % - zIdx  — Axial indices to evaluate (optional)

            if nargin < 2, zIdx = (1:absfilm(1).NZ).'; end

            X = 1-absfilm.W(zIdx,:)./absfilm.mix.NEARWALL.W(zIdx,:);
        end

        function stableThick = STABLETHICK(absfilm, zIdx)
            % Minimum stable base film thickness [m]
            %
            % Computes the minimum thickness a stable film will exist at 
            % before fracturing into rivulets

            model = absfilm.inputSet.model;

            if nargin < 2, zIdx = (1:absfilm(1).NZ).'; end

            switch model.MINSTABLETHICK
                case InputEnums.MINSTABLETHICK.CHUN
                    hflux = absfilm.HFLUX(zIdx,:);      % [W/m^2] wall heat flux 
                    hfg = absfilm.fluid.HFG;                               % [J/kg] latent heat of vaporization
                    v_f = 1/ (absfilm.fluid.RHOF);                         % [m^3/kg] specific volume of liquid phase
                    v_g = 1/ (absfilm.fluid.RHOG);                         % [m^3 /kg] specific volume of vapor phase
                    v_fg = v_g - v_f ;                                     % [m^3 /kg] specific volume difference between vapor and liquid phases
                    mu_f = absfilm.fluid.MUF;                              % [Pa*s] viscosity of liquid phase
                    mu_g = absfilm.fluid.MUG;                              % [Pa*s] viscosity of vapor phase
                    sigma = absfilm.fluid.SIGMA;                           % [N/m] Surface tension 
                    G = absfilm.mix.liquid.MFLUX(zIdx);                    % [kg/m^2] film mass flux 

                    stableThick = (hflux./(hfg*G)).^0.35 * v_fg*mu_f^2/sigma * 10^(8.8*(mu_g/mu_f)^0.617); %minimum stable film thickness

            end  
        end

        function htc = HTC(absfilm, zIdx)
            %HTC Heat transfer coefficient (W/m^2-K) 
            %
            % Computes the heat transfer coefficient for wet and dry films

            if nargin < 2, zIdx = (1:absfilm(1).NZ).'; end

            model = absfilm.inputSet.model;
            nwall = absfilm.inputSet.geometry.NWALL;                       % Number of walls

            PR_film = absfilm.fluid.PRANDTLF;                              % [-] Prandtl number for liquid film
            k_film = absfilm.fluid.KF;                                     % [W/m-k] thermal conductivity of saturated liquid
            k_vapor = absfilm.fluid.KG;                                    % [W/m-k] thermal conductivity of saturated vapor
            nu_vap = absfilm.mix.vapor.NU;                                 % [-] Nusselt number of the vapor phase
            D_h = absfilm.inputSet.geometry.HDIAM;                         % [m] hydraulic diameter of flow channel

            yplus = model.BASEYPLUS;                                       % [-] y-plus value base film is assumed to end at
            baseThick = absfilm.YPLUSTHICKTAUW(yplus, zIdx);               % [m] thickness of the basefilm
            %thermThick = baseThick/(PR_film^(1/3));                       % [m] thermal boundary layer thickness
            thick = absfilm.THICK(zIdx);                                   % [m] film thickness
            minstableThick = absfilm.STABLETHICK(zIdx);                    % [m] minimum stable film thickness

            % Expand nu_vap from [zIdx x 1] to [zIdx x nwall] so element-wise ops work
            nu_vap_wall = repmat(nu_vap, 1, nwall);                        % [zIdx x nwall]
            
            % Pre-allocate htc
            htc = zeros(length(zIdx), nwall);                              % [zIdx x nwall]
            
            % Boolean masks for each regime (element-wise, covers every [z, wall] pair)
            mask_dry    = (thick < minstableThick) | (thick == 0);         % film too thin or absent
            mask_base   = ~mask_dry & (baseThick <= thick);              % base film thinner than film
            mask_thin   = ~mask_dry & (thick < baseThick);                 % film thinner than base film
            
            % Apply HTC formula for each regime
            htc(mask_dry)  = nu_vap_wall(mask_dry)  .* (k_vapor / D_h);    % [W/m^2-K] dry/unstable
            htc(mask_base) = k_film ./ baseThick(mask_base);               % [W/m^2-K] base film limits
            htc(mask_thin) = k_film ./ thick(mask_thin);                   % [W/m^2-K] film limits
        end

        function htc_wake = HTCWAKE(absfilm, dry, zIdx)
            %HTC Heat transfer coefficient (W/m^2-K) 
            %
            % Computes the heat transfer coefficient post obstruction
            % Differs from HTC(absfilm,zIdx) in that dryout is determined
            % via wake dryout model, not via film thickness. This allows for non-constant width dry wakes
            %
            %INPUTS:
            % dry [zIdx x nwall] boolean of whether this node is dry or not

            if nargin < 3, zIdx = (1:absfilm(1).NZ).'; end

            model = absfilm.inputSet.model;
            nwall = absfilm.inputSet.geometry.NWALL;                       % Number of walls

            k_film = absfilm.fluid.KF;                                     % [W/m-k] thermal conductivity of saturated liquid
            k_vapor = absfilm.fluid.KG;                                    % [W/m-k] thermal conductivity of saturated vapor
            nu_vap = absfilm.mix.vapor.NU;                                 % [-] Nusselt number of the vapor phase
            D_h = absfilm.inputSet.geometry.HDIAM;                         % [m] hydraulic diameter of flow channel

            yplus = model.BASEYPLUS;                                       % [-] y-plus value base film is assumed to end at
            baseThick = absfilm.YPLUSTHICKTAUW(yplus, zIdx);               % [m] thickness of the basefilm
            thick = absfilm.THICK(zIdx);                                   % [m] film thickness
            

            % Expand nu_vap from [zIdx x 1] to [zIdx x nwall] so element-wise ops work
            nu_vap_wall = repmat(nu_vap, 1, nwall);                        % [zIdx x nwall]
            
            % Pre-allocate htc
            htc_wake = zeros(length(zIdx), nwall);                         % [zIdx x nwall]
            
            % Boolean masks for each regime (element-wise, covers every [z, wall] pair)
            mask_dry    = dry;                                             % film is dry, this is determined via the dry boolean
            mask_base   = ~mask_dry & (baseThick <= thick);                % base film thinner than film
            mask_thin   = ~mask_dry & (thick < baseThick);                 % film thinner than base film
            
            % Apply HTC formula for each regime
            htc_wake(mask_dry)  = nu_vap_wall(mask_dry)  .* (k_vapor / D_h);    % [W/m^2-K] dry/unstable
            htc_wake(mask_base) = k_film ./ baseThick(mask_base);               % [W/m^2-K] base film limits
            htc_wake(mask_thin) = k_film ./ thick(mask_thin);                   % [W/m^2-K] film limits
        end
    end



    methods(Access = private)

        function Cw = CW_TURB_CALC(absfilm,zIdx,C)
            %CW_TURB_CALC Private method to calculate the turbulent wall
            % friction factor [-]

            nwall = absfilm.inputSet.geometry.NWALL;                       % Number of walls
            Cw = repmat(C,length(zIdx),nwall);                             % [-]
        end

        function Cw = CW_LAM_CALC(absfilm,zIdx,C)
            %CW_LAM_CALC Private method to calculate the laminar wall
            % friction factor [-]

            RE = max(absfilm.RE(zIdx),1E-6);
            Cw = max(16./RE,C);                                            % [-]
        end

        function Cv = CV_WALLISTHICK_CALC(absfilm, thick, C)
            %CV_WALLISTHICK_CALC Private method to calculate the interfacial
            % shear factor [-]
            % 
            % Based on the :attr:`InputEnums.VAPORFRIC` = `WALLISTHICK` model.

            area = absfilm.inputSet.geometry.AREA;                         % [m^2] Cross-section area
            perim = absfilm.inputSet.geometry.PERIM;                       % [m]   Perimeter(s)
            Cv = C.*(1+(75/area).*sum(perim.*thick,2));                    % [-]
        end

        function entnum = ENTNUM(absfilm, zIdx)
            %ENTNUM Private method to calculate the entrainment number [-]
            %
            % Based on Okawa model assumptions.

            vapor = absfilm.mix.vapor;
            rhof  = absfilm.fluid.RHOF;                                    % [kg/m^3] Saturated liquid density
            rhog  = absfilm.fluid.RHOG;                                    % [kg/m^3] Saturated vapor density
            sig   = absfilm.fluid.SIGMA;                                   % [N/m] Surface tension
            area  = absfilm.inputSet.geometry.AREA;                        % [m^2] Coolant area
            perim = absfilm.inputSet.geometry.PERIM;                       % [m^2] Coolant area

            Wf = abs(absfilm.W(zIdx,:));

            % Wall friction factor (model consistent with entrainment correlation derivation)
            Cw = absfilm.CW_LAM_CALC(zIdx, 0.005);                         % [-] Wall friction factor, C=0.005

            % Film thickness (model consistent with entrainment correlation derivation)
            %delta0 = Wf./film.UEQUILS(mix,zIdx)./perim./rhof;              % Could use this simpler option instead if VAPORFRIC=WALLISTHICK could be selected specifically for this calculation
            slip = ones(size(Wf)); err=1;                                  % [-, -] Set initial guess and error for delta search
            for it = 1:100
                delta = (rhog/rhof).*slip.*Wf./max(1e-10,vapor.W(zIdx)).*area./perim; % [m] Film thickness(es)
                Cv = absfilm.CV_WALLISTHICK_CALC(delta, 0.005);            % [-] Interfacial friction factor, thick=delta, C=0.005
                newslip = sqrt(Cw./Cv.*(rhof/rhog));                       % [-] Slip formulation
                err = max(abs((newslip)./(slip)-1));                       % [-] Error
                slip = newslip;                                            % [-] Update slip
                if err < 0.01; break
                end
            end
            if err > 0.01
                disp('Okawa correlation : not converged')
            end

            entnum = Cv.*rhog.*absfilm.mix.JG(zIdx).^2.*delta./sig;        % [-] Entrainment number
        end

        function ment = OKAWAMENT(absfilm, zIdx, coefs)
            %OKAWAMENT Private method to calculate the entrainment mass flux [kg/m^2/s]
            %
            % Based on Okawa models.

            rhof   = absfilm.fluid.RHOF;                                   % [kg/m^3] Saturated liquid density
            rhog   = absfilm.fluid.RHOG;                                   % [kg/m^3] Saturated vapor density
            entnum = absfilm.ENTNUM(zIdx);                                 % [-] Entrainment number

            Refc = coefs(1); n = coefs(2);
            [ke, n2] = absfilm.OKAWACOEFS(entnum, coefs);

            ment = ke.*rhof.*entnum.^n2.*(rhof/rhog).^n;                   % [kg/m^2/s] Entrainment mass flux

            %ment(film.RE(zIdx)<=Refc) = 0;                                 % Set to 0 below critical film Reynolds
            deltaRe = 100;                                                 % [-] Set a Re window size across critical film Reynolds
            mult = min(1,max(0,(absfilm.RE(zIdx)-(Refc-deltaRe/2))./deltaRe)); % Set to 0 (linearly across Re window to avoid potential non-convergence)
            ment = mult.*ment;
        end

        function [ke, n2] = OKAWACOEFS(absfilm, entnum, coefs)
            %OKAWACOEFS Private method to determine the coefficients used in
            % :attr:`InputEnums.ENTRAINMENT` = `OKAWA` entrainment models
            % 
            % Based on the calculated entrainment number

            n2array = coefs(4:2:end);
            entbp = coefs(5:2:end);

            for k = 1:numel(entnum)
                kearray = coefs(3);
                for i = find(entbp-entnum(k)<=0)
                    kearray(i+1) = kearray(i)*entbp(i)^n2array(i)/entbp(i)^n2array(i+1);
                end
                ke(k) = kearray(end);
                n2(k) = n2array(length(kearray));
            end
            ke = reshape(ke,size(entnum)); n2 = reshape(n2,size(entnum));
        end

    end

    methods(Access = protected)

        function cpObj = copyElement(obj)
            %COPYELEMENT Override copyElement method to ensure correct references
            %with properties liquid and vapor during object copying

            % Make a shallow copy of all four properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);
        end

    end

end