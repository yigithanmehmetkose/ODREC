classdef ODREC
    %ODREC One Dimensional REgenerative Cooling
    % Conducts a regenerative cooling analysis. User can provide thrust
    % chamber profile as an external input. If a profile input is not
    % given, a bell contour nozzle and combustion chamber will be created
    % by the built-in function.
    % Function ODREC accepts invariant inputs which are as follows:
    % Function Thermal_Analysis accepts inputs that are variable
    % throughout an optimization which are as follows:

    properties (Access = public)
        I_sp            %specific impulse
        C_F             %thrust coefficient
        C_star          %characteristic velocity
        Tem_c           %coolant temperature
        p_c             %coolant pressure
        q               %wall heat flux
        T_w             %wall temperature
        Ma              %gas Mach number
        X               %axial coordinate
        Z               %radial coordinate
        X_g             %vapor fraction
        mdot_c          %coolant mass flowrate
        Cdot            %heat capacity
        htc             %heat transfer coefficient
    end

    properties (Access = private)
        %General Properties
        Fuel            %fuel name
        Oxidizer        %oxidizer name
        T_o_c           %oxidizer temperature in combustion chamber
        T_f_c           %fuel temperature in combustion chamber
        discretization  %discretization of each thrust chamber section
        eps             %surface roughness
        i_c             %coolant mass flowrate
        d_t             %throat diameter
        geom            %user defined geometry (if any)
        P_c_i           %initial coolant pressure
        T_c_i           %initial coolant temperature
        K_wall          %wall thermal conductivity
        exp_rat;        %nozzle expansion ratio
        L_star          %length of combustion chamber
        theta_n         %nozzle inlet angle in deg.
        theta_e         %nozzle exit angle in deg.
        %Coolant Thermophysical Properties
        T_l             %liquid temperature values
        mu_l            %liquid dynamic viscosity values
        K_l             %liquid thermal conductivity values
        Cp_l            %liquid specific heat values
        rho_l           %liquid density values
        P_sat           %pressure values on saturation line
        T_sat           %temperature values on saturation line
        h_f             %saturated liquid specific enthalpies
        h_g             %saturated vapor specific enthalpies
        mu_f            %saturated liquid dynamic viscosities
        mu_g            %saturated vapor specific enthalpies
        Cp_f            %saturated liquid specific heats
        Cp_g            %saturated vapor specific heats
        K_f             %saturated liquid thermal conductivities
        K_g             %saturated vapor thermal conductivities
        rho_f           %saturated liquid densities
        rho_g           %saturated vapor densities
        T_v             %vapor temperatures
        mu_v            %vapor dynamic viscosities
        K_v             %vapor thermal conductivities
        Cp_v            %vapor specific heats
        rho_v           %vapor densities
        P_1             %first pressure value in thermophysical property tables
        P_l             %last pressure value in thermophysical property tables
        N_p             %number of thermophysical property tables
    end

    methods  (Access = public)
        %% Class construction
        function obj = ODREC(varargin)
            %ODREC Construct an instance of ODREC analysis class.

            if nargin == 4      %if user does not provide thrust chamber profile
                inputs = varargin{1};
                l_prop = varargin{2};
                s_prop = varargin{3};
                v_prop = varargin{4};
                obj.geom = [];
            elseif nargin == 5  %if user does provide thrust chamber profile
                inputs = varargin{1};
                l_prop = varargin{2};
                s_prop = varargin{3};
                v_prop = varargin{4};
                obj.geom = varargin{5};
            end

            %Property Extraction:
            obj.Fuel = inputs.fuel;
            obj.Oxidizer = inputs.oxidizer;
            obj.T_o_c = inputs.T_o_c;
            obj.T_f_c = inputs.T_f_c;
            obj.discretization = inputs.discretization;
            obj.eps = inputs.eps;
            obj.i_c = inputs.i_c;


            obj.P_c_i = inputs.P_c_i;
            obj.T_c_i = inputs.T_c_i;
            obj.K_wall = inputs.k_wall;
            obj.L_star = inputs.L_star;
            obj.theta_e = inputs.theta_e;
            obj.theta_n = inputs.theta_n;
            obj.P_1 = inputs.P_1;
            obj.P_l = inputs.P_l;
            obj.N_p = inputs.N_p;

            %Liquid Coolant Properties:
            obj.T_l = l_prop(:,1);           %liquid temperatures
            obj.rho_l = l_prop(:,2);         %liquid dynamic viscosities
            obj.Cp_l = l_prop(:,3);          %liquid Prandtl numbers
            obj.mu_l = l_prop(:,4);          %liquid thermal conductivities
            obj.K_l = l_prop(:,5);           %liquid specific heats

            %Saturated Coolant Properties:
            obj.T_sat = s_prop(:,1);         %saturation pressures
            obj.P_sat = s_prop(:,2);         %saturation temperatures
            obj.rho_f = s_prop(:,3);         %saturated liquid densities
            obj.h_f = s_prop(:,4);           %saturated liquid enthalpies
            obj.Cp_f = s_prop(:,5);          %saturated liquid dynamic viscosities
            obj.mu_f = s_prop(:,6);          %saturated liquid dynamic viscosities
            obj.K_f = s_prop(:,7);           %saturated liquid thermal conductivities
            obj.rho_g = s_prop(:,8);         %saturated vapor densities
            obj.h_g = s_prop(:,9);           %saturated vapor enthalpies
            obj.Cp_g = s_prop(:,10);         %saturated liquid dynamic viscosities
            obj.mu_g = s_prop(:,11);         %saturated vapor dynamic viscosities
            obj.K_g = s_prop(:,12);          %saturated vapor thermal conductivities

            %Vapor Coolant Properties
            obj.T_v(:,:) = v_prop(:,1,:);    %vapor temperatures
            obj.rho_v(:,:) = v_prop(:,2,:);  %vapor dynamic viscosities
            obj.Cp_v(:,:) = v_prop(:,3,:);   %vapor Prandtl numbers
            obj.mu_v(:,:) = v_prop(:,4,:);   %vapor thermal conductivities
            obj.K_v(:,:) = v_prop(:,5,:);    %vapor specific heats
        end
        %% Geometrical calculation
        function obj = Thermal_Analysis(obj,inputs)
            %Calculates the thermophysical properties of the combustion gas
            %using NASA CEA program. Finds the throat dimension, mass
            %flowrates, and rocket performance parameters. Carries a
            %1-dimensional thermal and hydraulic analysis. Phase change is
            %accounted.

            %Constants
            g = DimVar(9.81,'m/s^2');                   %gravitational acceleration
            R_u = DimVar(8134,'kg-m^2/s^2-kmol-K');     %universal gas constant

            %Parameter Arrangement
            N_t = sum(obj.discretization);              %number of discretized locations
            N_v = length(obj.T_v(:,1));
            N_l = length(obj.T_l);
            T_c(N_t) = obj.T_c_i;                       %coolant inlet temperature
            P_c(N_t) = obj.P_c_i;                       %coolant inlet pressure
            k_wall = obj.K_wall;                        %wall thermal conductivity
            P_v = linspace(obj.P_1,obj.P_l,obj.N_p);    %vapor pressures
            T_int = zeros(1,N_t);                       %interpolated temperature (for phase calculation)

            %Read Inputs
            N = inputs.N;
            w = inputs.w;
            h = inputs.h;
            y = inputs.y;

            %Engine Properties
            P_cc = inputs.Pc;
            F = inputs.F;
            OF = inputs.OF;
            P_e = inputs.Pe;

            %CEA Integration
            reactants = [
                CEA.Reactant(obj.Fuel,          ...
                'Type','Fuel',                  ...
                'T',DimVar(obj.T_f_c,'K'),      ...
                'Q',DimVar(1,'kg'))             ...
                CEA.Reactant(obj.Oxidizer,      ...
                'Type','ox',                    ...
                'T',DimVar(obj.T_o_c,'K'),      ...
                'Q',DimVar(1,'kg'))             ...
                ];

            run = CEA.Run(reactants,            ...
                'ProblemType','Rocket',         ...
                'Flow','eq',                    ...
                'Pc',DimVar(P_cc,'Pa'),         ...
                'O/F',OF,                       ...
                'Outputs',{'t','p','gamma','m','viscosity','prandtl','conductivity','rho'});

            %Extracting gas properties
            k_cc = run.gamma(1);                %specific heat ratio in combustion chamber
            k_t = run.gamma(2);                 %specific heat ratio in throat
            P_t = run.Pressure(2).Value;        %pressure in throat
            M = run.MolarMass(2).Value;         %molar mass
            T_t = run.Temperature(2).Value;     %throat temperature
            T_cc = run.Temperature(1).Value;    %combustion chamber temperature
            mu_cc = run.viscosity(1).Value;     %viscosity in combustion chamber
            mu_t = run.viscosity(2).Value;      %viscosity in throat
            Pr_cc = run.prandtl(1);             %Prandtl number at combustion chamber
            Pr_t = run.prandtl(2);              %Prandtl number at throat

            %Basic propulsion properties
            k = (k_cc+k_t)/2;
            R = R_u.Value/M;                                                                    %specific gas constant
            Ma_e = sqrt(2/(k-1)*((P_e/P_t)^((1-k)/k)*(1+(k-1)/2)-1));                           %exit Mach number
            A_t = F/(P_cc*k*sqrt((2/(k+1))^((k+1)/(k-1))*2/(k-1)*(1-(P_e/P_cc)^((k-1)/k))));    %throat area
            D_t = sqrt(4*A_t/pi);                                                               %throat diameter
            obj.d_t = D_t;
            obj.exp_rat = 1/Ma_e*((1+(k-1)/2*Ma_e^2)/(1+(k-1)/2))^((k+1)/(2*k-2));              %nozzle expansion ratio
            mdot = A_t*P_cc*k*sqrt((2/(k+1))^((k+1)/(k-1))/(k*R*T_t));                          %total propellant mass flowrate
            mdot_o = OF/(OF+1)*mdot;                                                            %oxidizer mass flowrate
            mdot_f = mdot/(OF+1);                                                               %fuel mass flowrate
            mdot_p = [mdot_f mdot_o];

            %Performance
            obj.C_star = A_t*P_cc/mdot;         %characteristic velocity
            obj.C_F = F/(obj.C_star*mdot);      %thrust coefficient
            obj.I_sp = F/(mdot*g);              %specific impulse
            obj.mdot_c = mdot_p(obj.i_c)/N;     %coolant mass flowrate

            %Hot Gas Properties
            Cp_cc = R*k_cc/(k_cc-1);                                            %specific heat at combustion chamber
            Cp_t = R*k_t/(k_t-1);                                               %specific heat at throat
            mu_hg = @(T) ((mu_cc-mu_t)*T + mu_t*T_cc - mu_cc*T_t)/(T_cc-T_t);   %dynamic viscosity of the gas
            Pr_hg = @(T) ((Pr_cc-Pr_t)*T + Pr_t*T_cc - Pr_cc*T_t)/(T_cc-T_t);   %Prandtl number of the gas
            Cp_hg = @(T) ((Cp_cc-Cp_t)*T + Cp_t*T_cc - Cp_cc*T_t)/(T_cc-T_t);   %specific heat of the gas

            %Supersonic Flow Area Ratio
            k = (k_cc+k_t)/2;                                                   %gas specific heat ratio
            g = @(M,r) 1/M*((1+(k-1)/2*M^2)/(1+(k-1)/2))^((k+1)/(2*k-2)) - r;   %compressible flow area ratio equation

            %Channel Parameter Arrangement
            if length(w) == 1
                w = w*ones(1,N_t);
            end

            if length(h) == 1
                h = h*ones(1,N_t);
            end

            if length(y) == 1
                y = y*ones(1,N_t);
            end

            if isempty(inputs.w_b) == 0
                w_b = inputs.w_b;
                if length(w_b) == 1
                    w_b = w_b*ones(1,N_t);
                end
            end

            %Contour
            if isempty(obj.geom) == 1   %if user does not provide thrust chamber profile
                [x,Y] = Thrust_Chamber_Profile(obj);
            else                        %if user does provide thrust chamber profile
                x = obj.geom(:,1);
                Y = obj.geom(:,2);
            end

            % Mach Number Calculation
            [r_t,i_t] = min(Y);
            for i=1:N_t
                r = (Y(i)/r_t)^2;
                sol = @(x) g(x,r);
                if i>=i_t
                    M(i) = fzero(sol,[1,4]);
                else
                    M(i) = fzero(sol,[0.01,1]);
                end
            end

            % Geometry Calculation
            for i=N_t:-1:2
                L(i) = sqrt((Y(i)-Y(i-1))^2+(x(i)-x(i-1))^2);                                   %length of each section
                if Y(i-1) > Y(i)
                    K(i) = ((Y(i)/Y(i-1))^2-1)^2;                                               %expansion loss coefficient
                else
                    K(i) = 0.5-0.167*Y(i-1)/Y(i)-0.125*(Y(i-1)/Y(i))^2-0.208*(Y(i-1)/Y(i))^3;   %contraction loss coefficient
                end
            end

            % Initial Phase Determination
            for i=1:length(obj.T_sat)-1
                if obj.P_c_i > obj.P_sat(i) && obj.P_c_i < obj.P_sat(i+1)
                    T_int(N_t) = (obj.P_sat(i+1)-obj.P_c_i)/(obj.P_sat(i+1)-obj.P_sat(i))*(obj.T_sat(i)-obj.T_sat(i+1)) + obj.T_sat(i+1);
                    break
                elseif obj.P_c_i > max(obj.P_sat)
                    T_int(N_t) = max(obj.T_sat);
                end
            end

            % 1D Thermal Analysis:
            for t = N_t:-1:1

                %Liquid Phase
                if T_c(t) < T_int(t)
                    vf(t) = 0;
                    if T_c(t) < min(obj.T_l)
                        mu_c = (obj.T_l(2)-T_c(t))/(obj.T_l(2)-obj.T_l(1))*(obj.mu_l(1)-obj.mu_l(2)) + obj.mu_l(2);
                        K_c = (obj.T_l(2)-T_c(t))/(obj.T_l(2)-obj.T_l(1))*(obj.K_l(1)-obj.K_l(2)) + obj.K_l(2);
                        Cp_c = (obj.T_l(2)-T_c(t))/(obj.T_l(2)-obj.T_l(1))*(obj.Cp_l(1)-obj.Cp_l(2)) + obj.Cp_l(2);
                        rho_c = (obj.T_l(2)-T_c(t))/(obj.T_l(2)-obj.T_l(1))*(obj.rho_l(1)-obj.rho_l(2)) + obj.rho_l(2);
                        Pr_c = Cp_c*mu_c/K_c;
                    elseif T_c(t) > max(obj.T_l)
                        mu_c = (obj.T_l(N_l-1)-T_c(t))/(obj.T_l(N_l-1)-obj.T_l(N_l))*(obj.mu_l(N_l)-obj.mu_l(N_l-1)) + obj.mu_l(N_l-1);
                        K_c = (obj.T_l(N_l-1)-T_c(t))/(obj.T_l(N_l-1)-obj.T_l(N_l))*(obj.K_l(N_l)-obj.K_l(N_l-1)) + obj.K_l(N_l-1);
                        Cp_c = (obj.T_l(N_l-1)-T_c(t))/(obj.T_l(N_l-1)-obj.T_l(N_l))*(obj.Cp_l(N_l)-obj.Cp_l(N_l-1)) + obj.Cp_l(N_l-1);
                        rho_c = (obj.T_l(N_l-1)-T_c(t))/(obj.T_l(N_l-1)-obj.T_l(N_l))*(obj.rho_l(N_l)-obj.rho_l(N_l-1)) + obj.rho_l(N_l-1);
                        Pr_c = Cp_c*mu_c/K_c;
                    else
                        for i=1:length(obj.T_l)-1
                            if T_c(t) >= obj.T_l(i) && T_c(t) < obj.T_l(i+1)
                                mu_c = (obj.T_l(i+1)-T_c(t))/(obj.T_l(i+1)-obj.T_l(i))*(obj.mu_l(i)-obj.mu_l(i+1)) + obj.mu_l(i+1);
                                K_c = (obj.T_l(i+1)-T_c(t))/(obj.T_l(i+1)-obj.T_l(i))*(obj.K_l(i)-obj.K_l(i+1)) + obj.K_l(i+1);
                                Cp_c = (obj.T_l(i+1)-T_c(t))/(obj.T_l(i+1)-obj.T_l(i))*(obj.Cp_l(i)-obj.Cp_l(i+1)) + obj.Cp_l(i+1);
                                rho_c = (obj.T_l(i+1)-T_c(t))/(obj.T_l(i+1)-obj.T_l(i))*(obj.rho_l(i)-obj.rho_l(i+1)) + obj.rho_l(i+1);
                                Pr_c = Cp_c*mu_c/K_c;
                                break
                            end
                        end
                    end

                    %Saturation Phase
                elseif (T_c(t) == T_int(t))
                    for i=1:length(obj.h_f)-1
                        if T_c(t) > obj.T_sat(i) && T_c(t) < obj.T_sat(i+1)
                            %Interpolated Properties
                            mu_f_int = (obj.T_sat(i+1)-T_c(t))/(obj.T_sat(i+1)-obj.T_sat(i))*(obj.mu_f(i)-obj.mu_f(i+1)) + obj.mu_f(i+1);
                            mu_g_int = (obj.T_sat(i+1)-T_c(t))/(obj.T_sat(i+1)-obj.T_sat(i))*(obj.mu_g(i)-obj.mu_g(i+1)) + obj.mu_g(i+1);
                            Cp_f_int = (obj.T_sat(i+1)-T_c(t))/(obj.T_sat(i+1)-obj.T_sat(i))*(obj.Cp_f(i)-obj.Cp_f(i+1)) + obj.Cp_f(i+1);
                            Cp_g_int = (obj.T_sat(i+1)-T_c(t))/(obj.T_sat(i+1)-obj.T_sat(i))*(obj.Cp_g(i)-obj.Cp_g(i+1)) + obj.Cp_g(i+1);
                            K_f_int = (obj.T_sat(i+1)-T_c(t))/(obj.T_sat(i+1)-obj.T_sat(i))*(obj.K_f(i)-obj.K_f(i+1)) + obj.K_f(i+1);
                            K_g_int = (obj.T_sat(i+1)-T_c(t))/(obj.T_sat(i+1)-obj.T_sat(i))*(obj.K_g(i)-obj.K_g(i+1)) + obj.K_g(i+1);
                            rho_f_int = (obj.T_sat(i+1)-T_c(t))/(obj.T_sat(i+1)-obj.T_sat(i))*(obj.rho_f(i)-obj.rho_f(i+1)) + obj.rho_f(i+1);
                            rho_g_int = (obj.T_sat(i+1)-T_c(t))/(obj.T_sat(i+1)-obj.T_sat(i))*(obj.rho_g(i)-obj.rho_g(i+1)) + obj.rho_g(i+1);
                            %Coolant Properties
                            vf(t) = 1-(h_g_int-h_c(t))/(h_g_int-h_f_int);
                            mu_c = vf(t)*(mu_f_int-mu_g_int) + mu_f_int;
                            Cp_c = vf(t)*(Cp_f_int-Cp_g_int) + Cp_f_int;
                            K_c = vf(t)*(K_f_int-K_g_int) + K_f_int;
                            rho_c = vf(t)*(rho_f_int-rho_g_int) + rho_f_int;
                            Pr_c = Cp_c*mu_c/K_c;
                            break
                        end
                    end

                    %Vapor Phase
                elseif T_c(t) > T_int(t)
                    vf(t) = 1;
                    %Property arrays interpolated/extrapolated for pressure
                    if P_c(t) < min(P_v)
                        T_v_int = (P_v(2)-P_c(t))/(P_v(2)-P_v(1))*(obj.T_v(:,1)-obj.T_v(:,2)) + obj.T_v(:,2);
                        mu_v_int = (P_v(2)-P_c(t))/(P_v(2)-P_v(1))*(obj.mu_v(:,1)-obj.mu_v(:,2)) + obj.mu_v(:,2);
                        K_v_int = (P_v(2)-P_c(t))/(P_v(2)-P_v(1))*(obj.K_v(:,1)-obj.K_v(:,2)) + obj.K_v(:,2);
                        Cp_v_int = (P_v(2)-P_c(t))/(P_v(2)-P_v(1))*(obj.Cp_v(:,1)-obj.Cp_v(:,2)) + obj.Cp_v(:,2);
                        rho_v_int = (P_v(2)-P_c(t))/(P_v(2)-P_v(1))*(obj.rho_v(:,1)-obj.rho_v(:,2)) + obj.rho_v(:,2);
                    elseif P_c(t) > max(P_v)
                        T_v_int = (P_v(obj.N_p-1)-P_c(t))/(P_v(obj.N_p-1)-P_v(obj.N_p))*(obj.T_v(:,obj.N_p)-obj.T_v(:,obj.N_p-1)) + obj.T_v(:,obj.N_p-1);
                        mu_v_int = (P_v(obj.N_p-1)-P_c(t))/(P_v(obj.N_p-1)-P_v(obj.N_p))*(obj.mu_v(:,obj.N_p)-obj.mu_v(:,obj.N_p-1)) + obj.mu_v(:,obj.N_p-1);
                        K_v_int = (P_v(obj.N_p-1)-P_c(t))/(P_v(obj.N_p-1)-P_v(obj.N_p))*(obj.K_v(:,obj.N_p)-obj.K_v(:,obj.N_p-1)) + obj.K_v(:,obj.N_p-1);
                        Cp_v_int = (P_v(obj.N_p-1)-P_c(t))/(P_v(obj.N_p-1)-P_v(obj.N_p))*(obj.Cp_v(:,obj.N_p)-obj.Cp_v(:,obj.N_p-1)) + obj.Cp_v(:,obj.N_p-1);
                        rho_v_int = (P_v(obj.N_p-1)-P_c(t))/(P_v(obj.N_p-1)-P_v(obj.N_p))*(obj.rho_v(:,obj.N_p)-obj.rho_v(:,obj.N_p-1)) + obj.rho_v(:,obj.N_p-1);
                    else
                        for i=1:length(P_v)-1
                            if P_c(t) > P_v(i) && P_c(t) < P_v(i+1)
                                T_v_int = (P_v(i+1)-P_c(t))/(P_v(i+1)-P_v(i))*(obj.T_v(:,i)-obj.T_v(:,i+1)) + obj.T_v(:,i+1);
                                mu_v_int = (P_v(i+1)-P_c(t))/(P_v(i+1)-P_v(i))*(obj.mu_v(:,i)-obj.mu_v(:,i+1)) + obj.mu_v(:,i+1);
                                K_v_int = (P_v(i+1)-P_c(t))/(P_v(i+1)-P_v(i))*(obj.K_v(:,i)-obj.K_v(:,i+1)) + obj.K_v(:,i+1);
                                Cp_v_int = (P_v(i+1)-P_c(t))/(P_v(i+1)-P_v(i))*(obj.Cp_v(:,i)-obj.Cp_v(:,i+1)) + obj.Cp_v(:,i+1);
                                rho_v_int = (P_v(i+1)-P_c(t))/(P_v(i+1)-P_v(i))*(obj.rho_v(:,i)-obj.rho_v(:,i+1)) + obj.rho_v(:,i+1);
                                break
                            end
                        end
                    end
                    %Property values interpolated/extrapolated for temperature
                    if T_c(t) < min(T_v_int)
                        mu_c = (T_v_int(2)-T_c(t))/(T_v_int(2)-T_v_int(1))*(mu_v_int(1)-mu_v_int(2)) + mu_v_int(2);
                        K_c = (T_v_int(2)-T_c(t))/(T_v_int(2)-T_v_int(1))*(K_v_int(1)-K_v_int(2)) + K_v_int(2);
                        Cp_c = (T_v_int(2)-T_c(t))/(T_v_int(2)-T_v_int(1))*(Cp_v_int(1)-Cp_v_int(2)) + Cp_v_int(2);
                        rho_c = (T_v_int(2)-T_c(t))/(T_v_int(2)-T_v_int(1))*(rho_v_int(1)-rho_v_int(2)) + rho_v_int(2);
                        Pr_c = Cp_c*mu_c/K_c;
                    elseif T_c(t) > max(T_v_int)
                        mu_c = (T_v_int(N_v-1)-T_c(t))/(T_v_int(N_v-1)-T_v_int(N_v))*(mu_v_int(N_v)-mu_v_int(N_v-1)) + mu_v_int(N_v-1);
                        K_c = (T_v_int(N_v-1)-T_c(t))/(T_v_int(N_v-1)-T_v_int(N_v))*(K_v_int(N_v)-K_v_int(N_v-1)) + K_v_int(N_v-1);
                        Cp_c = (T_v_int(N_v-1)-T_c(t))/(T_v_int(N_v-1)-T_v_int(N_v))*(Cp_v_int(N_v)-Cp_v_int(N_v-1)) + Cp_v_int(N_v-1);
                        rho_c = (T_v_int(N_v-1)-T_c(t))/(T_v_int(N_v-1)-T_v_int(N_v))*(rho_v_int(N_v)-rho_v_int(N_v-1)) + rho_v_int(N_v-1);
                        Pr_c = Cp_c*mu_c/K_c;
                    else
                        for i=1:length(T_v_int)
                            if T_c(t) > T_v_int(i) && T_c(t) < T_v_int(i+1)
                                mu_c = (T_v_int(i+1)-T_c(t))/(T_v_int(i+1)-T_v_int(i))*(mu_v_int(i)-mu_v_int(i+1)) + mu_v_int(i+1);
                                K_c = (T_v_int(i+1)-T_c(t))/(T_v_int(i+1)-T_v_int(i))*(K_v_int(i)-K_v_int(i+1)) + K_v_int(i+1);
                                Cp_c = (T_v_int(i+1)-T_c(t))/(T_v_int(i+1)-T_v_int(i))*(Cp_v_int(i)-Cp_v_int(i+1)) + Cp_v_int(i+1);
                                rho_c = (T_v_int(i+1)-T_c(t))/(T_v_int(i+1)-T_v_int(i))*(rho_v_int(i)-rho_v_int(i+1)) + rho_v_int(i+1);
                                Pr_c = Cp_c*mu_c/K_c;
                                break
                            end
                        end
                    end
                end

                %Geometrical Properties
                d_i = 2*Y(t);                       %inner diameter
                d_o = d_i + 2*y(t);                 %outer diameter

                if isempty(inputs.w_b) == 1
                    w_b(t) = (pi*d_o - N*w(t))/N;   %fin base width
                end

                A_c = h(t)*w(t);                    %channel area
                P = 2*(h(t)+w(t));                  %channel perimeter
                d_h = 4*A_c/P;                      %hydraulic diameter

                %Compressible Flow
                T_0 = T_t*(1+(k-1)/2);                                                  %gas stagnation temperature
                T_g = T_0/(1+(k-1)/2*M(t)^2);                                           %gas free stream temperature
                T_aw = T_0*(1+(Pr_hg(T_g))^(1/3)*(k-1)/2*M(t)^2)/(1+(k-1)/2*M(t)^2);    %adiabatic wall temperature
               
                %Wall Temperature Iteration
                T_w_new = 410;
                Tw = 400;

                while abs(T_w_new-Tw) > 0.1
                    Tw = T_w_new;

                    %Bartz Correlation
                    sigma = (0.5*Tw/T_0*(1+(k-1)/2*(M(t))^2)+0.5)^-0.68*(1+(k-1)/2*(M(t))^2)^-0.12;
                    h_hg = 0.026/(D_t^0.2)*(mu_hg(T_0))^0.2*Cp_hg(T_0)/(Pr_hg(T_0))^0.6*(P_cc/obj.C_star)^0.8*(D_t^2/d_i^2)^0.9*sigma;

                    %Heat Transfer Coefficients
                    Re_c = 4*obj.mdot_c/(pi*d_h*mu_c);              %coolant Reynolds number
                    Nu_c = 0.023*Re_c^0.8*Pr_c^0.4;                 %coolant Nusselt number
                    u_c = Nu_c*K_c/d_h;                             %coolant convection coefficient
                    m = sqrt(2*u_c*w_b(t)/(k_wall));                %fin parameter
                    e_f = tanh(m*h(t)/w_b(t))/(m*h(t)/w_b(t));      %fin efficiency
                    u_c_f = u_c*(w(t)+2*e_f*h(t))/(w(t)+w_b(t));    %corrected coolant convection coefficient

                    q_w(t) = (T_aw-T_c(t))/(1/u_c_f + 1/h_hg + y(t)/k_wall);
                    T_w_new = T_aw-q_w(t)/h_hg;
                end

                T_w_i(t) = T_w_new;                                 %wall inner temperature
                Q = q_w(t)*pi*d_i*L(t);                             %heat rate

                %Haaland Formula
                f = 1/(1.9*log10(((obj.eps/d_h)^3.7)^1.11 + 6.9/Re_c))^2;                  %Darcy friction coefficient

                %Coolant State Determination
                if t>1
                    P_c(t-1) = P_c(t) - obj.mdot_c^2/(2*rho_c*A_c^2)*(f*(L(t))/d_h+K(t));  %coolant pressure
                    T_c(t-1) = T_c(t) + Q/(obj.mdot_c*N*Cp_c);                             %coolant temperature

                    %Phase Determination
                    for i=1:length(obj.T_sat)-1
                        if P_c(t-1) > obj.P_sat(i) && P_c(t-1) < obj.P_sat(i+1)
                            T_int(t-1) = (obj.P_sat(i+1)-P_c(t-1))/(obj.P_sat(i+1)-obj.P_sat(i))*(obj.T_sat(i)-obj.T_sat(i+1)) + obj.T_sat(i+1);
                            break
                        elseif P_c(t-1) > max(obj.P_sat)
                            T_int(t-1) = max(obj.T_sat);
                        end
                    end

                    %Phase Change Beginning
                    if T_c(t-1)>T_int(t-1) && T_c(t)<T_int(t)
                        for i=1:length(obj.h_f)-1
                            if T_c(t-1) > obj.T_sat(i) && T_c(t-1) < obj.T_sat(i+1)
                                %Enthalpy Interpolation
                                h_f_int = (obj.T_sat(i+1)-T_c(t-1))/(obj.T_sat(i+1)-obj.T_sat(i))*(obj.h_f(i)-obj.h_f(i+1)) + obj.h_f(i+1);
                                h_g_int = (obj.T_sat(i+1)-T_c(t-1))/(obj.T_sat(i+1)-obj.T_sat(i))*(obj.h_g(i)-obj.h_g(i+1)) + obj.h_g(i+1);
                                h_c(t-1) = h_f_int;      % There is an approximation here: +(T_int(t-1)-T_c(t-1))*mdot_c*Cp_c
                                T_c(t-1) = T_int(t-1);
                                break
                            end
                        end
                    end

                    %Phase Changing Coolant
                    if T_c(t) == T_int(t)
                        h_c(t-1) = h_c(t) + Q/(obj.mdot_c*N);
                        T_c(t-1) = T_int(t-1);
                        %Phase Change Ending
                        if h_c(t-1) >= h_g_int
                            T_c(t-1) = T_c(t) + Q/(obj.mdot_c*N*Cp_c);
                        end
                    end
                end
                C_c(t) = obj.mdot_c*N*Cp_c;
                u(t) = u_c_f;
            end
            obj.Tem_c = T_c;
            obj.p_c = P_c;
            obj.q = q_w;
            obj.T_w = T_w_i;
            obj.Ma = M;
            obj.X_g = vf;
            obj.Cdot = C_c;
            obj.htc = u;
            obj.X = x;
            obj.Z = Y;
        end
    end
    % Thrust chamber profile generation
    methods  (Access = private)
        function [X,Y] = Thrust_Chamber_Profile(obj)
            %Generates a thrust chamber profile.
            %User may not use this function and give the thrust chamber
            %profile externally.

            %Preprocessing
            R_t = obj.d_t/2;                                %throat radius
            R_c = (1.5*R_t*sind(-135) + 2.5*R_t)/cosd(45);  %combustion chamber radius

            nozzle_exp = obj.exp_rat;                       %nozzle expansion ratio
            L_c = obj.L_star*R_t^2/R_c^2;                   %length of combustion chamber
            Theta_n = obj.theta_n;                          %nozzle inlet angle in deg.
            Theta_e = obj.theta_e;                          %nozzle exit angle in deg.
            N1 = obj.discretization(1);                     %cylindrical part number of divisions
            N2 = obj.discretization(2);                     %cylindrical part number of divisions
            N3 = obj.discretization(3);                     %cylindrical part number of divisions
            N4 = obj.discretization(4);                     %cylindrical part number of divisions
            N5 = obj.discretization(5);                     %cylindrical part number of divisions

            a = 1.5*R_t*cosd(-135) - (R_c*cosd(45) - R_t*(1.5*sind(-135)+2.5));

            %Initialization
            x1 = zeros(1,N1);
            y1 = zeros(1,N1);
            x2 = zeros(1,N2);
            y2 = zeros(1,N2);
            x3 = zeros(1,N3);
            y3 = zeros(1,N3);
            x4 = zeros(1,N4);
            y4 = zeros(1,N4);
            x5 = zeros(1,N5);
            y5 = zeros(1,N5);

            %Cylinder Combustion Chamber Part:
            for i=1:N1
                x1(i) = -L_c*(N1-i)/(N1-1) + a - R_c*sind(45);
                y1(i) = R_c;
            end

            %Round Conical Section Part:
            for i=1:N2
                x2(i) = a - R_c*sind(45) + R_c*sind((i-1)/(N2-1)*45);
                y2(i) = R_c*cosd((i-1)/(N2-1)*45);
            end

            %Throat Entrance Part:
            for i=1:N3
                x3(i) = 1.5*R_t*cosd(-135+(i-1)/(N3-1)*45);
                y3(i) = 1.5*R_t*sind(-135+(i-1)/(N3-1)*45) + 2.5*R_t;
            end

            %Throat Exit Part:
            for i=1:N4
                x4(i) = 0.382*R_t*cosd(-90+(i-1)/(N4-1)*Theta_n);
                y4(i) = 0.382*R_t*sind(-90+(i-1)/(N4-1)*Theta_n) + 1.382*R_t;
            end

            %Bell Curve Part:
            N_x = 0.382*R_t*cosd(Theta_n-90);
            N_y = 0.382*R_t*sind(Theta_n-90) + 1.382*R_t;
            E_x = 0.8*(((sqrt(nozzle_exp)-1)-1)*R_t)/tand(15);
            E_y = sqrt(nozzle_exp)*R_t;
            m1 = tand(Theta_n);
            m2 = tand(Theta_e);
            C1 = N_y-m1*N_x;
            C2 = E_y-m2*E_x;
            Q_x = (C2-C1)/(m1-m2);
            Q_y = (m1*C2-m2*C1)/(m1-m2);

            for i=1:N5
                t = (i-1)/(N5-1);
                x5(i) = (1-t)^2*N_x + 2*(1-t)*t*Q_x + t^2*E_x;
                y5(i) = (1-t)^2*N_y + 2*(1-t)*t*Q_y + t^2*E_y;
            end

            %Arrays:
            X = [x1 x2 x3 x4 x5];
            Y = [y1 y2 y3 y4 y5];
        end
    end
end