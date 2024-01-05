# ODREC: *O*ne *D*imensional *RE*generative *C*ooling
ODREC is a program for the preliminary design and 1-D regenerative cooling analysis of liquid propellant rocket engines, which is developed using the object-oriented programming in MATLAB.

Please refer to the study (https://doi.org/10.3390/app14010071) for the details on the method/formulation and capabilities of ODREC.

There are two prerequisities to use ODREC:
* DimVar tool (https://github.com/tgvoskuilen/MatlabTools)
* MATLAB version of NASA CEA (https://github.com/PurdueH2Lab/MatlabCEA)

To use ODREC, an instance of the class should be created:
> Analysis = ODREC(inputs,l_prop,s_prop,v_prop);

where *inputs* is a struct with the following fields:
* *fuel*: the name of the fuel
* *oxidizer*: the name of the oxidizer
* *T_o_c*: temperature of the oxidizer entering into the combustion chamber in K
* *T_f_c*: temperature of the fuel entering into the combustion chamber in K
* *discretization*: an array of 5 elements, each of which denotes the number of sections each thrust chamber part is divided into. The parts are cylindrical combustion chamber part, circularly converging combustion chamber part, small converging nozzle part, fist diverging nozzle part, bell shaped nozzle part. If the user externally gives the thrust chamber profile as an input rather than generating within the ODREC, then *discretization* should be a scalar denoting the total number of discretization of the thrust chamber, rather than an array.
* *eps*: roughness of the inner thrust chamber wall in m
* *i_c*: index of the coolant, 1: fuel, 2: oxidizer
* *y*: thrust chamber wall thickness at channel base in m
* *L_star*: characteristic length of the combustion chamber in m
* *theta_n*: deflection angle of bell nozzle in deg
* *theta_e*: exit angle of the bell nozzle in deg
* *P_c_i*: inlet pressure of the coolant in Pa
* *T_c_i*: inlet temperature of the coolant in K
* *k_wall*: wall thermal conductivity in W/mK
* *P_1*: pressure in first table of vapor properties in Pa
* *P_l*: pressure in last table of vapor properties in Pa
* *N_p*: number of tables of vapor properties
* *P_c*: combustion chamber pressure in Pa
* *F*: thrust in N
* *OF*: oxidizer to fuel mass flow rate ratio
* *Pe*: ambient pressure in Pa
* *N*: number of cooling channels
* *w*: width of the cooling channels
* *h*: height of the cooling channels

Some remarks about *inputs*:
* Only those fuel and oxidizer names that are accepted in NASA CEA website (https://cearun.grc.nasa.gov/) can be used.
* Only those *T_o_c* and *T_f_c* values that are accepted in NASA CEA website can be used. For most propellants, these values are the saturation temperatures at the sea level. Note that *T_o_c* and *T_f_c* values do not change the analysis much.
* *y*,*w*, and *h* can be given as constants or variables throughout the thrust chamber. If they are variables, they can just be given as an array thoughout the axial position.
* If desired, an additional input field of *w_b*, which is the rib or channel wall width, can be provided as a constant or a variable. If not given, ODREC will calculate it using the other inputs.
* ODREC generates a thrust chamber profile if it is not given externally. If desired, the user should provide it as follows:
> Analysis = ODREC(inputs,l_prop,s_prop,v_prop,geom);

where *geom* is a matrix storing the axial position versus radial position of the thrust chamber profile:
> geom(:,1)=X;

> geom(:,2)=Y;

where X is the axial, Y is the radial position.

Some remarks about the coolant properties:















