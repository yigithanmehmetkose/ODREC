# ODREC: *O*ne *D*imensional *RE*generative *C*ooling
ODREC is a program for the preliminary design and 1-D regenerative cooling analysis of liquid propellant rocket engines, which is developed using the object-oriented programming in MATLAB.

Please refer to the study (https://doi.org/10.3390/app14010071) for the details on the method/formulation and capabilities of ODREC.

There are two prerequisities to use ODREC:
* DimVar tool (https://github.com/tgvoskuilen/MatlabTools)
* MATLAB version of NASA CEA (https://github.com/PurdueH2Lab/MatlabCEA)

To use ODREC, an instance of the class should be created:
> Analysis = ODREC(inputs,l_prop,s_prop,v_prop);

where *inputs* is a struct with the following fields:
* 
