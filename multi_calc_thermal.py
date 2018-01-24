# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:46:18 2017

@author: SunJingyu
"""

from Constants import *
from Data import *

#for name in ['7201','7209']:
for name in experiment_type_dict['TC-IP-TGMF']+experiment_type_dict['TC-OP-TGMF']+experiment_type_dict['TC-IP-TGMF-TBC']:
    filename = 'umat_cal_in718_%s.cal' % name
    outfile = open(filename, 'w')
    print >>outfile,"""
# -*- coding: utf-8 -*-
from Constants import *
from Data import *
from workbench_convection import workbench
name = '%s'
experiment_log = ExperimentLog(ExperimentLogFile)

air_temperature = 20.0 + 273.15
pressure = 6.5e5
dynamic_viscosity_0 = 1.7894e-5
rin = 0.0065/2.0
d = rin*2

experiment_log.output(name)
regular = r'.*'
load_type = experiment_log.obtainItem(name,'load_type',regular)[0]
regular = r'\d+\.?\d*'
temperature_mode = experiment_log.obtainItem(name,'temperature_mode',regular)
d_out = float(experiment_log.obtainItem(name,'d_out',regular)[0])
gauge_length = float(experiment_log.obtainItem(name,'gauge_length',regular)[0])
axial_strain = float(experiment_log.obtainItem(name,'axial_strain',regular)[0])
angel_strain = float(experiment_log.obtainItem(name,'angel_strain',regular)[0])
equivalent_strain = float(experiment_log.obtainItem(name,'equivalent_strain',regular)[0])
period = float(experiment_log.obtainItem(name,'period',regular)[0])
axial_temperature_phase = float(experiment_log.obtainItem(name,'axial_temperature_phase',regular)[0])
life = int(experiment_log.obtainItem(name,'comments',regular)[0])

for volume_flow in [50]:
    velocity = volume_flow/60.0/1000.0/(np.pi*rin**2)
    density = pressure/287.058/air_temperature
    thermal_conductivity = -0.00037+0.000103*air_temperature-4.657E-08*air_temperature**2
    specific_heat = 1070.3-0.564*air_temperature+0.001507*air_temperature**2-0.000001102*air_temperature**3-0.000000014*air_temperature**4
    dynamic_viscosity = dynamic_viscosity_0 * (air_temperature/288.15)**(1.5)*(288.15+110.4)/(air_temperature+110.4)
    kinematic_viscosity = dynamic_viscosity/density
    Re = density * velocity * d / dynamic_viscosity
    Pr = specific_heat*dynamic_viscosity/thermal_conductivity
    #Gnielinski correlation
    f = (0.79*np.log(Re)-1.64)**(-2.0)
    Nu = (f/8.0)*(Re-1000.0)*Pr/(1.0+12.7*(f/8.0)**0.5*(Pr**(2.0/3.0)-1))
    #Dittus-Boelter equation
    #Nu = 0.023*Re**(4.0/5.0)*Pr**0.4       

    film_coefficient = thermal_conductivity/d*Nu/1000.0
    sink_temperature = air_temperature - 273.15
    temperature_list = []
    
    workbench(name,loading_cycles=None,copy=True,
              film_coefficient=film_coefficient,
              sink_temperature=sink_temperature,
              temperature_list=temperature_list)
""" % name
    outfile.close()

outfile = open('multi_calc.bat', 'w')
#for name in ['7201','7209']:
for name in experiment_type_dict['TC-IP-TGMF']+experiment_type_dict['TC-OP-TGMF']+experiment_type_dict['TC-IP-TGMF-TBC']:
    print >>outfile,'start python umat_cal_in718_%s.cal' % name
outfile.close()