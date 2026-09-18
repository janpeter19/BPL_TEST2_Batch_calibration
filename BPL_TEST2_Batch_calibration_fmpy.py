# Setup data TEST2_Batch_calibration_fmpy 
# Author: Jan Peter Axelsson
#------------------------------------------------------------------------------------------------------------------
# 2026-09-10 - Created
# 2026-09-18 - Decrease the framework to what is necessary and move matlotlib to the other setup-file
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
#  Framework
#------------------------------------------------------------------------------------------------------------------

# Setup framework
import platform
import locale
from fmpy import simulate_fmu
from fmpy import read_model_description

# Set the environment - for Linux a JSON-file in the FMU is read
if platform.system() == 'Linux': locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

#------------------------------------------------------------------------------------------------------------------
#  Setup application FMU
#------------------------------------------------------------------------------------------------------------------

# Provde the right FMU and load for different platforms in user dialogue:
if platform.system() == 'Windows':
   print('Windows - run FMU pre-compiled JModelica 2.14')
   fmu_model ='BPL_TEST2_Batch_windows_jm_cs.fmu'     
   flag_vendor = 'JM' 
   flag_type = 'CS'
elif platform.system() == 'Linux':
   print('Linux - run FMU pre-compiled OpenModelica')
   fmu_model ='BPL_TEST2_Batch_linux_om_me.fmu'  
#  fmu_model ='BPL_TEST2_Batch_linux_2404_om_me.fmu'  
   flag_vendor = 'OM' 
   flag_type = 'ME'
else:    
   print('There is no FMU for this platform')

# Provide various opts-profiles
if flag_type in ['CS', 'cs']:
   opts_std = {'NCP': 500}
   opts_data = {'NCP': 12}
   opts_fast = {'NCP': 12}
elif flag_type in ['ME', 'me']:
   opts_std = {'NCP': 500}
   opts_data = {'NCP': 12}
   opts_fast = {'NCP': 12}
else:    
   print('There is no FMU for this platform')

# Extract model_description from fmu_model
model_description = read_model_description(fmu_model)

# Provide various MSL and BPL versions
if flag_vendor in ['JM', 'jm']:
#   MSL_usage = model.get('MSL.usage')[0]
   constants = [v for v in model_description.modelVariables if v.causality == 'local'] 
   MSL_usage = [x[1] for x in [(constants[k].name, constants[k].start) for k in range(len(constants))] if 'MSL.usage' in x[0]][0]   
   MSL_version = [x[1] for x in [(constants[k].name, constants[k].start) for k in range(len(constants))] if 'MSL.version' in x[0]][0]
   BPL_version = [x[1] for x in [(constants[k].name, constants[k].start) for k in range(len(constants))] if 'BPL.version' in x[0]][0] 
elif flag_vendor in ['OM', 'om']:
   MSL_usage = '4.1.0 - used components: none' 
   MSL_version = '4.1.0'
   BPL_version = 'Bioprocess Library version 2.3.2' 
else:    
   print('There is no FMU for this platform')
   
#------------------------------------------------------------------------------------------------------------------

# Simulation time
simulationTime = 5.0

# Dictionary of time discrete states
timeDiscreteStates = {} 

# Define a minimal compoent list of the model as a starting point for describe('parts')
component_list_minimum = ['bioreactor', 'bioreactor.culture']

# Provide process diagram on disk
fmu_process_diagram ='BPL_TEST2_Batch_process_diagram_om.png'

#------------------------------------------------------------------------------------------------------------------
#  Specific application constructs: stateValue, parValue, parLocation, parCheck,parValue diagrams, ax, lines
#------------------------------------------------------------------------------------------------------------------

# Create dictionaries parValue[] and parLocation[]
parValue = {}
parValue['V_start'] = 1.0
parValue['VX_start'] = 1.0
parValue['VS_start'] = 10.0

parValue['Y'] = 0.5
parValue['qSmax'] = 1.0
parValue['Ks'] = 0.1

parLocation = {}
parLocation['V_start'] = 'bioreactor.V_start'
parLocation['VX_start'] = 'bioreactor.m_start[1]' 
parLocation['VS_start'] = 'bioreactor.m_start[2]' 

parLocation['Y'] = 'bioreactor.culture.Y'
parLocation['qSmax'] = 'bioreactor.culture.qSmax'
parLocation['Ks'] = 'bioreactor.culture.Ks'

# Extra only for describe()
keyVariables = []
parLocation['mu'] = 'bioreactor.culture.mu'; keyVariables.append(parLocation['mu'])

# Parameter value check - especially for hysteresis to avoid runtime error
parCheck = []
parCheck.append("parValue['Y'] > 0")
parCheck.append("parValue['qSmax'] > 0")
parCheck.append("parValue['Ks'] > 0")
parCheck.append("parValue['V_start'] > 0")
parCheck.append("parValue['VX_start'] >= 0")
parCheck.append("parValue['VS_start'] >= 0")

# Create list of diagrams to be plotted by simu()
diagrams = []

# Create an empty list axes to be defined in newplot() and plotted by simu() or show()
ax = []

# Create list of pens for the diagrams
lines = ['-','--',':','-.']
