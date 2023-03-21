# Figure - Simulation of batch reactor 
#          with functions added to facilitate explorative simulation work
#          Here we use FMPy for simulation of the FMU
#
# Can start with: python -m fmpy.gui
#
# GNU General Public License v3.0
# Copyright (c) 2022, Jan Peter Axelsson, All rights reserved.
#------------------------------------------------------------------------------------------------------------------
# 2022-10-04 - First try to setup FMU-explore with FMPy for BPL_TEST2_Batch - first try works! (started with 0.9.1)
# 2022-10-05 - Import read_model_description and make some use of that
# 2022-10-11 - Include other variables like mu as outputs from simulation and system_info() more complete
# 2022-10-12 - Improved simu() with automatic extraction of variables from diagrams
# 2022-10-13 - Improved disp() to be as usual for pyfmi and also system_info() completed with FMU type
# 2022-10-13 - Improved the function describe() and now only 'parts' not implemented, also 'mu' does not work
# 2022-10-14 - Corrected describe() by separate handling of parameters and variables
# 2022-10-14 - Included in describe() now also 'liquidphase' and 'parts'
# 2022-10-15 - Improved simu() with argument output_interval special for fmpy and related to pyfmi ncp
# 2023-02-27 - Update FMU-explore to 0.9.6 in one leap and added list key_variables for logging
# 2023-02-28 - Place the list key_variables better
# 2023-03-09 - Prepare for Google Colab use
# 2023-03-21 - Clean-up and use standard FMU notation
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
#  Framework
#------------------------------------------------------------------------------------------------------------------

# Setup framework
import sys
import platform
import locale
import numpy as np 
import matplotlib.pyplot as plt 

from fmpy import simulate_fmu
from fmpy import read_model_description
import fmpy as fmpy

#from pyfmi import load_fmu
#from pyfmi.fmi import FMUException

from itertools import cycle
from importlib_metadata import version   # included in future Python 3.8

# Set the environment - for Linux a JSON-file in the FMU is read
if platform.system() == 'Linux': locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

#------------------------------------------------------------------------------------------------------------------
#  Setup application FMU
#------------------------------------------------------------------------------------------------------------------

# Provde the right FMU and load for different platforms in user dialogue:
global fmu_model, model, opts
if platform.system() == 'Windows':
   print('Windows - run FMU pre-compiled JModelica 2.14')
   fmu_model ='BPL_TEST2_Batch_windows_jm_cs.fmu'
   model_description = read_model_description(fmu_model)        
   flag_vendor = 'JM' 
   flag_type = 'CS'
elif platform.system() == 'Linux':
   print('Linux - run FMU pre-compiled OpenModelica 1.21.x')
   fmu_model ='BPL_TEST2_Batch_linux_om_me.fmu'  
   model_description = read_model_description(fmu_model)  
   flag_vendor = 'OM' 
   flag_type = 'ME'
else:    
   print('There is no FMU for this platform')

# Provide various opts-profiles
#if flag_type in ['CS', 'cs']:
#   opts_std = model.simulate_options()
#   opts_std['silent_mode'] = True
#   opts_std['ncp'] = 500 
#   opts_std['result_handling'] = 'binary'     
#elif flag_type in ['ME', 'me']:
#   opts_std = model.simulate_options()
#   opts_std["CVode_options"]["verbosity"] = 50 
#   opts_std['ncp'] = 500 
#   opts_std['result_handling'] = 'binary'  
#else:    
#   print('There is no FMU for this platform')

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
   MSL_usage = '3.2.3 - used components: none' 
   MSL_version = '3.2.3'
   BPL_version = 'Bioprocess Library version 2.1.1-beta' 
else:    
   print('There is no FMU for this platform')

# Simulation time
global simulationTime; simulationTime = 5.0

# Dictionary of time discrete states
timeDiscreteStates = {} 

# Define a minimal compoent list of the model as a starting point for describe('parts')
component_list_minimum = ['bioreactor', 'bioreactor.culture']

#------------------------------------------------------------------------------------------------------------------
#  Specific application constructs: stateDict, parDict, diagrams, newplot(), describe()
#------------------------------------------------------------------------------------------------------------------

# Create stateDict that later will be used to store final state and used for initialization in 'cont':
#stateDict = model.get_states_list()
global stateDict

# Create dictionaries parDict[] and parLocation[]
global parDict; parDict = {}
parDict['V_0'] = 1.0
parDict['VX_0'] = 1.0
parDict['VS_0'] = 10.0

parDict['Y'] = 0.5
parDict['qSmax'] = 1.0
parDict['Ks'] = 0.1

global parLocation; parLocation = {}
parLocation['V_0'] = 'bioreactor.V_0'
parLocation['VX_0'] = 'bioreactor.m_0[1]' 
parLocation['VS_0'] = 'bioreactor.m_0[2]' 

parLocation['Y'] = 'bioreactor.culture.Y'
parLocation['qSmax'] = 'bioreactor.culture.qSmax'
parLocation['Ks'] = 'bioreactor.culture.Ks'

# Extra only for describe()
global key_variables; key_variables = []
parLocation['mu'] = 'bioreactor.culture.mu'; key_variables.append(parLocation['mu'])

# Parameter value check - especially for hysteresis to avoid runtime error
global parCheck; parCheck = []
parCheck.append("parDict['Y'] > 0")
parCheck.append("parDict['qSmax'] > 0")
parCheck.append("parDict['Ks'] > 0")
parCheck.append("parDict['V_0'] > 0")
parCheck.append("parDict['VX_0'] >= 0")
parCheck.append("parDict['VS_0'] >= 0")

# Create list of diagrams to be plotted by simu()
global diagrams
diagrams = []

# Define standard diagrams
def newplot(title='Batch cultivation', plotType='TimeSeries'):
   """ Standard plot window
        title = ''
       two possible diagrams
        diagram = 'TimeSeries' default
        diagram = 'PhasePlane' """
    
   # Reset pens
   setLines()
   
   # Transfer of global axes to simu()      
   global ax, ax1, ax2    
    
   # Plot diagram 
   if plotType == 'TimeSeries':

      plt.figure()
      ax1 = plt.subplot(2,1,1)
      ax2 = plt.subplot(2,1,2)
    
      ax1.set_title(title)
      ax1.grid()
      ax1.set_ylabel('X and S [g/L]')
           
      ax2.grid()
      ax2.set_ylabel('mu [1/h]')
      ax2.set_xlabel('Time [h]') 
      
      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax1.plot(sim_res['time'],sim_res['bioreactor.c[1]'],color='b',linestyle=linetype)")
      diagrams.append("ax2.plot(sim_res['time'],sim_res['bioreactor.c[2]'],color='b',linestyle=linetype)")   

   elif plotType == 'TimeSeries2':

      plt.figure()
      ax1 = plt.subplot(2,1,1)
      ax2 = plt.subplot(2,1,2)
    
      ax1.set_title(title)
      ax1.grid()
      ax1.set_ylabel('X and S [g/L]')
           
      ax2.grid()
      ax2.set_ylabel('mu [1/h]')
      ax2.set_xlabel('Time [h]') 
      
      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax1.plot(sim_res['time'],sim_res['bioreactor.c[1]'],color='r',linestyle=linetype)")
      diagrams.append("ax1.plot(sim_res['time'],sim_res['bioreactor.c[2]'],color='b',linestyle=linetype)")   
      diagrams.append("ax2.plot(sim_res['time'],sim_res['bioreactor.culture.mu'],color='r',linestyle=linetype)") 

   elif plotType == 'Demo_1':
   
      plt.figure()
      ax1 = plt.subplot(2,1,1)
      ax2 = plt.subplot(2,1,2)

      ax1.set_title(title)
      ax1.grid()
      ax1.set_ylabel('S [g/L]')
      
      ax2.grid()
      ax2.set_ylabel('X [g/L]')
      ax2.set_xlabel('Time [h]') 
      
      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax1.plot(sim_res['time'],sim_res['bioreactor.c[2]'],color='b',linestyle=linetype)")
      diagrams.append("ax2.plot(sim_res['time'],sim_res['bioreactor.c[1]'],color='r',linestyle=linetype)")   
      
   elif plotType == 'Demo_2':
   
      plt.figure()
      ax1 = plt.subplot(2,1,1)
      ax2 = plt.subplot(2,1,2)

      ax1.set_title(title)
      ax1.grid()
      ax1.set_ylabel('S [g/L]')
      
      ax2.grid()
      ax2.set_ylabel('X [g/L]')
      ax2.set_xlabel('Time [h]') 
      
      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax1.plot(sim_res['time'],sim_res['bioreactor.c[2]'],'b*')")
      diagrams.append("ax2.plot(sim_res['time'],sim_res['bioreactor.c[1]'],'r*')") 

   elif plotType == 'PhasePlane':
       
      plt.figure()
      ax = plt.subplot(1,1,1)
    
      ax.set_title(title)
      ax.grid()
      ax.set_ylabel('S [g/L]')
      ax.set_xlabel('X [g/L]')

      # List of commands to be executed by simu() after a simulation         
      diagrams.clear()
      diagrams.append("ax.plot(sim_res['bioreactor.m[1]'],sim_res['bioreactor.m[2]'],color='b',linestyle=linetype)")
             
   else:
      print("Plot window type not correct")

# Define describtions partly coded here and partly taken from the FMU
def describe(name, decimals=3):
   """Look up description of culture, media, as well as parameters and variables in the model code"""
        
   if name == 'culture':
      print('Simplified text book model - only substrate S and cell concentration X')      
 
   elif name in ['broth', 'liquidphase', 'media']: 
      """Describe medium used"""
      
      X = model_get('liquidphase.X') 
      X_description = model_get_variable_description('liquidphase.X') 
      X_mw = model_get('liquidphase.mw[1]')
         
      S = model_get('liquidphase.S') 
      S_description = model_get_variable_description('liquidphase.S')
      S_mw = model_get('liquidphase.mw[2]')
         
      print('Reactor broth substances included in the model')
      print()
      print(X_description, '    index = ', X, 'molecular weight = ', X_mw, 'Da')
      print(S_description, 'index = ', S, 'molecular weight = ', S_mw, 'Da')
  
   elif name in ['parts']:
      describe_parts(component_list_minimum)

   elif name in ['MSL']:
      describe_MSL()

   else:
      describe_general(name, decimals)
      
#------------------------------------------------------------------------------------------------------------------
#  General code 
FMU_explore = 'FMU-explore for FMPy version 0.9.6'
#------------------------------------------------------------------------------------------------------------------

# Define function par() for parameter update
def par(parDict=parDict, parCheck=parCheck, parLocation=parLocation, *x, **x_kwarg):
   """ Set parameter values if available in the predefined dictionaryt parDict. """
   x_kwarg.update(*x)
   x_temp = {}
   for key in x_kwarg.keys():
      if key in parDict.keys():
         x_temp.update({key: x_kwarg[key]})
      else:
         print('Error:', key, '- seems not an accessible parameter - check the spelling')
   parDict.update(x_temp)
   
   parErrors = [requirement for requirement in parCheck if not(eval(requirement))]
   if not parErrors == []:
      print('Error - the following requirements do not hold:')
      for index, item in enumerate(parErrors): print(item)

# Define function init() for initial values update
def init(parDict=parDict, *x, **x_kwarg):
   """ Set initial values and the name should contain string '_0' to be accepted.
       The function can handle general parameter string location names if entered as a dictionary. """
   x_kwarg.update(*x)
   x_init={}
   for key in x_kwarg.keys():
      if '_0' in key: 
         x_init.update({key: x_kwarg[key]})
      else:
         print('Error:', key, '- seems not an initial value, use par() instead - check the spelling')
   parDict.update(x_init)

# Define fuctions similar to pyfmi model.get(), model.get_variable_descirption(), model.get_variable_unit()
def model_get(parLoc, model_description=model_description):
   """ Function corresponds to pyfmi model.get() but returns just a value and not a list"""
   par_var = model_description.modelVariables
   for k in range(len(par_var)):
      if par_var[k].name == parLoc:
         if par_var[k].variability in ['constant', 'fixed']:        
            value = float(par_var[k].start)        
         elif par_var[k].variability == 'continuous':
            try:
               timeSeries = sim_res[par_var[k].name]
               value = timeSeries[-1]
            except (AttributeError, ValueError):
               value = None
               print('Variable not logged')
         else:
            value = None
   return value

def model_get_variable_description(parLoc, model_description=model_description):
   """ Function corresponds to pyfmi model.get_variable_description() but returns just a value and not a list"""
   par_var = model_description.modelVariables
#   value = [x[1] for x in [(par_var[k].name, par_var[k].description) for k in range(len(par_var))] if parLoc in x[0]]
   value = [x.description for x in par_var if parLoc in x.name]   
   return value[0]
   
def model_get_variable_unit(parLoc, model_description=model_description):
   """ Function corresponds to pyfmi model.get_variable_unit() but returns just a value and not a list"""
   par_var = model_description.modelVariables
#   value = [x[1] for x in [(par_var[k].name, par_var[k].unit) for k in range(len(par_var))] if parLoc in x[0]]
   value = [x.unit for x in par_var if parLoc in x.name]
   return value[0]
      
# Define function disp() for display of initial values and parameters
def disp(name='', decimals=3, mode='short'):
   """ Display intial values and parameters in the model that include "name" and is in parLocation list.
       Note, it does not take the value from the dictionary par but from the model. """
   
   def dict_reverser(d):
      seen = set()
      return {v: k for k, v in d.items() if v not in seen or seen.add(v)}
   
   if mode in ['short']:
      k = 0
      for Location in [parLocation[k] for k in parDict.keys()]:
         if name in Location:
            if type(model_get(Location)) != np.bool_:
               print(dict_reverser(parLocation)[Location] , ':', np.round(model_get(Location),decimals))
            else:
               print(dict_reverser(parLocation)[Location] , ':', model_get(Location))               
         else:
            k = k+1
      if k == len(parLocation):
         for parName in parDict.keys():
            if name in parName:
               if type(model_get(Location)) != np.bool_:
                  print(parName,':', np.round(model_get(parLocation[parName]),decimals))
               else: 
                  print(parName,':', model_get(parLocation[parName])[0])

   if mode in ['long','location']:
      k = 0
      for Location in [parLocation[k] for k in parDict.keys()]:
         if name in Location:
            if type(model_get(Location)) != np.bool_:       
               print(Location,':', dict_reverser(parLocation)[Location] , ':', np.round(model_get(Location),decimals))
         else:
            k = k+1
      if k == len(parLocation):
         for parName in parDict.keys():
            if name in parName:
               if type(model_get(Location)) != np.bool_:
                  print(parLocation[parName], ':', dict_reverser(parLocation)[Location], ':', parName,':', 
                     np.round(model_get(parLocation[parName]),decimals))

# Line types
def setLines(lines=['-','--',':','-.']):
   """Set list of linetypes used in plots"""
   global linecycler
   linecycler = cycle(lines)

# Show plots from sim_res, just that
def show(diagrams=diagrams):
   """Show diagrams chosen by newplot()"""
   # Plot pen
   linetype = next(linecycler)    
   # Plot diagrams 
   for command in diagrams: eval(command)

# Define simulation
def simu(simulationTime=simulationTime, diagrams=diagrams, output_interval=None):
   global sim_res
   
   def extract_variables(diagrams):
       output = []
       variables = [v for v in model_description.modelVariables if v.causality == 'local']
       for j in range(len(diagrams)):
           for k in range(len(variables)):
               if variables[k].name in diagrams[j]:
                   output.append(variables[k].name)
       return output
   
   sim_res = simulate_fmu(
      filename = fmu_model,
      validate = False,
      start_time = 0,
      stop_time = simulationTime,
      output_interval = output_interval,
#      solver = 'CVode',
#      step_size = 1e-2,
      record_events = True,
      start_values = {parLocation[k]:parDict[k] for k in parDict.keys()},
      fmi_call_logger = None,
      output = list(set(extract_variables(diagrams) + key_variables))
   )
   # Plot diagrams from simulation
   linetype = next(linecycler)    
   for command in diagrams: eval(command)
      
# Describe model parts of the combined system
def describe_parts(component_list=[]):
   """List all parts of the model""" 
       
   def model_component(variable_name):
      i = 0
      name = ''
      finished = False
      if not variable_name[0] == '_':
         while not finished:
            name = name + variable_name[i]
            if i == len(variable_name)-1:
                finished = True 
            elif variable_name[i+1] in ['.', '(']: 
                finished = True
            else: 
                i=i+1
      if name in ['der', 'temp_1', 'temp_2', 'temp_3', 'temp_4', 'temp_5', 'temp_6', 'temp_7']: name = ''
      return name
    
#   variables = list(model.get_model_variables().keys())
   variables = [v.name for v in model_description.modelVariables]
        
   for i in range(len(variables)):
      component = model_component(variables[i])
      if (component not in component_list) \
      & (component not in ['','BPL', 'Customer', 'today[1]', 'today[2]', 'today[3]', 'temp_2', 'temp_3']):
         component_list.append(component)
      
   print(sorted(component_list, key=str.casefold))

# Describe MSL   
def describe_MSL(flag_vendor=flag_vendor):
   """List MSL version and components used"""
   print('MSL:', MSL_usage)
 
# Describe parameters and variables in the Modelica code
def describe_general(name, decimals):
  
   if name == 'time':
      description = 'Time'
      unit = 'h'
      print(description,'[',unit,']')
      
   elif name in parLocation.keys():
      description = model_get_variable_description(parLocation[name])
      value = model_get(parLocation[name])
      try:
         unit = model_get_variable_unit(parLocation[name])
      except FMUException:
         unit =''
      if unit =='':
         if type(value) != np.bool_:
            print(description, ':', np.round(value, decimals))
         else:
            print(description, ':', value)            
      else:
        print(description, ':', np.round(value, decimals), '[',unit,']')
                  
   else:
      description = model_get_variable_description(name)
      value = model_get(name)
      try:
         unit = model_get_variable_unit(name)
      except FMUException:
         unit =''
      if unit =='':
         if type(value) != np.bool_:
            print(description, ':', np.round(value, decimals))
         else:
            print(description, ':', value)     
      else:
         print(description, ':', np.round(value, decimals), '[',unit,']')
         
# Describe framework
def BPL_info():
   print()
   print('Model for bioreactor has been setup. Key commands:')
   print(' - par()       - change of parameters and initial values')
   print(' - init()      - change initial values only')
   print(' - simu()      - simulate and plot')
   print(' - newplot()   - make a new plot')
   print(' - show()      - show plot from previous simulation')
   print(' - disp()      - display parameters and initial values from the last simulation')
   print(' - describe()  - describe culture, broth, parameters, variables with values/units')
   print()
   print('Note that both disp() and describe() takes values from the last simulation')
   print()
   print('Brief information about a command by help(), eg help(simu)') 
   print('Key system information is listed with the command system_info()')

def system_info():
   """Print system information"""
#   FMU_type = model.__class__.__name__
   constants = [v for v in model_description.modelVariables if v.causality == 'local']
   
   print()
   print('System information')
   print(' -OS:', platform.system())
   print(' -Python:', platform.python_version())
   try:
       scipy_ver = scipy.__version__
       print(' -Scipy:',scipy_ver)
   except NameError:
       print(' -Scipy: not installed in the notebook')
   print(' -FMPy:', version('fmpy'))
   print(' -FMU by:', read_model_description(fmu_model).generationTool)
   print(' -FMI:', read_model_description(fmu_model).fmiVersion)
   if model_description.modelExchange is None:
      print(' -Type: CS')
   else:
      print(' -Type: ME')
   print(' -Name:', read_model_description(fmu_model).modelName)
   print(' -Generated:', read_model_description(fmu_model).generationDateAndTime)
   print(' -MSL:', MSL_version)    
   print(' -Description:', BPL_version)   
   print(' -Interaction:', FMU_explore)

#------------------------------------------------------------------------------------------------------------------
#  Startup
#------------------------------------------------------------------------------------------------------------------

BPL_info()
