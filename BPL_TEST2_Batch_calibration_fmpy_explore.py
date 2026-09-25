# Setup application functions BPL_TEST2_Batch_calibration, dependent on previous import from fmu_explore
# Author: Jan Peter Axelsson
# ------------------------------------------------------------------------------------------------------------------
# 2026-09-10 - Created
# 2026-09-18 - Brought in matplotlib.pylot since used here
# 2026-09-18 - Brought numpy back here since used in several notebooks
# 2026-09-25 - Change indentaiton from 3 spaces to 4 using black
# ------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------
#  Framework
# ------------------------------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------------------------------------------
#  Specific application constructs:  newplot(), describe()
# ------------------------------------------------------------------------------------------------------------------


# Define standard diagrams
def newplot(title="Batch cultivation", plotType="TimeSeries"):
    """Standard plot window
     title = ''
    two possible diagrams
     diagram = 'TimeSeries' default
     diagram = 'PhasePlane'"""

    # Reset pens
    resetPen()

    # Plot diagram
    if plotType == "TimeSeries":

        ax1 = plt.subplot(2, 1, 1)
        ax2 = plt.subplot(2, 1, 2)

        ax.clear()
        ax.append(ax1)
        ax.append(ax2)

        ax[0].set_title(title)
        ax[0].grid()
        ax[0].set_ylabel("X and S [g/L]")

        ax[1].grid()
        ax[1].set_ylabel("mu [1/h]")
        ax[1].set_xlabel("Time [h]")

        # List of commands to be executed by simu() after a simulation
        diagrams.clear()
        diagrams.append(
            "ax[0].plot(t,sim_res['bioreactor.c[1]'],color='r',linestyle=linetype)"
        )
        diagrams.append(
            "ax[0].plot(t,sim_res['bioreactor.c[2]'],color='b',linestyle=linetype)"
        )
        diagrams.append("ax[0].legend(['X','S'])")
        diagrams.append(
            "ax[1].plot(t,sim_res['bioreactor.culture.q[1]'],color='r',linestyle=linetype)"
        )

    elif plotType == "TimeSeries2":

        ax1 = plt.subplot(2, 1, 1)
        ax2 = plt.subplot(2, 1, 2)

        ax.clear()
        ax.append(ax1)
        ax.append(ax2)

        ax[0].set_title(title)
        ax[0].grid()
        ax[0].set_ylabel("S [g/L]")

        ax[1].grid()
        ax[1].set_ylabel("X [g/L]")
        ax[1].set_xlabel("Time [h]")

        # List of commands to be executed by simu() after a simulation
        diagrams.clear()
        diagrams.append(
            "ax[0].plot(t,sim_res['bioreactor.c[2]'],color='b',linestyle=linetype)"
        )
        diagrams.append(
            "ax[1].plot(t,sim_res['bioreactor.c[1]'],color='b',linestyle=linetype)"
        )

    elif plotType == "Demo_1":

        ax1 = plt.subplot(2, 1, 1)
        ax2 = plt.subplot(2, 1, 2)

        ax.clear()
        ax.append(ax1)
        ax.append(ax2)

        ax[0].set_title(title)
        ax[0].grid()
        ax[0].set_ylabel("S [g/L]")

        ax[1].grid()
        ax[1].set_ylabel("X [g/L]")
        ax[1].set_xlabel("Time [h]")

        # List of commands to be executed by simu() after a simulation
        diagrams.clear()
        diagrams.append(
            "ax[0].plot(sim_res['time'],sim_res['bioreactor.c[2]'],color='b',linestyle=linetype)"
        )
        diagrams.append(
            "ax[1].plot(sim_res['time'],sim_res['bioreactor.c[1]'],color='r',linestyle=linetype)"
        )

    elif plotType == "Demo_2":

        ax1 = plt.subplot(2, 1, 1)
        ax2 = plt.subplot(2, 1, 2)

        ax.clear()
        ax.append(ax1)
        ax.append(ax2)

        ax[0].set_title(title)
        ax[0].grid()
        ax[0].set_ylabel("S [g/L]")

        ax[1].grid()
        ax[1].set_ylabel("X [g/L]")
        ax[1].set_xlabel("Time [h]")

        # List of commands to be executed by simu() after a simulation
        diagrams.clear()
        diagrams.append("ax[0].plot(sim_res['time'],sim_res['bioreactor.c[2]'],'b*')")
        diagrams.append("ax[1].plot(sim_res['time'],sim_res['bioreactor.c[1]'],'r*')")

    elif plotType == "PhasePlane":

        ax1 = plt.subplot(1, 1, 1)

        ax.clear()
        ax.append(ax1)

        ax[0].set_title(title)
        ax[0].grid()
        ax[0].set_ylabel("S [g/L]")
        ax[0].set_xlabel("X [g/L]")

        # List of commands to be executed by simu() after a simulation
        diagrams.clear()
        diagrams.append(
            "ax.plot(sim_res['bioreactor.m[1]'],sim_res['bioreactor.m[2]'],color='b',linestyle=linetype)"
        )

    else:
        print("Plot window type not correct")


# Define describtions partly coded here and partly taken from the FMU
def describe(name, decimals=3):
    """Look up description of culture, media, as well as parameters and variables in the model code"""

    if name == "culture":
        print("Simplified text book model - only substrate S and cell concentration X")

    elif name in ["broth", "liquidphase", "media"]:
        """Describe medium used"""

        X = model_get("liquidphase.X")
        X_description = model_get_variable_description("liquidphase.X")
        X_mw = model_get("liquidphase.mw[1]")

        S = model_get("liquidphase.S")
        S_description = model_get_variable_description("liquidphase.S")
        S_mw = model_get("liquidphase.mw[2]")

        print("Reactor broth substances included in the model")
        print()
        print(X_description, "    index = ", X, "molecular weight = ", X_mw, "Da")
        print(S_description, "index = ", S, "molecular weight = ", S_mw, "Da")

    elif name in ["parts"]:
        describe_parts(component_list_minimum)

    elif name in ["MSL"]:
        describe_MSL()

    else:
        describe_general(name, decimals)


# ------------------------------------------------------------------------------------------------------------------
#  Startup
# ------------------------------------------------------------------------------------------------------------------

FMU_explore_info()
