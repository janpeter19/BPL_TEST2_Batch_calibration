# BPL_TEST2_Batch_calibration

Material of standard text book model of batch cultivation where the model is calibrated using time series data for substrate and cell concentration. It is possible to upload your own batch data and calibrate the model for your data!

A bioprocess model need to be calibrated by estimating key process parameters. The model here has only three parameters $Y$, $qSmax$ and $Ks$, but the last parameter has little influence. The intial values $V_0$, $VX_0$ and $VS_0$ are given.

We use here methods that only need evaluation of the loss function for each tested parameter setting, thus does not use derivative information. These methods are collected in Python Scipy-package in the function scipy.optimize.minimize() and it is easy to change method and adjust parameters of the algorithm. The advantage with these methods is that you can easilty address more complex models with har non-linearities and more parameters. The drawback is that optimization is more time consuming than methods that use derivative information. 

In Figure 1 below we see the results of fitting the model to simulated data and works perfectly well. Here parameters were given with min- and max-bounds. There are fewer methods that use bounds on parameters. Here we choose to work with SLSQP since it gave fasd and good convergence.

![](Fig_1_Batch_simu_data.png)

In Figure 2 we see the contour plot of the loss function, simplified by keeping the Ks-parameter value constant. Clearly there is a unique minima for the range of parameters considered. Further the contours has an ellipsoid-shape. The more narrow the form is the harder the parameter estimation problem is. In this case not too difficult. 

![](Fig_2_Batch_loss_function_contour.png)


A Jupyter notebook that go through the parameter estimation and generate the figures above you click here
[start BPL notebook](https://colab.research.google.com/github/janpeter19/BPL_TEST2_Batch_calibration/blob/main/BPL_TEST2_Batch_calibration_colab.ipynb).
In the menu choose Runtime/Run all. The installation takes just a few minutes. The subsequent execution of all the simulations take just about a minute or so.

Note, that you can easily change from working with simulated data to your own measurement data stored in an excel-file on your own computer. Just follow the instruction in the notebook!

License information:
* The binary-file with extension FMU is shared under the permissive MIT-license
* The other files are shared under the GPL 3.0 license