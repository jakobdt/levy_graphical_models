# levy_graphical_models
R code related to the paper *Graphical models for Lévy processes*.

The file <code>simulation.R</code> contains code for simulating a multivariate process with Hüsler-Reiss type jumps.
As it is the process will be 12-dimensional and the underlying graphical structure is easily plotted.
The file includes code for numerical approximation of certain quantities which are needed for the simulation.
Note that running these bits of code may take several minutes.

The file <code>tree_learning.R</code> contains the necessary code for learning an underlying tree structure given increments of the process.
At the bottom of the file there is an example which relies on the simulation functionality in <code>simulation.R</code>.

The file <code>recovery_probability.R</code> contains code for calculating the probability of estimating the correct tree.
This is based on Monte Carlo simulation and uses code from <code>simulation.R</code> and <code>tree_learning.R</code>.
Note that running the Monte Carlo simulation with the currently specified parameters may take a couple of hours.
If one is playing around with this it is adviced to take e.g. K = 10 to quickly get some results.
Once everything works one can go back to e.g. K = 1000.

The file <code>data_analysis.R</code> performs the analysis of stock price data which is described in the paper *Graphical models for Lévy processes*.
The code uses functionality from the file <code>tree_learning.R</code>.
The data is contained in the file <code>data_all.csv</code>. There is price data for 103 US companies for the period 2010-01-04 to 2015-12-31.
**Note that this file does not contain the actual prices! Instead the file contains the log-prices.**
