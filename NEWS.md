# PPTcirc 0.2.0

This is a minor release with various small features and bug fixes.

* `dsimpostppt()` bug fix: no longer simulates from mu by default (only simulates from mu when hm = 1).Additionally, all the results from simulating from alpha and mu are saved after the burn in period and thinning chain. 

* `postppt.plot()` gains new parameters: `xlim`, `tol` and `sep`. This will help the user to modify the limits in x-axis for lineal plots and choose the tol and separation of point in circular plots. 
