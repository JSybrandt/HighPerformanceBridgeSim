# Simulation

This is the typical simulation code with minor changes.
Firstly, we wrap SimulateDay and Bridge\_Pre\_Allocation
in functions in order to run them on palmetto.
Secondly, we change the name of "beta" (which is reserved 
on the palmetto version of matlab) to "bbeta".
