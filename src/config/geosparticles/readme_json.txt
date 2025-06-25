Explanation of parameters available in the JSON data files

Mandatory parameters:

rhop0: Dry rho

rh: List of RH values to calculate the humidified values at

rhDep: RH-dependence dictionary
Format: {"type": FORMAT, "params": PARAMETERS}

psd: Particle size distribution dictionary 
Format: {"type": FORMAT, "params": PARAMETERS}
    
ri: Refractive index definition
Format: {"format": FORMAT, "path": LIST OF PATHS}

Optional parameters:

maxrh: Treat rh values above this as if they were this value

hydrophobic: true/false, create special hydrophobic size bin
