import netCDF4
import os

def doConversion(infn, outfn, pfx):
  
  fn = os.path.join(pfx, infn)
  fn2 = os.path.join(pfx, outfn)
  f = netCDF4.Dataset(fn, 'r+')
  g = netCDF4.Dataset(fn2, 'w')

  #First import the netcdf4 library

  # Read en existing NetCDF file and create a new one
  # f is going to be the existing NetCDF file from where we want to import data
  # and g is going to be the new file.
  
  # To copy the global attributes of the netCDF file
  
  for attname in f.ncattrs():
      setattr(g,attname,getattr(f,attname))
  
  # To copy the dimension of the netCDF file

  dimsizes = {'radius': 2}
  
  for dimname,dim in list(f.dimensions.items()):
    # if you want to make changes in the dimensions of the new file
    # you should add your own conditions here before the creation of the dimension.
    dimsize = len(dim)
    if dimname in dimsizes:
      dimsize = dimsizes[dimname]
    g.createDimension(dimname,dimsize)
    # To copy the variables of the netCDF file

  for varname,ncvar in list(f.variables.items()):
    # if you want to make changes in the variables of the new file
    # you should add your own conditions here before the creation of the variable.
    var = g.createVariable(varname,ncvar.dtype,ncvar.dimensions)
    #Proceed to copy the variable attributes
    for attname in ncvar.ncattrs():
       setattr(var,attname,getattr(ncvar,attname))
    #Finally copy the variable data to the new created variable
    if var.shape != ncvar.shape:
      if len(var.shape) == 1:
        for vari in range(len(var.shape)):
          var[vari] = ncvar[0]
      elif len(var.shape) == 2:
        for ri in range(var.shape[0]):
          for rhi in range(var.shape[1]):
            if ri == 0: # hydrophobic case
              var[0,rhi] = ncvar[0,0] # set zero RH case 
            else: # regular case
              var[ri,rhi] = ncvar[0,rhi] # set zero RH case 
      elif len(var.shape) == 3:
        for ri in range(var.shape[0]):
          for rhi in range(var.shape[1]):
            if ri == 0: # hydrophobic case
              var[0,rhi,:] = ncvar[0,0,:] # set zero RH case 
            else: # regular case
              var[ri,rhi,:] = ncvar[0,rhi,:] # set zero RH case 
      elif len(var.shape) == 4:
        if var.shape[0] == 6: # this is pback that we have to handle separately
          for ri in range(var.shape[1]):
            for rhi in range(var.shape[2]):
              if ri == 0: # hydrophobic case
                var[:,0,rhi,:] = ncvar[:,0,0,:] # set zero RH case 
              else: # regular case
                var[:,ri,rhi,:] = ncvar[:,0,rhi,:] # set zero RH case 

        else: # scattering matrix elements -- should be len 5 with nPol as first index, but is not yet
          for ri in range(var.shape[0]):
            for rhi in range(var.shape[1]):
              if ri == 0: # hydrophobic case
                var[0,rhi,:,:] = ncvar[0,0,:,:] # set zero RH case 
              else: # regular case
                var[ri,rhi,:,:] = ncvar[0,rhi,:,:] # set zero RH case 
      else:
        print("unhandled number of dimensions")
        print(ncvar)
        sys.exit()

    else:
      var[:] = ncvar[:]


  f.close()
  g.close()
