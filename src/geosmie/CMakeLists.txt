#
# Cmake rules for geosmie
#

  esma_set_this()
  file(GLOB scripts *.py)
  install(PROGRAMS ${scripts} DESTINATION bin/${this})
  file(GLOB scripts *.j)
  install(PROGRAMS ${scripts} DESTINATION bin/${this})

  install(DIRECTORY data DESTINATION bin/${this})
  install(DIRECTORY namelists DESTINATION bin/${this})
