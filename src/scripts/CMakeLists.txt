#
# Cmake fules for scipts
#

# Install scripts
# ------------------
  file(GLOB scripts *sh)
  install(PROGRAMS ${scripts} DESTINATION bin)

  file(GLOB scripts *.j)
  install(PROGRAMS ${scripts} DESTINATION bin/${this})
