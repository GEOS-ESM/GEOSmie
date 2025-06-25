#
# Cmake fules for scipts
#

# Install scripts
# ------------------
  file(GLOB scripts *sh)
  install(PROGRAMS ${scripts} DESTINATION bin)

  file(GLOB scripts *.py)
  install(PROGRAMS ${scripts} DESTINATION bin)
