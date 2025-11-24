#
# Cmake rules for utils
#
# Python executables
# ------------------
  esma_set_this()
  file(GLOB scripts *.py)
  install(PROGRAMS ${scripts} DESTINATION bin/${this})
