#
# Cmake rules for gsf
#

  esma_set_this()
  file(GLOB scripts *.py)
  install(PROGRAMS ${scripts} DESTINATION bin/${this})

ecbuild_add_executable (
  TARGET spher_expan.x
  SOURCES spher_expan.f params.h
  NOINSTALL
  )

install (
    TARGETS spher_expan.x
    DESTINATION bin/${this}
    )
