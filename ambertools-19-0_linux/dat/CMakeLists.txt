
# we use this block so that we don't spam the console of whoever runs the install

install(DIRECTORY . USE_SOURCE_PERMISSIONS DESTINATION ${DATADIR} COMPONENT Data)

if(NOT (HOST_WINDOWS OR TARGET_WINDOWS))
	installtime_create_symlink(${CMAKE_INSTALL_POSTFIX}dat/leap/cmd/leaprc.protein.ff12polL ${CMAKE_INSTALL_POSTFIX}dat/leap/cmd/leaprc.protein.ff12pol Data)
endif()
