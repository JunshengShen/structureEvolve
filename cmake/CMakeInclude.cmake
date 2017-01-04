# ------------------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------------------
# Add new project into solution
function( Add_project projectName projectDir) # 1+ arguments: ARGV0 = ${projectName} / ARGV1 = ${projectDir}

if( ${ARGC} LESS 1 )
	message( SEND_ERROR "function Add_project: invalid number of arguments! (need project name and some additional common files)" )
else()
	####Gen DLL projectName	
	FILE(GLOB Project_${projectName}_SRC_SOURCES  App/${projectDir}/*.cpp App/${projectName}/*.def)
	FILE(GLOB Project_${projectName}_SRC_HEADERS  App/${projectDir}/*.h)
	
	set(_args ${ARGN})  
	
	add_library(${projectName} SHARED ${Project_${projectName}_SRC_SOURCES} 
	                        ${Project_${projectName}_SRC_HEADERS}                               				
				${_args}			
				      )
        ##Copy dll to 
        ADD_CUSTOM_COMMAND(TARGET ${projectName}
		  POST_BUILD          
	          COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${projectName}> ${destDir}
	)
	
endif()

endfunction()
