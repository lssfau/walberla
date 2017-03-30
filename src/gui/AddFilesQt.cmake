# Scans the given directory and returns all subfolders
macro (list_subdirectories  RETURN_VAL curdir)

  FILE ( GLOB _subdirs RELATIVE ${curdir} "${curdir}/*" )
  foreach(dir ${_subdirs})
    if(IS_DIRECTORY ${curdir}/${dir})
        LIST( APPEND ${RETURN_VAL} ${curdir}/${dir} ) 
    endif()
  endforeach(dir)
   
endmacro(list_subdirectories)


# Traverses given list of subfolders and 
# writes all files needed for compilation in OUT_FILES variable
# Processed Files:
#   - cpp files (are omitted if NO_CPPS is specified)
#   - scans all headers for Q_OBJECT macro, and adds moc files if necessary
#      (are omitted if NO_MOCS is specified) 
#   - adds rules to process .ui and .qrc files
#
# Example Usage:
#    add_files_qt( outFiles source1 sourceDir2 NO_CPPS NO_MOCS)
# 
# The position of the flag arguments NO_CPPS and NO_MOCS is not important so
#    add_files_qt( outFiles  NO_CPPS NO_MOCS source1 sourceDir2)
#    is also valid ( and equivalent)
macro ( add_files_qt OUT_FILES )
    
    set(argList ${ARGN}) # strangely argList behaves different than ARGN
    set( ADD_CPPS ON)
    set( ADD_MOCS ON)

    list(FIND argList  "NO_CPPS" noCppIndex)
    if( noCppIndex GREATER -1)
        set( ADD_CPPS OFF)
        list(REMOVE_ITEM argList "NO_CPPS" )
    endif()
    
    list(FIND argList  "NO_MOCS" noMocIndex)
    if( noMocIndex GREATER -1)
        set( ADD_MOCS OFF)
        list(REMOVE_ITEM argList "NO_MOCS" )
    endif()


    set (SOURCE_GROUP_PREFIX "")

    list(FIND argList  "SOURCE_GROUP" srcGroupIndex)
    if( srcGroupIndex GREATER -1)
        math(EXPR incrIndex  "${srcGroupIndex} + 1" )
        list(GET argList ${incrIndex} SOURCE_GROUP_PREFIX)
        list(REMOVE_AT argList ${srcGroupIndex} )
        list(REMOVE_AT argList ${srcGroupIndex} )
    endif()

    
  	foreach( DIRECTORY ${argList} ) 	
		set( ALL_HEADER )
		set( MOC_HEADER )	
		set( CPP )
		set( UI )
		set( RES )
		set( MOC_OUTPUT_FILES )
		set( UI_OUTPUT_FILES )
		set( RES_OUTPUT_FILES )
	
		file( GLOB ALL_HEADER . ${DIRECTORY}/*.h )
		file( GLOB UI         . ${DIRECTORY}/*.ui )
		file( GLOB RES        . ${DIRECTORY}/*.qrc )
		file( GLOB CPP        . ${DIRECTORY}/*.cpp )

		# Find files which need moc processing 
		# (are recognized on Q_OBJECT macro inside file)
		foreach( _current_HEADER ${ALL_HEADER} )
			GET_FILENAME_COMPONENT(_abs_HEADER ${_current_HEADER} ABSOLUTE)
			FILE( READ ${_abs_HEADER} _contents)
			STRING( REGEX MATCHALL "Q_OBJECT" _match  "${_contents}" )

			IF( _match)	
				LIST( APPEND MOC_HEADER ${_current_HEADER} )
			ENDIF (_match)
		endforeach( _current_HEADER )
	
	
		QT4_WRAP_CPP(MOC_OUTPUT_FILES ${MOC_HEADER} OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED -DBOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION )
		QT4_WRAP_UI(UI_OUTPUT_FILES ${UI} )
		QT4_ADD_RESOURCES(RES_OUTPUT_FILES ${RES})

		LIST( APPEND ${OUT_FILES} ${UI_OUTPUT_FILES} ${RES_OUTPUT_FILES} ${ALL_HEADER} )	
	    if( ADD_MOCS )
          	LIST( APPEND ${OUT_FILES} ${MOC_OUTPUT_FILES}  )	    
        endif()
	    if( ADD_CPPS )
          	LIST( APPEND ${OUT_FILES} ${CPP}  )	    
        endif()	
    

	endforeach( DIRECTORY )

endmacro ()
