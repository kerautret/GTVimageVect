
SET(SRCs
  tv-triangulation-color
)

  
FOREACH(FILE ${SRCs})
  add_executable(${FILE} ${FILE} ${SRC_COMP} ${SOURCE_HEADER_FILES})
  target_link_libraries( ${FILE}  ${CAIRO_LIBRAIRIES} ${DGTAL_LIBRARIES} ${GTVDep} )
ENDFOREACH(FILE)
