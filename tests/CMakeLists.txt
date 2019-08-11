SET(GTVV_TESTS
   testBasicVectoImage
)


SET(SOURCE_TESTS_HEADER_FILES
   ${PROJECT_SOURCE_DIR}/src/BasicVectoImageExporter.h
)


FOREACH(FILE ${GTVV_TESTS})
  add_executable(${FILE} ${FILE} ${SRC_COMP} ${SOURCE_TESTS_HEADER_FILES})
  target_link_libraries( ${FILE}  ${CAIRO_LIBRAIRIES} ${DGTAL_LIBRARIES} ${GTVDep}  )
ENDFOREACH(FILE)








 





