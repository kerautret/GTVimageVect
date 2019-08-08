SET(GTVV_TESTS
   testBasicVectoImage
)




FOREACH(FILE ${GTVV_TESTS})
  add_executable(${FILE} ${FILE} ${SRC_COMP})
  target_link_libraries( ${FILE}  ${CAIRO_LIBRAIRIES} ${DGTAL_LIBRARIES} ${GTVDep}  )
ENDFOREACH(FILE)








 





