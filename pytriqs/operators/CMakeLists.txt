
SET(PYTHON_SOURCES 
  ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
)

install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/operators)

triqs_python_extension(operators operators)

# All subdirs
add_all_subdirectories_with_cmakelist()
