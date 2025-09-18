###################################################################################################
#
# Module sqlite
#
###################################################################################################

add_library( walberla_sqlite )
add_library( walberla::sqlite ALIAS walberla_sqlite )
target_link_libraries( walberla_sqlite PUBLIC walberla::core PRIVATE sqlite3 )
target_sources( walberla_sqlite
        PUBLIC
        SQLite.h
)

target_sources( walberla_sqlite
        PRIVATE
        SQLite.cpp
)

###################################################################################################
