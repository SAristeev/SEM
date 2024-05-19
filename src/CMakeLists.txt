include_directories("include")

add_executable( SEM main.cpp )

target_link_libraries(SEM PRIVATE MKL::MKL)
target_link_libraries(SEM PRIVATE TBB::tbb)
target_link_libraries(SEM PRIVATE ${VTK_LIBRARIES})