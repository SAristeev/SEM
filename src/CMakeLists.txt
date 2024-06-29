add_executable(SEM main.cpp fc.cpp export2vtk.cpp solver.cpp parametric_hex.cpp gll.cpp debug_helper.cpp fem.cpp algebra.cpp)

target_include_directories(SEM PRIVATE "include")

target_link_libraries(SEM PRIVATE MKL::MKL)
target_link_libraries(SEM PRIVATE TBB::tbb)
target_link_libraries(SEM PRIVATE ${VTK_LIBRARIES})