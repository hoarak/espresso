file(GLOB Observables_SRC *.cpp)
add_library(Observables SHARED ${Observables_SRC})
add_dependencies(Observables EspressoConfig)
set_target_properties(Observables PROPERTIES MACOSX_RPATH TRUE)
