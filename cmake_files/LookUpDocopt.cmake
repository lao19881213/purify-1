# Looks up [Docopt](http://docopt.github.io/)
#
# - GIT_REPOSITORY: defaults to https://github.com/docopt/docopt.cpp.git)
# - GIT_TAG: defaults to v0.6.2 (latest tag)
# - BUILD_TYPE: defaults to Release
#
if(Docopt_ARGUMENTS)
  cmake_parse_arguments(Docopt "" "GIT_REPOSITORY;GIT_TAG;BUILD_TYPE" "" ${Docopt_ARGUMENTS})
endif()
if(NOT Docopt_GIT_REPOSITORY)
    set(Docopt_GIT_REPOSITORY https://github.com/docopt/docopt.cpp.git)
endif()
if(NOT Docopt_GIT_TAG)
    set(Docopt_GIT_TAG v0.6.2)
endif()
if(NOT Docopt_BUILD_TYPE)
    set(Docopt_BUILD_TYPE Release)
endif()


ExternalProject_Add(
    Lookup-Docopt
    PREFIX ${EXTERNAL_ROOT}
    GIT_REPOSITORY ${Docopt_GIT_REPOSITORY}
    GIT_TAG ${Docopt_GIT_TAG}
    CMAKE_ARGS
      -DBUILD_SHARED_LIBS=OFF
      -DCMAKE_INSTALL_PREFIX=${EXTERNAL_ROOT}
      -DCMAKE_BUILD_TYPE=${Sopt_BUILD_TYPE}
    INSTALL_DIR ${EXTERNAL_ROOT}
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)
add_recursive_cmake_step(Lookup-Docopt DEPENDEES install)
