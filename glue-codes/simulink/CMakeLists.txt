#
# Copyright 2023 National Renewable Energy Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

# List of module libraries to be linked into the MEX shared library. 
# This list will be linked twice to resolve undefined symbols due to linking order.
set(MEX_LIBS
  $<TARGET_FILE:openfast_prelib>
  $<TARGET_FILE:basicaerolib>
  $<TARGET_FILE:aerodynlib>
  $<TARGET_FILE:beamdynlib>
  $<TARGET_FILE:elastodynlib>
  $<TARGET_FILE:extptfm_mckflib>
  $<TARGET_FILE:feamlib>
  $<TARGET_FILE:hydrodynlib>
  $<TARGET_FILE:icedynlib>
  $<TARGET_FILE:icefloelib>
  $<TARGET_FILE:mappplib>
  $<TARGET_FILE:moordynlib>
  $<TARGET_FILE:orcaflexlib>
  $<TARGET_FILE:servodynlib-matlab> # MATLAB Specific
  $<TARGET_FILE:subdynlib>
  $<TARGET_FILE:extinflowlib>
  $<TARGET_FILE:extloadslib>
  $<TARGET_FILE:ifwlib>
  $<TARGET_FILE:extinflowtypeslib>
  $<TARGET_FILE:versioninfolib>
  $<TARGET_FILE:nwtclibs-matlab>    # MATLAB Specific
  $<TARGET_FILE:seastlib>
  $<TARGET_FILE:aerodisklib>
  $<TARGET_FILE:sedlib>
)

# Build the matlab shared library (mex) using the current toolchain.
# Recommpiles FAST_* files with COMPILE_SIMULINK defined.
# Links directly to library files, bypassing dependencies so MATLAB specific
# libraries can be substituted (nwtclibs and servodynlib). openfast_postlib
# is used as a dependency to ensure that all library files exist before
# the FAST_SFunc target is built.
matlab_add_mex(
  NAME FAST_SFunc
  SRC 
    src/FAST_SFunc.c
    ${PROJECT_SOURCE_DIR}/modules/openfast-library/src/FAST_Subs.f90
    ${PROJECT_SOURCE_DIR}/modules/openfast-library/src/FAST_Lin.f90
    ${PROJECT_SOURCE_DIR}/modules/openfast-library/src/FAST_Mods.f90
    ${PROJECT_SOURCE_DIR}/modules/openfast-library/src/FAST_Solver.f90
    ${PROJECT_SOURCE_DIR}/modules/openfast-library/src/FAST_Library.f90
  LINK_TO 
    ${MEX_LIBS}
    ${MEX_LIBS} # DO NOT REMOVE (needed to ensure no unresolved symbols)
    ${LAPACK_LIBRARIES} 
    ${CMAKE_DL_LIBS} 
    ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}
)
add_dependencies(FAST_SFunc openfast_postlib nwtclibs-matlab servodynlib-matlab)
target_compile_definitions(FAST_SFunc PUBLIC COMPILE_SIMULINK)
set_target_properties(FAST_SFunc PROPERTIES 
  Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/matlab)
target_include_directories(FAST_SFunc PUBLIC 
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/modules/openfast-library/src>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/modules/externalinflow/src>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/modules/extloads/src>
)
if(APPLE OR UNIX)
   target_compile_definitions(FAST_SFunc PRIVATE IMPLICIT_DLLEXPORT)
endif()

install(TARGETS FAST_SFunc
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION bin
  ARCHIVE DESTINATION lib
)
