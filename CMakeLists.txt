project ( ravi-simplex )
cmake_minimum_required ( VERSION 2.8 )

option(USE_LUA53 "Controls whether the library will be for Lua 5.3 or Ravi, default is OFF" OFF)

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake" )

include(${CMAKE_INSTALL_PREFIX}/cmake/lua-config.cmake)

message(STATUS "LUA_INCLUDE_DIR: ${LUA_INCLUDE_DIR}")
message(STATUS "LUA_LIBRARIES  : ${LUA_LIBRARIES}")

INSTALL(DIRECTORY lua/ DESTINATION "${LUADIR}")
