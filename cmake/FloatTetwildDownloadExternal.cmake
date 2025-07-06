################################################################################

include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(FLOAT_TETWILD_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(FLOAT_TETWILD_EXTRA_OPTIONS "")
endif()

# Shortcut function
function(float_tetwild_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${FLOAT_TETWILD_EXTERNAL}/${name}
        DOWNLOAD_DIR ${FLOAT_TETWILD_EXTERNAL}/.cache/${name}
        QUIET
        ${FLOAT_TETWILD_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################

## libigl
function(float_tetwild_download_libigl)
    float_tetwild_download_project(libigl
        GIT_REPOSITORY https://github.com/libigl/libigl.git
        #GIT_TAG        45cfc79fede992ea3923ded9de3c21d1c4faced1
        GIT_TAG        v2.6.0
    )
endfunction()

## Json
function(float_tetwild_download_json)
    float_tetwild_download_project(json
        GIT_REPOSITORY https://github.com/jdumas/json
        GIT_TAG        0901d33bf6e7dfe6f70fd9d142c8f5c6695c6c5b
    )
endfunction()

## Catch2
function(float_tetwild_download_catch2)
    float_tetwild_download_project(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.1.1
    )
endfunction()

## CLI11
function(float_tetwild_download_cli11)
    float_tetwild_download_project(cli11
            GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
            GIT_TAG        v2.5.0
    )
endfunction()

## tbb
set(TBB_TEST OFF CACHE BOOL "" FORCE)
function(float_tetwild_download_tbb)
    float_tetwild_download_project(tbb
    GIT_REPOSITORY https://github.com/uxlfoundation/oneTBB.git
    GIT_TAG v2022.2.0-rc1
    #    GIT_REPOSITORY https://github.com/wjakob/tbb.git
	#GIT_TAG 	 faab196203247a170a5efbda06b4f417dbba03f2
    )
endfunction()

## Sanitizers
function(float_tetwild_download_sanitizers)
    float_tetwild_download_project(sanitizers-cmake
        GIT_REPOSITORY https://github.com/arsenm/sanitizers-cmake.git
        GIT_TAG        6947cff3a9c9305eb9c16135dd81da3feb4bf87f
    )
endfunction()

## fmt
function(float_tetwild_download_fmt)
    float_tetwild_download_project(fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG        9.0.0
    )
endfunction()

## spdlog
function(float_tetwild_download_spdlog)
    float_tetwild_download_project(spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG        v1.11.0
    )
endfunction()

## Geogram LGPL
function(float_tetwild_download_geogram)
    float_tetwild_download_project(geogram
        GIT_REPOSITORY   https://github.com/BrunoLevy/geogram.git
        GIT_TAG          v1.9.6
        #GIT_REPOSITORY https://github.com/Yixin-Hu/geogram
        #GIT_TAG        b613750341a6cdd31ae8df80ecfc26ac7ca1a6ad
    )
endfunction()

## aabbcc
function(float_tetwild_download_aabbcc)
    float_tetwild_download_project(aabbcc
            GIT_REPOSITORY https://github.com/lohedges/aabbcc.git
            GIT_TAG        0c85e61362d384d70c71946826bfed0fb24a74ba
            )
endfunction()

### winding number
#function(float_tetwild_download_windingnumber)
#    float_tetwild_download_project(windingnumber
#            GIT_REPOSITORY https://github.com/alecjacobson/WindingNumber.git
#            GIT_TAG        bde8780ec848fa71c1294a0af50347e968b19493
#            )
#endfunction()

## exact envelope
function(float_tetwild_download_exact_envelope)
    float_tetwild_download_project(exact_envelope
            GIT_REPOSITORY https://github.com/wangbolun300/fast-envelope
            GIT_TAG        520ee04b6c69a802db31d1fd3a3e6e382d10ef98
            )
endfunction()
