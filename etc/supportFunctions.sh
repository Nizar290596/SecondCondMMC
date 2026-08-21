#======================================================================
# Support Functions
#======================================================================

# Function to check the cmake version
__cmake_version_above()
{
    CMAKE_VERSION_NUMBER=$(cmake --version | head -n 1 | cut -d " " -f 3)
    CMAKE_VERSION_MAJOR=$(echo ${CMAKE_VERSION_NUMBER} | cut -d . -f 1 )
    CMAKE_VERSION_Minor=$(echo ${CMAKE_VERSION_NUMBER} | cut -d . -f 2 )

    if [ $CMAKE_VERSION_MAJOR -lt 3 ]; then
        return 0
    fi
    if [ $CMAKE_VERSION_Minor -lt 12 ]; then
        return 0
    fi
    return 1
}

__die()
{
    RED='\033[0;31m'
    echo -e "${RED}$*"
    exit 1
}


__banner()
{
    echo "==============================================================================="
    echo $*
    echo "==============================================================================="
}



__bannerByFilePath()
{
    # Print out a banner to the current tutorial executed
    # Input: Provide path to tutorial directory
    # 
    # The path will then be deconstructed to give a descriptive banner 

    inputPath=$1

    currDirName=$caseName
    currPath=$inputPath

    dirPathToCase=()

    while [[ $currDirName != "tutorials" ]]; do
        currDirName=$(basename $currPath)
        # Append to list
        dirPathToCase+=($currDirName)
        currPath=$(dirname $currPath)
        if [[ $currPath == "/" ]]; then
            break;
        fi
    done

    # Create string to print
    bannerStr="."
    for (( idx=${#dirPathToCase[@]}-1 ; idx>=0 ; idx-- )) ; do
        temp="${dirPathToCase[idx]}"
        bannerStr="${bannerStr}/$temp"
    done


    __banner "Run tutorial: $bannerStr"
}

__warningBanner()
{
    BYellow='\033[1;33m'
    Color_Off='\033[0m'
    >&2 echo -e "${BYellow}==============================================================================="
    >&2 echo -e "Warning: $*"
    >&2 echo -e "=============================================================================== ${Color_Off}"
}


__abort()
{
        cat <<EOF
***************
*** ABORTED ***
***************
An error occurred. Exiting...
EOF
        exit 1
}

__checkOpenFOAMEnvironment()
{
    if [[ -z ${WM_PROJECT_DIR} ]]; then
        __banner Checking OpenFOAM environment
        __die OpenFOAM environment is not sourced
    fi
}


__checkOpenFOAMVersion()
{
    if [[ ${MMC_OPENFOAM_VERSION} != ${WM_PROJECT_VERSION} ]]; then
        __warningBanner Mismatch of OpenFOAM version
        echo "mmcFoam is written for OpenFOAM ${MMC_OPENFOAM_VERSION}"
        echo "You are trying to compile it with OpenFOAM ${WM_PROJECT_VERSION}"
    fi

    # Check for OpenFOAM patches
    __checkOpenFOAMPatches
}


__checkOpenFOAMPatches()
{
    # Check if OpenFOAM version includes mandatory patches
    for COMMIT_ID in ${MMC_OPENFOAM_PATCHES[@]}; do
        $(git -C ${WM_PROJECT_DIR} merge-base --is-ancestor $COMMIT_ID HEAD)
        FOUND_COMMIT=$?
        if [[ 0 -ne ${FOUND_COMMIT} ]]; then
            __warningBanner Missing OpenFOAM Patches
            echo "Install following OpenFOAM patch: ${COMMIT_ID}"
            echo "Navigate to your OpenFOAM directory, e.g., 'cd \${WM_PROJECT_DIR}' "\
                 "execute following commands:"
            echo "   git fetch origin"
            echo "   git merge ${COMMIT_ID}"
            echo "Now recompile OpenFOAM"
            echo "   ./Allwmake"
        fi
    done
}


__printMMCBanner()
{
    echo "==============================================================================="
    echo "                                       8888888888                              "
    echo "                                       888                                     "
    echo "                                       888                                     "
    echo "  88888b.d88b.  88888b.d88b.   .d8888b 8888888  .d88b.   8888b.  88888b.d88b.  "
    echo "  888 \"888 \"88b 888 \"888 \"88b d88P\"    888     d88\"\"88b     \"88b 888 \"888 \"88b "
    echo "  888  888  888 888  888  888 888      888     888  888 .d888888 888  888  888 "
    echo "  888  888  888 888  888  888 Y88b.    888     Y88..88P 888  888 888  888  888 "
    echo "  888  888  888 888  888  888  \"Y8888P 888      \"Y88P\"  \"Y888888 888  888  888 "
    echo "==============================================================================="
    echo "                                                                               "
}

__checkCPU()
{
    nProcs=$(nproc)
    nProcs=$(echo $nProcs | sed -e 's/ //g' )
    nProcs=$(($nProcs))

    if [[ $1 -gt $nProcs ]]; then
        echo "0"
        return 0
    fi
    
    echo "1"
    return 0
}


__checkGitSubmodule()
{
    if [[ ! -e ${MMC_PROJECT_USER_DIR}/thirdParty/movingAverage/.git ]]; then
    echo "${MMC_PROJECT_USER_DIR}/thirdParty/movingAverage/.git"
        >&2 echo "Git submodule of the moving average library is not available! Use:"
        >&2 echo "    git submodule update --init"
    fi
    if [[ ! -e ${MMC_PROJECT_USER_DIR}/tests/Catch2/.git ]]; then
        >&2 echo "Git submodule of the Catch2 library is not available! Use:"
        >&2 echo "    git submodule update --init"
    fi
}
