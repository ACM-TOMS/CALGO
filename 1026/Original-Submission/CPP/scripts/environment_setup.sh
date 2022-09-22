#! /bin/bash

SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")
CALS_REPO_PATH="${SCRIPTPATH}/.."

pushd . > /dev/null
echo -n "Checking ${CALS_REPO_PATH} for extern folder... "
cd ${CALS_REPO_PATH}/extern || (echo " " && echo "Can't find 'extern' folder. Exiting." && popd && exit)

echo "Found."

if [ ! -d "cmake" ]; then
  echo -n "Previous installation not found. Installing cmake to `pwd` ... "
  mkdir cmake
  wget https://github.com/Kitware/CMake/releases/download/v3.17.5/cmake-3.17.5-Linux-x86_64.tar.gz -qO- | tar xzf - -C cmake
  echo "Done"
fi

CMAKE_BIN_PATH=`pwd`/cmake/cmake-3.17.5-Linux-x86_64/bin
echo
echo "Copy the following line to your ${HOME}/.bashrc or ${HOME}/.zshrc to be able to find CMake 3.17.5 in a new terminal session:"
echo "export PATH=${CMAKE_BIN_PATH}:\${PATH}"
echo
echo -n "Adding ${CMAKE_BIN_PATH} to PATH variable... "
export PATH=${CMAKE_BIN_PATH}:${PATH}
echo "Done."

popd > /dev/null




