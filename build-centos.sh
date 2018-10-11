#!/bin/bash

rm -rf cross-$1
mkdir cross-$1
make distclean

if [ $1 == "lnx64" ]; then

    g++ *.cpp -o cross-$1/cross -s -O2 -std=c++11 -static
    qmake-qt4 gui
    make
    strip cross-gui
    mv cross-gui cross-$1/

elif [ $1 == "win32" ]; then

    i686-w64-mingw32-g++ *.cpp -o cross-$1/cross.exe -s -O2 -std=c++11 -static
    i686-w64-mingw32-qmake-qt4 gui
    make release
    i686-w64-mingw32-strip release/cross-gui.exe
    mv release/cross-gui.exe cross-$1/

elif [ $1 == "win64" ]; then

    x86_64-w64-mingw32-g++ *.cpp -o cross-$1/cross.exe -s -O2 -std=c++11 -static
    x86_64-w64-mingw32-qmake-qt4 gui
    make release
    x86_64-w64-mingw32-strip release/cross-gui.exe
    mv release/cross-gui.exe cross-$1/

fi

tar zcf cross-$1.tar.gz cross-$1
