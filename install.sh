#!/bin/bash

VER=v1.2.dev

PKG=cross-$VER-$1

TOP=$(pwd)

rm -rf $PKG
mkdir $PKG

make clean

if [[ $1 == "glnx64" ]]; then

    make || exit 1
    cp cross $TOP/$PKG/

    cd src
    make distclean
    qmake-qt5 gui || exit 1
    make || exit 1
    # sudo apt install libqt5widgets5
    # sudo yum install qt5-qtbase-gui
    strip -s cross-gui
    cp cross-gui $TOP/$PKG/

elif [[ $1 == "win32" ]]; then

    make win32 || exit 1
    cp cross.exe $TOP/$PKG/

    cd src
    make distclean
    mingw32-qmake-qt4 "CONFIG += static" gui || exit 1
    make release || exit 1
    mingw-strip -s release/cross-gui.exe
    cp release/cross-gui.exe $TOP/$PKG/

elif [[ $1 == "win64" ]]; then

    make win64 || exit 1
    cp cross.exe $TOP/$PKG/

    cd src
    make distclean
    mingw64-qmake-qt4 "CONFIG += static" gui || exit 1
    make release || exit 1
    mingw-strip -s release/cross-gui.exe
    cp release/cross-gui.exe $TOP/$PKG/

elif [[ $1 == "macos" ]]; then

    export PATH="/usr/local/opt/llvm@7/bin:$PATH"
    make macos || exit 1
    cp cross $TOP/$PKG/

    export PATH="/usr/local/opt/qt/bin:$PATH"
    export LDFLAGS="-L/usr/local/opt/qt/lib"
    export CPPFLAGS="-I/usr/local/opt/qt/include"

    cd src
    make distclean
    qmake gui || exit 1
    make || exit 1
    macdeployqt cross-gui.app
    mv cross-gui.app $TOP/$PKG/

else

    exit 1

fi

if [[ $1 == win* ]]; then
    zip -qr $PKG.zip $PKG
else
    tar zcf $PKG.tar.gz $PKG
fi
