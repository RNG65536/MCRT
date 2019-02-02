rmdir /s /q _build
mkdir _build
cd _build
cmake .. -G "Visual Studio 15 Win64" -Dembree_DIR="C:/Deps/embree/lib/cmake/embree"
REM cmake .. -G "Visual Studio 15 Win64"
cd ..