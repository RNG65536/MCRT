rmdir /s /q _build
mkdir _build
cd _build
REM cmake .. -G "Visual Studio 15 Win64" -Dembree_DIR="C:/Deps/embree/lib/cmake/embree"
cmake .. -G "Visual Studio 15 Win64"
cd ..