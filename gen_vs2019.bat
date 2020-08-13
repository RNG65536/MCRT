rmdir /s /q _build
mkdir _build
cd _build
cmake .. -G "Visual Studio 16" -A "x64" -DCMAKE_PREFIX_PATH="C:/Deps/embree/lib/cmake/embree"
cd ..