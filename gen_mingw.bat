rmdir /s /q _build
mkdir _build
cd _build
cmake .. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
cd ..