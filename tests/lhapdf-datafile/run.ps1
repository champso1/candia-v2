rm -Recurse .\testpdf
cmake --build ../.. --config Debug --target lhapdf-datafile
.\Debug\lhapdf-datafile.exe
