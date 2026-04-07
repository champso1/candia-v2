param(
	[Parameter(Mandatory)]
	[string]$buildtype
)

$testpdf_folder = ".\testpdf"
if (Test-Path $testpdf_folder) {
	Remove-Item -Recurse -Force $testpdf_folder
}
cmake --build ../.. --config $buildtype --target lhapdf-datafile
& ".\$buildtype\lhapdf-datafile.exe"
