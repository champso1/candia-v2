$DeleteOriginals = $false
if ($args -contains "--remove-files") {
    $DeleteOriginals = $true
    Write-Host "Will delete original files." -ForegroundColor Cyan
}

$DestDir = "backups"
if (-not (Test-Path $DestDir)) {
    New-Item -ItemType Directory -Path $DestDir | Out-Null
}

$Timestamp = Get-Date -Format "yyyy-MM-dd_HH-mm-ss"
$FileName = "back-$Timestamp.zip"
$FilePath = Join-Path $DestDir $FileName
$Targets = "tables", "diffs", "log", "latex", "data"

Write-Host "Compressing into $FilePath..."

try {
    # Check if targets exist to avoid errors
    $ValidTargets = $Targets | Where-Object { Test-Path $_ }

    if ($ValidTargets.Count -eq 0) {
        throw "None of the target folders/files were found."
    }

    # Compress files
    Compress-Archive -Path $ValidTargets -DestinationPath $FilePath -ErrorAction Stop
    Write-Host "Success!" -ForegroundColor Green

    # 4. Optional cleanup
    if ($DeleteOriginals) {
        Write-Host "Removing original files..." -ForegroundColor Yellow
        Remove-Item -Path $ValidTargets -Recurse -Force
        Write-Host "Done."
    }
}
catch {
    Write-Host "Failed to archive files: $($_.Exception.Message)" -ForegroundColor Red
    exit 1
}
