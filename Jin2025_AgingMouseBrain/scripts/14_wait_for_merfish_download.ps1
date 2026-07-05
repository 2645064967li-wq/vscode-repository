$ErrorActionPreference = Stop
$root = Split-Path -Parent $PSScriptRoot
$logDir = Join-Path $root results\spatial_region_feasibility
$log = Join-Path $logDir download_and_extract.log
while (Get-Process curl -ErrorAction SilentlyContinue) {
    Start-Sleep -Seconds 20
}
Set-Location $root
python scripts\13_extract_target_merfish_expression.py *>&1 |
    Out-File -LiteralPath $log -Encoding utf8
