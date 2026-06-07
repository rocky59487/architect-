# One-click verification gate for the FrameSolver engine (FrameCore-only, post-cleanup).
#   [1/3] standalone FrameCore fixtures (analytic golden oracles)
#   [2/3] UE headless automation (FrameCore.*)
#   [3/3] OpenSees offline cross-validation (skipped if openseespy absent)
# Prints a combined PASS/FAIL summary and sets the exit code (0 = all green).
#
# Usage:  powershell -ExecutionPolicy Bypass -File E:\project\ArchSim\Scripts\run_gate.ps1
param([switch]$RequireOpenSees)   # CI: fail (not skip) when openseespy is absent
$ErrorActionPreference = 'Continue'

$Root   = 'E:\project\ArchSim'
$Engine = 'E:\project\UE_5.7'
$UProj  = Join-Path $Root 'ArchSim.uproject'
$UeCmd  = Join-Path $Engine 'Engine\Binaries\Win64\UnrealEditor-Cmd.exe'
$Log    = Join-Path $Root 'Saved\Logs\ArchSim.log'

Write-Host '======================================================'
Write-Host ' FrameSolver verification gate'
Write-Host '======================================================'

# ---- [1/3] standalone gate ----
Write-Host ''
Write-Host '[1/3] standalone FrameCore gate (build.bat)...'
& (Join-Path $Root 'Plugins\FrameSolver\Standalone\build.bat') | Tee-Object -Variable StandaloneOut | Out-Null
$StandaloneRC = $LASTEXITCODE
$StandaloneLine = ($StandaloneOut | Select-String -Pattern 'ALL PASS|FAILURES' | Select-Object -Last 1)
Write-Host ("       standalone: {0} (exit {1})" -f $StandaloneLine, $StandaloneRC)

# ---- [2/3] UE headless automation ----
Write-Host ''
Write-Host '[2/3] UE headless automation...'
$ExecCmds = 'Automation RunTests FrameCore; Quit'
& $UeCmd $UProj "-ExecCmds=$ExecCmds" -unattended -nullrhi -nopause -nosplash -log | Out-Null
$UeRC = $LASTEXITCODE

# ASCII-only parse (Windows PowerShell 5.1 reads .ps1 as ANSI, which would corrupt
# localized result strings): the automation runner prints "TEST COMPLETE. EXIT CODE: N"
# (N = failed-test count, 0 = all green) and one "Test Completed." per test.
$UeExit = 1; $Total = 0
if (Test-Path $Log) {
    $ec = Select-String -Path $Log -Pattern 'TEST COMPLETE\. EXIT CODE: (\d+)' | Select-Object -Last 1
    if ($ec) { $UeExit = [int]$ec.Matches[0].Groups[1].Value }
    $tot = (Select-String -Path $Log -Pattern 'Test Completed\.' -AllMatches |
            ForEach-Object { $_.Matches.Count } | Measure-Object -Sum).Sum
    if ($tot) { $Total = [int]$tot }
}
Write-Host ("       UE automation: {0} tests run, exit code {1} (process exit {2})" -f $Total, $UeExit, $UeRC)

# ---- [3/3] OpenSees offline cross-validation (#14; skipped if openseespy absent) ----
Write-Host ''
Write-Host '[3/3] OpenSees offline cross-validation...'
$OsRC = 0; $OsState = 'skipped'
& python (Join-Path $Root 'Tools\opensees_compare.py') | Tee-Object -Variable OsOut | Out-Null
$OsRC = $LASTEXITCODE
if ($OsRC -eq 0)     { $OsState = 'PASS' }
elseif ($OsRC -eq 2) { $OsState = 'skipped (openseespy not installed)' }
else                 { $OsState = 'FAIL' }
Write-Host ("       OpenSees compare: {0} (exit {1})" -f $OsState, $OsRC)

# ---- verdict ----
Write-Host ''
Write-Host '======================================================'
# OpenSees leg: a real diff (1) always fails. Absence (2) is a soft skip locally, but a hard
# fail under -RequireOpenSees (CI), so the cross-validation can never silently vanish.
if ($OsRC -eq 2) {
    if ($RequireOpenSees) { Write-Host '       OpenSees REQUIRED but openseespy absent -> FAIL' -ForegroundColor Red }
    else { Write-Host '       (OpenSees skipped: openseespy absent; pass -RequireOpenSees to enforce in CI)' -ForegroundColor Yellow }
}
$OsOk = if ($RequireOpenSees) { $OsRC -eq 0 } else { $OsRC -ne 1 }
$GateOk = ($StandaloneRC -eq 0) -and ($UeExit -eq 0) -and ($Total -gt 0) -and $OsOk
if ($GateOk) {
    Write-Host (" GATE: PASS  (standalone OK, UE {0} tests green, OpenSees {1})" -f $Total, $OsState) -ForegroundColor Green
    exit 0
} else {
    Write-Host (" GATE: FAIL  (standalone exit {0}, UE exit {1}, {2} tests, OpenSees {3})" -f $StandaloneRC, $UeExit, $Total, $OsState) -ForegroundColor Red
    exit 1
}
