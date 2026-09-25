# Diagnostico del rescoring (19 sep 2026): quien esta vivo, quien mantiene el
# PC despierto, y a que hora apaga la tarea de Fredy.
$procs = Get-CimInstance Win32_Process -Filter "Name='pythonw.exe' OR Name='python.exe'"
$vig = $procs | Where-Object { $_.CommandLine -like '*vigilar_rescoring_focalizada*' }
$des = $procs | Where-Object { $_.CommandLine -like '*mantener_despierto*' }
$run = $procs | Where-Object { $_.CommandLine -like '*runner_mmgbsa*' }

Write-Output ("Vigilantes vivos : {0}" -f @($vig).Count)
foreach ($p in $vig) { Write-Output ("   PID {0} arrancado {1}" -f $p.ProcessId, $p.CreationDate) }
Write-Output ("Mantener despierto: {0}" -f @($des).Count)
foreach ($p in $des) { Write-Output ("   PID {0} arrancado {1}" -f $p.ProcessId, $p.CreationDate) }
Write-Output ("Runner en Windows : {0}  (el de verdad corre dentro de WSL)" -f @($run).Count)

Write-Output ""
Write-Output "Tarea de apagado ApagarPC_Fredy:"
$t = Get-ScheduledTask -TaskName 'ApagarPC_Fredy' -ErrorAction SilentlyContinue
if ($t) {
    Write-Output ("   Estado: {0}" -f $t.State)
    $t.Triggers | ForEach-Object { Write-Output ("   Disparador: {0}" -f $_.StartBoundary) }
    Write-Output ("   Habilitada: {0}" -f $t.Settings.Enabled)
} else {
    Write-Output "   No existe."
}
