import subprocess
out = subprocess.run(['powershell.exe','-NoProfile','-Command','Get-CimInstance Win32_Process -Filter "Name like \'python%\'" | Select-Object ProcessId,CommandLine | Format-List'], capture_output=True, text=True)
print(out.stdout[:4000])
