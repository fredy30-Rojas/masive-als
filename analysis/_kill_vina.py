import subprocess, re
out = subprocess.run(['powershell.exe','-NoProfile','-Command','Get-CimInstance Win32_Process | Select-Object ProcessId,CommandLine | Format-List'], capture_output=True, text=True)
procs = out.stdout.split('\n\n')
killed = 0
for p in procs:
    if 'vina' in p.lower() or 'validar_senuelos' in p.lower():
        m = re.search(r'ProcessId\s*:\s*(\d+)', p)
        if m:
            pid = m.group(1)
            r = subprocess.run(['taskkill','/PID',pid,'/F','/T'], capture_output=True, text=True)
            print('killed', pid, 'rc=', r.returncode)
            killed += 1
print('total killed:', killed)
