' Arranca el vigilante del rescoring focalizada sin ninguna ventana.
' El vigilante mira el CSV cada 5 minutos y avisa por Telegram y por voz UNA
' sola vez: cuando terminan las 179 moleculas, o cuando el runner se para.
' Registro: rescoring_local\vigilante_rescoring.log
' Estado:   rescoring_local\vigilante_rescoring_estado.json
Set fso = CreateObject("Scripting.FileSystemObject")
Set sh = CreateObject("WScript.Shell")
base = "C:\Users\Fredy\masive-als\analysis\rescoring_local\"
py = "C:\Users\Fredy\AppData\Local\Python\pythoncore-3.14-64\pythonw.exe"
If Not fso.FileExists(py) Then py = "pythonw.exe"
sh.CurrentDirectory = base
sh.Run """" & py & """ """ & base & "vigilar_rescoring_focalizada.py""", 0, False
