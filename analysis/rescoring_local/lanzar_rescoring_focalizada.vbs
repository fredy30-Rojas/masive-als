' Arranca el rescoring MM-GBSA de la lista focalizada sin ninguna ventana.
' El False del final es a proposito: no espera a que acabe. El rescoring son
' unas 40 horas, asi que este lanzador solo lo suelta y se va.
' Registro: rescoring_local\rescoring_focalizada_stdout.log
Set sh = CreateObject("WScript.Shell")
sh.Run """C:\Users\Fredy\masive-als\analysis\rescoring_local\rescoring_focalizada.bat""", 0, False
