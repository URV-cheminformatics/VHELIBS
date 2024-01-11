#!/usr/bin/sh
java -Dpython.import.site=false -jar one-jar/lib/jython.jar -c'import sys; sys.path.append("one-jar/src/pylib/"); sys.path.append("one-jar/lib/Jmol.jar"); import visualitzador; visualitzador.main([])'
