CLASSPATH=one-jar/src/pylib/:one-jar/lib/Jmol.jar
java -jar one-jar/lib/jython.jar -S -i -c'import sys; sys.path.append("one-jar/src/pylib/"); sys.path.append("one-jar/lib/Jmol.jar"); sys.argv.remove("-c"); import visualitzador; sv = visualitzador.main()' $@
