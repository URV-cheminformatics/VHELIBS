OLDDIR=`pwd`
rm *py.class
java -jar one-jar/lib/jython.jar -c'import sys; sys.path.append("one-jar/src/pylib/"); sys.path.append("one-jar/lib/Jmol.jar"); import visualitzador'
cd one-jar
ant clean
ant
cd $OLDDIR
