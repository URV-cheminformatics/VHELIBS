OLDDIR=`pwd`
rm *py.class
CLASSPATH=one-jar/lib/Jmol-slim.jar jython -c 'import visualitzador'
cd one-jar
ant clean
ant
cd $OLDDIR
