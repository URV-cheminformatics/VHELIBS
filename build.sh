OLDDIR=`pwd`
rm *py.class
CLASSPATH=one-jar/lib/Jmol.jar jython -c 'import visualitzador'
cd one-jar
ant clean
ant
cd $OLDDIR
