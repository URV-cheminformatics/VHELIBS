OLDDIR=`pwd`
rm *.class
javac PdbAtomJava.java
java -jar one-jar/lib/jython.jar -c'import sys; sys.path.append("one-jar/src/pylib/"); sys.path.append("one-jar/lib/Jmol.jar"); import visualitzador'
cp -vr *.class *.png swingutils one-jar/src/
cp -v Main.java  one-jar/java-src/
cd one-jar
ant clean
ant
cd $OLDDIR

python setup.py build_ext --inplace
