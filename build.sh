OLDDIR=`pwd`
rm *.class
javac PdbAtomJava.java
version=$(java -Dpython.import.site=false -jar one-jar/lib/jython.jar -c'import sys; sys.path.append("one-jar/src/pylib/"); sys.path.append("one-jar/lib/Jmol.jar"); import visualitzador; sys.stdout.write(visualitzador.VHELIBS_VERSION)')
cp -vr *.class *.png swingutils pdbx one-jar/src/
for d in pdbx swingutils; do
    rm one-jar/src/$d/*.class
    rm -r one-jar/src/$d/__pycache__
    rm one-jar/src/$d/*/*.class
    rm -r one-jar/src/$d/*/__pycache__
done
cp -v Main.java  one-jar/java-src/
cd one-jar
ant clean
ant
mv VHELIBS.jar VHELIBS-$version.jar
cd "$OLDDIR"

#python setup.py build_ext --inplace
