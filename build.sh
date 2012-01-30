cp *.py one-jar/src/
cp Main.java one-jar/java-src/
cd one-jar
ant clean
ant
rm src/*
rm java-src/*
cd ..
