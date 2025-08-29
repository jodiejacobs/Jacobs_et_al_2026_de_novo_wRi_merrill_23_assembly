#!/bin/sh



#  Path to Canu.

bin="/private/home/jomojaco/mambaforge/envs/improve_assembly/bin"

#  Report paths.

echo ""
echo "Found perl:"
echo "  " `which perl`
echo "  " `perl --version | grep version`
echo ""
echo "Found java:"
echo "  " `which /private/home/jomojaco/mambaforge/envs/improve_assembly/lib/jvm/bin/java`
echo "  " `/private/home/jomojaco/mambaforge/envs/improve_assembly/lib/jvm/bin/java -showversion 2>&1 | head -n 1`
echo ""
echo "Found canu:"
echo "  " $bin/canu
echo "  " `$bin/canu -version`
echo ""


#  Environment for any object storage.

export CANU_OBJECT_STORE_CLIENT=
export CANU_OBJECT_STORE_CLIENT_UA=
export CANU_OBJECT_STORE_CLIENT_DA=
export CANU_OBJECT_STORE_NAMESPACE=
export CANU_OBJECT_STORE_PROJECT=




rm -f canu.out
ln -s canu-scripts/canu.01.out canu.out

/usr/bin/env perl \
$bin/canu -p 'canu' 'genomeSize=1.3m' -raw -nanopore '/private/groups/russelllab/jodie/merrill_23_wRi_genome/flye/merge_bootcamp_data/assembly.fasta' canuIteration=0
