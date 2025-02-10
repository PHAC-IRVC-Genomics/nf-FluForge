#!/bin/bash
# make-tarball.sh
# A shell script for building a tarball of vadr models, by:
#  - cloning the corresponding branch of the vadr model git repo
#  - switching to the desired branch
#  - save current short git commit hash and other info 
#    to build-info.txt
#  - removing .git info
#  - unpacking gzipped .cm* and .hmm* files
#  - removing unwanted files (gunzip.sh, gzip.sh, make-tarball.sh)
#  - tarring and gzipping entire dir
# 
# usage: 
# make-tarball.sh <branch_name (e.g. \"develop\"> <directory_name (e.g. \"vadr-models-sarscov2-1.2-1dev0\")>"
#
# The following line will make the script fail if any commands fail
set -e

BBREPO="https://nawrockie@bitbucket.org/nawrockie/vadr-models-flu.git"

########################
# Validate correct usage
########################
# make sure correct number of cmdline arguments were used, exit if not
if [ "$#" -ne 2 ]; then
   echo "Usage: $0 <branch_name (e.g. \"develop\"> <directory_name (e.g. \"vadr-models-flu-1.5.1-1dev1\")>"
   exit 1
fi
BRANCH=$1;
DIRNAME=$2;

# if the directory already exists, exit
if [ -d $DIRNAME ]; then 
    echo "directory $DIRNAME already exists, remove it first"
    exit 1
fi

#  - cloning the corresponding branch of the vadr model git repo
git clone $BBREPO $DIRNAME
cd $DIRNAME

#  - switching to the desired branch
git checkout $BRANCH

#  - save current short git commit hash and other info 
#    to build-info.txt
echo -n "current git commit hash: " >  build-info.txt
git rev-parse --short HEAD >> build-info.txt
echo "git repo:    " $BBREPO >> build-info.txt
echo "git branch:  " $BRANCH >> build-info.txt
echo "dir name:    " $DIRNAME >> build-info.txt
echo "packaged on: `date`" >> build-info.txt

#  - removing .git info
rm -rf .git

#  - unpacking gzipped .cm* and .hmm* files
##sh ./gunzip.sh

#  - removing unwanted files (gunzip.sh, gzip.sh, make-tarball.sh)
##rm make-tarball.sh
##rm gzip.sh
##rm gunzip.sh
rm 00NOTES.txt

#  - tarring and gzipping entire dir
cd ..
tar -cvf $DIRNAME.tar $DIRNAME
gzip $DIRNAME.tar

