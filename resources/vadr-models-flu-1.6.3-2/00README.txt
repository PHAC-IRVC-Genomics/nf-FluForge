March 2024
vadr-models-flu-1.6.3-2

https://github.com/ncbi/vadr

VADR documentation:
https://github.com/ncbi/vadr/blob/master/README.md

Model download site:
https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/flu

See 00RELEASE-NOTES.txt for details on changes between model versions. 
See 00NOTES.txt for additional information on the models.

------------

Recommended command for flu annotation using vadr v1.6.3
(as of December 18, 2023) but still under testing and subject to 
change.

v-annotate.pl \
--split --cpu 4 -r \
--atgonly --xnocomp --nomisc \
--alt_fail extrant5,extrant3 \
--mkey flu \
--mdir <PATH-TO-THIS-MODEL-DIR> \
<fastafile> <outputdir>

The '--split --cpu 4' options will run v-annotate.pl multi-threaded on
4 threads. To change to '<n>' threads use '--split --cpu <n>', but
make sure you have <n> * 4G total RAM available. 

To run single threaded remove the '--split --cpu 4' options.

--

Contact eric.nawrocki@nih.gov for help.
