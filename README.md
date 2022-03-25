# Fourth Year Capstone Project Data Handling Scripts
Collection of data handling scripts used for a capstone research project.

## Benchmarking
The JSON files generated for the benchmarking process by ExpansionHunter are archived in a .tar.gz file.
You'll need to extract the files and move them to the appropriate directory

~~~
tar -xvzf exphunterbenchmarkdata.tar.gz --directory benchmark_data/
~~~

## HGDP RE analysis
The JSON files containing the inferred repeat expansion reads for the analysed genomes are archived.
To extract them:
~~~
tar -xvzf c9orf72_exphunter.tar.gz fmr1_exphunter.tar.gz nonrarevariants_exphunter.tar.gz
~~~
