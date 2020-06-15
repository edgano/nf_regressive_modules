# NF Regressive Modules

![Regressive Modules CI](https://github.com/edgano/nf_regressive_modules/workflows/nf_reg_modules%20CI/badge.svg)
![Tcoffee CI](https://github.com/edgano/nf_regressive_modules/workflows/tcoffee%20CI/badge.svg)

This repository aims to group all the new modes of the Regressive algorithm, presented in the manuscript "Fast and accurate large multiple sequence alignments using root-to-leave regressive computation".


#### For details on how to use the Regressive Multiple Sequence Alignment method, see the [T-Coffee documentation](https://tcoffee.readthedocs.io/en/latest/tcoffee_quickstart_regressive.html).


### Credits
This workflow was written by Edgar Garriga ([edgano](https://github.com/edgano)) and 
Jose Espinosa([JoseEspinosa](https://github.com/JoseEspinosa)) at the [Center for Genomic Regulation (CRG)](http://www.crg.eu).

The authors who contributed to the analysis and manuscript are:

* Edgar Garriga Nogales
* Jose Espinosa Carrasco
* Paolo Di Tommaso
* Cedric Notredame

### Pipeline
The pipeline for generating trees, alignments and performing the evaluations is built using 
[Nextflow](https://www.nextflow.io), a workflow tool to run tasks across 
multiple compute infrastructures in a very portable manner. It comes with a docker container 
making installation trivial and results highly reproducible.

### Pipeline Quick Start
Make sure you have either docker/singularity installed or the required dependencies listed 
in the last section.

Install the Nextflow runtime by running the following command:

    $ curl -fsSL get.nextflow.io | bash


When done, you can launch the pipeline execution by entering the command shown below:

    $ nextflow run edgano/nf_regressive_modules
    

By default the pipeline is executed against the provided example dataset. 
Check the *Pipeline parameters*  section below to see how enter your data on the program 
command line.     
  

### Containers

All the methods above are available in multiple [Docker](http://www.docker.com) images on DockerHub and the image is tested to be compatible with the [Singularity](http://singularity.lbl.gov/).

### Pipeline parameters

#### `--seqs` 
   
* Specifies the location of the input *fasta* file(s).
* Multiple files can be specified using the usual wildcards (*, ?), in this case make sure to surround the parameter string
  value by single quote characters (see the example below)

Example: 

    $ nextflow run edgano/nf_regressive_modules --seqs '/home/seqs/*.fasta'

This will handle each fasta file as a seperate sample.


#### `--refs` 

* Specifies the location of the reference *aligned fasta* file(s).


#### `--trees` 

* Specifies the location of input tree file(s).


#### `--align_methods` 

* Specifies which alignment methods should be used.
* Options include: "CLUSTALO,MAFFT-FFTNS1,MAFFT-SPARSECORE,MAFFT-GINSI,PROBCONS,UPP"


#### `--tree_methods` 

* Specifies which guide-tree / clustering methods should be used.
* Options include: "MBED,MAFFT-PARTTREE,FAMSA-SLINK,MAFFT-DPPARTTREE0"
* Default locations is `results/trees` directory.

#### `--buckets` 

* List of bucket sizes or maximum size of the subMSAs in the regressive proceedure.
* Default value is "1000" sequences.


#### `--regressive_align` 

* Flag to generate standard regressive MSAs.
* See `modules/templates/regressive_align` for the specific commands executed.


#### `--progressive_align` 

* Flag to perform MSAs in a progressive way.
* Progressive MSA is alignment where the guide-tree is provided as input.
* See `modules/templates/progressive_align` for the specific commands executed.


#### `--slave_align` 

* Flag to generate regressive MSAs with the capability to decide the guide-tree for the child subMSAs.
* Regressive MSA is alignment where the `-slave_tree_methods` should be provided.
* See `modules/templates/slave_align` for the specific commands executed.


#### `--dynamic_align` 

* Flag to generate regressive MSAs running different alignments methods concerning the deepth in the guide-tree.
* Regressive MSA is alignment where the `-dynamicSize` should be provided and a database `-db` (running PDB as default, included in the docker image)
* See `modules/templates/dynamic_align` for the specific commands executed.


#### `--pool_align` 

* Flag to generate regressive MSAs running with the capability to merge subMSAs with small number of sequences.
* See `modules/templates/pool_align` for the specific commands executed.


#### `--evaluate` 

* Flag to perform evaluation of the alignments, it includes TC, SoP and Col scores.
* Requires reference sequences to be provided with the `--refs` parameter.
* Default locations is `results/score` directory.

#### `--homoplasy` 

* Flag to perform generate the homoplasy information of the alignments.
* Default locations is `results/homoplasy` directory.

#### `--gapCount` 

* Flag to count the number of gaps on a progressive alignments.
* Default locations is `results/gaps` directory.

#### `--metrics` 

* Flag to generate the metrics (cpu/memory) of the run.
* Default locations is `results/metrics` directory.

#### `--easel` 

* Flag to get structural information (%identity) of the alignments.
* Default locations is `results/easel` directory.

#### `--output`

* Location of the results.
* Default locations is `results` directory.
