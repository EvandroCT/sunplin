## Synopsis

This software aims to provide massive parallel processing to generate multiple variants of a tree with phylogenetic uncertainty and calculate a patristic distance matrix for each version generated. This is done using CUDA, wich enable the use of NVIDIA GPU to speedup these operations wich are inherently time-consuming.

## Installation

Before installing this software, the **Cuda Toolkit  (at least 7.5)** must be already installed on your system. If haven't yet, instructions can be found [here](http://docs.NVIDIA.com/cuda/index.html).

At the root of the source directory, run the following command:

>$ make

If no errors occur, an executable file called 'sunplin' should be created.

## Usage Example

**Important: It's necessary to have a functioning NVIDIA GPU with at least a Compute Capability 3.0. You can find out this information about your device [here](https://developer.NVIDIA.com/cuda-gpus).**

Supose you have a newick file called 'newick.tree' at the same directory as the binary file created at compilation containing a phylogenetic tree with uncertainty (with absent species); also,  another file called 'put.list' wich has in each line a PUT (Phylogenetic Uncertainty Taxon) followed by a MDCC (Most Derived Consensus Clade),  separated by space.

The following command will generate 1000 diferent versions of the original tree (newick.tree) each one ramdomly expanded with the taxa contained in each line within the list of PUTs (put.list):

>$ ./sunplin 1000


By running the following command you will be able to pick up nonstandard input files:

>$ ./sunplin 1000 /path/to/newick\_file /path/to/list\_of\_puts\_file

Either commands above will prompt a message asking whether you wish to store the generated trees into a standard file called 'versions.tree'. Be aware that this operation could take several time.

## Motivation

Phylogenetic comparative analyses usually rely on a single consensus phylogenetic tree in order to
study evolutionary processes. However, most phylogenetic trees are incomplete with regard to species sampling, which may critically compromise analyses. Some approaches have been proposed to integrate non-molecular phylogenetic information into incomplete molecular phylogenies. An expanded tree approach consists of adding missing species to random locations within their clade. The information contained in the topology of the resulting expanded trees can be captured by the pairwise phylogenetic distance (patristic distance)  between species and stored in a matrix for further statistical analysis. Thus, the random expansion and processing of multiple phylogenetic trees can be used to estimate the phylogenetic uncertainty through a simulation procedure. Because of the computational burden required, unless this procedure is efficiently implemented, the analyses are of limited applicability.


## Tests

Using a NVIDIA GeForce GTX TITAN X graphics card and the standard example input files provided by this repository, we were capable of seamlessly running a simulation with up to 3000 replications. 

This software has been conceived with academic purpose, thus, tests was made within a very restricted range of environments, therefore, more tests have be done. All the comunity is welcome to do your own tests and provide us your feedback.

## Contributors

Evandro Taquary

Thiago Santos

Wellington Martins

## License

This software is available under GNU GPLv3 license.
