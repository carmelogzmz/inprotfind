# inprotfind    

A workflow for insect protein identification and phylogenetic analyses.

- Documentation
- GitHub
- PyPI
- License (BSD)


## Features

- **Cross-Platform**: Supports Linux (AMD and Aarch) and MacOS (Intel and Silicon)
- **Versatile**: Compatible with Python 3.6 to 3.12
- **Flexible**: Can be used both within Python scripts and directly from the terminal through shell scripting




## Prerequisites

**inprotfind** requires _mmseqs2_, _mafft_ and _fasttree_ installed in the same environment to work properly. The three tools are easily installable using [_conda_](https://docs.anaconda.com/miniconda/) and it is the way I recommend.

#### [_mmseqs2_](https://github.com/soedinglab/MMseqs2?tab=readme-ov-file)

```
conda install bioconda::mmseqs2
```

#### [_mafft_](https://mafft.cbrc.jp/alignment/server/index.html)

```
conda install bioconda::mafft
```

#### [_FastTree_](http://www.microbesonline.org/fasttree/)

```
conda install bioconda::fasttree
```


## Installing inprotfind

Install **inprotfind** with pip

```bash
  pip install inprotfind
```
    
## Using inprotfind

**inprotfind** includes a series of functions that allow you to compare a query protein sequence against a database containing over four million insect protein sequences (using the _find_matches_ function) to identify the potential type of protein and its possible species, genus, or family. Additionally, after aligning the most similar sequences (_align_sequences_ function), it can place the query sequence in a phylogenetic context alongside the most similar sequences from the database (_build_tree_ function).

### main functions

The following three functions allow the user to go through the complete workflow, from the query protein sequence to the tree where the query sequence is located in a phylogenetic context:


#### find\_matches(job\_name, query\_path, db\_name="arthropods\_NCBI")

It requires three arguments:

- `job_name` requires a name, chosen entirely by the user, to identify the search or task to be performed (e.g., Task01, BeeProtein, ProjectArboVirus, etc.). It is recommend that the user choose an informative name, to easy identification. The function creates a folder (in the working directory) with the same name where all the files created during the execution are stored. 
- `query_path` requires the path where the file with the query protein sequence (a **fasta file**).
- `db_name` requires the name of the target database to use for comparing the query protein sequence. If no name is passed. The function uses the default dabased (*"arthropods_NCBI"*), included in **inprotfind**. In case a name is passed by the user, the function search the database with the same name in the working directory within a folder named _databases_. NOTE: At the moment, the implementation to use private databases is not finished, so I recommend to leave the default value for `db_name` as otherwise the function will not work properly.

As an example, if the user wants to name the task as _HoverflyProtein01_, and the fasta file with the query protein sequence is in */home/{USER}/HoverflyProject/queries/query\_Hoverfly01.fasta* (and using the default database) The code will be:

_Python:_

```python
import inprotfind as ipf

ipf.find_matches("HoverflyProtein01", "home/USER/HoverflyProject/queries/query_Hoverfly01.fasta")
```

_Shell:_

```bash
inprotfind find_matches --job_name HoverflyProtein01 --query_path home/USER/HoverflyProject/queries/query_Hoverfly01.fasta
```

Once the function is executed, tt creates a folder named "HoverflyProtein01" (in the working directory) where all the files created during the execution are stored. The main result files are _best\_matches\_all.m8_, with the metadata of the 100 more similar sequences from the target database, and _best\_matches.m8_ with just the metadata of just the 30 first sequences included in _best\_matches\_all.m8_. 

#### align\_sequences(job\_name)

It only requires the `job_name` argument to select the files on which to apply the function. 

Following the same example from the previous function:

_Python:_
```python
import inprotfind as ipf

ipf.align_sequences("HoverflyProtein01")
```

_Shell:_
```bash
inprotfind align_sequences --job_name HoverflyProtein01
```

Once the function is executed, it collects the sequences which accession code is in _best\_matches.m8_, and add the query protein sequence in a multifasta file (_filtered\_sequences.fasta_), and then aligns the sequences using MAFFT (_aligned\_sequences.fasta_). Both _filtered\_sequences.fasta_ and _aligned\_sequences.fasta_ are stored in the "HoverflyProtein01" folder.

#### build\_tree(job\_name, tree_type = "simple")

It requires two arguments:

- `job_name` is required to find the _aligned\_sequences.fasta_ file. 
- `tree_type` can be `"simple"`, `"interactive"` or `"ascii"`. `"simple"`is selected by default, and plot a simple and static phylogenetic tree in a pop-up window (or in the Plots panel if using an IDE with this feature). `"interactive"`use the _ete3_ library to open the ete3 software and plot an interactive phylogenetic tree with different visualization options. `"ascii"` shows the phylogenetic tree directly drawn in the python console or in the terminal.

Finishing with the example, if the user selects `"interactive"`:

_Python:_

```python
import inprotfind as ipf

ipf.build_tree("HoverflyProtein01", "interactive")
```

_Shell:_

```bash
inprotfind build_tree --job_name HoverflyProtein01 --tree_type interactive
```

Once the function is exectuted, it creates a file named _tree.nwk_ (in the folder "HoverflyProtein01"). Then open the _ete3_ interface and plot the interactive phylogenetic tree.

### complementary functions

There are few secondary functions that run within the main functions to secure that all the dependencies and tools are correctly installed. 

- **verifying_mmseqs2()** checks if _MMseqs2_ is installed in the working environment.
- **verifying_mafft()** checks if _MAFFT_ is installed in the working environment.
- **verifying_fasttree()** checks if _FastTree_ is installed in the working environment.

If any of these tools is not installed, the execution stops and ask the user to install it.

## Workflow example

**inprotfind** includes ten query protein sequences as examples for the user to learn how to properly use the library. To use the examples, it is enough to use the word "example" followed by a number from 1 to 10 as argument for the `query_path`in the _find\_matches_ function. We now will use the example number 3 to run the complete workflow and get some results. We will name the task as "QueryExample3" and use the default database "arthropods_NCBI":

#### 1. Searching similar sequences in the target database

```python
import inprotfind as ipf

ipf.find_matches(job_name = "QueryExample3", query_path = "example3", db_name = "arthropods_NCBI")
```

 _Notice that when using one of the ten examples, we are not passing a path as argument to `query_path`but just a string saying "example1", "example2"... "example10". In this example we use "example3". The function recognizes that we are asking to work with an example and it will recover it from the directory where the library is installed._

Let's have a look to the content of the "QueryExample3" folder. It should look like this:

```
.
└── QueryExample3/
    ├── queryDB/
    │   ├── queryDB
    │   ├── queryDB.dbtype
    │   ├── queryDB.index
    │   ├── queryDB.lookup
    │   ├── queryDB.source
    │   ├── queryDB_h
    │   ├── queryDB_h.dbtype
    │   └── queryDB_h.index
    ├── resultDB/
    │   ├── resultDB
    │   ├── resultDB.dbtype
    │   └── resultDB.index
    ├── best_matches.m8
    ├── best_matches_all.m8
    └── db_name.txt
```

It contains the folders named _queryDB_ and _resultDB_ with the query sequence in mmseqs2 format and a database in mmseqs2 format with the 100 more similar sequences from the target database. The db_name.txt file just stores the name of the target database (this file is read for _align\_sequences_ function to know which database was used). Lastly, it also contains the files _best\_matches_all.m8_, with the metadata of the 100 similar sequences, and _best\_matches.m8_ with the first 30 more similar sequences. Let's see the first 10 lines of _best\_matches.m8_:

```
query_example03 XP_044591841.1  0.954   1306    59      0       5       1280    317     1622    0.0     2620    Cotesia glomerata       GCF_020080835.1 123269956       protein madd-4 isoform X1
query_example03 XP_014300707.1  0.845   1300    198     0       3       1280    316     1615    0.0     2269    Microplitis demolitor   GCF_000572035.2 LOC103578210    uncharacterized protein LOC103578210
query_example03 XP_034935882.1  0.727   1335    348     0       5       1280    289     1623    0.0     1954    Chelonus insularis      GCF_013357705.1 LOC118065002    protein madd-4
query_example03 XP_015112712.1  0.68    1275    408     0       5       1279    299     1573    0.0     1723    Diachasma alloeum       GCF_001412515.2 LOC107038239    protein madd-4 isoform X1
query_example03 XP_043590523.1  0.582   1276    531     0       5       1280    307     1578    0.0     1429    Bombus pyrosoma GCF_014825855.1 122571180       protein madd-4-like isoform X1
query_example03 XP_012247090.1  0.582   1276    531     0       5       1280    143     1413    0.0     1427    Bombus impatiens        GCF_000188095.3 LOC100745833    LOC100745833
query_example03 XP_033296944.1  0.581   1276    532     0       5       1280    304     1574    0.0     1425    Bombus bifarius GCF_011952205.1 LOC117204053    LOC117204053
query_example03 XP_043804322.1  0.582   1275    532     0       5       1279    302     1575    0.0     1425    Apis laboriosa  GCF_014066325.1 122721106       protein madd-4 isoform X1
query_example03 XP_031365701.1  0.581   1275    532     0       5       1279    331     1600    0.0     1425    Apis dorsata    GCF_000469605.1 LOC102680940    protein madd-4 isoform X1
query_example03 XP_033192796.1  0.58    1276    533     0       5       1280    304     1574    0.0     1423    Bombus vancouverensis nearcticus        GCF_011952275.1 LOC117158237    protein madd-4-like
```

The file follows a similar tabular BLAST file (`m8`) with few columns more with some metadate. The colums are:

- query id
- target id
- percentage of identical matches
- alignment length
- number of mismatches
- number of gap openings
- query start
- query end
- target start
- target end
- e-value
- bit score
- organism
- genome id
- gen id
- protein description

We see that the first match has a 95.4% of identical matches. The first protein is described as a madd-4 isoform X1 belonging to Cotesia glomerata. While we cannot say 100% sure that our query sequence belongs to this species, seeing the ten first species (all of then bees) we could say that is a bee, and likely from the _Cotesia_ genus (SPOILER: is _Cotesia rubecula_).

#### 2. Aligning sequences with MAFFT

```python
import inprotfind as ipf

ipf.align_sequences(job_name = "QueryExample3")
```
Let's see again the content of the "QueryExample3" folder:

```
.
└── QueryExample4/
    ├── filteredDB/
    │   ├── filteredDB
    │   ├── filteredDB.dbtype
    │   ├── filteredDB.index
    │   ├── filteredDB.lookup
    │   ├── filteredDB.source
    │   ├── filteredDB_h
    │   ├── filteredDB_h.dbtype
    │   └── filteredDB_h.index
    ├── queryDB/
    │   ├── queryDB
    │   ├── queryDB.dbtype
    │   ├── queryDB.index
    │   ├── queryDB.lookup
    │   ├── queryDB.source
    │   ├── queryDB_h
    │   ├── queryDB_h.dbtype
    │   └── queryDB_h.index
    ├── resultDB/
    │   ├── resultDB
    │   ├── resultDB.dbtype
    │   └── resultDB.index
    ├── aligned_sequences.fasta
    ├── best_matches.m8
    ├── best_matches_all.m8
    ├── db_name.txt
    └── filtered_sequences.fasta
```
It contains a new folder with the database in mmseqs2 format of the selected sequences (while resultDB contains the 100 best matches, filteredDB contains the 30 best). The file _filtered\_sequences.fasta_ store the protein sequences of the query and the 30 best matches, and the file _aligned\_sequences.fasta_ shows the alignment of this 31 sequences (30 best + query).

#### 3. Building the phylogenetic tree

```python
import inprotfind as ipf

ipf.build_tree(job_name = "QueryExample3", tree_type = "simple")
```
This function creates the _tree.nwk_ file from which the phylogenetic tree is drawn. We will not see again the directory tree as the only addition is _tree.nwk_ to the list of files. When the function finishes, the tree is drawn (in a pop-up window or in the Plots tab if using and IDE with this feature). It should look something like this:






Here, we can see that the query sequences is located in the tree together with the 30 more similar sequences.


#### 4. Checking the results

**inprotfind** includes a complementary function (_show\_example\_result_) to compare the results got and the official information about the protein. To do that, we just need to execute the function passing it the number of the example (from 1 to 10). As we use the example 3, we execute the following:

```python
import inprotfind as ipf

ipf.show_example_result(example=3)
```

_this function can be also execute with shell scripting by: `inprotfind show_example_result --example 3`_

Let's have a look to the output shows in the terminal:

```python
        example         organism           id                                                          annotation
query_example03 Cotesia_rubecula Crub007500.1 ADAMTS-like protein 3 OS=Homo sapiens OX=9606 GN=ADAMTSL3 PE=1 SV=4
```

We see here that the protein belongs to _Cotesia rubecula_, and according to the annotation, it is an ADAMTS-like protein. In the _best_matches.m8_ file we show that the best match corresponded to a protein madd-4 isoform X1 belonging to _Cotesia glomerata_.


## Summary

We have managed to progress from an unknown protein sequence to an approximation of the protein's nature, its specific origin, and its context in the phylogenetic tree. inprotfind emerges as a very user-friendly tool for the rapid characterization of insect proteins.
## Authors

Carmelo Gómez Martínez
- GitHub: [@carmelogzmz](https://www.github.com/carmelogzmz)
