# Mass Directed Acyclic Word Graph (DAWG) (Python)
A word graph made for the sole purpose of identify mass spectra data

## Background
### What is mass spectrometry data?
If you hadn't heard of mass spectrometry, I doubt you will get much out of this datastructure, but feel free to stick around and learn!

Mass spectrometry data is data that is generated from a mass spectrometer (duh). More specifically, this data structure is meant for identify peptides (short sequences of amino acids) from the list of numbers (floats) from a mass spectrometer run. Amino acids are described as a sequence masses. 

The primary use of this data structure, to be even MORE specific, is by using singly and doubly charged masses to identify a peptide. The output, at its most basic form, is a list of floating point numbers. It might look something like the following

```
[99.023, 140.743, 209.887, 288.115, 402.778]
```

This sequence of floating point numbers describe a particular sequence of amino acids. The way this works is as follows: 

Say our sequence of amino acids is:  `MALW` (don't shoot me I know the masses don't make sense). The above sequence of amino acids describe any arrangement of ion types and charges. For the sake of this example, we will limit our scope to what're called `b` and `y` ions with possible charges of `1` and `2`. The `b` ions describe amino acids in a left-to-right fashion, and `y` right-to-left. So for our example sequence, we can break it down to something like:

```
  b1  b2  b3
M | A | L | W
  y3  y2  y1
```

So `M` is described by both a `b` ion and a `y` ion. The difference is that `b` is the mass of `M` and `y` is really the mass of `ALW`. But you can see how these ions complement eachother. 

To dive a bit deeper, charges, as the name suggests, are the actual charges of the molecules and amino acid chains from the mass spectrometer. The most common of which are `1` and `2`. So we could have up to (in this case) 4 different combinations of ions (`b+, b++, y+, y++`) that describe each junction. 

### What is a DAWG?
A Directed Acyclic Word Graph (DAWG) is an extension of a [prefix tree](https://en.wikipedia.org/wiki/Trie) Prefix trees are a way to store a lot of sequencial data (typically strings) in an easy to search and somewhat compressed form. DAWGs are a [graph](https://en.wikipedia.org/wiki/Graph_(abstract_data_type)) extension of these. They do the same sort of thing, but in a more compressed manner. Nodes with the same value are merged together to reduce the space usage. The wikipedia page for these data structures can be found [here](https://en.wikipedia.org/wiki/Deterministic_acyclic_finite_state_automaton). 

### Great, but wheres the connection?
Mass spectrometry data, at its core, is sequential data. There are two mainstream ways to identify peptides (amino acid sequences) from mass spec data: (1) a database serach, (2) de novo sequencing. The application this data structure is geared towards is database searches. 

Protein databases are typically saved in `.fasta` files, which contain a little bit of header data and then the sequence. Entries look like the following: 

```
>sp|someIdentifier|someName some extra info
ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWX
YZABCDEFGHIJKLMNOPQRSTUVWXYZ...
```

At a very high level, database searches take a spectrum (the list of numbers) like the one shown above and go over each entry and each position in the database, generate a theoretical spectrum from some subsequence of the protein, and compare the two to see how well they match. 

So for our database, if we wanted to try to identify the sequence, we would essentially go through each position of each protein and generate a spectrum. So we would start with `ABCD...`, then move onto `BCDE...` and so on and so forth. 

Going through a database like this without some insane optimizations for each spectrum would take ages. So thats where the DAWG comes in. 

If we limit our searches to some upper bound of peptide length (say 20), then we can fairly quickly iterate through the entire database, find all sequences who's length is < 21, generate theoretical spectra for these sequences, and put them in a data structure to quickly search. 

#### Example
---
Our set of sequences would be (for now, just `b` ions, so left to right):
```
ABCDEFGHIJKLMNOPQRST
BCDEFGHIJKLMNOPQRSTU
ABCDEFGHIJKLMNOPQRST
...
XYZ
YZ
Z
```
Notice that after `Z` becomes the end letter, our sequences shrink on the left side. Since we are only looking at `b` ions in this example, we need to get all of the left starting sequences. If we had stuck to ONLY doing peptides of length 20, we would never find `Z`, `YZ`, `XYZ`, etc. in our searches. 

Once we have our set of sequences, we can generate theoretical spectra for this finite set. 

__NOTE:__ at this point in time, only singly charged and doubly charged masses are allowed into the DAWG.

This would look something like

```
Sequence: ABCD...
singly charged masses: [100.123, 167.9877, 299.773, ...]
doubly charged masses: [50.0625, 83.4588, 149.9932, ...]
```
We would then add these two charged seqeuences to the DAWG tagged at each step with the corresponding subsequence. So our DAWG at this point would look like

```
root --> {100.123, 50.0625, [A]} --> {167.9877, 83.4588, [AB]} 
                                                ||
                                                ||
                                                \/
                                     {299.773, 149.9932, [ABC]} 
```

Now if we had another sequence `AXC` enter the set with the mass sequences
```
singly: [100.123, 189.9877, 299.773, ...]
doubly: [50.0625, 94.998, 149.9932, ...]
```
we could add it to the graph, and our graph would then look like

```
root --> {100.123, 50.0625, [A]} --> {167.9877, 83.4588, [AB]} 
                ||                              ||
                ||                              ||
                \/                              \/
        {189.9877, 94.998, [AX]} ---> {299.773, 149.9932, [ABC, AXY]} 
```
Our graph is now connected, unlike a tree would have been (AXY would have become its own branch). Its important to note that we now have a bit of ambiguity in which path we took to get to ABC, AXY, but the important step is we've reduced the set of peptides that we care about from what seemed infinite to 2, so some post analysis of these sequences could narrow this down.

---

## Mass DAWG
This module is written in C++ 11 for its speed and ease of use to write. Python bindings are what are presented in this package
### Installation and usage

#### Install package from PyPI
```bash
$> pip install mass_dawg
```
#### Install from source
```bash
$>git clone https://github.com/zmcgrath96/mass_DAWG.git
```

__UPDATE__: As of August 4th 2020 (pip release `1.1.0`), sequences can be inserted out of order.

### Example 

#### Example main.py file
```python
from mass_dawg import PyMassDawg

if __name__ == '__main__':

    # create the mass dawg
    md = PyMassDawg()

    # add our first sequence
    singly = [200.2, 400.4, 600.6, 800.8]
    doubly = [100.1, 200.2, 300.3, 400.4]
    sequence = 'ABCD'

    # add it to the graph
    print('inserting ABCD into graph')
    md.insert(singly, doubly, sequence)

    # show the graph
    md.show()

    # add a new series to it. Graph should reconverge on "D"
    singly = [200.2, 400.4, 600.6, 800.8]
    doubly = [100.1, 200.4, 300.6, 400.4]
    sequence = "AXYD"

    # add the new one to the graph
    print('inserting AXYD into graph')
    md.insert(singly, doubly, sequence)

    # finishing the graph does a final round
    # of merging nodes with the same value
    # for this example, it'll turn the nodes with masses (400.4, 800.8)
    # into 1 node with kmers "ABCD", "AXYD"
    md.finish()

    md.show()

    # search for a seqeunce missing 2 masses
    print('searching for a sequence with gaps')
    searching = [600.6, 800.8]
    results = md.fuzzy_search(searching, 2, 10)
    print(f'Number of results: {len(results)}')

    for result in results:
        print(result)
```

#### Example output

```
inserting ABCD into graph
root
  |---> kmers: {A} 	 masses: 100.099998, 200.199997
    |---> kmers: {AB} 	 masses: 200.199997, 400.399994
      |---> kmers: {ABC} 	 masses: 300.299988, 600.599976
        |---> kmers: {ABCD} 	 masses: 400.399994, 800.799988
inserting AXYD into graph
root
  |---> kmers: {A} 	 masses: 100.099998, 200.199997
    |---> kmers: {AB} 	 masses: 200.199997, 400.399994
      |---> kmers: {ABC} 	 masses: 300.299988, 600.599976
        |---> kmers: {ABCD, AXYD} 	 masses: 400.399994, 800.799988
    |---> kmers: {AX} 	 masses: 200.399994, 400.399994
      |---> kmers: {AXY} 	 masses: 300.600006, 600.599976
        |---> kmers: {ABCD, AXYD} 	 masses: 400.399994, 800.799988

searching for a sequence with gaps
Number of results: 2
ABCD
AXYD
```

### Exposed MassDawg functions (API)
* __show()__: Print the graph to the console as a tree (merged nodes have their kmers put into a list)
* __insert(singly_sequence: list, doubly_sequence: list, kmer: str) -> None__: Insert a pair of singly charged and doubly charged masses into the dawg associated withthe kmer (all 3 parameters MUST be the same length)
* __fuzzy_search(sequence: list, gap_allowance: int, ppm_tol: int) -> None__: Search the graph for a sequence of floats allowing for up to gapAllowance missed masses. ppmTol is the allowed tolerance for a mass to fall within (ppm = parts per million)
* __finish()__: Go through the graph one final time to merge all remaining nodes that have not been checked for duplicates
