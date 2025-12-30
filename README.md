# MVQB

This repository implements the parallel maximum (α,β)-quasi biclique algorithm **MVQB** and the enumeration algorithm proposed in our paper. If you are using the code, please cite our paper.

## Compile the code

```
$ g++ -fopenmp -O3 MVQB.cpp -o mvqb
$ g++ -fopenmp -O3 MVQB_Max.cpp -o mvqb_max
$ g++ -fopenmp -O3 MVQB_Enum.cpp -o mvqb_enum
```

It generates executables "mvqb_max" and "mvqb_enum", which correspond to our parallel MVQB algorithm (Section 3.4) and enumeration algorithm (Section 3.5).

## Get datasets

The datasets used in our paper can be found from the links below.

`http://konect.cc`

## Data format

To provide a clearer understanding of the data format, we present a simple illustrate example as shown below.

Note that the three integers in the first line means the number of vertices _n_, the size of _L_ and the number of edges _m_. For the next _m_ lines, the first number is the vertice _i_, and the next numbers are the vertices connect to _i_.

```
11 5 16
0 5 6 7
1 6 7 8 10
2 7 8 9
3 8 9 5
4 9 5 6
5 0 3 4
6 0 1 4
7 0 1 2
8 1 2 3
9 2 3 4
10 1
```

## Run the code

The usage procedure for our program is as follows

You can use `-h ` to display the help menu.

```
$ ./mvqb_max -h
```

The algorithms support parallel execution. You can control the number of threads using OMP_NUM_THREADS.

An example of computing the maximum (α,β)-quasi biclique for our provided example graph '’Example.g‘ is as follows

```
$ OMP_NUM_THREADS=8
./mvqb_max -f "Example.g" -u 2 -v 2 -a 0.7 -b 0.7
```

An example of enumerating/counting all (α,β)-quasi bicliques is as follows

```
$ OMP_NUM_THREADS=8
./mvqb_enum -f "Example.g" -u 2 -v 2 -a 0.7 -b 0.7
```

## Result Analysis

If you are using our provided example above, then you will see the following outputs in mvqb.txt. (other files might be a little different)

```
---------------Example.g a,b,lb_L,lb_R: 0.7 0.7 2 2 --------------
4
Running Time: 0 sec
```

From the obtained empirical results, we can extract the size of maximum (α,β)-quasi biclique is 4.

## Algorithm files

MVQB_Max.cpp: our parallel maximum (α,β)-quasi biclique algorithm in the paper (Section 3.4).

MVQB_Enum.cpp: our parallel enumeration/counting algorithm in the paper (Section 3.5).

MVQB.cpp: The original sequential algorithm.

ENUM.cpp: The baseline-ENUM algorithm in the paper.

MVQB_c.cpp: MVQB without the graph coloring-based preprocessing technique (i.e. MVQB\\C).

MVQB_ubd.cpp: MVQB without using the degree-based upper bound (i.e. MVQB\\UBd).

MVQB_ubm.cpp: MVQB without missing degree-based upper bound (i.e. MVQB\\UBm).

MVQB_ubdm.cpp: MVQB without using either of the upper bounding rules (i.e. MVQB\\UBdm).

MVQB_rrd.cpp: MVQB without the degree-based reduction rule (i.e MVQB\\RRd).

MVQB_rrc.cpp: MVQB without the common neighbor-based reduction rule (i.e MVQB\\RRc).

MVQB_rrdc.cpp: MVQB without using either of the reduction rules (i.e MVQB\\RRdc).

MVQB-br.cpp: MVQB uses the same random branching vertex selection approach as the baseline ENUM.

MVQB-bm.cpp: MVQB selects the vertex with the minimum degree in the candidate set.