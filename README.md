# DistanceSCAN
Code Contributor: Kaixin Liu

If you have any questions, feel free to contact me. My email is lkx17@mails.tsinghua.edu.cn.

Please cite our paper if you choose to use our code. 

```
@inproceedings{DBLP:***,
  author    = {Kaixin Liu and
               Sibo Wang and
               Yong Zhang and
               Chunxiao Xing},
  title     = {An Efficient Algorithm for Distance-based Structural Graph Clustering},
  journal   = {PACMMOD},
  volume    = {1},
  number    = {45},
  pages     = {**--**},
  year      = {2023},
  doi       = {10.1145/3588725},
}
```
## Tested Environment
- Ubuntu
- C++ 14
- GCC 4.8
- Boost
- cmake

## Compile
```sh
$ cmake .
$ make
```

## Parameters
### Construct Sketches

```sh
./Distance_SCAN_SIGMOD --operation construct_sketches --algo <algorithm> [options]
```

- algo: the type of all-distances bottom-k sketches
    - distancescan_pst: storing all-distances sketches in persist search trees.
    - distancescan: storing all-distances sketches by max-heaps.
- options
    - --prefix \<prefix\>
    - --dataset \<dataset\>
    - -d \<the max distance threshold d\>
    - -k \<the number of samples in bottom-k sketches\>

- Example

```sh
$ ./Distance_SCAN_SIGMOD --operation construct_sketches --dataset ego-facebook --algo distancescan -k 16 -d 0.4 
```

### Query
```sh
./Distance_SCAN_SIGMOD --operation construct_sketches --algo <algorithm> [options]
```


- algo: the algorithm you prefer to run.
    - scan: SCAN for the distance-SCAN problem.
    - pscan: pscan for the distance-SCAN problem.
    - exact: the exact algorithm for the distance-SCAN problem.
    - distancescan: our method.
- options
    - --prefix \<prefix\>
    - --dataset \<dataset\>
    - -d \<the distance threshold d\>
    - -m \<the max distance threshold d of sketches\>
    - -k \<the number of samples in bottom-k sketches\>
    - -u \<the threshold of structurally similar neighbors $\mu$\>
    - -e \<the similarity threshold $\epsilon$ \>

- Example

```sh
$ ./Distance_SCAN_SIGMOD --operation query --dataset ego-facebook --algo distancescan -k 16 -u 5 -e 0.2 -d 0.3 -m 0.4 
```

## Data
You can download from https://snap.stanford.edu/data/,  http://law.di.unimi.it/datasets.php and https://www.aminer.cn/data/?nav=openData#Topic-coauthor.





