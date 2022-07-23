# osmtree

Source code for ICDE 2023 paper submission: "Indexing for Near Sorted Data"

## About 
The repository contains the source code for B+tree and SWARE implementations.
In the current version of the code, both implementations are generic, but the application files that test these index data structures only support integer data type. Also, the application files use the same value for both key and value of each entry. Future extensions of the code will support larger data types. 

Both data structures require a buffer pool allocation while running the code, which can be extended to the required amount if needed to run completely in memory. The buffer pool allocation is given in terms of number of blocks where each block is 4KB. For example, if you use an allocation of 1M blocks, then you are allocating 1M\*4KB = 4GB of memory for the tree data structure.

The SA B+-tree also requires the number of entries its in-memory buffer can hold while executing the program, along with the fill factor to maintain while bulkloading. Each entry is a key-value pair. 
## How to run

### Generating ingestion workload
Use the sortedness data generator from this repo: https://github.com/BU-DiSC/bods to generate ingestion keys (can specify payload size=0 to generate only keys). As mentioned above, the application files use the same value for both key and value of each entry (K,V pair). Note the path to the generated workload. 

### B+ Tree
1. compile using the command "make bplustree"
2. Run the executable file using the following format: 
```c
./test_base_index <ingestion_workload_path> <output_file_name> <buffer_pool_allocation> <K> <L> <#. queries>
```
For example, you would use:
```c
./test_base_index createdata_1000000-elems_100000-K_100000-L_1seed1632764083.dat sample.txt 1000000 100000 100000 200000
```
Here, we are ingesting a workload of 1M entries/keys with K=L=100,000. We are using a buffer pool cache for 1M blocks and are executing 200,000 point queries. The output latencies for both ingestion and point queries are written to "sample.txt".

### SA B+ tree
1. compile using the command "make satree"
2. Run the executable file using the following format: 
```c
./test_satree <ingestion_workload_path> <output_file_name> <buffer_pool_allocation> <K> <L> <#. entries> <swareBuffer allocation> <fill factor %>  <#. queries>
```
For example, you would use:
```c
/test_satree createdata_1000000-elems_100000-K_100000-L_1seed1632764083.dat swaresample.txt 1000000 100000 100000 1000000 10000 95 200000
```
Here, we are ingesting a workload of 1M entries/keys with K=L=100,000. We are using a buffer pool cache for 1M blocks and are executing 200,000 point queries. The output latencies for both ingestion and point queries are written to "swaresample.txt". The in-memory buffer will hold 10,000 entries (1% of 1M) and we are maintaining a fill factor of 95%. 
