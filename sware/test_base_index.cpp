#include <iostream>
#include "inc/betree.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <algorithm>
#include <chrono>

using namespace std;

typedef int key_type;


int main(int argc, char *argv[])
{

    cout << "**************************** NEW ANALYSIS ****************************" << endl;

    string input_file = "";
    string output_file = "";

    int num_nodes_for_buffer_pool, l;

    key_type k;
    key_type num_queries;
    if (argc > 6)
    {
        input_file = argv[1];
        output_file = argv[2];
        num_nodes_for_buffer_pool = atoi(argv[3]);
        l = atoi(argv[5]);
        k = atoi(argv[4]);
        num_queries = atol(argv[6]);
    }
    else
    {
        cout << "Insufficient arguments to the program!" << endl;
        cout << "Use: ./test_base_index <input_workload_file> <output_file> <num_nodes_in_buffer_pool> <k> <l> <num_queries>";
        return 0;
    }

    cout << "Input file name = " << input_file << endl;
    cout << "Output file name = " << output_file << endl;
    BeTree<key_type, key_type> tree("manager", "./tree_dat", 4096, num_nodes_for_buffer_pool);
    cout << "Caps = " << tree.getBlocksInMemoryCap() << endl;
    long cache_capacity = tree.getBlocksInMemoryCap();

    long int size = 0;
    key_type *data;

    ifstream infile(input_file, ios::in | ios::binary);
    if (!infile)
    {
        cout << "Cannot open file!" << endl;
        return 0;
    }

    ofstream outfile(output_file, ios::out | ios::app);
    if (!outfile)
    {
        cout << "Cannot open output file!" << endl;
        return 0;
    }

    FILE *file = fopen(input_file.c_str(), "rb");
    if (file == NULL)
        return 0;

    fseek(file, 0, SEEK_END);
    size = ftell(file);
    fclose(file);

    cout << "size = " << size << endl;

    data = new key_type[size / sizeof(key_type)];
    infile.read((char *)data, size);
    infile.close();

    key_type n = size / sizeof(key_type);
    int j = 0;

    cout << "Size of one data element =" << sizeof(data[0]) << endl;
    cout<<"Total inserts = "<<n<<endl;
    cout << "\t\t********** INSERTS **********" << endl;
    for (key_type i = 0; i < n; i++)
    {
        tree.insert(data[i] + 1, data[i] + 1);
    }

#ifdef TIMER
    cout << "Time taken for inserting" << n << " elements = " << tree.timer.insert_time << " nanoseconds" << endl;
    long double ins_time = tree.timer.insert_time / (1e+6);
    long double bulk_time = 0;
#endif

  
    key_type key = 0;
    key_type nops = num_queries;
  
    cout << "\n\t\t********** POINT-QUERY **********" << endl;

    for (int i = 0; i < nops; i++)
    {
        uint query_index = (rand() % n) + 1;
        if (!tree.query(query_index))
        {
            key++;
        }
    }

    cout << "Not found Keys = " << key << endl;
#ifdef TIMER
    cout << "Time taken for point-querying " << nops << " elements = " << tree.timer.point_query_time << " nanoseconds" << endl;
    long double pt_time = tree.timer.point_query_time / (1e+6);
#endif
  

    cout << "================== Final STATS ====================" << endl;
    tree.fanout();
    int av_fanout, max_fanout, min_fanout, med_fanout, num_int, num_leafs;
    av_fanout = tree.traits.average_fanout;
    max_fanout = tree.traits.max_fanout;
    min_fanout = tree.traits.min_fanout;
    med_fanout = tree.traits.median_fanout;
    num_int = tree.traits.num_internal_nodes;
    num_leafs = tree.traits.num_leaf_nodes;
    cout << "Average Fanout = " << tree.traits.average_fanout << endl;
    cout << "Max Fanout = " << tree.traits.max_fanout << endl;
    cout << "Min Fanout = " << tree.traits.min_fanout << endl;
    cout << "Median Fanout = " << tree.traits.median_fanout << endl;
    cout << "Num Internal Nodes, Leaf Nodes = " << tree.traits.num_internal_nodes << "," << tree.traits.num_leaf_nodes << endl;
    cout << "Internal Splits = " << tree.traits.internal_splits << endl;

    int depth = tree.depth();
    cout << "Depth of tree = " << depth << endl;


    #ifdef TIMER
    outfile<<k<<","<<l<<","<<n<<","<<tree.timer.insert_time<<","<<num_queries<<","<<tree.timer.point_query_time<<endl;
    #endif

    outfile.close();
    // delete[] elements;
    delete[] data;

    return 0;
}