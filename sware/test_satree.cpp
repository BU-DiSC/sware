#include <iostream>
#include <algorithm>
#include <random>
#include "inc/satree.h"
#include "inc/swareBuffer.h"
#include <unordered_set>

using namespace std;

typedef int key_type;


inline void showProgress(const uint32_t &n, const uint32_t &count);
inline void showProgress(const uint32_t &workload_size, const uint32_t &counter)
{
    if (counter / (workload_size / 100) >= 1)
    {
        for (int i = 0; i < 104; i++)
        {
            std::cout << "\b";
            fflush(stdout);
        }
    }
    for (int i = 0; i < counter / (workload_size / 100); i++)
    {
        std::cout << "=";
        fflush(stdout);
    }
    std::cout << std::setfill(' ') << std::setw(101 - counter / (workload_size / 100));
    std::cout << counter * 100 / workload_size << "%";
    fflush(stdout);
    if (counter == workload_size)
    {
        std::cout << "\n";
        return;
    }
}

int main(int argc, char *argv[])
{

    string input_file = "";
    string output_file = "";
    int l, num_nodes_for_buffer_pool, num_entries_buffer_can_hold;
    int fill_factor_p;

    key_type k, n;

    key_type num_queries;

    if (argc > 9)
    {
        input_file = argv[1];
        output_file = argv[2];
        num_nodes_for_buffer_pool = atoi(argv[3]);
        l = atoi(argv[5]);
        k = atoi(argv[4]);
        n = atoi(argv[6]);

            num_entries_buffer_can_hold = atoi(argv[7]);
        fill_factor_p = atoi(argv[8]);
        if (fill_factor_p < 1 || fill_factor_p > 100)
        {
            cout << "Fill Factor has to be between 1% and 100%" << endl;
            return 0;
        }
        num_queries = atol(argv[9]);
    }
    else
    {
        cout << "Invalid no. of arguments" << endl;
        cout << "Use: ./test_osmtree <input_workload_file> <output_file> <num_nodes_in_buffer_pool_cache> <k> <l> <num_inserts> <num_entries_osmBuffer_can_hold> <fill_factor_percentage> <num_queries>"<<endl;
        return 0;
    }

    cout << "Input file name = " << input_file << endl;
    cout << "Output file name = " << output_file << endl;
    double fill_factor = fill_factor_p / 100.0;
    OsmTree<key_type, key_type> tree("manager", "./tree_dat", 4096, num_nodes_for_buffer_pool, num_entries_buffer_can_hold, fill_factor);

    cout << "Caps = " << tree.getBlocksInMemoryCap() << endl;

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
    // fclose(input_file);
    infile.close();

    key_type num = size / sizeof(key_type);
    cout<<"Num inserts = "<<num<<endl;
    int j = 0;

    cout << "Size of one data element =" << sizeof(data[0]) << endl;

    cout << "\n\t\t********** Inserts **********" << endl;
    // insert into buffer

    key_type progress_counter = 0;
    key_type workload_size = n;
    for (key_type i = 0; i < n; i++, progress_counter++)
    {
        tree.osmInsert(data[i] + 1, data[i] + 1);
        if (n < 100)
            n = 100;
        if (progress_counter % (n / 100) == 0)
        {
            // showProgress(n, progress_counter);
        }
    }

    int cap = tree.getOsmBufCap();
    int t = tree.getOsmBufSize();
    int emp = 0;
    progress_counter = 0;

    int nops = num_queries;

    OsmBufferCounters counters = tree.getBufferCounters();
    cout << "Top inserts = " << top_inserts << endl;
    cout << "Num flushes = " << num_flushes << endl;
    cout << "Total number of sort attempts = " << counters.total_sort << endl;
    cout << "Number of quick sort = " << counters.num_quick << endl;
    cout << "Num adaptive sort = " << counters.num_adaptive << endl;
    cout << "Number of times Adaptive Sort failed = " << counters.num_adaptive_fail << endl;
#ifdef OSMTIMER
    cout << "Time in nanoseconds for insert = " << tree.osmTimer.insert_time << endl;
#endif

    cout << "\n\t\t********** Point QUERY **********" << endl;
    key_type x = 0;
    progress_counter = 0;
    for (int i = 0; i < nops; i++, progress_counter++)
    {
        uint query_index = (rand() % n) + 1;
        bool flag = tree.osmQuery(query_index);
        if (!flag)
        {
            x++;
        }
        if (n < 100)
            n = 100;
        if (progress_counter % (nops / 100) == 0)
        {
            // showProgress(nops, progress_counter);
        }
    }
    cout << "Not found keys = " << x << endl;

#ifdef OSMTIMER
    cout << "Time in nanoseconds for query = " << tree.osmTimer.point_query_time << endl;
#endif

    counters = tree.getBufferCounters();

    cout << "----------------- \t Important Query Metrics \t -----------------" << endl;

    cout << "#. of Point Lookups = " << nops << endl;
    cout << "#. of OSM Fence Pointer Scans = " << counters.osm_fence_queries << endl;
    cout << "#. times of OSM Fence Pointer returned yes (#. of times query went into buffer) = " << counters.osm_fence_positive << endl;
    cout << "#. of times Unsorted Section Zonemap queried = " << counters.unsorted_zonemap_queries << endl;
    cout << "#. of times unsorted section zonemap answered yes = " << counters.unsorted_zonemap_positive << endl;
    cout << "#. of Sequential Scans = " << counters.num_seq_scan << endl;
    cout << "#. of Global Bloom filter Queries = " << counters.global_bf_queried << endl;
    cout << "#. times Global Bloom Filter returned yes = " << counters.global_bf_positive << endl;
    cout << "#. of Zones Queried = " << counters.num_zones_queried << endl;
    cout << "#. of Zones returned yes = " << counters.num_zones_positive << endl;
    cout << "#. of times Sublevel BFs queried = " << counters.sublevel_bf_queried << endl;
    cout << "#. of times Sublevel BFs returned yes = " << counters.sublevel_bf_positive << endl;
    cout << "#. of Pages Scanned through sequential search = " << counters.seq_zone_queries << endl;
    cout << "#. of pages with positive result = " << counters.seq_zone_positive << endl;
    cout << "#. of Keys found in seq scan of unsorted section = " << counters.num_seq_found << endl;
    cout << "#. of interpolation search queries (post sorting subblocks) in unsorted section = " << counters.subblock_interpolation_queries << endl;
    cout << "#. of keys found through interpolation search in unsorted section = " << counters.subblock_interpolation_positive << endl;
    cout << "#. of keys found in the unsorted section of the buffer (section includes subsorted blocks for queries) = " << counters.num_unsorted_positive << endl;
    cout << endl;

    cout << "#. of times sorted section zonemap queried = " << counters.sorted_zonemap_queries << endl;
    cout << "#. of times sorted section zonemap answered yes = " << counters.sorted_zonemap_positive << endl;
    cout << "#. of times interpolation search performed = " << counters.num_bin_scan << endl;
    cout << "#. of times interpolation search yielded positive result (number of keys found in sorted section) =" << counters.num_bin_found << endl;
    cout << endl;

    cout << "#. of times Tree's Zonemap was queried = " << tree_zone_queries << endl;
    cout << "#. of times tree's zonemap returned yes = " << tree_zone_positive << endl;
    cout << "#. of times Tree was queried = " << num_tree_scan << endl;
    cout << "#. of times tree scan returned positive (#. of keys found in tree) = " << num_tree_found << endl
         << endl;

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

    #ifdef OSMTIMER
    outfile<<k<<","<<l<<","<<n<<","<<tree.osmTimer.insert_time<<","<<num_queries<<","<<tree.osmTimer.point_query_time<<endl;
    #endif
    outfile.close();
    delete[] data;
    return 0;
}