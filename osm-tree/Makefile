CC=g++
CFLAGS=-g -std=c++11
DATA_TYPE_FLAG=-DDATA_INT
BASE_INDEX_FLAG=-DBPLUS
BASE_INDEX_FILE=inc/betree.h

BASE_INDEX_TARGET=test_base_index

CREATE_DIR=@mkdir -p tree_dat/ 

all: bplustree osmtree

bplustree: 
	$(CREATE_DIR)
	$(CC) $(CFLAGS) $(BASE_INDEX_FILE) $(BASE_INDEX_TARGET).cpp -o $(BASE_INDEX_TARGET) $(BASE_INDEX_FLAG) -DTIMER

bepsilontree:
	$(CREATE_DIR)
	$(CC) $(CFLAGS) $(BASE_INDEX_FILE) $(BASE_INDEX_TARGET).cpp -o $(BASE_INDEX_TARGET)_bepsilon -DTIMER

osmtree: 
	$(CREATE_DIR)
	$(CC) $(CFLAGS) inc/osmtree.h test_osmtree.cpp -o test_osmtree -DOSMTIMER $(BASE_INDEX_FLAG) -DSPLIT80

clean:
	rm -rf $(BASE_INDEX_TARGET) $(BASE_INDEX_TARGET)_bepsilon test_osmtree