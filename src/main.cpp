#include <cstring>
#include <iterator>
#include <climits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <math.h>
#include <stdlib.h>
#include <thread>
#include <algorithm>

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/align.h>
#include <seqan/seq_io.h>

#include "MurmurHash3.h"
#include "CommandLineOptions.h"

using namespace std;
using namespace seqan;

/* Short Reads*/
StringSet<CharString> shortReadIds;
StringSet<CharString> shortReadSeqs;

/* Long Reads*/
StringSet<CharString> longReadIds;
StringSet<CharString> longReadSeqs;

vector< string> shortReads; // [read_id]
vector< string> longReads;  // [read_id]

vector< vector< vector<uint> > > shortMurmurhashFingerprints; // [hash_id][read_id][kmer_id]
vector< vector< vector<uint> > > longMurmurhashFingerprints;  // [hash_id][read_id][kmer_id]

vector < map < uint, vector < pair<int, int> > > > shortMinmers; // [hash_id]<hash_value -> [<read_id,piece_id>]>
vector < map < uint, vector < pair<int, int> > > > longMinmers;  // [hash_id]<hash_value -> [<read_id,piece_id>]>

vector < map < int, vector < int > > > filteredPairsF1;
vector < map < int, vector < int > > > filteredPairsF2;;

template <typename Tfilename, typename Tids, typename Tseqs>
void loadReads(Tfilename & fileLocation, Tids & ids, Tseqs & seqs, vector<string> & sequencingReads) {
	SeqFileIn seqFileIn;
	if (!open(seqFileIn, toCString(fileLocation))) {
		std::cerr << "ERROR: Could not open the file.\n";
		return;
	}
	try {
		readRecords(ids, seqs, seqFileIn);
	} catch (Exception const & e) {
		std::cout << "ERROR: " << e.what() << std::endl;
		return;
	}
	
	cout << length(seqs) << " sequencing reads are loaded from " << fileLocation << "." << endl;
	sequencingReads.clear();
	for (unsigned i = 0; i < length(seqs); ++i)
		sequencingReads.push_back(std::string(toCString(seqs[i])));
}

void hashKMersOfReadWithMurmurHash( string &read, vector<uint> &hashedKMers, int seed, int k){

	//output of murmurhash3 function
	uint32_t murmurHashOut;

	//for each kmer of the current read
	for( int i = 0; i < (read.size() - k + 1); i++){
		MurmurHash3_x86_32( read.c_str() + i, k, seed, &murmurHashOut);
		hashedKMers[i] = murmurHashOut;
	}
}

void murmurHashFingerprintsWithSeed (	vector < string > *reads,
					vector < vector < vector < uint > > > *murmurHashFingerprints, 
					int startIndex, int endIndex, int k){
	int readCount = reads->size();
	
	//for each read for current hash function
	for (int seed = startIndex; seed <= endIndex; ++seed)
		for (int readIndex = 0; readIndex < readCount; ++readIndex)
			hashKMersOfReadWithMurmurHash((*reads)[readIndex], (*murmurHashFingerprints)[seed][readIndex], seed, k);
}

//This functions creates h-many hash functions using MurmurHash3 to generate
//fingerprints for each read. Therefore it returns the vector that consists
//of the reads for which h-many fingerprints created.
//You can reach any fingerprint by this indexing order:
//hManyMurmurHashFingerPrints[hthHashFunction][nthRead][kthFingerPrint]
//parameter h defines how many fingerprint function will be created.
void createHMurmurHashFingerPrints(	vector< string > & reads,
					vector< vector< vector<uint> > > & hManyMurmurHashFingerprints, 
					int k, int hashCount, int numberOfThreads){
	if (numberOfThreads < 1) 
		numberOfThreads = 1;
	int readCount = reads.size();
	int batchSize = (hashCount + numberOfThreads - 1) / numberOfThreads;

	//creating h-many hash functions
	hManyMurmurHashFingerprints.clear();
	hManyMurmurHashFingerprints.resize(hashCount);

	for (int h = 0; h < hashCount; ++h){

		//total number of reads whose kmers to be hashed
		hManyMurmurHashFingerprints[h].resize(readCount);

		for( int readIndex = 0; readIndex < readCount; ++readIndex)
			//total number of kmers count for the current read
			hManyMurmurHashFingerprints[h][readIndex].resize(reads[readIndex].size() - k + 1);
	}
	
	vector< std::thread> threads;
	//for each hash function
	for (int t = 0; t < numberOfThreads; ++t) {
		int startIndex = t * batchSize;
		int endIndex = (t + 1) * batchSize - 1;
		if (endIndex >= hashCount) endIndex = hashCount - 1;
		threads.push_back(thread(murmurHashFingerprintsWithSeed, &reads, &hManyMurmurHashFingerprints, startIndex, endIndex, k));
	}
	for (auto& th : threads) 
		th.join();
}

//this functions takes vector of fingerprints of kmers of each read. It finds minimum mer in each set of
//fingerprints of reads. Then creates sketch consisting of h many minmer for each read. Those minmers are
//mapped into dictionary. Where a value (minmer) in the dictionary gives all the read indices that has
//the value as its min-mer. Therefore, it will be efficient to find all the reads that share same minmer.
//If specified pieceLength is shorter than the actual lengths of reads, minmer maps are created for each
//piece. The output vector holding those maps are the pieces of the read.
void createMinmers( 	vector < vector < vector < uint > > > & fingerprints,
			vector < map < uint, vector <  pair < int, int > > > > & minmers, 
			int pieceLength, int iterationLength, int kmerLength, int hashCount){
 
	minmers.clear();
	minmers.resize(hashCount);

	for (int h = 0; h < hashCount; ++h) {
		if ((h % 10) == 0) cout << "h="<<h<<" finished." << endl;
		minmers[h] = map <uint, vector < pair < int, int > > >();

		for (int r = 0; r < fingerprints[h].size(); ++r) {
			int readLength = fingerprints[h][r].size() + kmerLength - 1;
			int pieceCount = 1 + (readLength - pieceLength) / iterationLength;
			int kmerCountPerPiece = pieceLength - kmerLength + 1;

			for (int p = 0; p < pieceCount; ++p) {
				uint curMin = UINT_MAX;
				for (int k = p * iterationLength; k < (p * iterationLength + kmerCountPerPiece); ++k)
					if (fingerprints[h][r][k] < curMin)
						curMin = fingerprints[h][r][k];
				if (minmers[h].find(curMin) == minmers[h].end())
					minmers[h][curMin] = vector < pair<int, int> > (1, make_pair(r, p));
				else
					minmers[h][curMin].push_back(make_pair(r, p));
			}
		}
	}
}

void firstFilterElimination (	vector < map < uint, vector < pair < int, int > > > > & shortMinmers,
				vector < map < uint, vector < pair < int, int > > > > & longMinmers,
				vector < map < int, vector < int > > > & filteredPairs,
				int shortReadCount, int longReadCount, int hashCount, double threshold) {
	filteredPairs.clear();
	filteredPairs.resize(shortReadCount);
	for (int i = 0; i < shortReadCount; ++i)
		filteredPairs[i] = map < int, vector < int > >();

	map <pair<pair<int, int>, pair<int, int> >, int> countTable; // <<shortReadID, shortReadPieceID>, <longReadID, longReadPieceID>> -> count;
	countTable.clear();

	map < uint, vector < pair < int, int > > >::iterator itShort;
	for (int h = 0; h < hashCount; h++) {
		if ((h % 10) == 0) cout << "h="<<h<<" finished. countTable.size()=" << countTable.size() << endl;
		// iterate through every minmer value of shortRead dataset and find them among longRead minmers 
		for (itShort = shortMinmers[h].begin(); itShort != shortMinmers[h].end(); ++itShort) {
			uint minmerValue = itShort->first;
			vector < pair < int, int > > & shortList = itShort->second;
			
			if (longMinmers[h].find(minmerValue) != longMinmers[h].end()) {
				vector < pair < int, int > > & longList = longMinmers[h][minmerValue];

				for (int i = 0; i < shortList.size(); ++i)
					for (int j = 0; j < longList.size(); ++j)
						countTable[make_pair(shortList[i], longList[j])]++;
			}
		}
	}
	
	map < pair < pair < int, int>, pair < int, int> >, int >::iterator itCountTable;
	for (itCountTable = countTable.begin(); itCountTable != countTable.end(); ++itCountTable) {
		int count = itCountTable->second;
		int shortReadID = itCountTable->first.first.first;
		int longReadID = itCountTable->first.second.first;
		int longReadPieceID = itCountTable->first.second.second;
		if ( (double) count / hashCount >= threshold) {
			if (filteredPairs[shortReadID].find(longReadID) == filteredPairs[shortReadID].end())
				filteredPairs[shortReadID][longReadID] = vector < int >(1, longReadPieceID);
			else
				filteredPairs[shortReadID][longReadID].push_back(longReadPieceID);
		}
	}	
}

//TODO: rewrite this method
void secondFilterElimination( vector< map<int, vector<int> > > &pairsPassedFirstFilter,
                             vector< string> &shortReads, vector< string> &longReads,
                             vector< map<int, vector<int> > > &pairsPassedSecondFilter,
                             int pieceLength, int iterationLength, int kmerlength, float threshold){
    
    int numberOfFirstReadPieces;
    int numberOfKmersForEachReadPiece;
    pairsPassedSecondFilter.resize(pairsPassedFirstFilter.size());
    
    
    
    for( int firstRead = 0; firstRead < shortReads.size(); firstRead++){
        
        vector< uint32_t> firstReadFinprints;
        for( int sIndex = 0; sIndex < shortReads[firstRead].size(); sIndex++){
            
            string kmer = shortReads[firstRead].substr(sIndex, kmerlength);
            uint32_t shortReadMurmurOut;
            MurmurHash3_x86_32( kmer.c_str(), kmer.size(), 0, &shortReadMurmurOut);
            
            firstReadFinprints.push_back( shortReadMurmurOut);
        }
        
        numberOfFirstReadPieces = (((int)shortReads[firstRead].size()-pieceLength)/iterationLength) + 1;
        
        //Set of kmer length = length of the string - kmerlength + 1
        numberOfKmersForEachReadPiece =  pieceLength - kmerlength + 1;
        
        //will iterate once if piecelength equals to short read size. so, not worry about this iteration.
        for( int curShortReadPiece = 0; curShortReadPiece < numberOfFirstReadPieces; curShortReadPiece++){
            
            //key = second read
            for( map<int, vector<int> >::iterator curKey = pairsPassedFirstFilter[firstRead].begin(); curKey != pairsPassedFirstFilter[firstRead].end(); ++curKey){
                
                //for each piece passed first filter for that pair (curshortread - second read)
                for( int curLongReadPieceIndex = 0; curLongReadPieceIndex < pairsPassedFirstFilter[firstRead][curKey->first].size(); curLongReadPieceIndex++){
                    
                    int w = 0;
                    
                    //hasMatch indicates whether the current short read kmer on that index has a match
                    //with a long read piece's kmer earlier
                    vector< bool> hasMatch;
                    hasMatch.resize( firstReadFinprints.size(), 0);
                    
                    //take the piece in the long read that passed the first filter with the current short read
                    int curLongReadPiece = pairsPassedFirstFilter[firstRead][curKey->first][curLongReadPieceIndex];
                    for( int sIndex = curLongReadPiece*iterationLength; sIndex < (curLongReadPiece*iterationLength) + numberOfKmersForEachReadPiece; sIndex++){
                        
                        string longReadKMer = longReads[curKey->first].substr(sIndex, kmerlength);
                        uint32_t longReadMurmurOut;
                        MurmurHash3_x86_32( longReadKMer.c_str(), longReadKMer.size(), 0, &longReadMurmurOut);
                        
                        for( int curShortMer = 0; curShortMer < firstReadFinprints.size(); curShortMer++)
                            if( !hasMatch[curShortMer] && firstReadFinprints[curShortMer] == longReadMurmurOut){
                                w++;
                                hasMatch[curShortMer] = true;
                            }
                    }
                    
                    if( (float)w/(float)numberOfKmersForEachReadPiece >= threshold){
                        
                        pairsPassedSecondFilter[firstRead][curKey->first].push_back(pairsPassedFirstFilter[firstRead][curKey->first][curLongReadPieceIndex]);
                        break;
                    }
                }
            }
        }
    }
}

int main(int argc, const char * argv[]) {
	CommandLineOptions options;
	seqan::ArgumentParser::ParseResult parseRes = parseCommandLine(options, argc, argv);
	
	if (parseRes != seqan::ArgumentParser::PARSE_OK)
		return parseRes == seqan::ArgumentParser::PARSE_ERROR;
	
	loadReads (options.shortInputFile, shortReadIds, shortReadSeqs, shortReads);
	loadReads (options.longInputFile, longReadIds, longReadSeqs, longReads); 

	clock_t hashShortStart = clock();
	cout << "Creating fingerprints..." << endl;
	createHMurmurHashFingerPrints(shortReads, shortMurmurhashFingerprints, options.filter1KmerSize, options.hashCount, options.numberOfThreads);
	printf("Time taken for creating fingerprints (short reads): %.2fs\n", (double)(clock() - hashShortStart)/CLOCKS_PER_SEC);
	clock_t hashLongStart = clock();
	createHMurmurHashFingerPrints(longReads, longMurmurhashFingerprints, options.filter1KmerSize, options.hashCount, options.numberOfThreads);
	printf("Time taken for creating fingerprints (long reads): %.2fs\n", (double)(clock() - hashLongStart)/CLOCKS_PER_SEC);
	printf("Time taken for creating fingerprints (short + long reads): %.2fs\n", (double)(clock() - hashShortStart)/CLOCKS_PER_SEC);

	cout << "[INFO] Creating minmers..." << endl;
	clock_t minmerStart = clock();
	// long read length = short read length while creating sketches
	createMinmers(shortMurmurhashFingerprints, shortMinmers, options.pieceLength, options.iterationLength, options.filter1KmerSize, options.hashCount);
	createMinmers(longMurmurhashFingerprints, longMinmers, options.pieceLength, options.iterationLength, options.filter1KmerSize, options.hashCount);
	printf("[INFO] Time taken for creating minmers: %.2fs\n", (double)(clock() - minmerStart)/CLOCKS_PER_SEC);

	shortMurmurhashFingerprints.clear();
	longMurmurhashFingerprints.clear();	

	cout << "[INFO] Applying first filter..." << endl;
	clock_t fStart = clock();
	firstFilterElimination (shortMinmers, longMinmers, filteredPairsF1, shortReads.size(), longReads.size(), options.hashCount, options.filter1JaccardThreshold);
	printf("[INFO] Time taken for first filter: %.2fs\n", (double)(clock() - fStart)/CLOCKS_PER_SEC);
	//TODO: true mappings tend to pass from first filter as 2 or 3 consecutive pieces

	long long counter = 0;
        for (int i = 0; i < shortReads.size(); ++i) {
                counter += filteredPairsF1[i].size();
        }
        cout << "[INFO] after f1: " << counter << endl;

	cout << "[INFO] Applying second filter..." << endl;
	clock_t f2Start = clock();
	secondFilterElimination(filteredPairsF1, shortReads, longReads, filteredPairsF2, 
				options.pieceLength, options.iterationLength, options.filter2KmerSize, options.filter2JaccardThreshold);
	printf("[INFO] Time taken for second filter: %.2fs\n", (double)(clock() - f2Start)/CLOCKS_PER_SEC);
	counter  = 0;
	for (int i = 0; i < shortReads.size(); ++i)
		counter += filteredPairsF2[i].size();
	cout << "[INFO] after f2: " << counter << endl;

	return 0;
}
