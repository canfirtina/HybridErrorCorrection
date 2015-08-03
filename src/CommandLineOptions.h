#ifndef COMMAND_LINE_OPTIONS_H_
#define COMMAND_LINE_OPTIONS_H_

#include <iostream>
#include <seqan/arg_parse.h>

struct CommandLineOptions
{
	seqan::CharString longInputFile;
	seqan::CharString shortInputFile;
	unsigned filter1KmerSize;
	unsigned filter2KmerSize;
	unsigned hashCount;
	double filter1JaccardThreshold;
	double filter2JaccardThreshold;
	unsigned pieceLength;
	unsigned iterationLength;
	unsigned numberOfThreads;

	CommandLineOptions() :
		filter1KmerSize(10), filter2KmerSize(6), hashCount(100), filter1JaccardThreshold(0.001), filter2JaccardThreshold(0.1), 
		pieceLength(76), iterationLength(38), numberOfThreads(16)
	{}
};

seqan::ArgumentParser::ParseResult
parseCommandLine(CommandLineOptions & options, int argc, char const ** argv)
{
	using namespace std;
	seqan::ArgumentParser parser("error_correction_filters");

	addOption(parser, seqan::ArgParseOption("LI", "longInputFile", "fast{a,q} file which contains long reads.", seqan::ArgParseArgument::INPUT_FILE, "FILE"));
	addOption(parser, seqan::ArgParseOption("SI", "shortInputFile", "fast{a,q} file which contains short reads.", seqan::ArgParseArgument::INPUT_FILE, "FILE"));
	
	addOption(parser, seqan::ArgParseOption("k1", "filter1KmerSize", "k-mer size for 1st filter.", seqan::ArgParseArgument::INTEGER, "INT"));
	addOption(parser, seqan::ArgParseOption("k2", "filter2KmerSize", "k-mer size for 2nd filter.", seqan::ArgParseArgument::INTEGER, "INT"));
	addOption(parser, seqan::ArgParseOption("hashc", "hashCount", "Number of hash functions for 1st filter.", seqan::ArgParseArgument::INTEGER, "INT"));

	addOption(parser, seqan::ArgParseOption("j1", "filter1JaccardThreshold", "Jaccard threshold for 1st filter.", seqan::ArgParseArgument::DOUBLE, "DOUBLE"));
	addOption(parser, seqan::ArgParseOption("j2", "filter2JaccardThreshold", "Jaccard threshold for 2nd filter.", seqan::ArgParseArgument::DOUBLE, "DOUBLE"));

	addOption(parser, seqan::ArgParseOption("p", "pieceLength", "Piece length.", seqan::ArgParseArgument::INTEGER, "INT"));
	addOption(parser, seqan::ArgParseOption("i", "iterationLength", "Iteration length.", seqan::ArgParseArgument::INTEGER, "INT"));
	addOption(parser, seqan::ArgParseOption("nt", "numberOfThreads", "Number of threads to be used.", seqan::ArgParseArgument::INTEGER));

	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;
	
	getOptionValue(options.longInputFile, parser, "longInputFile");
	getOptionValue(options.shortInputFile, parser, "shortInputFile");
	getOptionValue(options.filter1KmerSize, parser, "filter1KmerSize");
	getOptionValue(options.filter2KmerSize, parser, "filter2KmerSize");
	getOptionValue(options.hashCount, parser, "hashCount");
	getOptionValue(options.filter1JaccardThreshold, parser, "filter1JaccardThreshold");
	getOptionValue(options.filter2JaccardThreshold, parser, "filter2JaccardThreshold");
	getOptionValue(options.pieceLength, parser, "pieceLength");
	getOptionValue(options.iterationLength, parser, "iterationLength");
	getOptionValue(options.numberOfThreads, parser, "numberOfThreads");
	
	return seqan::ArgumentParser::PARSE_OK;
}

#endif
