#ifndef FMINDEX_H
#define FMINDEX_H 1

#include "wat_array.h"
#include <istream>
#include <map>
#include <ostream>
#include <stdint.h>
#include <vector>

/** An FM index. */
class FMIndex
{
  public:
	int load(std::istream& is);
	int save(std::ostream& os);
	int read(const char *fname, std::vector<uint8_t> &s);
	int readQstring(const char *fname,
			std::vector<std::vector<uint8_t> > &qs);
	int buildFmIndex(const char *fnmae, int _percent);
	void search(const std::vector<uint8_t> &qs,
			std::pair<uint64_t, uint64_t> &res);
	void searchHamming(const std::vector<uint8_t> &qs,
			int dist,
			std::vector<std::pair<uint64_t, uint64_t> > &res);
	void searchEdit(const std::vector<uint8_t> &qs, int dist,
			std::vector<std::pair<uint64_t, uint64_t> > &res);
	size_t locate(uint64_t i);

  private:
	void calculateStatistics(const std::vector<uint8_t> &s);
	int buildBWT(const std::vector<uint8_t> &s,
			const std::vector<uint32_t> &sa,
			std::vector<uint64_t> &bwt);
	int buildSA(const std::vector<uint8_t> &s,
			std::vector<uint32_t> &sa);
	int buildSampledSA(const std::vector<uint8_t> &s,
			const std::vector<uint32_t> &sa);
	int buildWaveletTree(const std::vector<uint64_t> &bwt);
	void searchHamming(const std::vector<uint8_t> &qs,
			int i, size_t sp, size_t ep, int d, int dist,
			std::vector<std::pair<uint64_t, uint64_t> > &res);
	void searchEdit(const std::vector<uint8_t> &qs, int count,
			int sp, int ep, int dist,
			std::vector<std::pair<uint64_t, uint64_t> > &res);

	int percent;
	uint8_t alphaSize;
	std::vector<uint32_t> cf;
	std::map<unsigned char, uint8_t> mapping;
	std::map<uint8_t, unsigned char> rmapping;
	std::vector<uint32_t> sampledSA;
	std::vector<std::vector<int> > cols;
	wat_array::WatArray wa;
};

#endif
