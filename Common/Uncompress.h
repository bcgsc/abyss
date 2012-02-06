#ifndef UNCOMPRESS_H
#define UNCOMPRESS_H 1

bool uncompress_init();

namespace {
const bool uncompressInitialized = uncompress_init();
bool getUncompressInitialized() __attribute__((unused));
bool getUncompressInitialized() { return uncompressInitialized; }
}

#endif
