#ifndef UNCOMPRESS_H
#define UNCOMPRESS_H 1

bool uncompress_init();

/** Call the initalizer to force the Uncompress.o object file to be
 * linked in.
 */
const bool init_uncompress = uncompress_init();

#endif
