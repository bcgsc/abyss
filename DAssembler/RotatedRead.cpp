#include <iostream>
#include "RotatedRead.h"
#include "Sequence.h"

using namespace std;

//Constructor
RotatedRead::RotatedRead(const string& orig_seq, unsigned count)
	: seq(orig_seq), count(count)
{
    // Add all rotations of the read to the object
    string appended_seq = orig_seq + '$';
    used = false;
    for(unsigned int i = 0; i < appended_seq.size(); ++i)
    {
        string rotated_read = appended_seq.substr(i,(appended_seq.size()-i)) +
         appended_seq.substr(0,i);
        rotations.push_back(rotated_read);
    }
    
}
