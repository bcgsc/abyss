#ifndef ROTATION_H
#define ROTATION_H 1

#include <string>
#include <vector>

using namespace std;

class Rotation {
  
  public:
      Rotation(const string& orig_seq)
		  : seq(orig_seq), overlap(0) { }
      
      bool operator <(const Rotation& x) const
      {
            return seq < x.seq;
      }
      bool operator ==(const Rotation& x) const
      {
            return seq == x.seq;
      }
       
      string seq;
      unsigned overlap;  
};


#endif //ROTATION_H
