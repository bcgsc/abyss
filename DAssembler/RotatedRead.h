#ifndef ROTATEDREAD_H
#define ROTATEDREAD_H 1

#include <string>
#include <vector>

class RotatedRead{
  public:   
      RotatedRead(const std::string& orig_seq, unsigned count = 1);

      bool operator <(const RotatedRead& x) const
      {
          return seq < x.seq;
      }
      bool operator ==(const RotatedRead& x) const
      {
          return seq == x.seq;
      }
      
      std::string seq;
      std::vector<std::string> rotations;
      unsigned count;
      bool used;
};

#endif //ROTATEDREAD_H
