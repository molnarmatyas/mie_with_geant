#ifndef mieNormaNumericReader_h
#define mieNormaNumericReader_h 1

#include <random>
#include <vector>

class mieNormaNumericReader 
{
  public:
    mieNormaNumericReader();
    ~mieNormaNumericReader();
    void setXSect(const char *);
    std::mt19937 fGen;
    std::discrete_distribution<int> fDist;
    double generate(std::vector<double> &);
  
    std::vector<double> fMieXSect;
    std::vector<double> fTheta;
    std::vector<double> fWeights;
  private:
};

#endif