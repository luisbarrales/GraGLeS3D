#ifndef OUTOFBOUNDSEXCEPTION_h
#define OUTOFBOUNDSEXCEPTION_h

#include<iostream>
#include <sstream>
using namespace std;

class outOfBoundsException {
    int d,i;
public:
    outOfBoundsException(int d, int i);
    string what();
    int get_i();
};

#endif