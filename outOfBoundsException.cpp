#include"outOfBoundsException.h"

outOfBoundsException::outOfBoundsException(int d, int i) : d(d), i(i) {};
string outOfBoundsException::what() {
    ostringstream oss;
    oss << d << "," << i;
    return oss.str();
}
int outOfBoundsException::get_i() { return i; }
