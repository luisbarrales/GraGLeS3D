#ifndef _applic_h_
#define _applic_h_

#define DEBUG


#define _STR(x) _VAL(x)
#define _VAL(x) #x
#define ASSERT(cond) if(!(cond)) exitus("failed assertion:"__FILE__"line"_STR(__LINE__)":"#cond)
#define TOLFP(x) ((1.0)-(0.99999/((double) x))) //tolerance
#define ERRTXT(text) (text" : file: "__FILE__" line:"_STR(__LINE__)"\n")

#define NOW true
#define NOTYET false

#define NEVER true


#ifndef DEBUG
        #define QUICKASSERT(cond)       ((void)0)
#else
        #define QUICKASSERT(cond)       ASSERT(cond)
#endif

void exitus (const char *s);
#endif
