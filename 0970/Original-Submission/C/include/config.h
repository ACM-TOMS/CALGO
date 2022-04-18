
#if defined(__cplusplus)
extern "C" {
#endif

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS 1
#endif

#ifndef _CONFIG_H_
#define	_CONFIG_H_

// New defines
//#define VERIFY_RESULTS 1
//#define SPEED 1
//#define FILE_OUTPUT 1
#define KS 1

// Original stuff

//#define	WINDOWS32
//#define	PROTOTYPES
//#define	LITTLE_ENDIAN
//#define	LOWHI

/*
 * AUTO DEFINES (DON'T TOUCH!)
 */

#ifndef	CSTRTD
typedef char *CSTRTD;
#endif
#ifndef	BSTRTD
typedef unsigned char *BSTRTD;
#endif

#ifndef	BYTE
typedef unsigned char BYTE;
#endif
#ifndef	UINT
typedef unsigned int UINT;
#endif
#ifndef	USHORT
typedef unsigned short USHORT;
#endif
#ifndef	ULONG
typedef unsigned long ULONG;
#endif
#ifndef	DIGIT
typedef USHORT DIGIT;	/* 16-bit word */
#endif
#ifndef	DBLWORD
typedef ULONG DBLWORD;  /* 32-bit word */
#endif

#ifndef	WORD64
typedef ULONG WORD64[2];  /* 64-bit word */
#endif

#endif /* _CONFIG_H_ */

#if defined(__cplusplus)
}
#endif
